// =============================================================================================================================
//> This N-view triangulation code is built on top of the publicly released certifiable solver for multiple views triangulation
//  https://github.com/C-H-Chien/FastNViewTriangulation arises from the paper:
//  Garcia-Salguero, Mercedes, and Javier Gonzalez-Jimenez. "Certifiable solver for real-time N-view triangulation." 
//  IEEE Robotics and Automation Letters 8, no. 4 (2023): 1999-2005.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =============================================================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

#include "mex.h"
#include "../include/NViewsTypes.h"
#include "../include/NViewsUtils.h"
#include "../include/NViewsClass.h"
#include "../include/definitions.h"
#include "../utils/generatePointCloud.h"

#include <eigen3/Eigen/Dense>
#include <Eigen/Eigenvalues>  

using namespace NViewsTrian;

//> structured data for feature tracks
struct Feature_Track {

    unsigned Length;                                //> feature track length = number of features
    std::vector<Eigen::Vector3d> Locations;         //> feature point locations in pixels
    std::vector<Eigen::Matrix3d> Abs_Rots;          //> relative rotations
    std::vector<Eigen::Vector3d> Abs_Transls;       //> relative translations
    std::vector<Eigen::Vector3d> Optimal_Locations; //> corrected feature point locations in pixels

}; //> End of struct Feature_Track

void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {
    
    // if(nl != 5)                                 mexErrMsgTxt("Five inputs are required!");
    if(mxGetClassID(pr[0]) != mxDOUBLE_CLASS)   mexErrMsgTxt("Class of feature point locations must be double.");

    //> Get data from MATLAB
    Feature_Track feature_track_;
    double* Observation_Locations   = (double*) mxGetData(pr[0]);
    feature_track_.Length           = (unsigned) mxGetN(pr[0]);
    double* intrinsic_matrix        = (double*) mxGetData(pr[1]);
    double* rotations               = (double*) mxGetData(pr[2]);
    double* translations            = (double*) mxGetData(pr[3]);
    bool mexFunction_debug          = (bool) mxGetScalar(pr[4]);

    if (mexFunction_debug) 
        mexPrintf("Feature Track Length = %d frames\n", feature_track_.Length);
    
    //> Convert data to the need of Fast NView Triangulation
    Eigen::Matrix3d K;
    K << intrinsic_matrix[0], intrinsic_matrix[1], intrinsic_matrix[2], \
         intrinsic_matrix[3], intrinsic_matrix[4], intrinsic_matrix[5], \
         intrinsic_matrix[6], intrinsic_matrix[7], intrinsic_matrix[8];
    Eigen::Matrix3d inverse_K = K.inverse();
    Eigen::MatrixXd obs_i(3, feature_track_.Length);
    for (int i = 0; i < feature_track_.Length; i++) {
        //> feature locations (point observations)
        Eigen::Vector3d point = Eigen::Vector3d::Zero();
        point << Observation_Locations[(i)*2 + 0], Observation_Locations[(i)*2 + 1], 1.0;
        obs_i.col(i) = point;
        feature_track_.Locations.push_back(point);

        //> absolute rotation matrices
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        R << rotations[(i)*9 + 0], rotations[(i)*9 + 1], rotations[(i)*9 + 2], \
             rotations[(i)*9 + 3], rotations[(i)*9 + 4], rotations[(i)*9 + 5], \
             rotations[(i)*9 + 6], rotations[(i)*9 + 7], rotations[(i)*9 + 8];
        feature_track_.Abs_Rots.push_back(R);

        //> absolute translation vectors
        Eigen::Vector3d t = Eigen::Vector3d::Zero();
        t << translations[(i)*3 + 0], translations[(i)*3 + 1], translations[(i)*3 + 2];
        feature_track_.Abs_Transls.push_back(t);
    }

    //> Construct pairwise epipolar constraints
    int M_cameras = feature_track_.Length;
    Eigen::MatrixXd idx_matrix;
    int n_comb = generateM2Comb(M_cameras, idx_matrix);

    //> Construct pairwise relative poses and feature correspondences
    std::vector< PairObj > Feature_Track_Connections;
    for (int jj = 0; jj < n_comb; jj++) {
        PairObj feature_corr_i; 
        int id1 = idx_matrix(0, jj); 
        int id2 = idx_matrix(1, jj); 

        //> Compute relative poses
        Matrix3 R1 = feature_track_.Abs_Rots[id1];
        Matrix3 R2 = feature_track_.Abs_Rots[id2]; 
        Vector3 t1 = feature_track_.Abs_Transls[id1]; 
        Vector3 t2 = feature_track_.Abs_Transls[id2]; 
        Matrix4 P1 = Matrix4::Identity(); 
        Matrix4 P2 = Matrix4::Identity(); 
        P1.block<3,3>(0,0) = R1; 
        P1.block<3,1>(0,3) = t1;
        P2.block<3,3>(0,0) = R2; 
        P2.block<3,1>(0,3) = t2;
                
        Matrix4 Prel = P2 * P1.inverse(); 
        Matrix3 Rrel = Prel.block<3,3>(0,0); 
        Vector3 trel = Prel.block<3,1>(0,3);                 
        trel.normalize();
                
        //> Compute essential matrices
        Matrix3 Ess = Matrix3::Identity(); 
        Matrix3 Tx = Matrix3::Zero(); 
        Tx << 0, -trel(2), trel(1), trel(2), 0, -trel(0), -trel(1), trel(0), 0; 
        Ess = Tx * Rrel;

        //> Organize feature connection information
        feature_corr_i.id1 = id1; 
        feature_corr_i.id2 = id2; 
        feature_corr_i.F = Ess; 
        feature_corr_i.p1 = inverse_K * feature_track_.Locations[id1]; //str_out.obs[0].col(id1); 
        feature_corr_i.p2 = inverse_K * feature_track_.Locations[id2];
        Feature_Track_Connections.push_back(feature_corr_i);

        if (mexFunction_debug) {
            std::cout << "Indices:" << std::endl << feature_corr_i.id1 << ", " << feature_corr_i.id2 << std::endl;
            std::cout << "Essential matrix:" << std::endl << feature_corr_i.F << std::endl;
            std::cout << "Feature point pair:" << std::endl << "[" << feature_corr_i.p1(0) << ", " << feature_corr_i.p1(1) << ", " << feature_corr_i.p1(2) << "]" << std::endl;
            std::cout << "[" << feature_corr_i.p2(0) << ", " << feature_corr_i.p2(1) << ", " << feature_corr_i.p2(2) << "]" << std::endl;
        }
    }

    //> Run N-view triangulation!
    NViewsClass corr_N_view; 
    //> (1) Create constraint matrices
    corr_N_view.createProblemMatrices(Feature_Track_Connections, M_cameras); 
       
    //> (2) Run correction and check whether the corrected features are globally optimal solutions from N-view triangulation
    NViewsOptions options_corr; 
    options_corr.save_val_constr = false;
    options_corr.debug           = false; 
    NViewsResult Feature_Corrections = corr_N_view.correctObservations(options_corr);
    bool certified_global_optimum = corr_N_view.certified_global_optimum;
    std::cout << "Is solution a global optimum? " << (certified_global_optimum ? std::string("Yes") : std::string("No")) << std::endl;
 
    //> (3) Make corrections to the observed feature track point locations
    std::vector<Matrix4> proj_s; 
    std::vector<Vector3> Corrected_Features;  
        
    for (int jc=0; jc < M_cameras;jc++) {
        //> Projection matrix
        Matrix3 R = feature_track_.Abs_Rots[jc];  
        Vector3 t = feature_track_.Abs_Transls[jc]; 
        Matrix4 P1 = Matrix4::Identity(); 
        P1.block<3,3>(0,0) = R; 
        P1.block<3,1>(0,3) = t;
        proj_s.push_back(P1); 
                
        //> update observation by the correction from N-view triangulation 
        Vector3 pt = inverse_K * feature_track_.Locations[jc];  
        Vector3 delta_ref; 
        delta_ref << Feature_Corrections.sol_final( jc*2), Feature_Corrections.sol_final(jc*2 + 1), 0; 
                
        Corrected_Features.push_back(pt + delta_ref); 
    }

    //> (4) Triangulate corrected points to a common 3D point via linear triangulation and reproject to the images
    Vector3 Gamma_Corrected;    //> 3D point under the world coordinate
    Eigen::VectorXd Depths_Corrected; 
    double Linearity_Err = triangulateNPoint(proj_s, Corrected_Features, Gamma_Corrected, Depths_Corrected);
    for (int i = 0; i < feature_track_.Length; i++) {
        feature_track_.Optimal_Locations.push_back(K * Corrected_Features[i]);
    }
    std::vector<double> Reprojection_Errors = reproject_to_images( proj_s, feature_track_.Locations, K, Gamma_Corrected );

    //> Output reprojection errors to MATLAB
    if ( nl > 0 ) {
        //> 1) Corrected feature track locations
        const mwSize feature_location_dim[2] = { mwSize(feature_track_.Length), mwSize(2) };
        pl[0] = mxCreateNumericArray(2, feature_location_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_corrected_features = (double*) mxGetData(pl[0]);
        for (int i = 0; i < Corrected_Features.size(); i++) {
            (output_corrected_features)[(i)*2 + 0] = Corrected_Features[i](0);
            (output_corrected_features)[(i)*2 + 1] = Corrected_Features[i](1);
        }
        //> 2) reprojection errors
        const mwSize reproject_err_dim[1] = { mwSize(feature_track_.Length) };
        pl[1] = mxCreateNumericArray(1, reproject_err_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_reproj_errors = (double*) mxGetData(pl[1]);
        for (int i = 0; i < Corrected_Features.size(); i++) {
            (output_reproj_errors)[i] = Reprojection_Errors[i];
        }
        mexPrintf("[CPP] Complete N-view Triangulation.\n");
    }
    else {
        mexErrMsgTxt("No output on the MATLAB side is given!");
    }
}

