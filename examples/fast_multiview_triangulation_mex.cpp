// =============================================================================================================================
//> This N-view triangulation code is built on top of the publicly released certifiable solver for multiple views triangulation
//  https://github.com/C-H-Chien/FastNViewTriangulation arises from the paper:
//  Garcia-Salguero, Mercedes, and Javier Gonzalez-Jimenez. "Certifiable solver for real-time N-view triangulation." 
//  IEEE Robotics and Automation Letters 8, no. 4 (2023): 1999-2005.
//
//> Build this code in the matlab command window:
//  $ mex examples/fast_multiview_triangulation_mex.cpp src/NViewsCertifier.cpp src/NViewsClass.cpp src/NViewsUtils.cpp ...
//    utils/generatePointCloud.cpp -I/path/to/eigen3
//  [Note] the path to the eigen3 library depends on the directory where you install eigen3, e.g., -I/usr/include/eigen3/
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

//> Macros
#define VERBOSE (false)

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
        if (mexFunction_debug) 
            std::cout << "Obs" << i << " = (" << Observation_Locations[(i)*2 + 0] << ", " << Observation_Locations[(i)*2 + 1] << ")" << std::endl;
        obs_i.col(i) = point;
        feature_track_.Locations.push_back(point);

        //> absolute rotation matrices
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        R << rotations[(i)*9 + 0], rotations[(i)*9 + 1], rotations[(i)*9 + 2], \
             rotations[(i)*9 + 3], rotations[(i)*9 + 4], rotations[(i)*9 + 5], \
             rotations[(i)*9 + 6], rotations[(i)*9 + 7], rotations[(i)*9 + 8];
        feature_track_.Abs_Rots.push_back(R);
        if (mexFunction_debug) {
            std::cout << "R" << i << " =" << std::endl;
            std::cout << R << std::endl;
        }

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
                
        Matrix4 Prel = P2.inverse() * P1; 
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
    options_corr.max_iters            = 50;
    NViewsResult Feature_Corrections = corr_N_view.correctObservations(options_corr);
    bool certified_global_optimum = corr_N_view.certified_global_optimum;
#if VERBOSE
    std::cout << "Is solution a global optimum? " << (certified_global_optimum ? std::string("Yes") : std::string("No")) << std::endl;
#endif
 
    //> (3) Make corrections to the observed feature track point locations
    std::vector<Matrix4> proj_s; 
    std::vector<Vector3> Corrected_Features_in_Metric;
    std::vector<Vector3> Observed_Features_in_Metric;
        
    for (int jc=0; jc < M_cameras;jc++) {
        //> Projection matrix
        Matrix3 R = feature_track_.Abs_Rots[jc];  
        Vector3 t = feature_track_.Abs_Transls[jc]; 
        Matrix4 P1 = Matrix4::Identity(); 
        //> Projection matrix: projecting 3D points under a common world coordinate to the camera coordinate
        P1.block<3,3>(0,0) = R.inverse(); 
        P1.block<3,1>(0,3) = -R.inverse() * t;
        proj_s.push_back(P1); 
                
        //> update observation by the correction from N-view triangulation 
        Vector3 pt = inverse_K * feature_track_.Locations[jc];  
        Vector3 delta_ref; 
        delta_ref << Feature_Corrections.sol_final( jc*2), Feature_Corrections.sol_final(jc*2 + 1), 0;

        Corrected_Features_in_Metric.push_back( pt + delta_ref );
        Observed_Features_in_Metric.push_back( pt );
        feature_track_.Optimal_Locations.push_back( K * Corrected_Features_in_Metric[jc] ); 
    }

    //> (4) Triangulate corrected points to a common 3D point via linear triangulation and reproject to the images
    Vector3 Gamma_Corrected;            //> 3D corrected point under the world coordinate
    Eigen::VectorXd Depths_Corrected;   //> (?)
    double Linearity_Err = triangulateNPoint(proj_s, Corrected_Features_in_Metric, Gamma_Corrected, Depths_Corrected);
    std::vector<double> Reprojection_Errors = reproject_to_images( proj_s, feature_track_.Locations, K, Gamma_Corrected, mexFunction_debug );

    //> Used for uncertainty measurements
    Vector3 Gamma_Observed;             //> 3D observed point under the world coordinate
    Eigen::VectorXd Depths_Observed;    //> (?)
    triangulateNPoint(proj_s, Observed_Features_in_Metric, Gamma_Observed, Depths_Observed);

    //> Output reprojection errors to MATLAB
    if ( nl > 0 ) {
        //> 1) Corrected feature track locations
        const mwSize feature_location_dim[2] = { mwSize(2), mwSize(feature_track_.Length) };
        pl[0] = mxCreateNumericArray(2, feature_location_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_corrected_features = (double*) mxGetData(pl[0]);
        for (int i = 0; i < Corrected_Features_in_Metric.size(); i++) {
            (output_corrected_features)[(i)*2 + 0] = feature_track_.Optimal_Locations[i](0);
            (output_corrected_features)[(i)*2 + 1] = feature_track_.Optimal_Locations[i](1);
            if (mexFunction_debug) 
                std::cout << "Corrected Features p" << i << " = (" << feature_track_.Optimal_Locations[i](0) << ", " << feature_track_.Optimal_Locations[i](1) << ")" << std::endl;
        }
        //> 2) reprojection errors
        const mwSize reproject_err_dim[1] = { mwSize(feature_track_.Length) };
        pl[1] = mxCreateNumericArray(1, reproject_err_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_reproj_errors = (double*) mxGetData(pl[1]);
        for (int i = 0; i < Corrected_Features_in_Metric.size(); i++) {
            (output_reproj_errors)[i] = Reprojection_Errors[i];
        }
        //> 3) whether the solution is a global optimum
        const mwSize is_sol_global_optimal_dim[1] = { mwSize(1) };
        pl[2] = mxCreateNumericArray(1, is_sol_global_optimal_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_is_sol_global_optimal = (double*) mxGetData(pl[2]);
        output_is_sol_global_optimal[0] = (double)(certified_global_optimum);

        // pl[2] = mxCreateDoubleScalar( (double)certified_global_optimum );
        //> 4) triangulated 3D point, both corrected and observed
        const mwSize Gamma_dim[2] = { mwSize(3), mwSize(2) };
        pl[3] = mxCreateNumericArray(2, Gamma_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_Gamma = (double*) mxGetData(pl[3]);
        for (int i = 0; i < 3; i++) {
            output_Gamma[i]     = Gamma_Corrected[i];
            output_Gamma[i + 3] = Gamma_Observed[i];
        }
#if VERBOSE
        mexPrintf("[CPP] Complete N-view Triangulation.\n");
#endif
    }
    else {
        mexErrMsgTxt("No output on the MATLAB side is given!");
    }
}

