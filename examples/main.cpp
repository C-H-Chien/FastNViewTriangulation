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

#include "../include/NViewsTypes.h"
#include "../include/NViewsUtils.h"
#include "../include/NViewsClass.h"
#include "../include/definitions.h"
#include "../utils/generatePointCloud.h"

#include <Eigen/Dense>
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
    bool is_sol_global_optimal;
    std::vector<double> Reprojection_Errors;
    Eigen::Vector3d Gamma;

}; //> End of struct Feature_Track

void Multiview_Triangulation( Feature_Track& feature_track_, Eigen::Matrix3d K ) {

    Eigen::Matrix3d inverse_K = K.inverse();
    
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
        Matrix3 Rrel = R2 * R1.transpose();
        Vector3 trel = t2 - R2*(R1.transpose())*t1;
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
    }

    //> Run N-view triangulation!
    NViewsClass corr_N_view; 
    //> (1) Create constraint matrices
    corr_N_view.createProblemMatrices(Feature_Track_Connections, M_cameras); 
       
    //> (2) Run correction and check whether the corrected features are globally optimal solutions from N-view triangulation
    NViewsOptions options_corr; 
    options_corr.save_val_constr = false;
    options_corr.debug           = false; 
    options_corr.max_iters       = 50;
    NViewsResult Feature_Corrections = corr_N_view.correctObservations(options_corr);
    bool certified_global_optimum = corr_N_view.certified_global_optimum;
    feature_track_.is_sol_global_optimal = certified_global_optimum;
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
        // P1.block<3,3>(0,0) = R.inverse(); 
        // P1.block<3,1>(0,3) = -R.inverse() * t;
        P1.block<3,3>(0,0) = R; 
        P1.block<3,1>(0,3) = t;
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
    std::vector<double> Reprojection_Errors = reproject_to_images( proj_s, feature_track_.Locations, K, Gamma_Corrected, false );
    feature_track_.Reprojection_Errors = Reprojection_Errors;
    feature_track_.Gamma = Gamma_Corrected;
}

int main(int argc, char** argv) {

    //> define input data to multiview triangulation
    Feature_Track feature_track_;

    //> Number of covisible views
    feature_track_.Length = 4;

    //> Edge points across these covisible views
    Eigen::Vector3d edge_point;
    edge_point << 565.125, 475.835, 1.0;    //> image 6
    feature_track_.Locations.push_back(edge_point);
    edge_point << 540.16, 514.655, 1.0;     //> image 8
    feature_track_.Locations.push_back(edge_point);
    edge_point << 555.296, 465.492, 1.0;    //> image 0
    feature_track_.Locations.push_back(edge_point);
    edge_point << 293.883, 381.497, 1.0;    //> image 2
    feature_track_.Locations.push_back(edge_point);
    
    //> Intrinsic matrix
    double fx = 1111.11136542426;
    double fy = 1111.11136542426;
    double cx = 399.500000000000;
    double cy = 399.500000000000;
    Eigen::Matrix3d K;
    K << fx, 0, cx, 0, fy, cy, 0, 0, 1;

    //> Camera rotations and translations
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t = Eigen::Vector3d::Zero();
    R << 0.429735013300083,  -0.902955132327889,   0.000000066933705, \
        -0.548736945415301,  -0.261155204677633,   0.794157240654608, \
         0.717088259484231,   0.341277161677635,   0.607712377310032;
    t << 0.236609517076345, 0.007867526316630, -3.992988151863210;
    feature_track_.Abs_Rots.push_back(R);
    feature_track_.Abs_Transls.push_back(t);
    R << 0.640732122716973,  -0.767764577353699,  -0.000000195076584, \
        -0.616277981699726,  -0.514310272903074,   0.596394527025542, \
         0.457890695902617,   0.382128999349512,   0.802691436807719;
    t << 0.063516276608768, 0.267096387307519, -3.990567076807020;
    feature_track_.Abs_Rots.push_back(R);
    feature_track_.Abs_Transls.push_back(t);
    R << 0.577608355645856,  -0.816314158594680,   0.000000004914734, \
        -0.214887229677298,  -0.152050105145358,   0.964730136681733, \
         0.787522780472976,   0.557236123528490,   0.263240870427739;
    t << 0.119353058702026, -0.298896407235438, -3.987030846688150;
    feature_track_.Abs_Rots.push_back(R);
    feature_track_.Abs_Transls.push_back(t);
    R << 0.683624427989991,   0.729834098129015,   0.000000062252118, \
        -0.318453509830794,   0.298290465243610,   0.899783417421079, \
        -0.656692540009721,   0.615113911740625,  -0.436336873514503;
    t << -0.706728946734671, -0.439810638182904, -3.912429209864560;
    feature_track_.Abs_Rots.push_back(R);
    feature_track_.Abs_Transls.push_back(t);

    //> Run multiview triangulation
    Multiview_Triangulation( feature_track_, K );

    //> print out multiview triangulation results
    //> (1) corrected point locations
    std::cout << "Corrected point locations:" << std::endl;
    for (int  i = 0; i < feature_track_.Length; i++) {
        Eigen::Vector3d corrected_location = feature_track_.Optimal_Locations[i];
        std::cout << "Point " << i << " = (" << corrected_location(0) << ", " << corrected_location(1) << ", " << corrected_location(2) << ")" << std::endl;
    }
    //> (2) corrected 3D point
    std::cout << "3D point = (" << feature_track_.Gamma(0) << ", " << feature_track_.Gamma(1) << "," << feature_track_.Gamma(2) << ")" << std::endl;
}

