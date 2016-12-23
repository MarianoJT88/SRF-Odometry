//====================================================
//  Project: Laser odometry
//  Authors: Mariano Jaimez Tarifa, Javier G. Monroy
//           MAPIR group, University of Malaga, Spain
//  Date: January 2016
//====================================================


#include <mrpt/poses/CPose2D.h>
#include <mrpt/utils/CTicTac.h>
#include <Eigen/Dense>
#include <iostream>


//#define M_LOG2E 1.44269504088896340736 //log2(e)
//inline float log2(const float x){
//    return  log(x) * M_LOG2E;
//}


class RF2O_RefS {
public:

    //Scans and cartesian coordinates: 1 - New, 2 - Old, 3 - Ref
    Eigen::ArrayXf range_wf;
    std::vector<Eigen::ArrayXf> range_1, range_2, range_3;
    std::vector<Eigen::ArrayXf> range_12, range_13, range_warped;
    std::vector<Eigen::ArrayXf> xx_1, xx_2, xx_3, xx_12, xx_13, xx_warped;
    std::vector<Eigen::ArrayXf> yy_1, yy_2, yy_3, yy_12, yy_13, yy_warped;
    std::vector<Eigen::ArrayXf> range_3_warpedTo2, xx_3_warpedTo2, yy_3_warpedTo2;

    //Rigid transformations and velocities (twists: vx, vy, w)
    std::vector<Eigen::MatrixXf> transformations; //T13
    Eigen::Matrix3f overall_trans_prev; // T23
    Eigen::Vector3f kai_abs, kai_loc;
    Eigen::Vector3f kai_loc_old, kai_loc_level;

    //Solver
    Eigen::MatrixXf A,Aw;
    Eigen::VectorXf B,Bw;
    Eigen::Matrix3f cov_odo;
	
    //Aux variables
    Eigen::ArrayXf dtita_12, dtita_13;
    Eigen::ArrayXf dt_12, dt_13;
    Eigen::ArrayXf weights_12, weights_13;
    Eigen::Array<bool, Eigen::Dynamic, 1> null_12, null_13;
    Eigen::Array<bool, Eigen::Dynamic, 1> outliers;

	float fovh;
    unsigned int cols, cols_i;
	unsigned int width;
	unsigned int ctf_levels;
	unsigned int image_level, level;
	unsigned int num_valid_range;
	unsigned int iter_irls;
	float g_mask[5];
    bool no_ref_scan;
    bool new_ref_scan;
    unsigned int method_ref_scan; //0 - ours, 1 - trans and rot thres, 2 - MAD(res)


    //Laser poses (most recent and previous)
    mrpt::poses::CPose2D laser_pose;
    mrpt::poses::CPose2D laser_oldpose;
	bool test;
    unsigned int method; //0 - consecutive scan alignment, 1 - keyscan alignment, 2 - multi-scan (hybrid) alignment

    //To measure runtimes
    mrpt::utils::CTicTac	clock;
    float                   runtime;


    //Methods
    void initialize(unsigned int size, float FOV_rad, unsigned int odo_method);
    void createScanPyramid();
	void calculateCoord();
	void performWarping();
    void performBestWarping();
    void warpScan3To2();
    void calculateRangeDerivatives();
	void computeWeights();
    void solveSystemQuadResiduals3Scans();
    void solveSystemSmoothTruncQuad3Scans();
    void solveSystemSmoothTruncQuadOnly13();
    void solveSystemSmoothTruncQuadOnly12();
    void solveSystemMCauchy();
    void solveSystemMTukey();
    void solveSystemTruncatedQuad();

	void filterLevelSolution();
	void PoseUpdate();
    void updateReferenceScan();
	void odometryCalculation();
};


// Extra functionality - Put variables and functions as class member
//==================================================================

//ofstream	m_fres;
//void OpenResFile();
//void WriteTrajFile();

//void CLaserOdo::OpenResFile()
//{
//	m_fres.open(".../resdata.txt");
//}

//void CLaserOdo::WriteTrajFile()
//{
//  //Write here what you want to save/analyze
//
//	//auxpose = m_pose_f - transf;
//	//auxpose.getAsQuaternion(quat);
//
//	////printf("\n Quaternions (est): %f, %f, %f, %f, %f, %f, %f", m_pose_f[0], m_pose_f[1], m_pose_f[2], quat(2), quat(3), -quat(1), -quat(0));
//	//
//	//char aux[24];
//	//sprintf(aux,"%.04f", timestamp_obs);
//	//m_fres << aux << " ";
//	//m_fres << m_pose_f[0] << " ";
//	//m_fres << m_pose_f[1] << " ";
//	//m_fres << m_pose_f[2] << " ";
//	//m_fres << quat(2) << " ";
//	//m_fres << quat(3) << " ";
//	//m_fres << -quat(1) << " ";
//	//m_fres << -quat(0) << endl;
//}


//void CLaserOdo::computeNormals()
//{
//    ArrayXf normx(cols), normy(cols), norm_ang(cols);
//    normx.fill(0.f); normy.fill(0.f); norm_ang.fill(0.f);

//	const float incr_tita = fovh/float(cols_i-1);
//	for (unsigned int u=0; u<cols_i; u++)
//	{
//        if (null(u) == false)
//		{
//			const float tita = -0.5f*fovh + float(u)*incr_tita;
//			const float alfa = -atan2(2.f*dtita(u), 2.f*range[image_level](u)*incr_tita);
//			norm_ang(u) = tita + alfa;
//			if (norm_ang(u) < -M_PI)
//				norm_ang(u) += 2.f*M_PI;
//			else if (norm_ang(u) < 0.f)
//				norm_ang(u) += M_PI;
//			else if (norm_ang(u) > M_PI)
//				norm_ang(u) -= M_PI;

//			normx(u) = cos(tita + alfa);
//			normy(u) = sin(tita + alfa);
//		}
//	}
//}

