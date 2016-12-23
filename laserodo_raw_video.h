/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */


#include <mrpt/obs/CRawlog.h>
#include <mrpt/obs/CObservationOdometry.h>
#include <mrpt/opengl.h>
#include <mrpt/opengl/CPlanarLaserScan.h>
#include <mrpt/maps/COccupancyGridMap2D.h>
#include <mrpt/gui.h>
#include <mrpt/system/filesystem.h>
#include <mrpt/utils/round.h>

#include "laser_odometry_v1.h"
#include "laser_odometry_standard.h"
#include "laser_odometry_3scans.h"
#include "laser_odometry_refscans.h"
#include "polar_match.h"
#include "csm/csm_all.h"
//#include "csm/sm/csm/csm_all.h"


using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::opengl;
using namespace mrpt::maps;
using namespace mrpt::obs;
using namespace mrpt::gui;
using namespace mrpt::poses;
using namespace CSM;
using namespace std;



struct TRobotLaser {
	CObservation2DRangeScan m_scan;
    CObservation2DRangeScan m_scan_old;
	int						m_level;
	int						m_segments;
};


class CLaserodoInterface
{
public:

	CPose2D					new_pose;
	CPose2D					last_pose;
	COccupancyGridMap2D		map;
	TRobotLaser				laser;
	CDisplayWindow3D		window;
	COpenGLScenePtr			scene;
    CRawlog                 dataset;
    unsigned int dataset_id;
    unsigned int rawlog_count;
    bool dataset_finished;
    bool localized;
    float laser_min_range;
    float old_camera_angle;

    //RF2O
    RF2O   odo;
    RF2O_RefS odo_test;

    //Polar scan matcher
    CPose2D     new_psm_pose, old_psm_pose;
    PMScan      ls, ls_ref;
    bool        matchFailed;

    //Canonical scan matcher
    sm_params params;               //Structure containing the input parameters
    sm_result result;               //Output of the scan matching
    LDP laser_ref, laser_sens;      //the laser scans to be matched (in the PL-ICP data structure)
    int recover_from_error;
    CPose2D new_csm_pose, old_csm_pose;

    //Wheel odometry
    CPose2D new_gt_pose, old_gt_pose, pose_ini;

    //Results
    vector<CPose3D>	real_poses, est_poses, test_poses, psm_poses, csm_poses;
    float est_time, test_time, psm_time, csm_time;


    void initializeEverything()
    {
        //Read Rawlog
        string folder = "/usr/wiss/jaimez/Dropbox/LaserOdo - shared/Datasets/";
        //string folder = "/home/mariano/Dropbox/LaserOdo - shared/Datasets/";

        dataset_id = 4;
        string dset;

        switch (dataset_id) {
        case 1:
            dset = "noMovement_hokuyoSmall_2.rawlog";
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180); laser.m_segments = 181; laser.m_scan.maxRange = 5.6f;
            break;
        case 2:
            dset = "RobotStopped_MovingObjects_4.rawlog";
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180); laser.m_segments = 181; laser.m_scan.maxRange = 56.f;
            break;
        case 3:
            //dset = "RealExperiment_MapBuilding.rawlog";
            //dset = "exp_Giraff_ORIGINAL.rawlog";
            dset = "exp_Giraff_LASER1.rawlog";
            //laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180); laser.m_segments = 181; laser.m_scan.maxRange = 5.6f;
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(240); laser.m_segments = 682; laser.m_scan.maxRange = 5.6f;
            break;
        case 4:
            //dset = "Freiburg_Indoor/Freiburg_Indoor_Building_079.rawlog";
            dset = "Freiburg_Indoor/freiburg_gt.rawlog"; //It works!!!
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180); laser.m_segments = 360; laser.m_scan.maxRange = 80.f;
            break;
        case 5:
            dset = "MIT_CSAIL/mit_csail.rawlog";
            //dset = "MIT_CSAIL/mit_csail_gt.rawlog";
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180.5); laser.m_segments = 361; laser.m_scan.maxRange = 80.f;
            break;
        case 6:
            dset = "NDT_orebro_withGT/localization.rawlog";
            laser_min_range = 0.7f; laser.m_scan.aperture = DEG2RAD(270); laser.m_segments = 541; laser.m_scan.maxRange = 29.f;
            break;
        }

        string filename = folder + dset;

        rawlog_count = 0;
        dataset_finished = false;
        localized = false;
        if (!dataset.loadFromRawLogFile(filename))
            throw std::runtime_error("\nCouldn't open rawlog dataset file for input...");

        new_pose = CPose2D(0.f, 0.f, 0.05f);
        last_pose = new_pose;

        //Read one scan
        readScanRawlog();
        laser.m_scan_old = laser.m_scan;

        //Psm
        matchFailed = false;
        pm_readScan(laser.m_scan, &ls);
        pm_preprocessScan(&ls);


        //Csm
        laser_ref = cast_CObservation2DRangeScan_to_LDP(laser.m_scan);
        // For the first scan, set estimate = odometry
        copy_d(laser_ref->odometry, 3, laser_ref->estimate);
        pm_init();                       //Contains precomputed range bearings and their sines/cosines
        Init_Params_csm();

        //RF2O
        odo.initialize(laser.m_segments, laser.m_scan.aperture, 3);
        odo_test.initialize(laser.m_segments, laser.m_scan.aperture, 2);

//        odo.initialize(laser.m_segments, laser.m_scan.aperture, false);
//        odo_test.initialize(laser.m_segments, laser.m_scan.aperture, true);
        loadFirstScanRF2O();

        //Scene and poses  
        setRF2OPose(new_pose);

        new_psm_pose = new_pose; old_psm_pose = new_pose;
        new_csm_pose = new_pose; old_csm_pose = new_pose;
        new_gt_pose = new_pose; old_gt_pose = new_pose;
        est_time = 0.f; test_time = 0.f; psm_time = 0.f; csm_time = 0.f;

        initializeScene();
    }

	

    void readScanRawlog()
	{
        CObservationPtr alfa = dataset.getAsObservation(rawlog_count);
        old_gt_pose = new_gt_pose;

        while (!IS_CLASS(alfa, CObservation2DRangeScan))
        {
            rawlog_count++;
            if (dataset.size() <= rawlog_count)
            //if (4280 <= rawlog_count)
            {
                dataset_finished = true;
                return;
            }

            if (IS_CLASS(alfa, CObservationOdometry))
            {
                CObservationOdometryPtr obs_odo = CObservationOdometryPtr(alfa);

                if (!localized)
                {
                    pose_ini = obs_odo->odometry;
                    localized = true;
                }
                else
                    new_gt_pose = obs_odo->odometry - pose_ini;

            }

            alfa = dataset.getAsObservation(rawlog_count);
        }

        CObservation2DRangeScanPtr obs2D = CObservation2DRangeScanPtr(alfa);

        laser.m_scan_old = laser.m_scan;
        laser.m_scan = *obs2D;
        for (unsigned int i=0; i<laser.m_segments; i++)
            if ((laser.m_scan.scan[i] > 0.995f*laser.m_scan.maxRange)||(laser.m_scan.scan[i] < laser_min_range)) //(i < 70) || (i>471))
            {
                laser.m_scan.scan[i] = 0.f;
                laser.m_scan.validRange[i] = false;
            }
            else
            {
                laser.m_scan.validRange[i] = true;
            }

        rawlog_count++;

        if (dataset.size() <= rawlog_count)
            dataset_finished = true;

        printf("\n Scan %d read", rawlog_count);
//        printf("\n Scan aperture = %f, segments = %d", RAD2DEG(laser.m_scan.aperture), laser.m_segments);
//        for (unsigned int i=0; i<laser.m_segments; i++)
//            printf("\n range = %f", laser.m_scan.scan[i]);
	}

    void initializeScene()
    {
        CPose3D robotpose3d;
        robotpose3d.x(0);
        robotpose3d.y(0);
        robotpose3d.setYawPitchRoll(0,0,0);

        //The display window is created
        mrpt::global_settings::OCTREE_RENDER_MAX_POINTS_PER_NODE = 1000000;
        //window.setWindowTitle("Reactive Navigation. Robot motion simulation");
        window.getDefaultViewport()->setCustomBackgroundColor(TColorf(1.f,1.f,1.f));
        window.resize(720,540); //indow.resize(1600,1200);
        window.setPos(1400,0);
        if (dataset_id == 4)        window.setCameraZoom(20);
        else if (dataset_id == 5)   window.setCameraZoom(30);
        //window.setCameraZoom(30); //75 for MIT, 45 for Freiburg (for the images of the sub-maps)
        window.setCameraAzimuthDeg(-90);
        old_camera_angle = -90;
        window.setCameraElevationDeg(90);
        scene = window.get3DSceneAndLock();


        //Robots
//        CSetOfObjectsPtr robot_real = opengl::stock_objects::CornerXYZ();
//        robot_real->setScale(0.5f);
//        //robot_real->setColor(0.15,0.35,1);
//        robot_real->setName("laser");
////        CPolyhedronPtr neck = robot_real->getByClass<CPolyhedron>(2);
////        neck->setColor(0.5,0.5,0.5);
//        scene->insert( robot_real );

        //Position of the laser
        CSpherePtr laser = opengl::CSphere::Create(0.03, 20, 20);
        laser->setColor(1.f, 0.f, 0.f);
        scene->insert( laser );


        //Trajectory
        CSetOfLinesPtr traj_lines = opengl::CSetOfLines::Create();
        traj_lines->setColor(1,0,0);
        traj_lines->setLineWidth(4);
        scene->insert( traj_lines );


        window.unlockAccess3DScene();
        window.repaint();
    }

    void updateScene(unsigned int iter)
    {
        const unsigned int max_number_scans = 300; //Freiburg - 200

        //Update scene
        scene = window.get3DSceneAndLock();
        //CPose2D pose = odo_test.laser_pose, pose_old = odo_test.laser_oldpose;
        //CPose2D pose = new_psm_pose, pose_old = old_psm_pose;
        //CPose2D pose = new_csm_pose, pose_old = old_csm_pose;
        CPose2D pose = new_gt_pose, pose_old = old_gt_pose;


        //Camera follows the robot (and uses the gt orientation with smoothing)
        window.setCameraPointingToPoint(pose[0], pose[1], 0.f);
        const float w = 0.05f;
        float ang_dif = RAD2DEG(pose[2] - new_gt_pose[2]);
        if (ang_dif > 90)       ang_dif -= 360;
        else if (ang_dif < -90) ang_dif += 360;
        const float new_cam_angle = ang_dif-90;

        if (new_cam_angle - old_camera_angle > 90)         old_camera_angle += 360;
        else if (new_cam_angle - old_camera_angle < -90)   old_camera_angle -= 360;
        float filtered_cam_angle = w*new_cam_angle + (1.f-w)*old_camera_angle;

        printf("\n New angle = %f, old angle = %f, filtered = %f", new_cam_angle, old_camera_angle, filtered_cam_angle);

        old_camera_angle = filtered_cam_angle;



        window.setCameraAzimuthDeg(filtered_cam_angle);

        //Robot
//        CRenderizablePtr robot_real = scene->getByName("laser");
//        robot_real->setPose(pose);

        //Insert new laser
        const unsigned int repr_level = round(log2(round(float(odo.width)/float(odo.cols))));
        CPointCloudPtr gl_laser = opengl::CPointCloud::Create();
        gl_laser->setColor(0.f, 0.f, 1.f);
        gl_laser->setPointSize(3.f);
        gl_laser->setPose(pose);
        scene->insert( gl_laser );

        //Get color
//            Eigen::Vector3f ini_c; ini_c << 1.f,0.f,0.f;
//            Eigen::Vector3f end_c; end_c << 0.f,0.f,1.f;

//            const float alpha = float(dataset.size() - rawlog_count)/float(dataset.size());
//            Eigen::Vector3f c_now = alpha*end_c + (1.f - alpha)*ini_c;
//            float r,g,b;
//            utils::jet2rgb(alpha, r, g, b);


        //if ((iter+4) % 5 == 0)
        {
            for (unsigned int i=0; i<odo.cols; i++)
                if ((odo.range[repr_level](i) < 3.8f)&&(odo.range[repr_level](i) > 0.05f))
                {
                    //const float x_old = odo.xx[repr_level](i);
                    //const float y_old = odo.yy[repr_level](i);
                    //const float x = x_old*cos(pose[2]) - y_old*sin(pose[2]) + pose[0];
                    //const float y = x_old*sin(pose[2]) + y_old*cos(pose[2]) + pose[1];
                    //gl_laser->push_back(x, y, 0.f, 0.f, 0.f, 1.f);
                    gl_laser->insertPoint(odo.xx[repr_level](i), odo.yy[repr_level](i), 0.f);
                }
        }

        //Position of the laser
        CSpherePtr laser = scene->getByClass<CSphere>(0);
        laser->setPose(pose);


        //Trajectory
        CSetOfLinesPtr traj_lines = scene->getByClass<CSetOfLines>(0);
        traj_lines->appendLine(pose_old[0], pose_old[1], 0.01, pose[0], pose[1], 0.01);

        //Remove old observations and trajectory
        if (traj_lines->size() > max_number_scans)
        {
            //traj_lines->removeFirstLine();
            CPointCloudPtr old_laser = scene->getByClass<CPointCloud>(0);
            CRenderizablePtr obj = old_laser;
            scene->removeObject(obj);
        }


        window.unlockAccess3DScene();
        window.repaint();
    }

	void resetScene()
	{
		scene = window.get3DSceneAndLock();
        CPose3D robotpose3d = new_pose;
		CRenderizablePtr obj;

		odo.laser_pose = CPose2D(robotpose3d);
		odo.laser_oldpose = CPose2D(robotpose3d);
		odo.kai_abs.assign(0.f);
		odo.kai_loc.assign(0.f);
		odo.kai_loc_old.assign(0.f);

		odo_test.laser_pose = CPose2D(robotpose3d);
		odo_test.laser_oldpose = CPose2D(robotpose3d);
		odo_test.kai_abs.assign(0.f);
		odo_test.kai_loc.assign(0.f);
		odo_test.kai_loc_old.assign(0.f);

        new_psm_pose = CPose2D(robotpose3d);
        old_psm_pose = CPose2D(robotpose3d);

        new_csm_pose = CPose2D(robotpose3d);
        old_csm_pose = CPose2D(robotpose3d);

		//Robots
        obj = scene->getByName("laser");
        obj->setPose(robotpose3d);


		//Laser
		const unsigned int repr_level = round(log2(round(float(odo.width)/float(odo.cols))));
		CPointCloudColouredPtr gl_laser;
		gl_laser = scene->getByClass<CPointCloudColoured> (0);
		gl_laser->clear();
		for (unsigned int i=0; i<odo.cols; i++)
			gl_laser->push_back(odo.xx[repr_level](i), odo.yy[repr_level](i), 0.1, 1-sqrt(odo.weights(i)), sqrt(odo.weights(i)), 0);

        gl_laser->setPose(robotpose3d);

		//Trajectories
		opengl::CSetOfLinesPtr traj_lines_real;
		traj_lines_real = scene->getByClass<CSetOfLines> (0);
        traj_lines_real->clear();


		window.unlockAccess3DScene();
		std::string legend;
		legend.append("--------------------------------------------\n");
		legend.append("| n - Simulate one step \n");
		legend.append("| s - Start/stop continuous sim. \n");
		legend.append("| m - Move the target \n");
		legend.append("| c - Reset estimated pose \n");
		legend.append("| e - Exit \n");
		legend.append("--------------------------------------------\n");
		legend.append(format("\n        %.02fFPS", window.getRenderingFPS()));
		
		window.addTextMessage(5,180, legend,utils::TColorf(1,1,1),"Arial",13);
		window.repaint();
	}

    void runRF2O()
    {
        //Run the robust nonlinear version
        for (unsigned int i=0; i<odo.width; i++)
            odo.range_wf(i) = laser.m_scan.scan[i];

        odo.odometryCalculation();
        est_time += odo.runtime;

        //Run the robust nonlinear version
        for (unsigned int i=0; i<odo_test.width; i++)
            odo_test.range_wf(i) = laser.m_scan.scan[i];

        odo_test.odometryCalculation();
        test_time += odo_test.runtime;
    }

    void loadFirstScanRF2O()
    {
        //Nonlinear version
        for (unsigned int i=0; i<odo_test.width; i++)
            odo_test.range_wf(i) = laser.m_scan.scan[i];
        odo_test.createScanPyramid();
    }

    void setRF2OPose(const CPose2D &reset_pose)
    {
        odo_test.laser_pose = reset_pose;
        odo_test.laser_oldpose = reset_pose;
    }

	void computeErrors(unsigned int react_freq)
	{
        const unsigned int size_v = csm_poses.size();
        CPose3D incr_est, incr_test, incr_psm, incr_csm;
        CPose3D dif_est, dif_test, dif_psm, dif_csm;
        float aver_trans_error_est = 0.f, aver_rot_error_est = 0.f;
        float aver_trans_error_test = 0.f, aver_rot_error_test = 0.f;
        float aver_trans_error_psm = 0.f, aver_rot_error_psm = 0.f;
        float aver_trans_error_csm = 0.f, aver_rot_error_csm = 0.f;
		for (unsigned int i=0; i<size_v-react_freq; i++)
		{
			//Compute increments per second
			incr_est = est_poses.at(i+react_freq) - est_poses.at(i);
			incr_test = test_poses.at(i+react_freq) - test_poses.at(i);
            incr_psm = psm_poses.at(i+react_freq) - psm_poses.at(i);
            incr_csm = csm_poses.at(i+react_freq) - csm_poses.at(i);

			//Compute differences between the estimates and the real one
            dif_est = incr_est;
            dif_test = incr_test;
            dif_psm = incr_psm;
            dif_csm = incr_csm;

			//Add their contribution to the average errors
            aver_trans_error_est += square(dif_est[0]) + square(dif_est[1]);
            aver_trans_error_test += square(dif_test[0]) + square(dif_test[1]);
            aver_trans_error_psm += square(dif_psm[0]) + square(dif_psm[1]);
            aver_trans_error_csm += square(dif_csm[0]) + square(dif_csm[1]);
            aver_rot_error_est += square(57.3f*dif_est.yaw());
            aver_rot_error_test += square(57.3f*dif_test.yaw());
            aver_rot_error_psm += square(57.3f*dif_psm.yaw());
            aver_rot_error_csm += square(57.3f*dif_csm.yaw());
		}

        aver_trans_error_est = sqrtf(aver_trans_error_est/(size_v-react_freq));
        aver_trans_error_test = sqrtf(aver_trans_error_test/(size_v-react_freq));
        aver_trans_error_psm = sqrtf(aver_trans_error_psm/(size_v-react_freq));
        aver_trans_error_csm = sqrtf(aver_trans_error_csm/(size_v-react_freq));
        aver_rot_error_est = sqrtf(aver_rot_error_est/(size_v-react_freq));
        aver_rot_error_test = sqrtf(aver_rot_error_test/(size_v-react_freq));
        aver_rot_error_psm = sqrtf(aver_rot_error_psm/(size_v-react_freq));
        aver_rot_error_csm = sqrtf(aver_rot_error_csm/(size_v-react_freq));

		//Show the errors
        printf("\n Aver_trans_est = %f m, Aver_rot_est = %f grad", aver_trans_error_est, aver_rot_error_est);
        printf("\n Aver_trans_test = %f m, Aver_rot_test = %f grad", aver_trans_error_test, aver_rot_error_test);
        printf("\n Aver_trans_psm = %f m, Aver_rot_psm = %f grad", aver_trans_error_psm, aver_rot_error_psm);
        printf("\n Aver_trans_csm = %f m, Aver_rot_csm = %f grad", aver_trans_error_csm, aver_rot_error_csm);
        printf("\n Runtimes: rf2o_lin = %f, rf2o_rob = %f, psm = %f, csm = %f", est_time/size_v, test_time/size_v, psm_time/size_v, csm_time/size_v);
        fflush(stdout);
	}

    //One of the compared methods (PSM)
    void runPolarScanMatching()
    {
        //Read new scan
        ls_ref = ls;
        pm_readScan(laser.m_scan, &ls);

        CTicTac clock; clock.Tic();

        //preprocess the scans...
        pm_preprocessScan(&ls);

        //Initial seed
        if (matchFailed)
        {
            ls.rx = 0.f; ls_ref.rx = 0.f;
            ls.ry = 0.f; ls_ref.ry = 0.f;
            ls.th = 0.f; ls_ref.th = 0.f;
        }
        else
        {
            ls.rx = 100.f*new_psm_pose[1]; ls_ref.rx = 100.f*old_psm_pose[1];
            ls.ry = 100.f*new_psm_pose[0]; ls_ref.ry = 100.f*old_psm_pose[0];
            ls.th = new_psm_pose[2]; ls_ref.th = old_psm_pose[2];
        }


        matchFailed = false;
        try
        {
            pm_psm(&ls_ref,&ls);
        }catch(int err)
        {
            cerr << "Error caught: failed match." << endl;
            matchFailed = true;
        }

//        PM_TYPE err_idx = pm_error_index2(&ls_ref,&ls);
//        cout <<" err: "<<err_idx;
        //printf("\n Motion (x, y, phi) = (%f, %f, %f)", 0.01f*ls.rx, 0.01f*ls.ry, ls.th);

        const float psm_t = 1000.f*clock.Tac();
        psm_time += psm_t;

//        printf("\nPSM runtime = %f ms", psm_t);
//        fflush(stdout);

        //Saturate values to max permitted increments
        if (abs(ls.ry) > 50.f) ls.ry = 0.f;
        if (abs(ls.rx) > 50.f) ls.rx = 0.f;
        if (abs(ls.th) > DEG2RAD(60.f)) ls.th = 0.f;

        old_psm_pose = new_psm_pose;
        new_psm_pose = old_psm_pose + CPose2D(0.01f*ls.ry,0.01f*ls.rx,ls.th);
    }

    void Init_Params_csm()
    {
        params.max_angular_correction_deg = 90;             //"Maximum angular displacement between scans"
        params.max_linear_correction = 2.0;                 //"Maximum translation between scans (m)"
        params.max_iterations = 1000;                       //"When we had enough"

        params.epsilon_xy = 0.0001;                         //"A threshold for stopping (m)"
        params.epsilon_theta = 0.0001;                      //"A threshold for stopping (rad)"

        params.max_correspondence_dist = 2.0;
        params.sigma = 0.01;                                //Noise in the scan

        params.use_corr_tricks = 1;                         //"Use smart tricks for finding correspondences."
        params.restart = 1;                                 //"Restart: Restart if error is over threshold"
        params.restart_threshold_mean_error = 0.01;
        params.restart_dt= 0.01;
        params.restart_dtheta = 1.5 * 3.14 /180;

        params.clustering_threshold = 0.05;                 //"Max distance for staying in the same clustering"
        params.orientation_neighbourhood = 3;               //"Number of neighbour rays used to estimate the orientation."

        params.use_point_to_line_distance = 1;              //"If 0, it's vanilla ICP."
        params.do_alpha_test = 0;                           //"Discard correspondences based on the angles"
        params.do_alpha_test_thresholdDeg = 20.0;           //

        params.outliers_maxPerc = 0.95;
        params.outliers_adaptive_order =0.7;
        params.outliers_adaptive_mult=2.0;
        params.do_visibility_test = 0;

        params.outliers_remove_doubles = 1;                 //"no two points in laser_sens can have the same corr."
        params.do_compute_covariance = 0;                   //"If 1, computes the covariance of ICP using the method http://purl.org/censi/2006/icpcov ."
        params.debug_verify_tricks = 0;                     //"Checks that find_correspondences_tricks gives the right answer.

        params.laser[0] = 0.0;                              //Pose of sensor with respect to robot: used for computing the first estimate given the odometry.
        params.laser[1] = 0.0;
        params.laser[2] = 0.0;

        params.min_reading = 0.0;                           //Don't use readings less than min_reading (m)
        params.max_reading = 79.0;                        //Don't use readings longer than max_reading (m)

        params.use_ml_weights = 0;                          //"If 1, the field 'true_alpha' (or 'alpha') in the first scan is used to compute the incidence beta, and the factor (1/cos^2(beta)) used to weight the correspondence.");
        params.use_sigma_weights = 0;                       //"If 1, the field 'readings_sigma' in the second scan is used to weight the correspondence by 1/sigma^2"
        recover_from_error = 0;
        printf("[Init_Params]: Params structure sucessfully initialized.\n\n");
    }

    LDP cast_CObservation2DRangeScan_to_LDP(CObservation2DRangeScan &laser_obj)
    {
        int n;
        n = laser_obj.scan.size();		//Num of samples (rays) of the scan laser

        LDP ld = ld_alloc_new(n);       //create new structure initialized with "n" rays

        ld->min_theta = -laser_obj.aperture/2;			//rad
        ld->max_theta = laser_obj.aperture/2;			//rad
        float angle_ray = laser_obj.aperture/(n-1);    //rad

        double index_cero_angle = floor( (n-1)/2);

        //scan readings
        for (int i=0; i<n; i++)
        {
            if (laser_obj.validRange[i] == false)
            {
                ld->readings[i] = NAN;                  //*(ld->readings+i)
                ld->valid[i] = 0;
            }
            else
            {
                ld->readings[i] = static_cast<double>(laser_obj.scan[i]);   //m
                ld->valid[i] = 1;                       //int
            }

            ld->theta[i] = (i-index_cero_angle)*angle_ray;             //rad
            ld->readings_sigma[i] = laser_obj.stdError;
            ld->cluster[i] = -1;                        //No cluster selected
        }

        //To avoid failure due to precission, copy min and max theta (see ld_valid_fields)
        ld->theta[0] = ld->min_theta;
        ld->theta[n-1] = ld->max_theta;

        //Set initial guess (odometry) to 0 (mandatory for proper functionality)
        ld->odometry[0] = 0.0;
        ld->odometry[1] = 0.0;
        ld->odometry[2] = 0.0;

    /*
        mrpt::poses::CPose3D LaserPoseOnTheRobot;
        laser_obj.getSensorPose(LaserPoseOnTheRobot);
        //Set the initial pose
        laser_pose = LaserPoseOnTheRobot;
        laser_oldpose = LaserPoseOnTheRobot;

        jo_read_double_array(jo, "odometry", ld->odometry, 3, NAN);
        jo_read_double_array(jo, "estimate", ld->estimate, 3, NAN);
        jo_read_double_array(jo, "true_pose", ld->true_pose, 3, NAN);
    */

        double time_seconds = mrpt::system::timestampToDouble(laser_obj.timestamp);
        double fractpart, intpart;
        fractpart = modf(time_seconds, &intpart);
        ld->tv.tv_sec = intpart;
        ld->tv.tv_usec = fractpart;

        if(!ld_valid_fields(ld))
        {
            sm_error("[Init] Invalid laser data in first scan.\n");
            return ld;
        }

        //Structure is ready!
        return ld;
    }

    /* Runs the PL-ICP */
    bool runCanonicalScanMatching()
    {
        laser_sens = cast_CObservation2DRangeScan_to_LDP(laser.m_scan);
        CTicTac clock; clock.Tic();

        try
        {
            //set the scan data into the input structure
            params.laser_ref  = laser_ref;
            params.laser_sens = laser_sens;

            // Set first guess as the difference in odometry
            if (any_nan(params.laser_ref->odometry,3) || any_nan(params.laser_sens->odometry,3) )
            {
                sm_error("The 'odometry' field is set to NaN so I don't know how to get an initial guess. I usually use the difference in the odometry fields to obtain the initial guess.\n");
                sm_error("  laser_ref->odometry = %s \n",  friendly_pose(params.laser_ref->odometry) );
                sm_error("  laser_sens->odometry = %s \n", friendly_pose(params.laser_sens->odometry) );
                sm_error(" I will quit it here. \n");
                return false;
            }

            double odometry[3];
            pose_diff_d(laser_sens->odometry, laser_ref->odometry, odometry);
            double ominus_laser[3], temp[3];
            ominus_d(params.laser, ominus_laser);
            oplus_d(ominus_laser, odometry, temp);
            oplus_d(temp, params.laser, params.first_guess);

            // Do the actual work
            sm_icp(&params, &result);

        }
        catch(exception &e)
        {
            cerr << "[Matching] Returning false due to exception: " << endl;
            cerr << e.what() << endl;
            return false;
        }
        catch(...)
        {
            cerr << "[Matching] Returning false due to unknown exception: " << endl;
            return false;
        }

        if(!result.valid)
        {
            if(recover_from_error)
            {
                sm_info("One ICP matching failed. Because you passed  -recover_from_error, I will try to recover."
                " Note, however, that this might not be good in some cases. \n");
                sm_info("The recover is that the displacement is set to 0. No result stats is output. \n");

                /* For the first scan, set estimate = odometry */
                copy_d(laser_ref->estimate, 3, laser_sens->estimate);
                ld_free(laser_ref); laser_ref = laser_sens;
            }
            else
            {
                sm_error("One ICP matching failed. Because I process recursively, I will stop here.\n");
                sm_error("Use the option -recover_from_error if you want to try to recover.\n");
                ld_free(laser_ref);
                return 2;
            }
        }
        else
        {
            /* Add the result to the previous estimate */
            oplus_d(laser_ref->estimate, result.x, laser_sens->estimate);

            old_csm_pose = new_csm_pose;
            new_csm_pose = old_csm_pose + CPose2D(result.x[0],result.x[1],result.x[2]);

            //Free reference scan, and update with the sensed one.
            ld_free(laser_ref); laser_ref = laser_sens;
        }

        const float csm_t = 1000.f*clock.Tac();
        csm_time += csm_t;
//        printf("\nCSM runtime = %f ms \n", csm_t);
    }

    void saveResults()
    {
        ofstream	m_fres;
        char        aux[100];
        int         nFile = 0;
        bool        free_name = false;
        string      name;

        while (!free_name)
        {
            nFile++;
            sprintf(aux, "results_%03u.txt", nFile );
            name = aux;
            free_name = !system::fileExists(name);
        }
        m_fres.open(aux);

        //Save all poses
        for (unsigned int i=0; i<real_poses.size(); i++)
        {
            m_fres << real_poses[i][0] << " " << real_poses[i][1] << " " << real_poses[i][3] << " ";
            m_fres << est_poses[i][0] << " " << est_poses[i][1] << " " << est_poses[i][3] << " ";
            m_fres << test_poses[i][0] << " " << test_poses[i][1] << " " << test_poses[i][3] << " ";
            m_fres << csm_poses[i][0] << " " << csm_poses[i][1] << " " << csm_poses[i][3] << " ";
            m_fres << psm_poses[i][0] << " " << psm_poses[i][1] << " " << psm_poses[i][3] << " ";

            m_fres << "\n";
        }

        m_fres.close();

        printf("\n Results saved in file \n");
    }

};



