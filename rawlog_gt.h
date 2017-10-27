/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */


#include <mrpt/obs/CRawlog.h>
#include <mrpt/obs/CObservationOdometry.h>
#include <mrpt/opengl.h>
#include <mrpt/opengl/CPlanarLaserScan.h>
#include <mrpt/maps/COccupancyGridMap2D.h>
#include <mrpt/gui.h>
#include <mrpt/utils/round.h>


using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::opengl;
using namespace mrpt::maps;
using namespace mrpt::obs;
using namespace mrpt::gui;
using namespace mrpt::poses;
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


    //Results
    vector<CPose3D>	real_poses, est_poses, test_poses, psm_poses, csm_poses;
    float est_time, test_time, psm_time, csm_time;
    CPose3D pose_ini;


    void initializeEverything()
    {
        //Read Rawlog
        string folder = "C:/Users/jaimez/Dropbox/LaserOdo - shared/Datasets/";
        //string folder = "/home/mariano/Dropbox/LaserOdo - shared/Datasets/";

        dataset_id = 3;
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
            dset = "Freiburg_Indoor/freiburg_gt.rawlog";
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180); laser.m_segments = 360; laser.m_scan.maxRange = 80.f;
            break;
        case 5:
            dset = "MIT_CSAIL/mit_csail_gt.rawlog";
            laser_min_range = 0.05f; laser.m_scan.aperture = DEG2RAD(180.5); laser.m_segments = 361; laser.m_scan.maxRange = 80.f;
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

        est_time = 0.f; test_time = 0.f; psm_time = 0.f; csm_time = 0.f;

        initializeScene();
    }

	

    void readScanRawlog()
	{
        CObservationPtr alfa = dataset.getAsObservation(rawlog_count);

        while (!IS_CLASS(alfa, CObservation2DRangeScan))
        {
            rawlog_count++;
            if (dataset.size() <= rawlog_count)
            {
                dataset_finished = true;
                return;
            }

            if (IS_CLASS(alfa, CObservationOdometry) && (alfa->sensorLabel == "ODOMETRY"))
            {
                CObservationOdometryPtr obs_odo = CObservationOdometryPtr(alfa);
                last_pose = new_pose;

                if (!localized)
                {
                    pose_ini = obs_odo->odometry;
                    localized = true;
                }
                else
                {
                    new_pose = obs_odo->odometry; // - pose_ini;
                    updatePoseInMap();
                }

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
        window.resize(2200,1400);
        window.setPos(1400,0);
        window.setCameraZoom(45);
        window.setCameraAzimuthDeg(-90);
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


        //The laserscan is inserted
        CPointCloudColouredPtr gl_laser = CPointCloudColoured::Create();
        gl_laser->setPointSize(4.f);
        gl_laser->setPose(robotpose3d);
        gl_laser->enablePointSmooth(true);
        scene->insert( gl_laser );


        //Trajectories
        CSetOfLinesPtr traj_lines = opengl::CSetOfLines::Create();
        traj_lines->setColor(1,0,0);
        traj_lines->setLineWidth(4);
        scene->insert( traj_lines );


        window.unlockAccess3DScene();
        window.repaint();
    }

    void updateScene(unsigned int iter)
    {

        //Stop the program to move the scene before taking the snapshot
        const float threshold = float(dataset.size())/5.f;
        if (floor(float(rawlog_count)/threshold) != floor(float(rawlog_count-2)/threshold))
            system::os::getch();


        scene = window.get3DSceneAndLock();

        //Camera follows the robot
        window.setCameraPointingToPoint(new_pose[0], new_pose[1], 0.f);

        //Robot
//        CRenderizablePtr robot_real = scene->getByName("laser");
//        robot_real->setPose(new_pose);

        //Laser
        CPointCloudColouredPtr gl_laser;
        gl_laser = scene->getByClass<CPointCloudColoured> (0);

        //Get color
//            Eigen::Vector3f ini_c; ini_c << 1.f,0.f,0.f;
//            Eigen::Vector3f end_c; end_c << 0.f,0.f,1.f;

//            const float alpha = float(dataset.size() - rawlog_count)/float(dataset.size());
//            Eigen::Vector3f c_now = alpha*end_c + (1.f - alpha)*ini_c;
//            float r,g,b;
//            utils::jet2rgb(alpha, r, g, b);

//        if (floor(float(rawlog_count)/threshold) != floor(float(rawlog_count-2)/threshold))
            gl_laser->clear();

        //if ((iter+4) % 5 == 0)
        {
            for (unsigned int i=0; i<laser.m_segments; i++)
                if ((laser.m_scan.scan[i] < 15.8f)&&(laser.m_scan.scan[i] > 0.05f))
                {
                    const float tita = -0.5f*laser.m_scan.aperture + (float(i) + 0.5f)*laser.m_scan.aperture/float(laser.m_segments);
                    const float x = laser.m_scan.scan[i]*cos(tita);
                    const float y = laser.m_scan.scan[i]*sin(tita);
                    const float x_trans = x*cos(new_pose[2]) - y*sin(new_pose[2]) + new_pose[0];
                    const float y_trans = x*sin(new_pose[2]) + y*cos(new_pose[2]) + new_pose[1];
                    gl_laser->push_back(x_trans, y_trans, 0.f, 0.f, 0.f, 1.f);
                }
        }

        //Trajectories
        CSetOfLinesPtr traj_lines = scene->getByClass<CSetOfLines>(0);

        if (floor(float(rawlog_count)/threshold) != floor(float(rawlog_count-2)/threshold))
        {
            system::os::getch();
            traj_lines->clear();
        }

        traj_lines->appendLine(last_pose[0], last_pose[1], 0.01, new_pose[0], new_pose[1], 0.01);

        window.unlockAccess3DScene();
        window.repaint();
    }

    void updatePoseInMap()
    {

        scene = window.get3DSceneAndLock();

        //Camera follows the robot
        window.setCameraPointingToPoint(new_pose[0], new_pose[1], 0.f);

        //Trajectories
        CSetOfLinesPtr traj_lines = scene->getByClass<CSetOfLines>(0);

        traj_lines->appendLine(last_pose[0], last_pose[1], 0.01, new_pose[0], new_pose[1], 0.01);

        window.unlockAccess3DScene();
        window.repaint();
    }

    void clearScene()
    {

        //Update scene
        scene = window.get3DSceneAndLock();

        //Laser
        CPointCloudColouredPtr gl_laser;
        gl_laser = scene->getByClass<CPointCloudColoured> (0);
        gl_laser->clear();


        //Trajectories
        CSetOfLinesPtr traj_lines = scene->getByClass<CSetOfLines>(0);
        traj_lines->clear();

        window.unlockAccess3DScene();
        window.repaint();
    }

	void resetScene()
	{
		scene = window.get3DSceneAndLock();
        CPose3D robotpose3d = new_pose;
		CRenderizablePtr obj;

		//Robots
        obj = scene->getByName("laser");
        obj->setPose(robotpose3d);

		//Laser
//		const unsigned int repr_level = round(log2(round(float(odo.width)/float(odo.cols))));
//		CPointCloudColouredPtr gl_laser;
//		gl_laser = scene->getByClass<CPointCloudColoured> (0);
//		gl_laser->clear();
//		for (unsigned int i=0; i<odo.cols; i++)
//			gl_laser->push_back(odo.xx[repr_level](i), odo.yy[repr_level](i), 0.1, 1-sqrt(odo.weights(i)), sqrt(odo.weights(i)), 0);

//        gl_laser->setPose(robotpose3d);

		//Trajectories
		opengl::CSetOfLinesPtr traj_lines_real;
		traj_lines_real = scene->getByClass<CSetOfLines> (0);
        traj_lines_real->clear();


		window.unlockAccess3DScene();	
		window.repaint();
	}



	void computeErrors(unsigned int react_freq)
	{
//        const unsigned int size_v = csm_poses.size();
//        CPose3D incr_est, incr_test, incr_psm, incr_csm;
//        CPose3D dif_est, dif_test, dif_psm, dif_csm;
//        float aver_trans_error_est = 0.f, aver_rot_error_est = 0.f;
//        float aver_trans_error_test = 0.f, aver_rot_error_test = 0.f;
//        float aver_trans_error_psm = 0.f, aver_rot_error_psm = 0.f;
//        float aver_trans_error_csm = 0.f, aver_rot_error_csm = 0.f;
//		for (unsigned int i=0; i<size_v-react_freq; i++)
//		{
//			//Compute increments per second
//			incr_est = est_poses.at(i+react_freq) - est_poses.at(i);
//			incr_test = test_poses.at(i+react_freq) - test_poses.at(i);
//            incr_psm = psm_poses.at(i+react_freq) - psm_poses.at(i);
//            incr_csm = csm_poses.at(i+react_freq) - csm_poses.at(i);

//			//Compute differences between the estimates and the real one
//            dif_est = incr_est;
//            dif_test = incr_test;
//            dif_psm = incr_psm;
//            dif_csm = incr_csm;

//			//Add their contribution to the average errors
//            aver_trans_error_est += square(dif_est[0]) + square(dif_est[1]);
//            aver_trans_error_test += square(dif_test[0]) + square(dif_test[1]);
//            aver_trans_error_psm += square(dif_psm[0]) + square(dif_psm[1]);
//            aver_trans_error_csm += square(dif_csm[0]) + square(dif_csm[1]);
//            aver_rot_error_est += square(57.3f*dif_est.yaw());
//            aver_rot_error_test += square(57.3f*dif_test.yaw());
//            aver_rot_error_psm += square(57.3f*dif_psm.yaw());
//            aver_rot_error_csm += square(57.3f*dif_csm.yaw());
//		}

//        aver_trans_error_est = sqrtf(aver_trans_error_est/(size_v-react_freq));
//        aver_trans_error_test = sqrtf(aver_trans_error_test/(size_v-react_freq));
//        aver_trans_error_psm = sqrtf(aver_trans_error_psm/(size_v-react_freq));
//        aver_trans_error_csm = sqrtf(aver_trans_error_csm/(size_v-react_freq));
//        aver_rot_error_est = sqrtf(aver_rot_error_est/(size_v-react_freq));
//        aver_rot_error_test = sqrtf(aver_rot_error_test/(size_v-react_freq));
//        aver_rot_error_psm = sqrtf(aver_rot_error_psm/(size_v-react_freq));
//        aver_rot_error_csm = sqrtf(aver_rot_error_csm/(size_v-react_freq));

//		//Show the errors
//        printf("\n Aver_trans_est = %f m, Aver_rot_est = %f grad", aver_trans_error_est, aver_rot_error_est);
//        printf("\n Aver_trans_test = %f m, Aver_rot_test = %f grad", aver_trans_error_test, aver_rot_error_test);
//        printf("\n Aver_trans_psm = %f m, Aver_rot_psm = %f grad", aver_trans_error_psm, aver_rot_error_psm);
//        printf("\n Aver_trans_csm = %f m, Aver_rot_csm = %f grad", aver_trans_error_csm, aver_rot_error_csm);
//        printf("\n Runtimes: rf2o_lin = %f, rf2o_rob = %f, psm = %f, csm = %f", est_time/size_v, test_time/size_v, psm_time/size_v, csm_time/size_v);
//        fflush(stdout);
	}
};



