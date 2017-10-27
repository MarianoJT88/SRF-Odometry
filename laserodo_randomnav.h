/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */


#include <mrpt/nav/reactive/CReactiveNavigationSystem3D.h>
#include <mrpt/opengl.h>
#include <mrpt/opengl/CPlanarLaserScan.h>
#include <mrpt/maps/COccupancyGridMap2D.h>
#include <mrpt/utils/CRobotSimulator.h>
#include <mrpt/gui.h>
#include <mrpt/utils/round.h>
#include <mrpt/system/filesystem.h>
#include <mrpt/math/lightweight_geom_data.h>

#include "map.xpm"
#include "map_lines_rf2o.xpm"
#include "map_lab_rf2o.xpm"
#include "map_lab.xpm"
#include "laser_odometry_v1.h"
#include "laser_odometry_standard.h"
#include "laser_odometry_3scans.h"
#include "laser_odometry_refscans.h"


using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::nav;
using namespace mrpt::opengl;
using namespace mrpt::maps;
using namespace mrpt::obs;
using namespace mrpt::gui;
using namespace mrpt::poses;
using namespace mrpt::math;
using namespace std;



struct TRobotLaser {
	CObservation2DRangeScan m_scan;
    CObservation2DRangeScan m_scan_old;
	int						m_level;
	int						m_segments;
};


class CMyReactInterface : public CReactiveInterfaceImplementation
{
public:

	CPose2D					new_pose;
	CPose2D					last_pose;
	CPose2D					target;
	CRobotSimulator			robotSim;
	TRobotShape				robotShape;
	COccupancyGridMap2D		map;
	TRobotLaser				laser;
	CDisplayWindow3D		window;
	COpenGLScenePtr			scene;


    //RF2O
    RF2O_RefS    odo;
    RF2O_standard      odo_test; //Should be RF2O or RF2O_standard


    //Results
    vector<CPose3D>	real_poses, est_poses, test_poses;
    vector<float> aver_res, trunc_aver_res, median_res;
    float est_time, test_time;

	
    bool getCurrentPoseAndSpeeds( poses::CPose2D &curPose, float &curV, float &curW)
    {
        robotSim.getRealPose( curPose );
        curV = robotSim.getV();
        curW = robotSim.getW();
        return true;
    }

	bool changeSpeeds( float v, float w )
	{
		robotSim.movementCommand(v,w);
		return true;
	}

	bool senseObstacles( CSimplePointsMap 	&obstacles )
	{
		last_pose = new_pose;
        laser.m_scan_old = laser.m_scan;
		robotSim.getRealPose(new_pose);
		CPose3D robotpose3d(new_pose[0], new_pose[1], 0, new_pose[2], 0, 0);

		obstacles.clear();

		//Laser scans
		map.laserScanSimulator(laser.m_scan, new_pose, 0.5f, laser.m_segments, laser.m_scan.stdError, 1, 0);
		
		//Correct points with null observations (set them to 0)
		for (unsigned int i=0; i<laser.m_scan.scan.size(); i++)
            if (laser.m_scan.scan[i] > 0.995f*laser.m_scan.maxRange)
            {
				laser.m_scan.scan[i] = 0.f;
                laser.m_scan.validRange[i] = false;
            }
            else
            {
                laser.m_scan.validRange[i] = true;
            }
		obstacles.insertObservation(&laser.m_scan);

		return true;
	}

	void loadMaps( const utils::CConfigFileBase &ini )
	{
		CImage myImg;

        //Original synthetic map
        //float resolution = 0.025f; //0.02
        //myImg.loadFromXPM(map_xpm);

        //Synthetic map made of lines
//        float resolution = 0.02f;
//        myImg.loadFromXPM(map_lines_rf2o_xpm);

        //Lab map (new)
        float resolution = 0.045f;
        myImg.loadFromXPM(map_lab_rf2o_xpm);

        map.loadFromBitmap(myImg,resolution);


		std::cout << std::endl << "The map have been loaded successfully.";
	}

	void loadConfiguration( const utils::CConfigFileBase &ini )
	{
		unsigned int num_levels;
		std::vector <float> lasercoord, xaux, yaux;

		//Read laser params
		ini.read_vector("LASER_CONFIG","LASER_POSE", std::vector<float> (0), lasercoord , true);
		laser.m_scan.maxRange = ini.read_float("LASER_CONFIG","LASER_MAX_RANGE", 50, true);
		laser.m_scan.aperture = DEG2RAD(ini.read_float("LASER_CONFIG","LASER_APERTURE", 180, true));
		laser.m_scan.stdError = ini.read_float("LASER_CONFIG","LASER_STD_ERROR", 0.05, true);
		laser.m_scan.sensorPose.setFromValues(lasercoord[0],lasercoord[1],lasercoord[2],lasercoord[3],lasercoord[4],lasercoord[5]);
		laser.m_level = ini.read_int("LASER_CONFIG","LASER_LEVEL", 1, true);
		laser.m_segments = ini.read_int("LASER_CONFIG","LASER_SEGMENTS", 181, true);

		//Read config params which describe the robot shape
		num_levels = ini.read_int("ROBOT_CONFIG","HEIGHT_LEVELS", 1, true);
		robotShape.polygons.resize(num_levels);
		robotShape.heights.resize(num_levels);
		for (unsigned int i=1;i<=num_levels;i++)
		{
			robotShape.heights[i-1] = ini.read_float("ROBOT_CONFIG",format("LEVEL%d_HEIGHT",i), 1, true);
			ini.read_vector("ROBOT_CONFIG",format("LEVEL%d_VECTORX",i), std::vector<float> (0), xaux, false);
			ini.read_vector("ROBOT_CONFIG",format("LEVEL%d_VECTORY",i), std::vector<float> (0), yaux, false);
			ASSERT_(xaux.size() == yaux.size());
			for (unsigned int j=0;j<xaux.size();j++)
			{
				robotShape.polygons[i-1].AddVertex(xaux[j], yaux[j]);
			}
		}

		//Read other params associated with the robot model and its navigation
		float tau = ini.read_float("NAVIGATION_CONFIG","ROBOTMODEL_TAU", 0, true);
		float delay = ini.read_float("NAVIGATION_CONFIG","ROBOTMODEL_DELAY", 0, true);
        float x_ini, y_ini;
		robotSim.setDelayModelParams(tau, delay);
		robotSim.resetStatus();
		robotSim.setOdometryErrors(0);

        //Generate random initial pose within the map
        bool valid_ini_pose = false;
        while (!valid_ini_pose)
        {
            const float rand_num_x = float(std::rand()%1000)/1000.f;
            const float rand_num_y = float(std::rand()%1000)/1000.f;
            x_ini = map.getXMin() + rand_num_x*(map.getXMax() - map.getXMin());
            y_ini = map.getYMin() + rand_num_y*(map.getYMax() - map.getYMin());

            //Check the candidate initial pose
            const float range = 0.5f;
            valid_ini_pose = true;
            for (float u=-0.5f; u<=0.5f; u+= 0.1f)
                for (float v=-0.5f; v<=0.5f; v+= 0.1f)
                    if (map.getPos(x_ini + u*range, y_ini + v*range) < 0.7f)
                            valid_ini_pose = false;
        }

        robotSim.setRealPose(CPose2D(x_ini, y_ini, 0.f));
	}

    void initilizeEverything()
    {
        odo.initialize(laser.m_segments, laser.m_scan.aperture, false);
        odo_test.initialize(laser.m_segments, laser.m_scan.aperture, 3);
        initializeScene();
        loadFirstScanRF2O();
        setRF2OPose(new_pose);
        est_time = 0.f; test_time = 0.f;
    }

	void initializeScene()
	{
		CPose3D robotpose3d;
		robotpose3d.x(robotSim.getX());
		robotpose3d.y(robotSim.getY());
		robotpose3d.setYawPitchRoll(robotSim.getPHI(),0,0);

		//The display window is created
		mrpt::global_settings::OCTREE_RENDER_MAX_POINTS_PER_NODE = 10000;
		window.setWindowTitle("Reactive Navigation. Robot motion simulation");
        window.resize(1200,800);
        window.setPos(800,0);
        window.setCameraZoom(22);
        window.setCameraAzimuthDeg(-90);
        window.setCameraElevationDeg(90);
		scene = window.get3DSceneAndLock();

		//Map is inserted
		CSetOfObjectsPtr gl_grid = CSetOfObjects::Create();
		map.getAs3DObject(gl_grid);
		scene->insert(gl_grid);

		//A CornerXYZ object is inserted as an absolute frame of reference
		//CSetOfObjectsPtr ref = opengl::stock_objects::CornerXYZ();
		//ref->setLocation(0,0,0);
		//scene->insert( ref );

		//The target is inserted
		CDiskPtr target = opengl::CDisk::Create(0.4, 0.3);
		target->setLocation(0, 0, 0);
		target->setColor(0.2,0.3,0.9);
		scene->insert( target );

		//Robots
		//CPolyhedronPtr robot_real;
		//robot_real = opengl::CPolyhedron::CreateCustomPrism(robotShape.polygons[0], robotShape.heights[0]);
		//robot_real->setName("robot_real");
		//robot_real->setPose(robotpose3d);
  //      robot_real->setColor(0.f,0.f,0.f);
		//robot_real->setWireframe(true);
		//robot_real->setLineWidth(2);
		//scene->insert( robot_real );

  //      CPolyhedronPtr robot_est;
  //      robot_est = opengl::CPolyhedron::CreateCustomPrism(robotShape.polygons[0], robotShape.heights[0]);
  //      robot_est->setName("robot_est");
  //      robot_est->setPose(robotpose3d);
  //      robot_est->setColor(1,0.4,0);
  //      robot_est->setWireframe(true);
  //      robot_est->setLineWidth(2);
  //      scene->insert( robot_est );

		//CPolyhedronPtr robot_test;
		//robot_test = opengl::CPolyhedron::CreateCustomPrism(robotShape.polygons[0], robotShape.heights[0]);
		//robot_test->setName("robot_test");
		//robot_test->setPose(robotpose3d);
  //      robot_test->setColor(0.f,0.8f,0.f);
		//robot_test->setWireframe(true);
		//robot_test->setLineWidth(2);
		//scene->insert( robot_test );


        //Initialization of the scans. One scan is simulated
        //-------------------------------------------------------------------------
		CSimplePointsMap auxpoints;
		senseObstacles( auxpoints );
        laser.m_scan_old = laser.m_scan;
        //--------------------------------------------------------------------------

		//The laserscan is inserted
		//CPointCloudColouredPtr gl_laser = CPointCloudColoured::Create();
		//gl_laser->setPointSize(5.f);
		//gl_laser->setPose(robotpose3d);
		//scene->insert( gl_laser );

		CPlanarLaserScanPtr gl_scan = CPlanarLaserScan::Create();
		gl_scan->enableLine(true);
		gl_scan->enableSurface(true);
		gl_scan->enablePoints(true);
		gl_scan->setName(format("Laser"));
		gl_scan->setScan(laser.m_scan);
		scene->insert(gl_scan);

		//Trajectories
		CSetOfLinesPtr traj_lines_real = opengl::CSetOfLines::Create();
        traj_lines_real->setColor(0,0,0);
		traj_lines_real->setLineWidth(4);
		scene->insert( traj_lines_real );
        opengl::CSetOfLinesPtr traj_lines_est = opengl::CSetOfLines::Create();
        traj_lines_est->setColor(1,0.4,0);
        traj_lines_est->setLineWidth(4);
        scene->insert( traj_lines_est );
		opengl::CSetOfLinesPtr traj_lines_test = opengl::CSetOfLines::Create();
        traj_lines_test->setColor(0,0.8,0);
		traj_lines_test->setLineWidth(4);
		scene->insert( traj_lines_test );


		window.unlockAccess3DScene();
		window.repaint();
	}

	void updateScene()
	{
		scene = window.get3DSceneAndLock();
		CPose3D robotpose3d;
		CRenderizablePtr obj;
		robotpose3d.x(robotSim.getX());
		robotpose3d.y(robotSim.getY());
		robotpose3d.setYawPitchRoll(robotSim.getPHI(),0,0);

		//Robots
		//obj = scene->getByName("robot_real");
		//obj->setPose(robotpose3d);

  //      obj = scene->getByName("robot_est");
  //      obj->setPose(odo.laser_pose);

		//obj = scene->getByName("robot_test");
		//obj->setPose(odo_test.laser_pose);


		const unsigned int repr_level = round(log2(round(float(odo.width)/float(odo.cols))));

		//Laser
		CPose3D laserpose;
		laser.m_scan.getSensorPose(laserpose);
		//CPointCloudColouredPtr gl_laser;
		//gl_laser = scene->getByClass<CPointCloudColoured> (0);
		//gl_laser->clear();
		//gl_laser->setPose(robotpose3d);
  //      for (unsigned int i=0; i<odo.cols; i++)
  //      {
  //          if (odo.outliers(i) == true)
  //              gl_laser->push_back(odo_test.xx[repr_level](i), odo_test.yy[repr_level](i), 0.1, 0, 0, 1);
  //          else
  //              gl_laser->push_back(odo_test.xx[repr_level](i), odo_test.yy[repr_level](i), 0.1, 1-sqrt(odo_test.weights(i)), sqrt(odo_test.weights(i)), 0);
  //      }



		CPlanarLaserScanPtr gl_scan = scene->getByClass<CPlanarLaserScan>(0);
		gl_scan->setScan(laser.m_scan);
		gl_scan->setPose(robotpose3d);



		////Trajectories
		//opengl::CSetOfLinesPtr traj_lines_real;
		//traj_lines_real = scene->getByClass<CSetOfLines> (0);
		//traj_lines_real->appendLine(last_pose[0], last_pose[1], 0.2, new_pose[0], new_pose[1], 0.2);

  //      opengl::CSetOfLinesPtr traj_lines_est;
  //      traj_lines_est = scene->getByClass<CSetOfLines> (1);
  //      traj_lines_est->appendLine(odo.laser_oldpose[0], odo.laser_oldpose[1], 0.2, odo.laser_pose[0], odo.laser_pose[1], 0.2);

		//opengl::CSetOfLinesPtr traj_lines_test;
		//traj_lines_test = scene->getByClass<CSetOfLines> (2);
		//traj_lines_test->appendLine(odo_test.laser_oldpose[0], odo_test.laser_oldpose[1], 0.2, odo_test.laser_pose[0], odo_test.laser_pose[1], 0.2);


		window.unlockAccess3DScene();
		window.repaint();
	}

	void resetScene()
	{
		scene = window.get3DSceneAndLock();
		CPose3D robotpose3d;
		CRenderizablePtr obj;
		robotpose3d.x(robotSim.getX());
		robotpose3d.y(robotSim.getY());
		robotpose3d.setYawPitchRoll(robotSim.getPHI(),0,0);

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


		//Robots
		obj = scene->getByName("robot_real");
		obj->setPose(robotpose3d);

        obj = scene->getByName("robot_est");
        obj->setPose(odo.laser_pose);

		obj = scene->getByName("robot_test");
		obj->setPose(odo_test.laser_pose);


		//Laser
		const unsigned int repr_level = round(log2(round(float(odo.width)/float(odo.cols))));

		CPose3D laserpose;
		laser.m_scan.getSensorPose(laserpose);
		CPointCloudColouredPtr gl_laser;
		gl_laser = scene->getByClass<CPointCloudColoured> (0);
		gl_laser->clear();
		for (unsigned int i=0; i<odo.cols; i++)
            gl_laser->push_back(odo_test.xx[repr_level](i), odo_test.yy[repr_level](i), 0.1, 1-sqrt(odo_test.weights(i)), sqrt(odo_test.weights(i)), 0);

		gl_laser->setPose(robotpose3d + laserpose);

		//Trajectories
		opengl::CSetOfLinesPtr traj_lines_real;
		traj_lines_real = scene->getByClass<CSetOfLines> (0);
		traj_lines_real->clear();

        opengl::CSetOfLinesPtr traj_lines_est;
        traj_lines_est = scene->getByClass<CSetOfLines> (1);
        traj_lines_est->clear();

		opengl::CSetOfLinesPtr traj_lines_test;
		traj_lines_test = scene->getByClass<CSetOfLines> (2);
		traj_lines_test->clear();


		window.unlockAccess3DScene();
		window.repaint();
	}


	CAbstractReactiveNavigationSystem::TNavigationParams createNewTarget(float x,
					float y, float targetAllowedDistance, bool targetIsRelative = 0)
	{
		CPose2D robotpose;
		CAbstractReactiveNavigationSystem::TNavigationParams navparams;
		navparams.target = mrpt::math::TPoint2D(x,y);
		navparams.targetAllowedDistance = targetAllowedDistance;
		navparams.targetIsRelative = targetIsRelative;
		if (targetIsRelative == 0)
			target = CPose2D(x, y, 0);
		else
		{
			robotSim.getRealPose(robotpose);
			target = CPose2D(x ,y, 0) + robotpose;
		}
		return navparams;
	}

    void runRF2O()
    {
        //Run the first version
        for (unsigned int i=0; i<odo.width; i++)
            odo.range_wf(i) = laser.m_scan.scan[i];

        odo.odometryCalculation();
        est_time += odo.runtime;

//        float res1, res2, res3;
//        odo.computeAverageResiduals(res1, res2, res3);
//        aver_res.push_back(res1);
//        trunc_aver_res.push_back(res2);
//        median_res.push_back(res3);

        //Run the second version
        for (unsigned int i=0; i<odo_test.width; i++)
            odo_test.range_wf(i) = laser.m_scan.scan[i];

        odo_test.odometryCalculation();
        test_time += odo_test.runtime;
    }

    void loadFirstScanRF2O()
    {
        //First version
        for (unsigned int i=0; i<odo.width; i++)
            odo.range_wf(i) = laser.m_scan.scan[i];
        odo.createScanPyramid();

        //Second version
        for (unsigned int i=0; i<odo_test.width; i++)
            odo_test.range_wf(i) = laser.m_scan.scan[i];
        odo_test.createScanPyramid();
    }

    void setRF2OPose(const CPose2D &reset_pose)
    {
        odo.laser_pose = reset_pose;
        odo.laser_oldpose = reset_pose;

        odo_test.laser_pose = reset_pose;
        odo_test.laser_oldpose = reset_pose;
    }

	void computeErrors(unsigned int react_freq)
	{
		const unsigned int size_v = real_poses.size();
        CPose3D incr_real, incr_est, incr_test;
        CPose3D dif_est, dif_test;
        float aver_trans_error_est = 0.f, aver_rot_error_est = 0.f;
        float aver_trans_error_test = 0.f, aver_rot_error_test = 0.f;

		for (unsigned int i=0; i<size_v-react_freq; i++)
		{
			//Compute increments per second
			incr_real = real_poses.at(i+react_freq) - real_poses.at(i);
			incr_est = est_poses.at(i+react_freq) - est_poses.at(i);
			incr_test = test_poses.at(i+react_freq) - test_poses.at(i);

			//Compute differences between the estimates and the real one
			dif_est = incr_real - incr_est;
			dif_test = incr_real - incr_test;

			//Add their contribution to the average errors
			aver_trans_error_est += sqrtf(square(dif_est[0]) + square(dif_est[1]));
			aver_trans_error_test += sqrtf(square(dif_test[0]) + square(dif_test[1]));
			aver_rot_error_est += 57.3f*abs(dif_est.yaw());
			aver_rot_error_test += 57.3f*abs(dif_test.yaw());
		}

		//Show the errors
        printf("\n Aver_trans_est = %f m, Aver_rot_est = %f grad", aver_trans_error_est/(size_v-react_freq), aver_rot_error_est/(size_v-react_freq));
        printf("\n Aver_trans_test = %f m, Aver_rot_test = %f grad", aver_trans_error_test/(size_v-react_freq), aver_rot_error_test/(size_v-react_freq));
        printf("\n Runtimes: rf2o_lin = %f, rf2o_rob = %f", est_time/size_v, test_time/size_v);
        fflush(stdout);
	}

    void saveResults(unsigned int freq)
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


        //Save freq, real_pose (x,y,theta), est_pose (x,y,theta)
        for (unsigned int i=0; i<real_poses.size(); i++)
        {
            m_fres << freq << " ";
            m_fres << real_poses[i][0] << " " << real_poses[i][1] << " " << real_poses[i][3] << " ";
            m_fres << est_poses[i][0] << " " << est_poses[i][1] << " " << est_poses[i][3] << " ";
           //m_fres << aver_res[i] << " " << trunc_aver_res[i] << " " << median_res[i] << " ";
            m_fres << "\n";
        }

        m_fres.close();

        printf("\n Results saved in file \n");
    }
};



