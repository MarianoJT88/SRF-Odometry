/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */


#include <iostream>
#include <mrpt/system/threads.h> // sleep()
#include <mrpt/utils/CConfigFile.h>
#include <mrpt/utils/CConfigFileMemory.h>
#include "laserodo_randomnav.h"


const char *default_cfg_txt =
	"; ---------------------------------------------------------------\n"
	"; FILE: Reactive Parameters.txt\n"
	";\n"
	";  MJT @ JUIN-2013\n"
	"; ---------------------------------------------------------------\n\n\n"

	"[ROBOT_CONFIG]\n"

	"Name = MyRobot\n\n"

	"HEIGHT_LEVELS = 1 \n\n"	//Only one level works with this simulator + odometry!!!!!

	";Indicate the geometry of each level \n\n"

	";Type: Polyhedric 	(LEVELX_HEIGHT, LEVELX_VECTORX, LEVELX_VECTORY) \n\n"

	"LEVEL1_HEIGHT = 0.6 \n"
	"LEVEL1_VECTORX = -0.2 0.4 -0.2 \n"
	"LEVEL1_VECTORY = -0.3 0 0.3 \n\n"

	"[LASER_CONFIG] \n\n"

	";Indicate the laser parameters. This information must be consistent with that included before \n"
	";Laser pose is relative to the robot coordinate system. \n"
	";Information required: 	LASERX_POSE, LASERY_POSE, LASERX_MAX_RANGE, LASERX_APERTURE \n"
	";							LASERX_STD_ERROR, LASERX_LEVEL, LASERX_SEGMENTS \n\n"

    "LASER_POSE = 0 0 0.05 0 0 0 \n"
    "LASER_MAX_RANGE = 30 \n"
    "LASER_APERTURE = 270 \n"
    "LASER_STD_ERROR = 0.01 \n"
	"LASER_LEVEL = 1 \n"
    "LASER_SEGMENTS = 1080 \n\n"


	"[NAVIGATION_CONFIG] \n\n"

	"; 0: VFF,  1: ND \n"
	"HOLONOMIC_METHOD = 1 \n\n"


	";	Parameters for the navigation \n"
	"; ---------------------------------------------------- \n\n"

	"weights = 0.5 0.05 0.5 2.0 0.5 0.3 \n\n"

	"; 1: Free space \n"
	"; 2: Dist. in sectors \n"
	"; 3: Heading toward target \n"
	"; 4: Closer to target (euclidean) \n"
	"; 5: Hysteresis \n"
	"; 6: Security Distance \n\n"

	"DIST_TO_TARGET_FOR_SENDING_EVENT = 0.5	; Minimum distance to target for sending the end event. Set to 0 to send it just on navigation end \n\n"

	"VMAX_MPS = 0.70			; Speed limits - mps \n"
	"WMAX_DEGPS = 60			; dps \n"
    "SPEEDFILTER_TAU = 0		; The 'TAU' time constant of a first order lowpass filter for speed commands (s) \n"
	"ROBOTMODEL_DELAY = 0		; The delay until motor reaction (s) \n"
    "ROBOTMODEL_TAU = 0 		; The 'TAU' time constant of a first order robot model (s) \n"
	"MAX_DISTANCE_PTG = 2		; Marks the maximum distance regarded by the reactive navigator (m) \n"
	"GRID_RESOLUTION = 0.02 	; Resolutions used to build the collision_grid \n\n\n"



	";	PTGs	.All of them has the same fields to fill, but they don't use all of them. \n"
	";----------------------------------------------------------------------------------- \n"
	";	Types:	1 - Circular arcs \n"
	";			2 - alpha - A, Trajectories with asymptotical heading \n"
	";			3 - C|C,S, R = vmax/wmax, Trajectories to move backward and then forward \n"
	";			4 - C|C,s, like PTG 3, but if t > threshold -> v = w = 0 \n"
	";			5 - CS, Trajectories with a minimum turning radius \n"
	";			6 - alpha - SP, Trajectories built upon a spiral segment \n"
	";			7 - \n\n"


    "PTG_COUNT = 3			;Number of path models used \n\n"

	"PTG1_TYPE = 1 \n"
	"PTG1_NALFAS = 121 \n"
	"PTG1_VMAX = 0.5 \n"
	"PTG1_WMAX = 45 \n"
	"PTG1_K = 1 \n"
	"PTG1_AV = 57.3 \n"
	"PTG1_AW = 57.3 \n\n"

	"PTG2_TYPE = 2 \n"
	"PTG2_NALFAS = 121 \n"
	"PTG2_VMAX = 0.5 \n"
    "PTG2_WMAX = 45 \n"
	"PTG2_K = 1.0 \n"
	"PTG2_AV = 57.3 \n"
	"PTG2_AW = 57.3 \n\n"

	"PTG3_TYPE = 5 \n"
	"PTG3_NALFAS = 121 \n"
	"PTG3_VMAX = 0.5 \n"
	"PTG3_WMAX = 45 \n"
	"PTG3_K = 1.0 \n"
	"PTG3_AV = 57.3 \n"
	"PTG3_AW = 57.3 \n\n"



	";	Parameters for the 'Nearness diagram' Holonomic method \n"
	"; ------------------------------------------------------------ \n\n"

	"[ND_CONFIG] \n"
	"factorWeights = 1.0 2.0 0.5 1.0 \n"
	"; 1: Free space \n"
	"; 2: Dist. in sectors \n"
	"; 3: Closer to target (euclidean) \n"
	"; 4: Hysteresis \n"

	"WIDE_GAP_SIZE_PERCENT = 0.25			; The robot travels nearer to obstacles if this parameter is small. \n"
	"										; The smaller it is, the closer the selected direction is respect to \n"
	"										; the Target direction in TP-Space (under some conditions) \n"
	"MAX_SECTOR_DIST_FOR_D2_PERCENT = 0.25	; \n"
	"RISK_EVALUATION_SECTORS_PERCENT = 0.25	; \n"
	"RISK_EVALUATION_DISTANCE = 0.7			; Parameter used to decrease speed if obstacles are closer than this threshold \n"
	"										; in normalized ps-meters [0,1] \n"
	"TARGET_SLOW_APPROACHING_DISTANCE = 0.8	; Used to decrease speed gradually when the target is going to be reached \n"
    "TOO_CLOSE_OBSTACLE = 0.05				; In normalized ps-meters [0,1] \n\n\n"


	";	Parameters for the VFF Holonomic method \n"
	"; ------------------------------------------------------------ \n\n"

	"[VFF_CONFIG] \n\n"

	"TARGET_SLOW_APPROACHING_DISTANCE = 0.8	; Used to decrease speed gradually when the target is going to be reached \n"
	"TARGET_ATTRACTIVE_FORCE = 7.5			; Use it to control the relative weight of the target respect to the obstacles \n\n\n";



// ------------------------------------------------------
//						MAIN
// ------------------------------------------------------


int main()
{

    //Initial steps. Load configuration
    //------------------------------------------------------
    CMyReactInterface ReactInterface;
    CReactiveNavigationSystem3D rn3d (ReactInterface, false, false);

    utils::CConfigFileMemory configNavigation(default_cfg_txt);
    rn3d.loadConfigFile( configNavigation );
    ReactInterface.loadMaps( configNavigation );
    ReactInterface.loadConfiguration( configNavigation );


    //Initialize all methods
    ReactInterface.initilizeEverything();
    rn3d.initialize();


    //Main loop
    //-----------------------------------------------------
    bool use_reactive = true;
    bool stop = 0;
    bool working = false;
    bool one_step = false;
    int pushed_key = 0;
    int iter_count = 0;
    TPoint3D last_Target_Pos(0,0,0);
    unsigned int odo_freq = 5;
    unsigned int react_sim_per_est = 1;
    float sim_period = 1.f/float(odo_freq*react_sim_per_est);
    CTicTac target_clock;
    CPose2D old_incr = CPose2D(0,0,0);
    target_clock.Tic();


    while (!stop)
    {
        if (ReactInterface.window.keyHit())
            pushed_key = ReactInterface.window.getPushedKey();
        else
            pushed_key = 0;

        switch (pushed_key) {

        case 's':
            //Start/stop continuous estimation
            working = !working;
            break;

        case 'n':
            //Estimate next step
            one_step = true;
            break;

        case 'e':
            //Exit program
            stop = 1;
            rn3d.cancel();
            break;

        case 'r':
            //Compute statistics (deviations of the odometry from ground truth)
            ReactInterface.computeErrors(odo_freq);
            break;

        case 'c':
            //Reset pose estimation
            ReactInterface.resetScene();
            ReactInterface.real_poses.clear();
            ReactInterface.est_poses.clear(); ReactInterface.est_time = 0.f;
            ReactInterface.test_poses.clear(); ReactInterface.test_time = 0.f;
            ReactInterface.aver_res.clear(); ReactInterface.trunc_aver_res.clear(); ReactInterface.median_res.clear();
            break;
        }


        //						Run reactive + odometry
        //==============================================================================
        if ((working)||(one_step))
        {
            //Execute navigation
            if (use_reactive)
            {
                ReactInterface.robotSim.simulateInterval( sim_period );
                rn3d.navigationStep();

                //Modify the speed manually for the experiments
//                ReactInterface.robotSim.simulateInterval( 0.000001f ); //To update internal v,w and retrieve their last value according to the reactive command
//                //const float rand_scale = 0.75f*(1.f + sin(0.015f*iter_count));

//                const float vel_scaling = 3.f;
//                const float rand_scale = (odo_freq/vel_scaling)*0.75f*(1.f + sin(0.015f*iter_count));

//                const float new_v = rand_scale*ReactInterface.robotSim.getV();
//                const float new_w = rand_scale*ReactInterface.robotSim.getW();
                //printf("\n new_v = %f, new_w = %f", new_v, new_w);

                //ReactInterface.changeSpeeds(new_v, new_w);
            }

            //Create random increment of pose for the robot
            else
            {
                CPose2D robot_pose;
                ReactInterface.robotSim.getRealPose(robot_pose);

                //Create new valid pose within the map bounds:
                bool valid_target = false;
                const float max_displ = 7.f, max_rot = 6.f;
                while (!valid_target)
                {
                    const float rand_num_x = float(std::rand()%1000)/1000.f;
                    const float rand_num_y = float(std::rand()%1000)/1000.f;
                    const float rand_num_theta = float(std::rand()%1000)/1000.f;

                    const float alpha = 0.9f;
                    const float new_incr_x = (1.f - alpha)*max_displ*(0.5f - rand_num_x) + alpha*old_incr[0];
                    const float new_incr_y = (1.f - alpha)*max_displ*(0.5f - rand_num_y) + alpha*old_incr[1];
                    const float new_incr_theta = (1.f - alpha)*max_rot*(0.5f - rand_num_theta) + alpha*old_incr[2];

                    const float new_x = new_incr_x + robot_pose[0];
                    const float new_y = new_incr_y + robot_pose[1];
                    const float new_theta = new_incr_theta + robot_pose[2];

                    //Check the new valid target
                    const float range = 0.8f;
                    valid_target = true;
                    for (float u=-0.5f; u<=0.5f; u+= 0.1f)
                        for (float v=-0.5f; v<=0.5f; v+= 0.1f)
                            if (ReactInterface.map.getPos(new_x + u*range, new_y + v*range) < 0.7f)
                                    valid_target = false;

                    if (valid_target == true)
                    {
                        robot_pose = CPose2D(new_x, new_y, new_theta);
                        old_incr = CPose2D(new_incr_x, new_incr_y, new_incr_theta);
                    }

                    old_incr *= 0.99f;
                }

                ReactInterface.robotSim.setRealPose(robot_pose);
                //ReactInterface.new_pose = robot_pose;
                CSimplePointsMap auxpoints;
                ReactInterface.senseObstacles(auxpoints);
            }

            //---------------------------------------------------------------------------------------------
            //Save results and finish when the max_num_iters is reached
//            if (ReactInterface.real_poses.size() >= 5000)
//            {
//                ReactInterface.saveResults(odo_freq);
//                stop = true;
//            }
            //---------------------------------------------------------------------------------------------

            iter_count++;

            if (iter_count % react_sim_per_est == 0)
            {
                //Execute odometry
                ReactInterface.runRF2O();

                //Add the new poses
                ReactInterface.real_poses.push_back(CPose3D(ReactInterface.new_pose));
                ReactInterface.est_poses.push_back(CPose3D(ReactInterface.odo.laser_pose));
                ReactInterface.test_poses.push_back(CPose3D(ReactInterface.odo_test.laser_pose));
            }


            //Create new target if we have reached the previous one or we cannot reach it
            if ((rn3d.IDLE == rn3d.getCurrentState())||(target_clock.Tac() > 5.f)&&(use_reactive))
            {
                target_clock.Tic();
                CSimplePointsMap auxpoints;
                ReactInterface.senseObstacles( auxpoints );

                //Create new random target within the map bounds:
                bool valid_target = false;
                const float maxdist_next_target = 10.f;
                while (!valid_target)
                {
                    const float rand_num_x = float(std::rand()%1000)/1000.f;
                    const float rand_num_y = float(std::rand()%1000)/1000.f;
//                        last_Target_Pos.x = ReactInterface.map.getXMin() + rand_num_x*(ReactInterface.map.getXMax() - ReactInterface.map.getXMin());
//                        last_Target_Pos.y = ReactInterface.map.getYMin() + rand_num_y*(ReactInterface.map.getYMax() - ReactInterface.map.getYMin());
                    last_Target_Pos.x = maxdist_next_target*(0.5f - rand_num_x) + ReactInterface.new_pose[0];
                    last_Target_Pos.y = maxdist_next_target*(0.5f - rand_num_y) + ReactInterface.new_pose[1];

                    //Check whether the new target is valid
                    const float range = 0.8f;
                    valid_target = true;
                    for (float u=-0.5f; u<=0.5f; u+= 0.1f)
                        for (float v=-0.5f; v<=0.5f; v+= 0.1f)
                            if (ReactInterface.map.getPos(last_Target_Pos.x + u*range, last_Target_Pos.y + v*range) < 0.7f)
                                    valid_target = false;
                }

                const CAbstractReactiveNavigationSystem::TNavigationParams  nav_params = ReactInterface.createNewTarget(last_Target_Pos.x, last_Target_Pos.y, 0.3, 0);
                rn3d.navigate(&nav_params);
                ReactInterface.scene->getByClass<CDisk>(0)->setLocation(last_Target_Pos.x, last_Target_Pos.y, last_Target_Pos.z);
            }

            //Update the scene and slow down the simulation
            ReactInterface.updateScene();
            //mrpt::system::sleep(1);	//To slow it down

            if (one_step)
                one_step = false;
        }
    }

    return 0;
}

