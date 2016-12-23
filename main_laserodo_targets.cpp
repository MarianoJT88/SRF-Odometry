/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */

#include <iostream>
#include <mrpt/system/threads.h> // sleep()
#include <mrpt/utils/CConfigFile.h>
#include <mrpt/utils/CConfigFileMemory.h>
#include "laserodo_targets.h"


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

    "LASER_POSE = 0 0 0.4 0 0 0 \n"
    "LASER_MAX_RANGE = 5.6 \n"
	"LASER_APERTURE = 240 \n"
    "LASER_STD_ERROR = 0.01 \n" //0.01
	"LASER_LEVEL = 1 \n"
    "LASER_SEGMENTS = 682 \n\n"


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

    "X0 = -7                    ; Initial location (meters) \n"
    "Y0 = -4 \n"
    "PHI0 = 0					; In degrees \n"
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
	"PTG2_WMAX = 55 \n"
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
	"TOO_CLOSE_OBSTACLE = 0.03				; In normalized ps-meters [0,1] \n\n\n"


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

    //Initial steps. Load configuration from file or default
    //------------------------------------------------------
    CMyReactInterface ReactInterface;
    CReactiveNavigationSystem3D rn3d (ReactInterface, false, false);

    utils::CConfigFileMemory configNavigation(default_cfg_txt);
    rn3d.loadConfigFile( configNavigation );
    ReactInterface.loadMaps( configNavigation );
    ReactInterface.loadConfiguration( configNavigation );


    //Initialize all methods
    ReactInterface.initializeEverything();
    rn3d.initialize();

    //Set the targets
    ReactInterface.x_target[0] = -2.4f; ReactInterface.y_target[0] = -2.1f;
    ReactInterface.x_target[1] = -5.5f; ReactInterface.y_target[1] = 0.7f;
    ReactInterface.x_target[2] = 0.2f; ReactInterface.y_target[2] = 1.f;
    ReactInterface.x_target[3] = 2.4f; ReactInterface.y_target[3] = -1.4f;
    ReactInterface.x_target[4] = 5.9f; ReactInterface.y_target[4] = -3.6f;
    ReactInterface.x_target[5] = 6.4f; ReactInterface.y_target[5] = 2.1f;
    ReactInterface.x_target[6] = -7.f; ReactInterface.y_target[6] = 3.8f;


    //                          Main loop
    //-----------------------------------------------------------------

    bool stop = 0;
    bool working = false;
    bool one_step = false;
    bool moving_target = 0;
    int pushed_key = 0;
    int iter_count = 0;
    unsigned int odo_freq = 1;
    unsigned int react_sim_per_est = 10;
    float sim_period = 1.f/float(odo_freq*react_sim_per_est);	//Hokuyo max frequency is 10Hz
    unsigned int num_target = 0;


    while (!stop)
    {
        if (ReactInterface.window.keyHit())
            pushed_key = ReactInterface.window.getPushedKey();
        else
            pushed_key = 0;

        switch (pushed_key) {

        case 'm':
            //Move the target
            moving_target = 1;
            break;

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

        case 'f':
            //Save results in file
            ReactInterface.saveResults(odo_freq);
            break;

        case 'r':
            //Compute statistics (deviations of the odometry from ground truth)
            ReactInterface.computeErrors(odo_freq);
            break;

        case 'l':
            //Compute statistics (deviations of the odometry from ground truth) for different segment lengths
            ReactInterface.computeErrorsPerLength();
            break;

        case 'c':
            //Reset pose estimation
            ReactInterface.resetScene();
            ReactInterface.real_poses.clear();
            ReactInterface.est_poses.clear(); ReactInterface.est_time = 0.f;
            ReactInterface.test_poses.clear(); ReactInterface.test_time = 0.f;
            ReactInterface.psm_poses.clear(); ReactInterface.psm_time = 0.f;
            ReactInterface.csm_poses.clear(); ReactInterface.csm_time = 0.f;
            break;

        case 'x':
            //Save scans
            ReactInterface.saveScans();
            break;

        }

        //						Run reactive + odometry
        //==============================================================================
        if ((working)||(one_step))
        {
            //Execute navigation
            //mrpt::system::sleep(30);	//I don't know why but it is necessary...
            rn3d.navigationStep();
            ReactInterface.robotSim.simulateInterval( sim_period );

            iter_count++;

            if (iter_count % react_sim_per_est == 0)
            {
                //Execute odometry
                ReactInterface.runRF2O();

                //Execute PSM
                ReactInterface.runPolarScanMatching();

                //Execute CSM
                ReactInterface.runCanonicalScanMatching();

                //Add the new poses
                ReactInterface.real_poses.push_back(CPose3D(ReactInterface.new_pose));
                ReactInterface.est_poses.push_back(CPose3D(ReactInterface.odo.laser_pose));
                ReactInterface.test_poses.push_back(CPose3D(ReactInterface.odo_test.laser_pose));
                ReactInterface.psm_poses.push_back((CPose3D(ReactInterface.new_psm_pose)));
                ReactInterface.csm_poses.push_back((CPose3D(ReactInterface.new_csm_pose)));

            }


            //Sense new obstacles and update scene
            if ((rn3d.IDLE == rn3d.getCurrentState())||(rn3d.SUSPENDED == rn3d.getCurrentState()))
            {
                CSimplePointsMap auxpoints;
                ReactInterface.senseObstacles( auxpoints );

                if (num_target < 7)
                {
                    const CAbstractReactiveNavigationSystem::TNavigationParams  nav_params = ReactInterface.createNewTarget(ReactInterface.x_target[num_target], ReactInterface.y_target[num_target], 0.3, 0);
                    rn3d.navigate(&nav_params);

                    // Move an object to the position picked by the user:
                    ReactInterface.scene->getByClass<CDisk>(0)->setLocation(ReactInterface.x_target[num_target], ReactInterface.y_target[num_target], 0.f);
                    num_target++;
                }
                else
                    working = false;

            }
            ReactInterface.updateScene();

            if (one_step)
                one_step = false;
        }
    }

    return 0;

}

