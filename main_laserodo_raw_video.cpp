/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */

#include <iostream>
#include <mrpt/system/threads.h> // sleep()
#include <mrpt/utils/CConfigFile.h>
#include <mrpt/utils/CConfigFileMemory.h>
#include "laserodo_raw_video.h"



// ------------------------------------------------------
//						MAIN
// ------------------------------------------------------

int main()
{

    //Initial steps. Load configuration from file or default
    //------------------------------------------------------
    CLaserodoInterface odoInterface;

    //Initialize all methods
    odoInterface.initializeEverything();

    bool stop = 0;
    bool working = false;
    bool one_step = false;
    int pushed_key = 0;
    int iter_count = 0;
    unsigned int decimate = 1;
    const float dt = 0.025; //Freiburg 0.025, MIT 0.03
    utils::CTicTac clock;
    clock.Tic();


    while (!stop)
    {
        if (odoInterface.window.keyHit())
            pushed_key = odoInterface.window.getPushedKey();
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
            break;

        case 'r':
            //Compute statistics (deviations of the odometry from ground truth)
            odoInterface.computeErrors(10/decimate);
            break;

        case 'c':
            //Reset pose estimation
            odoInterface.resetScene();
            odoInterface.real_poses.clear();
            odoInterface.est_poses.clear(); odoInterface.est_time = 0.f;
            odoInterface.test_poses.clear(); odoInterface.test_time = 0.f;
            odoInterface.psm_poses.clear(); odoInterface.psm_time = 0.f;
            odoInterface.csm_poses.clear(); odoInterface.csm_time = 0.f;
            break;

        case 'x':
            //Save results
            odoInterface.saveResults();
            break;
        }


        //						Run reactive + odometry
        //==============================================================================
        if (((working)||(one_step))&&(!odoInterface.dataset_finished))
        {
            while (clock.Tac() < dt);
            clock.Tic();

            //Execute navigation
            odoInterface.readScanRawlog();
            iter_count++;

            if (iter_count % decimate == 0)
            {
                //Execute odometry
                odoInterface.runRF2O();

                //Execute PSM
                odoInterface.runPolarScanMatching();

                //Execute CSM
                odoInterface.runCanonicalScanMatching();

                //Add the new poses
                odoInterface.real_poses.push_back(CPose3D(odoInterface.new_gt_pose));
                odoInterface.est_poses.push_back(CPose3D(odoInterface.odo.laser_pose));
                odoInterface.test_poses.push_back(CPose3D(odoInterface.odo_test.laser_pose));
                odoInterface.psm_poses.push_back(CPose3D(odoInterface.new_psm_pose));
                odoInterface.csm_poses.push_back(CPose3D(odoInterface.new_csm_pose));
            }

            odoInterface.updateScene(iter_count);


            if (one_step)
                one_step = false;
        }
    }

    return 0;
}

