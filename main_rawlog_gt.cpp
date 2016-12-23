/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */

#include <iostream>
#include <mrpt/system/threads.h> // sleep()
#include <mrpt/utils/CConfigFile.h>
#include <mrpt/utils/CConfigFileMemory.h>
#include "rawlog_gt.h"



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
            //Clear scene
            odoInterface.clearScene();
            break;

        }


        //						Run reactive + odometry
        //==============================================================================
        if (((working)||(one_step))&&(!odoInterface.dataset_finished))
        {
            //Execute navigation
            odoInterface.readScanRawlog();
            iter_count++;

            odoInterface.updateScene(iter_count);
            mrpt::system::sleep(10);	//To slow it down

            if (one_step)
                one_step = false;
        }
    }

    return 0;
}

