/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */

#include "laser_odometry_3scans.h"


using namespace mrpt::utils;
using namespace Eigen;
using namespace std;


void RF2O_3S::initialize(unsigned int size, float FOV_rad, bool is_test)
{
    test = is_test;
    cols = size;
    width = size;
    fovh = FOV_rad;
    ctf_levels = ceilf(log2(cols) - 4.3f);
    iter_irls = 5;
    fps = 1.f;	//In Hz
	
    //Resize original range scan
    range_wf.resize(width);

    //Resize the transformation matrices
    transformations.resize(ctf_levels);
    for (unsigned int i = 0; i < ctf_levels; i++)
    {
        transformations[i].resize(3,3);
        transformations[i].setIdentity();
    }

	//Resize pyramid
	unsigned int s, cols_i;
    const unsigned int pyr_levels = round(log2(round(float(width)/float(cols)))) + ctf_levels;
    range_1.resize(pyr_levels); range_2.resize(pyr_levels); range_3.resize(pyr_levels);
    range_12.resize(pyr_levels); range_13.resize(pyr_levels);
    xx_1.resize(pyr_levels); xx_2.resize(pyr_levels); xx_3.resize(pyr_levels);
    xx_12.resize(pyr_levels); xx_13.resize(pyr_levels);
    yy_1.resize(pyr_levels); yy_2.resize(pyr_levels); yy_3.resize(pyr_levels);
    yy_12.resize(pyr_levels); yy_13.resize(pyr_levels);
    range_warped.resize(pyr_levels); xx_warped.resize(pyr_levels); yy_warped.resize(pyr_levels);
    range_3_warpedTo2.resize(pyr_levels); xx_3_warpedTo2.resize(pyr_levels); yy_3_warpedTo2.resize(pyr_levels);

	for (unsigned int i = 0; i<pyr_levels; i++)
    {
        s = pow(2.f,int(i));
        cols_i = ceil(float(width)/float(s));

        range_1[i].resize(cols_i); range_2[i].resize(cols_i); range_3[i].resize(cols_i);
        range_12[i].resize(cols_i); range_13[i].resize(cols_i);
        range_1[i].fill(0.f); range_2[i].fill(0.f); range_3[i].fill(0.f);
        xx_1[i].resize(cols_i); xx_2[i].resize(cols_i); xx_3[i].resize(cols_i);
        xx_12[i].resize(cols_i); xx_13[i].resize(cols_i);
        xx_1[i].fill(0.f); xx_2[i].fill(0.f); xx_3[i].fill(0.f);
        yy_1[i].resize(cols_i); yy_2[i].resize(cols_i); yy_3[i].resize(cols_i);
        yy_12[i].resize(cols_i); yy_13[i].resize(cols_i);
        yy_1[i].fill(0.f); yy_2[i].fill(0.f); yy_3[i].fill(0.f);

		if (cols_i <= cols)
		{
            range_warped[i].resize(cols_i); xx_warped[i].resize(cols_i); yy_warped[i].resize(cols_i);
            range_3_warpedTo2[i].resize(cols_i); xx_3_warpedTo2[i].resize(cols_i); yy_3_warpedTo2[i].resize(cols_i);
		}
    }

    //Resize aux variables
    dt_12.resize(cols); dt_13.resize(cols);
    dtita_12.resize(cols); dtita_13.resize(cols);
    weights_12.resize(cols); weights_13.resize(cols);
    null_12.resize(cols); null_13.resize(cols);
    null_12.fill(false); null_13.fill(false);
	cov_odo.assign(0.f);
    outliers.resize(cols);
    outliers.fill(false);


	//Compute gaussian mask
	g_mask[0] = 1.f/16.f; g_mask[1] = 0.25f; g_mask[2] = 6.f/16.f; g_mask[3] = g_mask[1]; g_mask[4] = g_mask[0];

    //Initialize "last velocity" as zero
	kai_abs.assign(0.f);
	kai_loc_old.assign(0.f);
    overall_trans_prev.setIdentity();
}


void RF2O_3S::createScanPyramid()
{
	const float max_range_dif = 0.3f;
	
    //Push scans back
    range_2.swap(range_3); xx_2.swap(xx_3); yy_2.swap(yy_3);
    range_1.swap(range_2); xx_1.swap(xx_2); yy_1.swap(yy_2);


    //The number of levels of the pyramid does not match the number of levels used
    //in the odometry computation (because we sometimes want to finish with lower resolutions)

    unsigned int pyr_levels = round(log2(round(float(width)/float(cols)))) + ctf_levels;

    //Generate levels
    for (unsigned int i = 0; i<pyr_levels; i++)
    {
        unsigned int s = pow(2.f,int(i));
        cols_i = ceil(float(width)/float(s));
		const unsigned int i_1 = i-1;

        //              First level -> Filter, not downsample
        //------------------------------------------------------------------------
        if (i == 0)
		{
			for (unsigned int u = 0; u < cols_i; u++)
            {	
				const float dcenter = range_wf(u);
					
				//Inner pixels
                if ((u>1)&&(u<cols_i-2))
                {		
					if (dcenter > 0.f)
					{	
                        float sum = 0.f, weight = 0.f;

						for (int l=-2; l<3; l++)
						{
							const float abs_dif = abs(range_wf(u+l)-dcenter);
							if (abs_dif < max_range_dif)
							{
								const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
								weight += aux_w;
								sum += aux_w*range_wf(u+l);
							}
						}
                        range_1[i](u) = sum/weight;
					}
					else
                        range_1[i](u) = 0.f;

                }

                //Boundary
                else
                {
                    if (dcenter > 0.f)
					{						
                        float sum = 0.f, weight = 0.f;

						for (int l=-2; l<3; l++)	
						{
							const int indu = u+l;
							if ((indu>=0)&&(indu<cols_i))
							{
								const float abs_dif = abs(range_wf(indu)-dcenter);										
								if (abs_dif < max_range_dif)
								{
									const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
									weight += aux_w;
									sum += aux_w*range_wf(indu);
								}
							}
						}
                        range_1[i](u) = sum/weight;
					}
					else
                        range_1[i](u) = 0.f;

                }
            }
		}

        //                              Downsampling
        //-----------------------------------------------------------------------------
        else
        {            
			for (unsigned int u = 0; u < cols_i; u++)
            {
                const int u2 = 2*u;		
                const float dcenter = range_1[i_1](u2);
					
				//Inner pixels
                if ((u>0)&&(u<cols_i-1))
                {		
					if (dcenter > 0.f)
					{	
                        float sum = 0.f, weight = 0.f;

						for (int l=-2; l<3; l++)
						{
                            const float abs_dif = abs(range_1[i_1](u2+l)-dcenter);
							if (abs_dif < max_range_dif)
							{
								const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
								weight += aux_w;
                                sum += aux_w*range_1[i_1](u2+l);
							}
						}
                        range_1[i](u) = sum/weight;
					}
					else
                        range_1[i](u) = 0.f;

                }

                //Boundary
                else
                {
                    if (dcenter > 0.f)
					{						
                        float sum = 0.f, weight = 0.f;
                        const unsigned int cols_i2 = range_1[i_1].rows();

						for (int l=-2; l<3; l++)	
						{
							const int indu = u2+l;

							if ((indu>=0)&&(indu<cols_i2))
							{
                                const float abs_dif = abs(range_1[i_1](indu)-dcenter);
								if (abs_dif < max_range_dif)
								{
									const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
									weight += aux_w;
                                    sum += aux_w*range_1[i_1](indu);
								}
							}
						}
                        range_1[i](u) = sum/weight;

					}
					else
                        range_1[i](u) = 0.f;
                }
            }
        }

        //Calculate coordinates "xy" of the points
        for (unsigned int u = 0; u < cols_i; u++) 
		{
            if (range_1[i](u) > 0.f)
			{
                const float tita = -0.5f*fovh + float(u)*fovh/float(cols_i-1);
                xx_1[i](u) = range_1[i](u)*cos(tita);
                yy_1[i](u) = range_1[i](u)*sin(tita);
			}
			else
			{
                xx_1[i](u) = 0.f;
                yy_1[i](u) = 0.f;
			}
        }
    }
}

void RF2O_3S::calculateCoord()
{		
    null_12.fill(false);
    null_13.fill(false);
    num_valid_range = 0;

    for (unsigned int u = 0; u < cols_i; u++)
	{
        //Coordinates 12
        if ((range_2[image_level](u) == 0.f) || (range_warped[image_level](u) == 0.f))
		{
            range_12[image_level](u) = 0.f;
            xx_12[image_level](u) = 0.f;
            yy_12[image_level](u) = 0.f;
            null_12(u) = true;
		}
		else
		{
            range_12[image_level](u) = 0.5f*(range_2[image_level](u) + range_warped[image_level](u));
            xx_12[image_level](u) = 0.5f*(xx_2[image_level](u) + xx_warped[image_level](u));
            yy_12[image_level](u) = 0.5f*(yy_2[image_level](u) + yy_warped[image_level](u));
            if ((u>0)&&(u<cols_i-1))
                num_valid_range++;
		}

        //Coordinates 13
        if ((range_3_warpedTo2[image_level](u) == 0.f) || (range_warped[image_level](u) == 0.f))
        {
            range_13[image_level](u) = 0.f;
            xx_13[image_level](u) = 0.f;
            yy_13[image_level](u) = 0.f;
            null_13(u) = true;
        }
        else
        {
            range_13[image_level](u) = 0.5f*(range_3_warpedTo2[image_level](u) + range_warped[image_level](u));
            xx_13[image_level](u) = 0.5f*(xx_3_warpedTo2[image_level](u) + xx_warped[image_level](u));
            yy_13[image_level](u) = 0.5f*(yy_3_warpedTo2[image_level](u) + yy_warped[image_level](u));
            if ((u>0)&&(u<cols_i-1))
                num_valid_range++;
        }
	}
}

void RF2O_3S::calculateRangeDerivatives()
{	
    //Compute distances between points
    Eigen::ArrayXf rtita_12(cols_i), rtita_13(cols_i);
    rtita_12.fill(1.f); rtita_13.fill(1.f);

	for (unsigned int u = 0; u < cols_i-1; u++)
    {
        const float dist_12 = square(xx_12[image_level](u+1) - xx_12[image_level](u))
                            + square(yy_12[image_level](u+1) - yy_12[image_level](u));

        const float dist_13 = square(xx_13[image_level](u+1) - xx_13[image_level](u))
                            + square(yy_13[image_level](u+1) - yy_13[image_level](u));

        if (dist_12  > 0.f)
            rtita_12(u) = sqrtf(dist_12);

        if (dist_13  > 0.f)
            rtita_13(u) = sqrtf(dist_13);
	}

    //Spatial derivatives
    for (unsigned int u = 1; u < cols_i-1; u++)
    {
        dtita_12(u) = (rtita_12(u-1)*(range_12[image_level](u+1)-range_12[image_level](u)) + rtita_12(u)*(range_12[image_level](u) - range_12[image_level](u-1)))/(rtita_12(u)+rtita_12(u-1));
        dtita_13(u) = (rtita_13(u-1)*(range_13[image_level](u+1)-range_13[image_level](u)) + rtita_13(u)*(range_13[image_level](u) - range_13[image_level](u-1)))/(rtita_13(u)+rtita_13(u-1));
    }

    dtita_12(0) = dtita_12(1);
    dtita_12(cols_i-1) = dtita_12(cols_i-2);

    dtita_13(0) = dtita_13(1);
    dtita_13(cols_i-1) = dtita_13(cols_i-2);

	//Temporal derivative
	for (unsigned int u = 0; u < cols_i; u++)
    {
        dt_12(u) = fps*(range_warped[image_level](u) - range_2[image_level](u));
        dt_13(u) = fps*(range_warped[image_level](u) - range_3_warpedTo2[image_level](u));
    }
}

void RF2O_3S::computeWeights()
{
	//The maximum weight size is reserved at the constructor
    weights_12.fill(0.f);
    weights_13.fill(0.f);
	
	//Parameters for error_linearization
	const float kdtita = 1.f;
	const float kdt = kdtita/square(fps);
	const float k2d = 0.2f;
    const float sensor_sigma = 4e-4f;
	
	for (unsigned int u = 1; u < cols_i-1; u++)
    {
        if (null_12(u) == false)
		{	
			//							Compute derivatives
			//-----------------------------------------------------------------------
            const float ini_dtita = range_2[image_level](u+1) - range_2[image_level](u-1);
            const float final_dtita = range_warped[image_level](u+1) - range_warped[image_level](u-1);

			const float dtitat = ini_dtita - final_dtita;
            const float dtita2 = dtita_12(u+1) - dtita_12(u-1);

            const float w_der = kdt*square(dt_12(u)) + kdtita*square(dtita_12(u)) + k2d*(abs(dtitat) + abs(dtita2)) + sensor_sigma;

            weights_12(u) = sqrtf(1.f/w_der);
		}

        if (null_13(u) == false)
        {
            //							Compute derivatives
            //-----------------------------------------------------------------------
            const float ini_dtita = range_3_warpedTo2[image_level](u+1) - range_3_warpedTo2[image_level](u-1);
            const float final_dtita = range_warped[image_level](u+1) - range_warped[image_level](u-1);

            const float dtitat = ini_dtita - final_dtita;
            const float dtita2 = dtita_13(u+1) - dtita_13(u-1);

            const float w_der = kdt*square(dt_13(u)) + kdtita*square(dtita_13(u)) + k2d*(abs(dtitat) + abs(dtita2)) + sensor_sigma;

            weights_13(u) = sqrtf(1.f/w_der);
        }
    }

    const float max_w = max(weights_12.maxCoeff(), weights_13.maxCoeff());
    const float inv_max_w = 1.f/max_w;
    weights_12 = inv_max_w*weights_12;
    weights_13 = inv_max_w*weights_13;
}


void RF2O_3S::solveSystemQuadResiduals3Scans()
{
    A.resize(num_valid_range,3);
    B.resize(num_valid_range);
    unsigned int cont = 0;
    const float kdtita = (cols_i-1)/fovh;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
    {
        if (null_12(u) == false)
        {
            // Precomputed expressions
            const float tw = weights_12(u);
            const float tita = -0.5f*fovh + u/kdtita;

            //Fill the matrix A
            A(cont, 0) = tw*(cos(tita) + dtita_12(u)*kdtita*sin(tita)/range_12[image_level](u));
            A(cont, 1) = tw*(sin(tita) - dtita_12(u)*kdtita*cos(tita)/range_12[image_level](u));
            A(cont, 2) = tw*(-yy_12[image_level](u)*cos(tita) + xx_12[image_level](u)*sin(tita) - dtita_12(u)*kdtita); //?????
            B(cont) = tw*(-dt_12(u));

            cont++;
        }

        if (null_13(u) == false)
        {
            // Precomputed expressions
            const float tw = weights_13(u);
            const float tita = -0.5f*fovh + u/kdtita;

            //Fill the matrix A
            A(cont, 0) = tw*(cos(tita) + dtita_13(u)*kdtita*sin(tita)/range_13[image_level](u));
            A(cont, 1) = tw*(sin(tita) - dtita_13(u)*kdtita*cos(tita)/range_13[image_level](u));
            A(cont, 2) = tw*(-yy_13[image_level](u)*cos(tita) + xx_13[image_level](u)*sin(tita) - dtita_13(u)*kdtita);
            B(cont) = tw*(-dt_13(u));

            cont++;
        }
    }

    //Solve the linear system of equations using a minimum least squares method
    MatrixXf AtA, AtB;
    AtA.multiply_AtA(A);
    AtB.multiply_AtB(A,B);
    kai_loc_level = AtA.ldlt().solve(AtB);

    //Covariance matrix calculation
    VectorXf res = A*kai_loc_level - B;
    cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();
}


//void RF2O_3S::solveSystemMCauchy()
//{
//	A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
//    B.resize(num_valid_range); Bw.resize(num_valid_range);
//	unsigned int cont = 0;
//	const float kdtita = float(cols_i-1)/fovh;

//	//Fill the matrix A and the vector B
//	//The order of the variables will be (vx, vy, wz)

//	for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//		{
//			// Precomputed expressions
//			const float tw = weights(u);
//			const float tita = -0.5*fovh + u/kdtita;

//			//Fill the matrix A
//			A(cont, 0) = tw*(cos(tita) + dtita(u)*kdtita*sin(tita)/range_inter[image_level](u));
//			A(cont, 1) = tw*(sin(tita) - dtita(u)*kdtita*cos(tita)/range_inter[image_level](u));
//			A(cont, 2) = tw*(-yy[image_level](u)*cos(tita) + xx[image_level](u)*sin(tita) - dtita(u)*kdtita);
//            B(cont) = tw*(-dt(u));

//			cont++;
//		}
	
//	//Solve the linear system of equations using a minimum least squares method
//	MatrixXf AtA, AtB;
//	AtA.multiply_AtA(A);
//	AtB.multiply_AtB(A,B);
//    kai_loc_level = AtA.ldlt().solve(AtB);
//    VectorXf res = A*kai_loc_level - B;
//	//cout << endl << "max res: " << res.maxCoeff();
//	//cout << endl << "min res: " << res.minCoeff();

//    //Compute the average dt and res
//    float aver_dt = 0.f, aver_res = 0.f; unsigned int ind = 0;
//    for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//        {
//            aver_dt += fabsf(dt(u));
//            aver_res += fabsf(res(ind++));
//        }
//    aver_dt /= cont; aver_res /= cont;
//    const float k = 10.f/aver_dt; //200

//    ////Compute the energy
//	//float energy = 0.f;
//	//for (unsigned int i=0; i<res.rows(); i++)
//	//	energy += log(1.f + square(k*res(i)));
//	//printf("\n\nEnergy(0) = %f", energy);

//    //Solve iterative reweighted least squares
//    //===================================================================
//    for (unsigned int i=1; i<=iter_irls; i++)
//    {
//        cont = 0;

//		for (unsigned int u = 1; u < cols_i-1; u++)
//            if (null(u) == false)
//			{
//                const float res_weight = sqrtf(1.f/(1.f + square(k*res(cont))));

//                //Fill the matrix Aw
//                Aw(cont,0) = res_weight*A(cont,0);
//                Aw(cont,1) = res_weight*A(cont,1);
//                Aw(cont,2) = res_weight*A(cont,2);
//                Bw(cont) = res_weight*B(cont);
//                cont++;
//            }

//        //Solve the linear system of equations using a minimum least squares method
//        AtA.multiply_AtA(Aw);
//        AtB.multiply_AtB(Aw,Bw);
//        kai_loc_level = AtA.ldlt().solve(AtB);
//        res = A*kai_loc_level - B;

//		////Compute the energy
//		//energy = 0.f;
//		//for (unsigned int j=0; j<res.rows(); j++)
//		//	energy += log(1.f + square(k*res(j)));
//		//printf("\nEnergy(%d) = %f", i, energy);
//    }

//    //Covariance calculation
//	cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();
//}

//void RF2O_3S::solveSystemMTukey()
//{
//    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
//    B.resize(num_valid_range); Bw.resize(num_valid_range);
//    unsigned int cont = 0;
//    const float kdtita = float(cols_i-1)/fovh;

//    //Fill the matrix A and the vector B
//    //The order of the variables will be (vx, vy, wz)

//    for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//        {
//            // Precomputed expressions
//            const float tw = weights(u);
//            const float tita = -0.5*fovh + u/kdtita;

//            //Fill the matrix A
//            A(cont, 0) = tw*(cos(tita) + dtita(u)*kdtita*sin(tita)/range_inter[image_level](u));
//            A(cont, 1) = tw*(sin(tita) - dtita(u)*kdtita*cos(tita)/range_inter[image_level](u));
//            A(cont, 2) = tw*(-yy[image_level](u)*cos(tita) + xx[image_level](u)*sin(tita) - dtita(u)*kdtita);
//            B(cont) = tw*(-dt(u));

//            cont++;
//        }

//    //Solve the linear system of equations using a minimum least squares method
//    MatrixXf AtA, AtB;
//    AtA.multiply_AtA(A);
//    AtB.multiply_AtB(A,B);
//    kai_loc_level = AtA.ldlt().solve(AtB);
//    VectorXf res = A*kai_loc_level - B;
//    //cout << endl << "max res: " << res.maxCoeff();
//    //cout << endl << "min res: " << res.minCoeff();

//    //Compute the median of res
//    vector<float> aux_vector;
//    for (unsigned int k = 0; k<res.rows(); k++)
//        aux_vector.push_back(res(k));
//    std::sort(aux_vector.begin(), aux_vector.end());
//    float res_median = aux_vector.at(res.rows()/2);

//    //Compute the median absolute deviation
//    aux_vector.clear();
//    for (unsigned int k = 0; k<res.rows(); k++)
//        aux_vector.push_back(abs(res(k) - res_median));
//    std::sort(aux_vector.begin(), aux_vector.end());
//    float mad = aux_vector.at(res.rows()/2);

//    //Find the m-estimator constant
//    float c = 5.f*mad;

//    ////Compute the energy
//    //float energy = 0.f;
//    //for (unsigned int i=0; i<res.rows(); i++)
//    //	energy += log(1.f + square(k*res(i)));
//    //printf("\n\nEnergy(0) = %f", energy);

//    //Solve iterative reweighted least squares
//    //===================================================================
//    for (unsigned int i=1; i<=iter_irls; i++)
//    {
//        cont = 0;

//        for (unsigned int u = 1; u < cols_i-1; u++)
//            if (null(u) == false)
//            {
//                float res_weight;
//                if (abs(res(cont)) <= c)    res_weight = square(1.f - square(res(cont)/c));
//                else                        res_weight = 0.f;

//                //Fill the matrix Aw
//                Aw(cont,0) = res_weight*A(cont,0);
//                Aw(cont,1) = res_weight*A(cont,1);
//                Aw(cont,2) = res_weight*A(cont,2);
//                Bw(cont) = res_weight*B(cont);
//                cont++;
//            }

//        //Solve the linear system of equations using a minimum least squares method
//        AtA.multiply_AtA(Aw);
//        AtB.multiply_AtB(Aw,Bw);
//        kai_loc_level = AtA.ldlt().solve(AtB);
//        res = A*kai_loc_level - B;

//        ////Compute the energy
//        //energy = 0.f;
//        //for (unsigned int j=0; j<res.rows(); j++)
//        //	energy += log(1.f + square(k*res(j)));
//        //printf("\nEnergy(%d) = %f", i, energy);

//        //Recompute c
//        //-------------------------------------------------
//        //Compute the median of res
//        aux_vector.clear();
//        for (unsigned int k = 0; k<res.rows(); k++)
//            aux_vector.push_back(res(k));
//        std::sort(aux_vector.begin(), aux_vector.end());
//        res_median = aux_vector.at(res.rows()/2);

//        //Compute the median absolute deviation
//        aux_vector.clear();
//        for (unsigned int k = 0; k<res.rows(); k++)
//            aux_vector.push_back(abs(res(k) - res_median));
//        std::sort(aux_vector.begin(), aux_vector.end());
//        mad = aux_vector.at(res.rows()/2);

//        //Find the m-estimator constant
//        c = 5.f*mad;

//    }

//    //Covariance calculation
//    cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();

//    //Update the outlier mask
//    cont = 0; outliers.fill(false); unsigned int num_outliers = 0;
//    for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//            if (abs(res(cont++)) > c)
//            {
//                outliers(u) = true;
//                num_outliers++;
//            }

//    printf("\n Num_outliers = %d", num_outliers);
//}

//void RF2O_3S::solveSystemTruncatedQuad()
//{
//    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
//    B.resize(num_valid_range); Bw.resize(num_valid_range);
//    unsigned int cont = 0;
//    const float kdtita = float(cols_i-1)/fovh;

//    //Fill the matrix A and the vector B
//    //The order of the variables will be (vx, vy, wz)

//    for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//        {
//            // Precomputed expressions
//            const float tw = weights(u);
//            const float tita = -0.5*fovh + u/kdtita;

//            //Fill the matrix A
//            A(cont, 0) = tw*(cos(tita) + dtita(u)*kdtita*sin(tita)/range_inter[image_level](u));
//            A(cont, 1) = tw*(sin(tita) - dtita(u)*kdtita*cos(tita)/range_inter[image_level](u));
//            A(cont, 2) = tw*(-yy[image_level](u)*cos(tita) + xx[image_level](u)*sin(tita) - dtita(u)*kdtita);
//            B(cont) = tw*(-dt(u));

//            cont++;
//        }

//    //Solve the linear system of equations using a minimum least squares method
//    MatrixXf AtA, AtB;
//    AtA.multiply_AtA(A);
//    AtB.multiply_AtB(A,B);
//    kai_loc_level = AtA.ldlt().solve(AtB);
//    VectorXf res = A*kai_loc_level - B;
//    //cout << endl << "max res: " << res.maxCoeff();
//    //cout << endl << "min res: " << res.minCoeff();

//    //Compute the median of res
//    vector<float> aux_vector;
//    for (unsigned int k = 0; k<res.rows(); k++)
//        aux_vector.push_back(res(k));
//    std::sort(aux_vector.begin(), aux_vector.end());
//    float res_median = aux_vector.at(res.rows()/2);

//    //Compute the median absolute deviation
//    aux_vector.clear();
//    for (unsigned int k = 0; k<res.rows(); k++)
//        aux_vector.push_back(abs(res(k) - res_median));
//    std::sort(aux_vector.begin(), aux_vector.end());
//    float mad = aux_vector.at(res.rows()/2);

//    //Find the m-estimator constant
//    float c = 5.f*mad;

//    ////Compute the energy
//    //float energy = 0.f;
//    //for (unsigned int i=0; i<res.rows(); i++)
//    //	energy += log(1.f + square(k*res(i)));
//    //printf("\n\nEnergy(0) = %f", energy);

//    //Solve iterative reweighted least squares
//    //===================================================================
//    for (unsigned int i=1; i<=iter_irls; i++)
//    {
//        cont = 0;

//        for (unsigned int u = 1; u < cols_i-1; u++)
//            if (null(u) == false)
//            {
//                float res_weight;
//                if (abs(res(cont)) <= c)    res_weight = 1.f;
//                else                        res_weight = 0.f;

//                //Fill the matrix Aw
//                Aw(cont,0) = res_weight*A(cont,0);
//                Aw(cont,1) = res_weight*A(cont,1);
//                Aw(cont,2) = res_weight*A(cont,2);
//                Bw(cont) = res_weight*B(cont);
//                cont++;
//            }

//        //Solve the linear system of equations using a minimum least squares method
//        AtA.multiply_AtA(Aw);
//        AtB.multiply_AtB(Aw,Bw);
//        kai_loc_level = AtA.ldlt().solve(AtB);
//        res = A*kai_loc_level - B;

//        ////Compute the energy
//        //energy = 0.f;
//        //for (unsigned int j=0; j<res.rows(); j++)
//        //	energy += log(1.f + square(k*res(j)));
//        //printf("\nEnergy(%d) = %f", i, energy);

//        //Recompute c
//        //-------------------------------------------------
//        //Compute the median of res
//        aux_vector.clear();
//        for (unsigned int k = 0; k<res.rows(); k++)
//            aux_vector.push_back(res(k));
//        std::sort(aux_vector.begin(), aux_vector.end());
//        res_median = aux_vector.at(res.rows()/2);

//        //Compute the median absolute deviation
//        aux_vector.clear();
//        for (unsigned int k = 0; k<res.rows(); k++)
//            aux_vector.push_back(abs(res(k) - res_median));
//        std::sort(aux_vector.begin(), aux_vector.end());
//        mad = aux_vector.at(res.rows()/2);

//        //Find the m-estimator constant
//        c = 5.f*mad;

//    }

//    //Covariance calculation
//    cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();

//    //Update the outlier mask
//    cont = 0; outliers.fill(false); unsigned int num_outliers = 0;
//    for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//            if (abs(res(cont++)) > c)
//            {
//                outliers(u) = true;
//                num_outliers++;
//            }

//    printf("\n Num_outliers = %d", num_outliers);
//}

void RF2O_3S::solveSystemSmoothTruncQuad3Scans()
{
    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
    B.resize(num_valid_range); Bw.resize(num_valid_range);
    unsigned int cont = 0;
    const float kdtita = float(cols_i-1)/fovh;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
    {
        const float tita = -0.5f*fovh + u/kdtita;
        const float cos_tita = cos(tita);
        const float sin_tita = sin(tita);

        if (null_12(u) == false)
        {
            // Precomputed expressions
            const float tw = weights_12(u);


            //Fill the matrix A
            A(cont, 0) = tw*(cos_tita + dtita_12(u)*kdtita*sin_tita/range_12[image_level](u));
            A(cont, 1) = tw*(sin_tita - dtita_12(u)*kdtita*cos_tita/range_12[image_level](u));
            A(cont, 2) = tw*(-yy_12[image_level](u)*cos_tita + xx_12[image_level](u)*sin_tita - dtita_12(u)*kdtita);
            B(cont) = tw*(-dt_12(u));

            cont++;
        }

        if (null_13(u) == false)
        {
            // Precomputed expressions
            const float tw = weights_13(u);

            //Fill the matrix A
            A(cont, 0) = tw*(cos_tita + dtita_13(u)*kdtita*sin_tita/range_13[image_level](u));
            A(cont, 1) = tw*(sin_tita - dtita_13(u)*kdtita*cos_tita/range_13[image_level](u));
            A(cont, 2) = tw*(-yy_13[image_level](u)*cos_tita + xx_13[image_level](u)*sin_tita - dtita_13(u)*kdtita);
            B(cont) = tw*(-dt_13(u));

            cont++;
        }
    }

    //Solve the linear system of equations using a minimum least squares method
    MatrixXf AtA, AtB;
    AtA.multiply_AtA(A);
    AtB.multiply_AtB(A,B);
    kai_loc_level = AtA.ldlt().solve(AtB);
    VectorXf res = A*kai_loc_level - B;
    //cout << endl << "max res: " << res.maxCoeff();
    //cout << endl << "min res: " << res.minCoeff();

    //Compute the median of res
    vector<float> aux_vector;
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(res(k));
    std::sort(aux_vector.begin(), aux_vector.end());
    const float res_median = aux_vector.at(res.rows()/2);

    //Compute the median absolute deviation
    aux_vector.clear();
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(abs(res(k) - res_median));
    std::sort(aux_vector.begin(), aux_vector.end());
    const float mad = aux_vector.at(res.rows()/2);

    //Find the m-estimator constant
    const float c = 4.f*mad;
    const float c_inv = 1.f/c;

    //Compute the energy
    float new_energy = 0.f, last_energy;
    for (unsigned int i=0; i<res.rows(); i++)
    {
        if (abs(res(i)) < c)     new_energy += 0.5f*square(res(i))*(1.f - 0.5f*square(res(i)/c));
        else                     new_energy += 0.25f*square(c);
    }
    //printf("\n\nEnergy(0) = %f", new_energy);
    last_energy = 2.f*new_energy;
    unsigned int iter = 1;

    //Solve iterative reweighted least squares
    //===================================================================
    while ((new_energy < 0.995f*last_energy)&&(iter < 10))
    {
        cont = 0;
        last_energy = new_energy;

        for (unsigned int u = 1; u < cols_i-1; u++)
        {
            if (null_12(u) == false)
            {
                float res_weight;
                if (abs(res(cont)) <= c)    res_weight = 1.f - square(res(cont)*c_inv);
                else                        res_weight = 0.f;

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
            }

            if (null_13(u) == false)
            {
                float res_weight;
                if (abs(res(cont)) <= c)    res_weight = 1.f - square(res(cont)*c_inv);
                else                        res_weight = 0.f;

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
            }
        }

        //Solve the linear system of equations using a minimum least squares method
        AtA.multiply_AtA(Aw);
        AtB.multiply_AtB(Aw,Bw);
        kai_loc_level = AtA.ldlt().solve(AtB);
        res = A*kai_loc_level - B;

        //Compute the energy
        new_energy = 0.f;
        for (unsigned int i=0; i<res.rows(); i++)
        {
            if (abs(res(i)) < c)    new_energy += 0.5f*square(res(i))*(1.f - 0.5f*square(res(i)/c));
            else                    new_energy += 0.25f*square(c);
        }
        //printf("\nEnergy(%d) = %f", iter, new_energy);
        iter++;

//        //Recompute c
//        //-------------------------------------------------
//        //Compute the median of res
//        aux_vector.clear();
//        for (unsigned int k = 0; k<res.rows(); k++)
//            aux_vector.push_back(res(k));
//        std::sort(aux_vector.begin(), aux_vector.end());
//        res_median = aux_vector.at(res.rows()/2);

//        //Compute the median absolute deviation
//        aux_vector.clear();
//        for (unsigned int k = 0; k<res.rows(); k++)
//            aux_vector.push_back(abs(res(k) - res_median));
//        std::sort(aux_vector.begin(), aux_vector.end());
//        mad = aux_vector.at(res.rows()/2);

//        //Find the m-estimator constant
//        c = 5.f*mad;
    }

    //Covariance calculation
    cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();

    //Update the outlier mask
//    cont = 0; outliers.fill(false); unsigned int num_outliers = 0;
//    for (unsigned int u = 1; u < cols_i-1; u++)
//        if (null(u) == false)
//            if (abs(res(cont++)) > c)
//            {
//                outliers(u) = true;
//                num_outliers++;
//            }

//    printf("\n Num_outliers = %d", num_outliers);
}


void RF2O_3S::performWarping()
{
    Matrix3f acu_trans;
    acu_trans.setIdentity();
    for (unsigned int i=0; i<=level; i++)
        acu_trans = transformations[i]*acu_trans;

    ArrayXf wacu(cols_i);
    wacu.fill(0.f);
    range_warped[image_level].fill(0.f);

    const float cols_lim = float(cols_i-1);
    const float kdtita = cols_lim/fovh;

    for (unsigned int j = 0; j<cols_i; j++)
    {
        if (range_1[image_level](j) > 0.f)
        {
            //Transform point to the warped reference frame
            const float x_w = acu_trans(0,0)*xx_1[image_level](j) + acu_trans(0,1)*yy_1[image_level](j) + acu_trans(0,2);
            const float y_w = acu_trans(1,0)*xx_1[image_level](j) + acu_trans(1,1)*yy_1[image_level](j) + acu_trans(1,2);
            const float tita_w = atan2(y_w, x_w);
            const float range_w = sqrt(x_w*x_w + y_w*y_w);

            //Calculate warping
            const float uwarp = kdtita*(tita_w + 0.5*fovh);

            //The warped pixel (which is not integer in general) contributes to all the surrounding ones
            if ((uwarp >= 0.f)&&(uwarp < cols_lim))
            {
                const int uwarp_l = uwarp;
                const int uwarp_r = uwarp_l + 1;
                const float delta_r = float(uwarp_r) - uwarp;
                const float delta_l = uwarp - float(uwarp_l);

                //Very close pixel
                if (abs(round(uwarp) - uwarp) < 0.05f)
                {
                    range_warped[image_level](round(uwarp)) += range_w;
                    wacu(round(uwarp)) += 1.f;
                }
                else
                {
                    const float w_r = square(delta_l);
                    range_warped[image_level](uwarp_r) += w_r*range_w;
                    wacu(uwarp_r) += w_r;

                    const float w_l = square(delta_r);
                    range_warped[image_level](uwarp_l) += w_l*range_w;
                    wacu(uwarp_l) += w_l;
                }
            }
        }
    }

    //Scale the averaged range and compute coordinates
    for (unsigned int u = 0; u<cols_i; u++)
    {
        if (wacu(u) > 0.f)
        {
            const float tita = -0.5f*fovh + float(u)/kdtita;
            range_warped[image_level](u) /= wacu(u);
            xx_warped[image_level](u) = range_warped[image_level](u)*cos(tita);
            yy_warped[image_level](u) = range_warped[image_level](u)*sin(tita);
        }
        else
        {
            range_warped[image_level](u) = 0.f;
            xx_warped[image_level](u) = 0.f;
            yy_warped[image_level](u) = 0.f;
        }
    }
}

void RF2O_3S::warpScan3To2()
{
    //Use the previous transformation to warp the scan 3 forward (to 2)
    Matrix3f acu_trans_inv = overall_trans_prev.inverse();

    //Create forward-warped scans for every level
    for (unsigned int i=0; i<ctf_levels; i++)
    {
        unsigned int s = pow(2.f,int(ctf_levels-(i+1)));
        cols_i = ceil(float(cols)/float(s));
        image_level = ctf_levels - i + round(log2(round(float(width)/float(cols)))) - 1;


        ArrayXf wacu(cols_i); wacu.fill(0.f);
        range_3_warpedTo2[image_level].fill(0.f);

        const float cols_lim = float(cols_i-1);
        const float kdtita = cols_lim/fovh;

        for (unsigned int j = 0; j<cols_i; j++)
        {
            if (range_3[image_level](j) > 0.f)
            {
                //Transform point to the warped reference frame
                const float x_w = acu_trans_inv(0,0)*xx_3[image_level](j) + acu_trans_inv(0,1)*yy_3[image_level](j) + acu_trans_inv(0,2);
                const float y_w = acu_trans_inv(1,0)*xx_3[image_level](j) + acu_trans_inv(1,1)*yy_3[image_level](j) + acu_trans_inv(1,2);
                const float tita_w = atan2(y_w, x_w);
                const float range_w = sqrt(x_w*x_w + y_w*y_w);

                //Calculate warping
                const float uwarp = kdtita*(tita_w + 0.5*fovh);

                //The warped pixel (which is not integer in general) contributes to all the surrounding ones
                if ((uwarp >= 0.f)&&(uwarp < cols_lim))
                {
                    const int uwarp_l = uwarp;
                    const int uwarp_r = uwarp_l + 1;
                    const float delta_r = float(uwarp_r) - uwarp;
                    const float delta_l = uwarp - float(uwarp_l);

                    //Very close pixel
                    if (abs(round(uwarp) - uwarp) < 0.05f)
                    {
                        range_3_warpedTo2[image_level](round(uwarp)) += range_w;
                        wacu(round(uwarp)) += 1.f;
                    }
                    else
                    {
                        const float w_r = square(delta_l);
                        range_3_warpedTo2[image_level](uwarp_r) += w_r*range_w;
                        wacu(uwarp_r) += w_r;

                        const float w_l = square(delta_r);
                        range_3_warpedTo2[image_level](uwarp_l) += w_l*range_w;
                        wacu(uwarp_l) += w_l;
                    }
                }
            }
        }

        //Scale the averaged range and compute coordinates
        for (unsigned int u = 0; u<cols_i; u++)
        {
            if (wacu(u) > 0.f)
            {
                const float tita = -0.5f*fovh + float(u)/kdtita;
                range_3_warpedTo2[image_level](u) /= wacu(u);
                xx_3_warpedTo2[image_level](u) = range_3_warpedTo2[image_level](u)*cos(tita);
                yy_3_warpedTo2[image_level](u) = range_3_warpedTo2[image_level](u)*sin(tita);
            }
            else
            {
                range_3_warpedTo2[image_level](u) = 0.f;
                xx_3_warpedTo2[image_level](u) = 0.f;
                yy_3_warpedTo2[image_level](u) = 0.f;
            }
        }
    }
}


void RF2O_3S::odometryCalculation()
{
    //==================================================================================
    //						DIFERENTIAL  ODOMETRY  MULTILEVEL
    //==================================================================================

    clock.Tic();
    createScanPyramid();
    warpScan3To2();

    //Coarse-to-fine scheme
    for (unsigned int i=0; i<ctf_levels; i++)
    {
        //Previous computations
        transformations[i].setIdentity();

        level = i;
        unsigned int s = pow(2.f,int(ctf_levels-(i+1)));
        cols_i = ceil(float(cols)/float(s));
        image_level = ctf_levels - i + round(log2(round(float(width)/float(cols)))) - 1;

        for (unsigned int k=0; k<5; k++)
        {
            //1. Perform warping
            if ((i == 0)&&(k == 0))
            {
                range_warped[image_level] = range_1[image_level];
                xx_warped[image_level] = xx_1[image_level];
                yy_warped[image_level] = yy_1[image_level];
            }
            else
                performWarping();


            //2. Calculate inter coords
            calculateCoord();

            //3. Compute derivatives
            calculateRangeDerivatives();

            //4. Compute weights
            computeWeights();

            //5. Solve odometry
            if (num_valid_range > 3)
            {
                //solveSystemQuadResiduals3Scans();
                solveSystemSmoothTruncQuad3Scans();
            }

            //6. Filter solution
            filterLevelSolution();

            if (kai_loc_level.norm() < 0.05f)
            {
                //printf("\n Number of non-linear iterations: %d", k+1);
                break;
            }
        }
    }

    runtime = 1000.f*clock.Tac();
    cout << endl << "Time odometry (ms): " << runtime;

    //Update poses
    PoseUpdate();
}

void RF2O_3S::filterLevelSolution()
{
    //		Calculate Eigenvalues and Eigenvectors
    //----------------------------------------------------------
    SelfAdjointEigenSolver<Matrix3f> eigensolver(cov_odo);
    if (eigensolver.info() != Success)
    {
        printf("\n Eigensolver couldn't find a solution. Pose is not updated");
        return;
    }
	
    //First, we have to describe both the new linear and angular speeds in the "eigenvector" basis
    //-------------------------------------------------------------------------------------------------
    Matrix3f Bii = eigensolver.eigenvectors();
    Vector3f kai_b = Bii.colPivHouseholderQr().solve(kai_loc_level);


    //Second, we have to describe both the old linear and angular speeds in the "eigenvector" basis too
    //-------------------------------------------------------------------------------------------------
    Vector3f kai_loc_sub;

    //Important: we have to substract the solutions from previous levels
    Matrix3f acu_trans = Matrix3f::Identity();
    for (unsigned int i=0; i<=level; i++)
        acu_trans = transformations[i]*acu_trans;

    kai_loc_sub(0) = -fps*acu_trans(0,2);
    kai_loc_sub(1) = -fps*acu_trans(1,2);
    if (acu_trans(0,0) > 1.f)
        kai_loc_sub(2) = 0.f;
    else
        kai_loc_sub(2) = -fps*acos(acu_trans(0,0))*sign(acu_trans(1,0));
    kai_loc_sub += kai_loc_old;

    Vector3f kai_b_old = Bii.colPivHouseholderQr().solve(kai_loc_sub);

    //Filter speed
    const float cf = 5e3f*expf(-int(level)), df = 0.02f*expf(-int(level));
    //const float cf = 50e3f*expf(-int(level)), df = 0.2f*expf(-int(level));

    Vector3f kai_b_fil;
    for (unsigned int i=0; i<3; i++)
    {
        kai_b_fil(i,0) = (kai_b(i,0) + (cf*eigensolver.eigenvalues()(i,0) + df)*kai_b_old(i,0))/(1.f + cf*eigensolver.eigenvalues()(i,0) + df);
        //kai_b_fil_f(i,0) = (1.f*kai_b(i,0) + 0.f*kai_b_old_f(i,0))/(1.0f + 0.f);
    }

    //Transform filtered speed to local reference frame and compute transformation
    Vector3f kai_loc_fil = Bii.inverse().colPivHouseholderQr().solve(kai_b_fil);

    //transformation
    const float incrx = kai_loc_fil(0)/fps;
    const float incry = kai_loc_fil(1)/fps;
    const float rot = kai_loc_fil(2)/fps;

    Matrix3f new_trans = Matrix3f::Identity();
    new_trans(0,0) = cos(rot);
    new_trans(0,1) = -sin(rot);
    new_trans(1,0) = sin(rot);
    new_trans(1,1) = cos(rot);

    Matrix2f V = Matrix2f::Identity();
    Vector2f incr; incr << incrx, incry;
    if (abs(rot) > 0.001f)
    {
        const float V1 = sin(rot)/rot;
        const float V2 = (1.f - cos(rot))/rot;
        V << V1, -V2, V2, V1;
    }

    new_trans(0,2) = (V*incr)(0);
    new_trans(1,2) = (V*incr)(1);

    transformations[level] = new_trans*transformations[level];
}

void RF2O_3S::PoseUpdate()
{
    //First, compute the overall transformation
    //---------------------------------------------------
    Matrix3f acu_trans = Matrix3f::Identity();
    for (unsigned int i=1; i<=ctf_levels; i++)
        acu_trans = transformations[i-1]*acu_trans;
    overall_trans_prev = acu_trans;


    //				Compute kai_loc and kai_abs
    //--------------------------------------------------------
    kai_loc(0) = fps*acu_trans(0,2);
    kai_loc(1) = fps*acu_trans(1,2);
    if (acu_trans(0,0) > 1.f)
        kai_loc(2) = 0.f;
    else
        kai_loc(2) = fps*acos(acu_trans(0,0))*sign(acu_trans(1,0));

    float phi = laser_pose.phi();

    kai_abs(0) = kai_loc(0)*cos(phi) - kai_loc(1)*sin(phi);
    kai_abs(1) = kai_loc(0)*sin(phi) + kai_loc(1)*cos(phi);
    kai_abs(2) = kai_loc(2);


    //						Update poses
    //-------------------------------------------------------
    laser_oldpose = laser_pose;
    mrpt::poses::CPose2D pose_aux_2D(acu_trans(0,2), acu_trans(1,2), kai_loc(2)/fps);
    laser_pose = laser_pose + pose_aux_2D;



    //                  Compute kai_loc_old
    //-------------------------------------------------------
    phi = laser_pose.phi();
    kai_loc_old(0) = kai_abs(0)*cos(phi) + kai_abs(1)*sin(phi);
    kai_loc_old(1) = -kai_abs(0)*sin(phi) + kai_abs(1)*cos(phi);
    kai_loc_old(2) = kai_abs(2);
}

