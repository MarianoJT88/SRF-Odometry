/* Project: Laser odometry
   Author: Mariano Jaimez Tarifa
   Date: January 2015 */

#include "laser_odometry_nosym.h"


using namespace mrpt::utils;
using namespace Eigen;
using namespace std;


void RF2O_standard::initialize(unsigned int size, float FOV_rad, unsigned int odo_ID)
{
    ID = odo_ID;
    cols = size;
    width = size;
    fovh = FOV_rad*size/(size-1); //Exact for simulation, but I don't know how the datasets are given...
    ctf_levels = ceilf(log2(cols) - 4.3f);
    iter_irls = 8;
    filter_velocity = true;
	
    //Resize original range scan
    range_wf.resize(width);

    //Resize the transformation matrix
    transformations.resize(ctf_levels);
    for (unsigned int i = 0; i < ctf_levels; i++)
        transformations[i].resize(3,3);

	//Resize pyramid
	unsigned int s, cols_i;
    const unsigned int pyr_levels = round(log2(round(float(width)/float(cols)))) + ctf_levels;
    range.resize(pyr_levels);
    range_old.resize(pyr_levels);
    xx.resize(pyr_levels);
    xx_old.resize(pyr_levels);
    yy.resize(pyr_levels);
    yy_old.resize(pyr_levels);
	range_warped.resize(pyr_levels);
	xx_warped.resize(pyr_levels);
	yy_warped.resize(pyr_levels);

	for (unsigned int i = 0; i<pyr_levels; i++)
    {
        s = pow(2.f,int(i));
        cols_i = ceil(float(width)/float(s));

        range[i].resize(cols_i); range_old[i].resize(cols_i);
        range[i].fill(0.f); range_old[i].fill(0.f);
        xx[i].resize(cols_i); xx_old[i].resize(cols_i);
        xx[i].fill(0.f); xx_old[i].fill(0.f);
        yy[i].resize(cols_i); yy_old[i].resize(cols_i);
        yy[i].fill(0.f); yy_old[i].fill(0.f);

		if (cols_i <= cols)
		{
            range_warped[i].resize(cols_i);
            xx_warped[i].resize(cols_i);
            yy_warped[i].resize(cols_i);
		}
    }

    //Resize aux variables
    dt.resize(cols);
    dtita.resize(cols);
    weights.resize(cols);
    null.resize(cols);
    null.fill(false);
	cov_odo.assign(0.f);
    outliers.resize(cols);
    outliers.fill(false);


	//Compute gaussian mask
	g_mask[0] = 1.f/16.f; g_mask[1] = 0.25f; g_mask[2] = 6.f/16.f; g_mask[3] = g_mask[1]; g_mask[4] = g_mask[0];

    //Initialize "last velocity" as zero
	kai_abs.assign(0.f);
	kai_loc_old.assign(0.f);
}


void RF2O_standard::createScanPyramid()
{
	const float max_range_dif = 0.3f;
	
	//Push the frames back
	range_old.swap(range);
	xx_old.swap(xx);
	yy_old.swap(yy);


    //The number of levels of the pyramid does not match the number of levels used
    //in the odometry computation (because we sometimes want to finish with lower resolutions)

    unsigned int pyr_levels = round(log2(round(float(width)/float(cols)))) + ctf_levels;

    //Generate levels
    for (unsigned int i = 0; i<pyr_levels; i++)
    {
        unsigned int s = pow(2.f,int(i));
        cols_i = ceil(float(width)/float(s));
		const unsigned int i_1 = i-1;

        //              First level -> Filter, not downsample (I have tested it and it works better with this filtering)
        //---------------------------------------------------------------------------------------------------------------
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
                        range[i](u) = sum/weight;
                    }
                    else
                        range[i](u) = 0.f;

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
                        range[i](u) = sum/weight;
                    }
                    else
                        range[i](u) = 0.f;

                }
            }
		}

        //                              Downsampling
        //-----------------------------------------------------------------------------
        else
        {            
            const unsigned int cols_prev_level = range[i_1].rows();

            //Odd number of elements in the previous level
            if ((cols_prev_level % 2) == 1)
                for (unsigned int u = 0; u < cols_i; u++)
                {
                    const int u2 = 2*u;
                    const float dcenter = range[i_1](u2);

                    //Inner pixels
                    if ((u>0)&&(u<cols_i-1))
                    {
                        if (dcenter > 0.f)
                        {
                            float sum = 0.f, weight = 0.f;

                            for (int l=-2; l<3; l++)
                            {
                                const float abs_dif = abs(range[i_1](u2+l)-dcenter);
                                if (abs_dif < max_range_dif)
                                {
                                    const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
                                    weight += aux_w;
                                    sum += aux_w*range[i_1](u2+l);
                                }
                            }
                            range[i](u) = sum/weight;
                        }
                        else
                            range[i](u) = 0.f;
                    }

                    //Boundary
                    else
                    {
                        if (dcenter > 0.f)
                        {
                            float sum = 0.f, weight = 0.f;
                            const unsigned int cols_i2 = range[i_1].rows();

                            for (int l=-2; l<3; l++)
                            {
                                const int indu = u2+l;

                                if ((indu>=0)&&(indu<cols_i2))
                                {
                                    const float abs_dif = abs(range[i_1](indu)-dcenter);
                                    if (abs_dif < max_range_dif)
                                    {
                                        const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
                                        weight += aux_w;
                                        sum += aux_w*range[i_1](indu);
                                    }
                                }
                            }
                            range[i](u) = sum/weight;

                        }
                        else
                            range[i](u) = 0.f;
                    }
                }

            //Even number of elements in the previous level
            else
            {
                for (unsigned int u = 0; u < cols_i; u++)
                {
                    const int u2 = 2*u;

                    if ((range[i_1](u2) == 0.f)&&(range[i_1](u2+1) == 0.f))
                        range[i](u) = 0.f;
                    else if (range[i_1](u2) == 0.f)
                        range[i](u) = range[i_1](u2+1);
                    else if (range[i_1](u2+1) == 0.f)
                        range[i](u) = range[i_1](u2);
                    else
                        range[i](u) = 0.5f*(range[i_1](u2) + range[i_1](u2+1));
                }
            }
        }

        //Calculate coordinates "xy" of the points
        for (unsigned int u = 0; u < cols_i; u++) 
		{
            if (range[i](u) > 0.f)
			{
                const float tita = -0.5f*fovh + (float(u) + 0.5f)*fovh/float(cols_i);
				xx[i](u) = range[i](u)*cos(tita);
				yy[i](u) = range[i](u)*sin(tita);
			}
			else
			{
				xx[i](u) = 0.f;
				yy[i](u) = 0.f;
			}
		}
    }
}

void RF2O_standard::calculateCoord()
{		
    null.fill(false);
    num_valid_range = 0;

    for (unsigned int u = 0; u < cols_i; u++)
	{
		if ((range_old[image_level](u) == 0.f) || (range_warped[image_level](u) == 0.f))
		{
            null(u) = true;
		}
		else
		{
            null(u) = false;
            if ((u>0)&&(u<cols_i-1))
                num_valid_range++;
		}
	}
}

void RF2O_standard::calculaterangeDerivativesSurface()
{	
    //Compute distances between points
    rtita.resize(cols_i);
    rtita.fill(1.f);

	for (unsigned int u = 0; u < cols_i-1; u++)
    {
        const float dist = square(xx_warped[image_level](u+1) - xx_warped[image_level](u))
                         + square(yy_warped[image_level](u+1) - yy_warped[image_level](u));
		if (dist  > 0.f)
            rtita(u) = sqrtf(dist);
	}

    //Spatial derivatives
    for (unsigned int u = 1; u < cols_i-1; u++)
        dtita(u) = (rtita(u-1)*(range_warped[image_level](u+1)-range_warped[image_level](u)) + rtita(u)*(range_warped[image_level](u) - range_warped[image_level](u-1)))/(rtita(u)+rtita(u-1));

	dtita(0) = dtita(1);
	dtita(cols_i-1) = dtita(cols_i-2);

	//Temporal derivative
	for (unsigned int u = 0; u < cols_i; u++)
        dt(u) = range_warped[image_level](u) - range_old[image_level](u);

}


void RF2O_standard::computeWeights()
{
	//The maximum weight size is reserved at the constructor
    weights.fill(0.f);
	
	//Parameters for error_linearization
    const float kd = 0.01f; // 1.f
    const float k2d = 2e-4f; //0.2f
    const float sensor_sigma = 4e-4f;
	
	for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
		{	
			//							Compute derivatives
			//-----------------------------------------------------------------------
            //const float ini_dtita = range_old[image_level](u+1) - range_old[image_level](u-1);
            //const float final_dtita = range_warped[image_level](u+1) - range_warped[image_level](u-1);
            //const float dtitat = ini_dtita - final_dtita;
			const float dtita2 = dtita(u+1) - dtita(u-1);

            const float w_der = kd*(square(dt(u)) + square(dtita(u))) + k2d*square(dtita2) + sensor_sigma;
            weights(u) = sqrtf(1.f/w_der);
		}

    const float inv_max = 1.f/weights.maxCoeff();
    //printf("\n Max weight = %f", 1.f/inv_max);
	weights = inv_max*weights;
}


void RF2O_standard::solveSystemQuadResiduals()
{
	A.resize(num_valid_range,3);
    B.resize(num_valid_range);
	unsigned int cont = 0;
    const float kdtita = (cols_i)/fovh;
    const float inv_kdtita = 1.f/kdtita;

	//Fill the matrix A and the vector B
	//The order of the variables will be (vx, vy, wz)

	for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
		{
			// Precomputed expressions
			const float tw = weights(u);
            const float tita = -0.5f*fovh + (float(u) + 0.5f)*inv_kdtita;
            const float cos_tita = cos(tita);
            const float sin_tita = sin(tita);

			//Fill the matrix A
            A(cont, 0) = tw*(cos_tita + dtita(u)*kdtita*sin_tita/range_warped[image_level](u));
            A(cont, 1) = tw*(sin_tita - dtita(u)*kdtita*cos_tita/range_warped[image_level](u));
            A(cont, 2) = tw*(-yy_warped[image_level](u)*cos_tita + xx_warped[image_level](u)*sin_tita - dtita(u)*kdtita);
            B(cont) = tw*(-dt(u));

			cont++;
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

void RF2O_standard::solveSystemQuadResidualsNoPreW()
{
    A.resize(num_valid_range,3);
    B.resize(num_valid_range);
    unsigned int cont = 0;
    const float kdtita = (cols_i)/fovh;
    const float inv_kdtita = 1.f/kdtita;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
        {
            // Precomputed expressions
            const float tita = -0.5f*fovh + (float(u) + 0.5f)*inv_kdtita;
            const float cos_tita = cos(tita);
            const float sin_tita = sin(tita);

            //Fill the matrix A
            A(cont, 0) = cos_tita + dtita(u)*kdtita*sin_tita/range_warped[image_level](u);
            A(cont, 1) = sin_tita - dtita(u)*kdtita*cos_tita/range_warped[image_level](u);
            A(cont, 2) = -yy_warped[image_level](u)*cos_tita + xx_warped[image_level](u)*sin_tita - dtita(u)*kdtita;
            B(cont) = -dt(u);

            cont++;
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


void RF2O_standard::solveSystemMCauchy()
{
	A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
    B.resize(num_valid_range); Bw.resize(num_valid_range);
	unsigned int cont = 0;
    const float kdtita = float(cols_i)/fovh;

	//Fill the matrix A and the vector B
	//The order of the variables will be (vx, vy, wz)

	for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
		{
			// Precomputed expressions
			const float tw = weights(u);
            const float tita = -0.5f*fovh + (float(u) + 0.5f)/kdtita;

			//Fill the matrix A
            A(cont, 0) = tw*(cos(tita) + dtita(u)*kdtita*sin(tita)/range_warped[image_level](u));
            A(cont, 1) = tw*(sin(tita) - dtita(u)*kdtita*cos(tita)/range_warped[image_level](u));
            A(cont, 2) = tw*(-yy_warped[image_level](u)*cos(tita) + xx_warped[image_level](u)*sin(tita) - dtita(u)*kdtita);
            B(cont) = tw*(-dt(u));

			cont++;
		}
	
	//Solve the linear system of equations using a minimum least squares method
	MatrixXf AtA, AtB;
	AtA.multiply_AtA(A);
	AtB.multiply_AtB(A,B);
    kai_loc_level = AtA.ldlt().solve(AtB);
    VectorXf res = A*kai_loc_level - B;
	//cout << endl << "max res: " << res.maxCoeff();
	//cout << endl << "min res: " << res.minCoeff();

    //Compute the average dt and res
    float aver_dt = 0.f, aver_res = 0.f; unsigned int ind = 0;
    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
        {
            aver_dt += fabsf(dt(u));
            aver_res += fabsf(res(ind++));
        }
    aver_dt /= cont; aver_res /= cont;
    const float k = 10.f/aver_dt; //200

    ////Compute the energy
    //float energy = 0.f;
    //for (unsigned int i=0; i<res.rows(); i++)
    //	energy += log(1.f + square(k*res(i)));
    //printf("\n\nEnergy(0) = %f", energy);

    //Solve iterative reweighted least squares
    //===================================================================
    for (unsigned int i=1; i<=iter_irls; i++)
    {
        cont = 0;

        for (unsigned int u = 1; u < cols_i-1; u++)
            if (null(u) == false)
            {
                const float res_weight = sqrtf(1.f/(1.f + square(k*res(cont))));

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
            }

        //Solve the linear system of equations using a minimum least squares method
        AtA.multiply_AtA(Aw);
        AtB.multiply_AtB(Aw,Bw);
        kai_loc_level = AtA.ldlt().solve(AtB);
        res = A*kai_loc_level - B;

        ////Compute the energy
        //energy = 0.f;
        //for (unsigned int j=0; j<res.rows(); j++)
        //	energy += log(1.f + square(k*res(j)));
        //printf("\nEnergy(%d) = %f", i, energy);
    }

    //Covariance calculation
	cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();
}

void RF2O_standard::solveSystemTruncatedQuad()
{
    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
    B.resize(num_valid_range); Bw.resize(num_valid_range);
    unsigned int cont = 0;
    const float kdtita = float(cols_i)/fovh;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
        {
            // Precomputed expressions
            const float tw = weights(u);
            const float tita = -0.5f*fovh + (float(u) + 0.5f)/kdtita;

            //Fill the matrix A
            A(cont, 0) = tw*(cos(tita) + dtita(u)*kdtita*sin(tita)/range_warped[image_level](u));
            A(cont, 1) = tw*(sin(tita) - dtita(u)*kdtita*cos(tita)/range_warped[image_level](u));
            A(cont, 2) = tw*(-yy[image_level](u)*cos(tita) + xx[image_level](u)*sin(tita) - dtita(u)*kdtita);
            B(cont) = tw*(-dt(u));

            cont++;
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
    float res_median = aux_vector.at(res.rows()/2);

    //Compute the median absolute deviation
    aux_vector.clear();
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(abs(res(k) - res_median));
    std::sort(aux_vector.begin(), aux_vector.end());
    float mad = aux_vector.at(res.rows()/2);

    //Find the m-estimator constant
    float c = 5.f*mad;

    ////Compute the energy
    //float energy = 0.f;
    //for (unsigned int i=0; i<res.rows(); i++)
    //	energy += log(1.f + square(k*res(i)));
    //printf("\n\nEnergy(0) = %f", energy);

    //Solve iterative reweighted least squares
    //===================================================================
    for (unsigned int i=1; i<=iter_irls; i++)
    {
        cont = 0;

        for (unsigned int u = 1; u < cols_i-1; u++)
            if (null(u) == false)
            {
                float res_weight;
                if (abs(res(cont)) <= c)    res_weight = 1.f;
                else                        res_weight = 0.f;

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
            }

        //Solve the linear system of equations using a minimum least squares method
        AtA.multiply_AtA(Aw);
        AtB.multiply_AtB(Aw,Bw);
        kai_loc_level = AtA.ldlt().solve(AtB);
        res = A*kai_loc_level - B;

        ////Compute the energy
        //energy = 0.f;
        //for (unsigned int j=0; j<res.rows(); j++)
        //	energy += log(1.f + square(k*res(j)));
        //printf("\nEnergy(%d) = %f", i, energy);

        //Recompute c
        //-------------------------------------------------
        //Compute the median of res
        aux_vector.clear();
        for (unsigned int k = 0; k<res.rows(); k++)
            aux_vector.push_back(res(k));
        std::sort(aux_vector.begin(), aux_vector.end());
        res_median = aux_vector.at(res.rows()/2);

        //Compute the median absolute deviation
        aux_vector.clear();
        for (unsigned int k = 0; k<res.rows(); k++)
            aux_vector.push_back(abs(res(k) - res_median));
        std::sort(aux_vector.begin(), aux_vector.end());
        mad = aux_vector.at(res.rows()/2);

        //Find the m-estimator constant
        c = 5.f*mad;

    }

    //Covariance calculation
    cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();

    //Update the outlier mask
    cont = 0; outliers.fill(false); unsigned int num_outliers = 0;
    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
            if (abs(res(cont++)) > c)
            {
                outliers(u) = true;
                num_outliers++;
            }

    printf("\n Num_outliers = %d", num_outliers);
}

void RF2O_standard::solveSystemSmoothTruncQuad()
{
    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
    B.resize(num_valid_range); Bw.resize(num_valid_range);
    unsigned int cont = 0;

    const float kdtita = float(cols_i)/fovh;
    const float inv_kdtita = 1.f/kdtita;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
        {
            // Precomputed expressions
            const float tw = weights(u);
            const float tita = -0.5f*fovh + (float(u) + 0.5f)*inv_kdtita;
            const float cos_tita = cos(tita);
            const float sin_tita = sin(tita);

            //Fill the matrix A
            A(cont, 0) = tw*(cos_tita + dtita(u)*kdtita*sin_tita/range_warped[image_level](u));
            A(cont, 1) = tw*(sin_tita - dtita(u)*kdtita*cos_tita/range_warped[image_level](u));
            A(cont, 2) = tw*(-yy_warped[image_level](u)*cos_tita + xx_warped[image_level](u)*sin_tita - dtita(u)*kdtita);
            B(cont) = tw*(-dt(u));

            cont++;
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
    float res_median = aux_vector.at(res.rows()/2);

    //Compute the median absolute deviation
    aux_vector.clear();
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(abs(res(k) - res_median));
    std::sort(aux_vector.begin(), aux_vector.end());
    float mad = aux_vector.at(res.rows()/2);

    //Find the m-estimator constant
    const float c = 4.f*mad; //This seems to be the best (4) - 5
    const float inv_c = 1.f/c;
    const float squared_c = square(c);

    //Compute the energy
    float new_energy = 0.f, last_energy;
    for (unsigned int i=0; i<res.rows(); i++)
    {
        if (abs(res(i)) < c)     new_energy += 0.5f*square(res(i))*(1.f - 0.5f*square(res(i)*inv_c));
        else                     new_energy += 0.25f*squared_c;
    }
    //printf("\n\nEnergy(0) = %f", new_energy);
    last_energy = 2.f*new_energy;
    unsigned int iter = 1;

    //Solve iteratively reweighted least squares
    //===================================================================
    while ((new_energy < 0.995f*last_energy)&&(iter < 10))
    {
        cont = 0;
        last_energy = new_energy;

        for (unsigned int u = 1; u < cols_i-1; u++)
            if (null(u) == false)
            {
                float res_weight;
                if (abs(res(cont)) <= c)    res_weight = (1.f - square(res(cont)*inv_c));
                else                        res_weight = 0.f;

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
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
            if (abs(res(i)) < c)    new_energy += 0.5f*square(res(i))*(1.f - 0.5f*square(res(i)*inv_c));
            else                    new_energy += 0.25f*squared_c;
        }
        //printf("\nEnergy(%d) = %f", iter, new_energy);
        iter++;


        //Recompute c
        //-------------------------------------------------
        //Compute the median of res
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
//        c = 4.f*mad;

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

void RF2O_standard::solveSystemSmoothTruncQuadNoPreW()
{
    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
    B.resize(num_valid_range); Bw.resize(num_valid_range);
    unsigned int cont = 0;
    const float kdtita = float(cols_i)/fovh;
    const float inv_kdtita = 1.f/kdtita;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
        {
            // Precomputed expressions
            const float tita = -0.5f*fovh + (float(u) + 0.5f)*inv_kdtita;
            const float cos_tita = cos(tita);
            const float sin_tita = sin(tita);

            //Fill the matrix A
            A(cont, 0) = cos_tita + dtita(u)*kdtita*sin_tita/range_warped[image_level](u);
            A(cont, 1) = sin_tita - dtita(u)*kdtita*cos_tita/range_warped[image_level](u);
            A(cont, 2) = -yy_warped[image_level](u)*cos_tita + xx_warped[image_level](u)*sin_tita - dtita(u)*kdtita;
            B(cont) = -dt(u);

            cont++;
        }

    //Solve the linear system of equations using a minimum least squares method
    MatrixXf AtA, AtB;
    AtA.multiply_AtA(A);
    AtB.multiply_AtB(A,B);
    kai_loc_level = AtA.ldlt().solve(AtB);
    VectorXf res = A*kai_loc_level - B;


    //Compute the median of res
    vector<float> aux_vector;
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(res(k));
    std::sort(aux_vector.begin(), aux_vector.end());
    float res_median = aux_vector.at(res.rows()/2);

    //Compute the median absolute deviation
    aux_vector.clear();
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(abs(res(k) - res_median));
    std::sort(aux_vector.begin(), aux_vector.end());
    float mad = aux_vector.at(res.rows()/2);

    //Find the m-estimator constant
    const float c = 4.f*mad; //This seems to be the best (4)
    const float inv_c = 1.f/c;

    //Compute the energy
    float new_energy = 0.f, last_energy;
    for (unsigned int i=0; i<res.rows(); i++)
    {
        if (abs(res(i)) < c)     new_energy += 0.5f*square(res(i))*(1.f - 0.5f*square(res(i)*inv_c));
        else                     new_energy += 0.25f*square(c);
    }
    //printf("\n\nEnergy(0) = %f", new_energy);
    last_energy = 2.f*new_energy + 1.f;
    unsigned int iter = 1;

    //Solve iteratively reweighted least squares
    //===================================================================
    while ((new_energy < 0.995f*last_energy)&&(iter < 10))
    {
        cont = 0;
        last_energy = new_energy;

        for (unsigned int u = 1; u < cols_i-1; u++)
            if (null(u) == false)
            {
                float res_weight;
                if (abs(res(cont)) <= c)    res_weight = (1.f - square(res(cont)*inv_c));
                else                        res_weight = 0.f;

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
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
            if (abs(res(i)) < c)    new_energy += 0.5f*square(res(i))*(1.f - 0.5f*square(res(i)*inv_c));
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
}

void RF2O_standard::solveSystemSmoothTruncQuadFromBeginning()
{
    A.resize(num_valid_range,3); Aw.resize(num_valid_range,3);
    B.resize(num_valid_range); Bw.resize(num_valid_range);
    unsigned int cont = 0;
    const float kdtita = float(cols_i)/fovh;

    //Fill the matrix A and the vector B
    //The order of the variables will be (vx, vy, wz)

    for (unsigned int u = 1; u < cols_i-1; u++)
        if (null(u) == false)
        {
            // Precomputed expressions
            const float tw = weights(u);
            const float tita = -0.5f*fovh + (float(u) + 0.5f)/kdtita;

            //Fill the matrix A
            A(cont, 0) = tw*(cos(tita) + dtita(u)*kdtita*sin(tita)/range_warped[image_level](u));
            A(cont, 1) = tw*(sin(tita) - dtita(u)*kdtita*cos(tita)/range_warped[image_level](u));
            A(cont, 2) = tw*(-yy_warped[image_level](u)*cos(tita) + xx_warped[image_level](u)*sin(tita) - dtita(u)*kdtita);
            B(cont) = tw*(-dt(u));

            cont++;
        }

    //Create variables and compute initial residuals
    MatrixXf AtA, AtB;
    VectorXf res = B;
    //cout << endl << "max res: " << res.maxCoeff();
    //cout << endl << "min res: " << res.minCoeff();

    //Compute the median of res
    vector<float> aux_vector;
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(res(k));
    std::sort(aux_vector.begin(), aux_vector.end());
    float res_median = aux_vector.at(res.rows()/2);

    //Compute the median absolute deviation
    aux_vector.clear();
    for (unsigned int k = 0; k<res.rows(); k++)
        aux_vector.push_back(abs(res(k) - res_median));
    std::sort(aux_vector.begin(), aux_vector.end());
    float mad = aux_vector.at(res.rows()/2);

    //Find the m-estimator constant
    float c = 5.f*mad;

    ////Compute the energy
    //float energy = 0.f;
    //for (unsigned int i=0; i<res.rows(); i++)
    //	energy += log(1.f + square(k*res(i)));
    //printf("\n\nEnergy(0) = %f", energy);

    //Solve iterative reweighted least squares
    //===================================================================
    for (unsigned int i=0; i<=iter_irls; i++)
    {
        cont = 0;

        for (unsigned int u = 1; u < cols_i-1; u++)
            if (null(u) == false)
            {
                float res_weight;
                if (abs(res(cont)) <= c)    res_weight = (1.f - square(res(cont)/c));
                else                        res_weight = 0.f;

                //Fill the matrix Aw
                Aw(cont,0) = res_weight*A(cont,0);
                Aw(cont,1) = res_weight*A(cont,1);
                Aw(cont,2) = res_weight*A(cont,2);
                Bw(cont) = res_weight*B(cont);
                cont++;
            }

        //Solve the linear system of equations using a minimum least squares method
        AtA.multiply_AtA(Aw);
        AtB.multiply_AtB(Aw,Bw);
        kai_loc_level = AtA.ldlt().solve(AtB);
        res = A*kai_loc_level - B;

        ////Compute the energy
        //energy = 0.f;
        //for (unsigned int j=0; j<res.rows(); j++)
        //	energy += log(1.f + square(k*res(j)));
        //printf("\nEnergy(%d) = %f", i, energy);

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
}


void RF2O_standard::performWarping()
{
	Matrix3f acu_trans; 
	acu_trans.setIdentity();
    for (unsigned int i=0; i<=level; i++)
        acu_trans = transformations[i]*acu_trans;

    ArrayXf wacu(cols_i);
    wacu.fill(0.f);
    range_warped[image_level].fill(0.f);

	const float cols_lim = float(cols_i-1);
    const float kdtita = cols_i/fovh;

	for (unsigned int j = 0; j<cols_i; j++)
	{				
		if (range[image_level](j) > 0.f)
		{
			//Transform point to the warped reference frame
			const float x_w = acu_trans(0,0)*xx[image_level](j) + acu_trans(0,1)*yy[image_level](j) + acu_trans(0,2);
			const float y_w = acu_trans(1,0)*xx[image_level](j) + acu_trans(1,1)*yy[image_level](j) + acu_trans(1,2);
			const float tita_w = atan2(y_w, x_w);
			const float range_w = sqrt(x_w*x_w + y_w*y_w);

			//Calculate warping
            const float uwarp = kdtita*(tita_w + 0.5*fovh) - 0.5f;

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
            const float tita = -0.5f*fovh + (float(u) + 0.5f)/kdtita;
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

void RF2O_standard::performBestWarping()
{
    Matrix3f acu_trans;
    acu_trans.setIdentity();
    for (unsigned int i=0; i<=level; i++)
        acu_trans = transformations[i]*acu_trans;

    ArrayXf x_trans(cols_i), y_trans(cols_i), u_trans(cols_i), range_trans(cols_i);
    x_trans.fill(0.f); y_trans.fill(0.f); range_trans.fill(0.f);
    range_warped[image_level].fill(0.f);

    const float kdtita = float(cols_i)/fovh;

    //Transform points to the reference frame of the old scan
    for (unsigned int u = 0; u<cols_i; u++)
    {
        if (range[image_level](u) != 0.f)
        {
            //Transform point to the warped reference frame
            x_trans(u) = acu_trans(0,0)*xx[image_level](u) + acu_trans(0,1)*yy[image_level](u) + acu_trans(0,2);
            y_trans(u) = acu_trans(1,0)*xx[image_level](u) + acu_trans(1,1)*yy[image_level](u) + acu_trans(1,2);
            range_trans(u) = sqrtf(square(x_trans(u)) + square(y_trans(u)));
            const float tita_trans = atan2(y_trans(u), x_trans(u));
            u_trans(u) = kdtita*(tita_trans + 0.5f*fovh) - 0.5f;
        }
    }

    //Check projection for each segment
    for (unsigned int u = 0; u<cols_i-1; u++)
    {
        if ((range_trans(u) == 0.f) || (range_trans(u+1) == 0.f))
            continue;
        else if (floorf(u_trans(u)) != floorf(u_trans(u+1)))
        {
           const float range_l = (floorf(u_trans(u)) > floorf(u_trans(u+1))) ? range_trans(u+1) : range_trans(u);
           const float range_r = (floorf(u_trans(u)) > floorf(u_trans(u+1))) ? range_trans(u) : range_trans(u+1);
           const float u_trans_l = (floorf(u_trans(u)) > floorf(u_trans(u+1))) ? u_trans(u+1) : u_trans(u);
           const float u_trans_r = (floorf(u_trans(u)) > floorf(u_trans(u+1))) ? u_trans(u) : u_trans(u+1);
           const int u_l = min(floorf(u_trans(u)), floorf(u_trans(u+1)));
           const int u_r = max(floorf(u_trans(u)), floorf(u_trans(u+1)));

           for (unsigned int u_segment=u_l+1; (u_segment<=u_r)&&(u_segment<cols_i)&&(u_segment>=0); u_segment++)
           {
               const float range_interp = ((u_segment - u_trans_l)*range_r + (u_trans_r - u_segment)*range_l)/(u_trans_r - u_trans_l);
               if ((range_warped[image_level](u_segment) == 0.f)||(range_interp < range_warped[image_level](u_segment)))
                   range_warped[image_level](u_segment) = range_interp;
           }
        }
    }

    //Compute coordinates
    for (unsigned int u = 0; u<cols_i; u++)
    {
        const float tita = -0.5f*fovh + (float(u) + 0.5f)/kdtita;
        xx_warped[image_level](u) = range_warped[image_level](u)*cos(tita);
        yy_warped[image_level](u) = range_warped[image_level](u)*sin(tita);
    }
}

void RF2O_standard::performFastWarping()
{
    //Warp the second image and count the amount of pixels projected to each of the pixels in the first image
    //Camera parameters (which also depend on the level resolution)
    Matrix3f acu_trans, acu_trans_inv;
    acu_trans.setIdentity();
    for (unsigned int i=0; i<=level; i++)
        acu_trans = transformations[i]*acu_trans;

    acu_trans_inv = acu_trans.inverse();
    range_warped[image_level].fill(0.f);

    const float kdtita = cols_i/fovh;

    for (unsigned int u = 0; u<cols_i; u++)
    {
        const float r = range_old[image_level](u);
        if (r > 0.f)
        {
            //Transform point to the warped reference frame **********************************************
            const float x_w = acu_trans_inv(0,0)*xx_old[image_level](u) + acu_trans_inv(0,1)*yy_old[image_level](u) + acu_trans_inv(0,2);
            const float y_w = acu_trans_inv(1,0)*xx_old[image_level](u) + acu_trans_inv(1,1)*yy_old[image_level](u) + acu_trans_inv(1,2);
            const float tita_w = atan2(y_w, x_w);

            //Calculate warping
            const float uwarp = kdtita*(tita_w + 0.5*fovh) - 0.5f;

            float r_2;
            interpolateRange(r_2, uwarp);

            if (r_2 > 0.f)
            {
                float r_w = sqrtf(square(x_w) + square(y_w));
                range_warped[image_level](u) = r_2 - (r_w-r);
            }
            else
                range_warped[image_level](u) = 0.f;

            //Transform point back to coordinates of the old scan -> xx_warped[image_level](u) and yy...
            const float tita_original = -0.5f*fovh + (u + 0.5f)*fovh/cols_i;
            xx_warped[image_level](u) = range_warped[image_level](u)*cos(tita_original);
            yy_warped[image_level](u) = range_warped[image_level](u)*sin(tita_original);
        }
    }
}

void RF2O_standard::interpolateRange(float &range_pixel, float uwarp)
{
    if ((uwarp <= 0.f)||(uwarp >= cols_i-1))
    {
        range_pixel = 0.f;
        return;
    }
    else
    {
        const unsigned int u_l = floorf(uwarp);
        const unsigned int u_r = ceilf(uwarp);
        const float range_l = range[image_level](u_l);
        const float range_r = range[image_level](u_r);

        if (range_l*range_r == 0)
            range_pixel = 0.f;
        else
            range_pixel = (uwarp - u_l)*range_r + (u_r - uwarp)*range_l;
    }
}

void RF2O_standard::odometryCalculation()
{
	//==================================================================================
	//						DIFERENTIAL  ODOMETRY  MULTILEVEL
	//==================================================================================

    clock.Tic();
    transf_acu_per_iteration.clear();
    transf_level.clear();
    acu_trans_overall.setIdentity();
    createScanPyramid();

    //Coarse-to-fine scheme
    for (unsigned int i=0; i<ctf_levels; i++)
    {
        //Previous computations
        transformations[i].setIdentity();

        level = i;
        unsigned int s = pow(2.f,int(ctf_levels-(i+1)));
        cols_i = ceil(float(cols)/float(s));
        image_level = ctf_levels - i + round(log2(round(float(width)/float(cols)))) - 1;

        const unsigned int nonlin_iters = 3;
        for (unsigned int k = 0; k<nonlin_iters; k++)
        {
            //1. Perform warping
            if ((i == 0)&&(k == 0))
            {
                range_warped[image_level] = range[image_level];
                xx_warped[image_level] = xx[image_level];
                yy_warped[image_level] = yy[image_level];
            }
            else
                performBestWarping();

            //1.5 - Check residuals
//            float r1, r2, r3;
//            computeAverageResiduals(r1, r2, r3);
//            if (k == 0) printf("\n");
//            printf("\n Level = %d, iter = %d, aver_r = %f, trunc_aver_r = %f, median_r = %f", level, k, r1, r2, r3);


            //2. Calculate inter coords
            calculateCoord();

            //3. Compute derivatives
            calculaterangeDerivativesSurface();

            //4. Compute weights
            computeWeights();

            //5. Solve odometry
            if (num_valid_range > 3)
            {
                if (ID == 0)
                    //solveSystemSmoothTruncQuad();
                    solveSystemQuadResidualsNoPreW();
                else if (ID == 1)
                    solveSystemQuadResiduals();
                else if (ID == 2)
                    solveSystemSmoothTruncQuadNoPreW();
                else
                    solveSystemSmoothTruncQuad();
            }

            //6. Filter solution
            filterLevelSolution();


            if (kai_loc_level.norm() < 0.05f)
            {
                //printf("\n Number of non-linear iterations: %d", k+1);
                break;
            }
            else if (k == nonlin_iters-1)
            {
                ;//printf("\n Number of non-linear iterations: %d", k+1);
            }
        }

        //Check residuals end level
//        float r1, r2, r3;
//        computeAverageResiduals(r1, r2, r3);
//        printf("\n Level = %d, iter = -, aver_r = %f, trunc_aver_r = %f, median_r = %f", level, r1, r2, r3);
    }

    runtime = 1000.f*clock.Tac();
    cout << endl << "Time odometry (ms): " << runtime;

    //Update poses
    PoseUpdate();
}

void RF2O_standard::filterLevelSolution()
{
    Vector3f kai2Pose = kai_loc_level;

    if (filter_velocity)
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

        kai_loc_sub(0) = -acu_trans(0,2);
        kai_loc_sub(1) = -acu_trans(1,2);
        if (acu_trans(0,0) > 1.f)
            kai_loc_sub(2) = 0.f;
        else
            kai_loc_sub(2) = -acos(acu_trans(0,0))*sign(acu_trans(1,0));
        kai_loc_sub += kai_loc_old;

        Vector3f kai_b_old = Bii.colPivHouseholderQr().solve(kai_loc_sub);

        //Filter speed
        //const float cf = 15e3f*expf(-int(level)), df = 0.05f*expf(-int(level));
        const float cf = 5e3f*expf(-int(level)), df = 0.02f*expf(-int(level));

        Vector3f kai_b_fil;
        for (unsigned int i=0; i<3; i++)
        {
            kai_b_fil(i,0) = (kai_b(i,0) + (cf*eigensolver.eigenvalues()(i,0) + df)*kai_b_old(i,0))/(1.f + cf*eigensolver.eigenvalues()(i,0) + df);
            //kai_b_fil_f(i,0) = (1.f*kai_b(i,0) + 0.f*kai_b_old_f(i,0))/(1.0f + 0.f);
        }

        //Transform filtered speed to local reference frame and compute transformation
        Vector3f kai_loc_fil = Bii.inverse().colPivHouseholderQr().solve(kai_b_fil);
        kai2Pose = kai_loc_fil;
    }

	//transformation
    const float incrx = kai2Pose(0);
    const float incry = kai2Pose(1);
    const float rot = kai2Pose(2);

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
    acu_trans_overall = new_trans*acu_trans_overall;

    //To keep track of every single iteration
    transf_level.push_back(level);
    transf_acu_per_iteration.push_back(acu_trans_overall);
}

void RF2O_standard::PoseUpdate()
{
	//First, compute the overall transformation
	//---------------------------------------------------
    Matrix3f acu_trans = Matrix3f::Identity();
	for (unsigned int i=1; i<=ctf_levels; i++)
		acu_trans = transformations[i-1]*acu_trans;


	//				Compute kai_loc and kai_abs
	//--------------------------------------------------------
    kai_loc(0) = acu_trans(0,2);
    kai_loc(1) = acu_trans(1,2);
	if (acu_trans(0,0) > 1.f)
		kai_loc(2) = 0.f;
	else
        kai_loc(2) = acos(acu_trans(0,0))*sign(acu_trans(1,0));

    float phi = laser_pose.phi();

	kai_abs(0) = kai_loc(0)*cos(phi) - kai_loc(1)*sin(phi);
	kai_abs(1) = kai_loc(0)*sin(phi) + kai_loc(1)*cos(phi);
	kai_abs(2) = kai_loc(2);


	//						Update poses
	//-------------------------------------------------------
	laser_oldpose = laser_pose;
    mrpt::poses::CPose2D pose_aux_2D(acu_trans(0,2), acu_trans(1,2), kai_loc(2));
	laser_pose = laser_pose + pose_aux_2D;



    //                  Compute kai_loc_old
	//-------------------------------------------------------
	phi = laser_pose.phi();
	kai_loc_old(0) = kai_abs(0)*cos(phi) + kai_abs(1)*sin(phi);
	kai_loc_old(1) = -kai_abs(0)*sin(phi) + kai_abs(1)*cos(phi);
	kai_loc_old(2) = kai_abs(2);
}

void RF2O_standard::computeAverageResiduals(float &res1, float &res2, float &res3)
{
    //First, warp R2 towards R1
    //level = ctf_levels - 1;
    //image_level = 0;
    performWarping();

    res1 = 0.f; res2 = 0.f; res3 = 0.f;
    const float tau = 0.1f;
    vector<float> vec_res;
    unsigned int cont = 0;

    for (unsigned int v=0; v<cols_i; v++)
    {
        if ((range_warped[image_level](v) != 0.f)&&(range_old[image_level](v) != 0.f))
        {
            const float res = abs(range_warped[image_level](v) - range_old[image_level](v));
            res1 += res;
            res2 += min(res, tau);
            vec_res.push_back(res);
            cont++;
        }
    }

    res1 /= cont;
    res2 /= cont;

    std::sort(vec_res.begin(), vec_res.end());
    res3 = vec_res.at(vec_res.size()/2);

    //printf("\n Aver res = %f, Trunc aver res = %f, median res = %f", res1, res2, res3);
}

