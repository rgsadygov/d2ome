// This is the main DLL file.

#include "stdafx.h"
#include "LBFGS_Caller.h"


using namespace System;
using namespace std;
using namespace System::IO;
using namespace System::Collections;
using namespace System::Collections::Generic;
using namespace LBFGS_Library_Call;

int One_Compartment(float *Incorporations, float param, 
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g);

int Two_Compartment(float *Incorporations, float param1, float param2,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2);

int Three_Compartment(float *Incorporations, float param1, float param2, float param3,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2, lbfgsfloatval_t *g3);

//int One_Compartment_exponential(array <float> ^intensities, float param, float a, float b,
int One_Compartment_exponential(float *Incorporations, float param,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g);

int Two_Compartment_exponential(float *Incorporations, float param1, float param2,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2);

float *All_Times, *y_exmp;

float a, b;   // these are coefficients to be used in general OneCompartment (one exponential) case

int nTimePoints;           // All_Times, *y_exmp, nTimePoints are intermediaries, between the manage code variables fResponse, fTime, nData and static function evaluate
/*
*  x - is an array that holds the parameters
*  n is the number of parameters
*/

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x,
				lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{

	lbfgsfloatval_t fx = 0.0;

	int k = 9;             //the number of time points 
		 
	fx = 0;

	if(n == 1)
	//for (i = 0; i < n; i += 2)         //call for one-compartment model
	{
		lbfgsfloatval_t t1 = 0;

		lbfgsfloatval_t t2 = 0; 
			
		One_Compartment_exponential(y_exmp,  x[0], &t1,  &t2);

		//One_Compartment(y_exmp,  x[0], &t1,  &t2);

		//g[i+1] = 0.;

		fx = t1; g[0] = t2;
	}
	else if(2 == n)     // two-parameter, two-compartment call
	{
		lbfgsfloatval_t t1 = 0, t2 = 0, t3 = 0.;

		Two_Compartment(y_exmp, x[0], x[1], 
		&t1, &t2, &t3);

		fx = t1; g[0] = t2; g[1] = t3;
	}
	else if(3 == n) // three-parameter, three-compartment call
	{

		lbfgsfloatval_t t1 = 0, t2 = 0, t3 = 0., t4 = 0.;

		Three_Compartment(y_exmp, x[0], x[1], x[2],
			&t1, &t2, &t3, &t4);

		fx = t1; g[0] = t2; g[1] = t3; g[2] = t4;


	}
	else if(4 == n)     // two-parameter, two-compartment call
	{
		lbfgsfloatval_t t1 = 0, t2 = 0, t3 = 0.;

		Two_Compartment_exponential(y_exmp, x[0], x[1], 
		&t1, &t2, &t3);

		fx = t1; g[0] = t2; g[1] = t3;
	}

	return fx;
}

/*
*
* this function just reports the progress
*
*/

static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls )
{
    /*printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n"); */  
    return 0;
}

 
/*
*   time is the values of time series,
*   Incorporations are the values of incorporatation levels
*   nData - number of time series data, = # of time points, # of Incorporations
*   n - is the number of parameters, it has to be specific/pecular because of the
*   interpretation
*   g_current is the gradient,
    fx_current is the current value of the function, which is the sum of squares of differences between the theoretical
	(1 - exp(-param*t)) and experimental, Incorporations, data.
*/
int One_Compartment(float *Incorporations, float param, 
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g)
{
	int j;

    lbfgsfloatval_t fx_current, g_current = 0.0;

	lbfgsfloatval_t t1 = 0;

	lbfgsfloatval_t t2 = 0; 

	for(j = 0; j < nTimePoints; j++)
	{
		t1 = 1.0 - exp(-param * All_Times[j]);

		t2 = t2 + (Incorporations[j] - t1) * (Incorporations[j] - t1);

		g_current = g_current - 2.0 * (Incorporations[j] - t1) * All_Times[j] * exp(-param * All_Times[j]);
	}

	fx_current = t2;

	*fx = fx_current; *g = g_current;

	return 0;
}


/*
*   time is the values of time series,
*   Incorporations are the values of incorporatation levels
*   nData - number of time series data, = # of time points, # of Incorporations
*   n - is the number of parameters, it has to be specific/pecular because of the
*   interpretation
*   g_current1 is the gradient with respect to param1,
*   g_current2 is the gradient wrt to param2.
*   fx_current is the current value of the object function, which is the sum of squares of differences between the theoretical
*	(1 - exp(-param*t)) and experimental, Incorportations, data.
*    the algorithm minimizes the value of the objective function, fx.
*
*/
int Two_Compartment(float *Incorporations, float param1, float param2,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2)
{
	int j;

    lbfgsfloatval_t fx_current, g_current1 = 0.0, g_current2 = 0.0;

	lbfgsfloatval_t t1 = 0;

	lbfgsfloatval_t t2 = 0; 

	for(j = 0; j < nTimePoints; j++)   //nTimes = nData number of time points in time-course experiment
	{
		/*t1 = 1.0 - (param2 * exp(-All_Times[j]*param1)- param1 * exp(-All_Times[j] * param2))/
		                  (param2 - param1);*/

		t1 = a/(1.+exp(-param2) ) + (a - a/(1.+exp(-param2) )) * exp( -All_Times[j] * exp(param1) );

		t2 = t2 + (Incorporations[j] - t1) * (Incorporations[j] - t1);

		g_current1 = g_current1 + 2.0 * (Incorporations[j] - t1) * (a - a / (1. + exp(-param2) )) * All_Times[j] * exp(param1)* exp(-All_Times[j] * exp(param1));

		g_current2 = g_current2 - 2.0 * (Incorporations[j] - t1) * a * exp(-param2) / pow((1. + exp(-param2)), 2.) * 
			(1 - exp(-All_Times[j] * exp(param1))  );




		/*g_current1 = g_current1 - 2 * t1 *  ( (All_Times[j]*(param2 - param1) - 1) * exp(-param1 * All_Times[j]) + exp(-param2 * All_Times[j]) )*param2/
												((param2 - param1) * (param2 - param1));

		g_current2 = g_current2 - 2 * t1 * (exp(-param1 * All_Times[j]) - (1 + (param2 - param1) * All_Times[j])*exp(-param2 * All_Times[j]))*param1/
												((param2 - param1)*(param2 - param1));*/

		//printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g_current1, param1,  exp(-param1 * Time[j]));
	}

	fx_current = t2;

	*fx = fx_current; *g1 = g_current1; *g2 = g_current2;


	//printf("TwoParam=> f: %10.5f g1: %10.5f g2: %10.5f\n", fx_current, g_current1, g_current2);

	return 0;
}


int Three_Compartment(float *Incorporations, float param1, float param2, float param3,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2, lbfgsfloatval_t *g3)
{
	int j;

	lbfgsfloatval_t fx_current, g_current1 = 0.0, g_current2 = 0.0, g_current3 = 0.0;

	lbfgsfloatval_t t1 = 0;

	lbfgsfloatval_t t2 = 0;

	for (j = 0; j < nTimePoints; j++)   //nTimes = nData number of time points in time-course experiment
	{
		
		/*t1 = 1.0 - (param2 * exp(-All_Times[j]*param1)- param1 * exp(-All_Times[j] * param2))/
		(param2 - param1);*/

		t1 = exp(param3) / (1. + exp(-param2)) + (exp(param3) - exp(param3) / (1. + exp(-param2))) * exp(-All_Times[j] * exp(param1));

		t2 = t2 + (Incorporations[j] - t1) * (Incorporations[j] - t1);

		g_current1 = g_current1 + 2.0 * (Incorporations[j] - t1) * (exp(param3) - exp(param3) / (1. + exp(-param2))) * All_Times[j] * exp(param1)* exp(-All_Times[j] * exp(param1));

		g_current2 = g_current2 - 2.0 * (Incorporations[j] - t1) * exp(param3) * exp(-param2) / pow((1. + exp(-param2)), 2.) *
			(1 - exp(-All_Times[j] * exp(param1)));


		g_current3 = g_current3 - 2.0 * (Incorporations[j] - t1) * t1;



		/*g_current1 = g_current1 - 2 * t1 *  ( (All_Times[j]*(param2 - param1) - 1) * exp(-param1 * All_Times[j]) + exp(-param2 * All_Times[j]) )*param2/
		((param2 - param1) * (param2 - param1));

		g_current2 = g_current2 - 2 * t1 * (exp(-param1 * All_Times[j]) - (1 + (param2 - param1) * All_Times[j])*exp(-param2 * All_Times[j]))*param1/
		((param2 - param1)*(param2 - param1));*/

		//printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g_current1, param1,  exp(-param1 * Time[j]));
	}

	fx_current = t2;

	*fx = fx_current; *g1 = g_current1; *g2 = g_current2; *g3 = g_current3;


	//printf("TwoParam=> f: %10.5f g1: %10.5f g2: %10.5f\n", fx_current, g_current1, g_current2);

	return 0;
}









/*
*   time is the values of time series,
*   Incorporations are the values of incorporatation levels
*   nData - number of time series data, = # of time points, # of Incorporations
*   n - is the number of parameters, it has to be specific/pecular because of the
*   interpretation
*   g_current is the gradient,
    fx_current is the current value of the function, which is the sum of squares of differences between the theoretical
	(1 - exp(-t*exp(param))) and experimental, Incorporations, data.

	(a + b *exp(-t*exp(param)) )  is the functional form

	where are minimizign the residual sum of squares

	0.5 * (Incorportaions - (a + b * exp*-t*exp(param))) )^2

	t1 = (a + b *exp(-t*exp(param)) ) 
*
*   a and b are defined as global variables, and set values in Optimize()
*/
//int One_Compartment_exponential(float *Incorporations, float param, 
int One_Compartment_exponential(float *Incorporations, float param, 
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g)
{
	int j;

    lbfgsfloatval_t fx_current, g_current = 0.0;
	
	lbfgsfloatval_t t1 = 0;
	
	lbfgsfloatval_t t2 = 0; 

	for(j = 0; j < nTimePoints; j++)
	{
	
		/*t1 = 1.0 - exp(-All_Times[j]*exp(param));
		
		t2 = t2 + (Incorporations[j] - t1) * (Incorporations[j] - t1);
		
		g_current = g_current - 2.0 * (Incorporations[j] - t1) * All_Times[j] * exp(-All_Times[j]*exp(param))*exp(param);*/

		//printf("Values %f, %d\n", Incorporations[j], nTimePoints);

		t1 = a + b * exp(-All_Times[j]*exp(param));

		t2 = t2 + (Incorporations[j] - t1) * (Incorporations[j] - t1);

		g_current = g_current + 2.0 * (Incorporations[j] - t1) * b * All_Times[j] * exp(param)* exp(-All_Times[j]*exp(param));
	}

	fx_current = t2;
	
	*fx = fx_current; *g = g_current;

	//printf("OneParam=> f: %10.5f g: %10.5f \n", fx_current, g_current);

	return 0;
}


/*
*   time is the values of time series,
*   Incorporations are the values of incorporatation levels
*   nData - number of time series data, = # of time points, # of Incorporations
*   n - is the number of parameters, it has to be specific/pecular because of the
*   interpretation
*   g_current1 is the gradient with respect to param1,
*   g_current2 is the gradient wrt to param2.
*   fx_current is the current value of the object function, which is the sum of squares of differences between the theoretical
*	(1 - exp(-param*t)) and experimental, Incorporations, data.
*    the algorithm minimizes the value of the objective function, fx.
*
*/
int Two_Compartment_exponential(float *Incorporations, float param1, float param2,
       lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2)
{
       int j;
 
       float expParam1, expParam2;
 
    lbfgsfloatval_t fx_current, g_current1 = 0.0, g_current2 = 0.0;
 
       lbfgsfloatval_t t1 = 0;
 
       lbfgsfloatval_t t2 = 0;
 
       expParam1 = exp(param1); expParam2 = exp(param2);
		/*
       if(false)
       {
             puts("OLD VERSION");
 
             for(j = 0; j < nTimePoints; j++)   //nTimes = nData number of time points in time-course experiment
             {
                    t1 = 1.0 - (param2 * exp(-All_Times[j]*param1)- param1 * exp(-All_Times[j] * param2))/
                                                 (param2 - param1);
 
                    t1 = Incorporations[j] - t1;
 
                    t2 = t2 + t1 * t1;
 
                    g_current1 = g_current1 - 2 * t1 *  ( (All_Times[j]*(param2 - param1) - 1) * exp(-param1 * All_Times[j]) + exp(-param2 * All_Times[j]) )*param2/
                                                                                      ((param2 - param1) * (param2 - param1));
 
                    g_current2 = g_current2 - 2 * t1 * (exp(-param1 * All_Times[j]) - (1 + (param2 - param1) * All_Times[j])*exp(-param2 * All_Times[j]))*param1/
                                                                                      ((param2 - param1)*(param2 - param1));
 
                    printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g_current1, param1,  exp(-param1 * All_Times[j]));
             }
       }
	   */
       //else
       //{
             for(j = 0; j < nTimePoints; j++)   //nTimes = nData number of time points in time-course experiment
             {
                    t1 = 1.0 - (exp(param2) * exp(-All_Times[j]*exp(param1)) - exp(param1) * exp(-All_Times[j] * exp(param2)) )/
                                                 (exp(param2) - exp(param1) );
 
                    t1 = Incorporations[j] - t1;
 
                    t2 = t2 + t1 * t1;
 
                    g_current1 = g_current1 - 2 * t1 * exp(param1) * ( (All_Times[j]*(exp(param2) - exp(param1)) - 1) * exp(-exp(param1) * All_Times[j]) + exp(-exp(param2) * All_Times[j]) ) *
                                                               exp(param2)/((exp(param2) - exp(param1)) * (exp(param2) - exp(param1)) );
 
                    g_current2 = g_current2 - 2 * t1 * exp(param2) * (exp(-exp(param1) * All_Times[j]) - (1 + (exp(param2) - exp(param1) ) * All_Times[j])*exp(-exp(param2) * All_Times[j])) *
                                                               exp(param1)/((exp(param2) - exp(param1) )*(exp(param2) - exp(param1)) );
 
                   // printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g_current1, exp(param1),  exp(-exp(param1) * All_Times[j]));
             }     
      // }
 
       fx_current = t2;
 
       *fx = fx_current; *g1 = g_current1; *g2 = g_current2;
 
       //printf("g = %10.5f %10.5f\n", (*g1), param1);
 
       return 0;
}




void  LBFGS::InitializeTime()
{
	All_Times = (float *) calloc(nData, sizeof(float));

	y_exmp =  (float *) calloc(nData, sizeof(float));

	if(NULL == All_Times || NULL == y_exmp)
	{
		printf("Cannot Allocate Memory for Time\n");

		exit (1);
	}

	for(int i=0; i < nData; i++)
	{
		All_Times[i] = fTime[i];
	}

	nTimePoints = nData;
}

int LBFGS::Optimize(float * fResponse, float a1, float a2, float *rkd, float *rks, float *fx1)
{
	int i, j, k, l, ret = 0;

	lbfgsfloatval_t fx;

    lbfgsfloatval_t *x = lbfgs_malloc(nParams);
    
	if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");

        return 1;
    }

	lbfgs_parameter_t param;


	//set a and b to be used in OneCompartment Model

	a = a1; b = a2;
		
	j = k = l = 0;

	for(i = 0; i < nData; i++)
	{
		y_exmp[i] = fResponse[i];
	}

	/* Initialize the variables. */
	for (i = 0; i < nParams; i += 2) 
	{
		//x[i] = -1.2;
		//*rkd = (rand()%100)*1e-2;
		//*rks = (rand()%100)*1e-1;
		//*rkd = (float) dist[0];
		//*rks = (float) dist[1];
		//fprintf(fp,"%f, %f",*rkd,*rks);

		x[i] = 0.1;//*rkd//0.1;//0.0413;
	    
		if(nParams > 1 && i+1< nParams)
		{
			x[i+1] = 1;//*rks;//1;//100;//8.2928e+11;//0.1;
		}
	}
	//x[0] = 0.1; x[1] = 0.11; x[2] = 0.12;

	/* Initialize the parameters for the L-BFGS optimization. */
	lbfgs_parameter_init(&param);
	/*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/

	/*
		Start the L-BFGS optimization; this will invoke the callback functions
		evaluate() and progress() when necessary.
		*/

	ret = lbfgs(nParams, x, &fx, evaluate, progress, NULL, &param);

	*fx1 = fx;
	
	for(i = 0; i < nParams; i++)
	{
		fParams[i] = x[i];
	}

	/* Report the result. */
	//printf("L-BFGS optimization terminated with status code = %d\n", ret);
	//printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);

	lbfgs_free(x);

	return ret;
}
/*
* release the memory allocated for fTime
*/
void  LBFGS::Release_Memory()
{
	free(fTime);

	free(All_Times);

	free(fParams);
}
