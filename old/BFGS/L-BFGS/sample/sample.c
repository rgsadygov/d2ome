#include <stdio.h>
#include <lbfgs.h>
#include <math.h>
#include <string.h>




int One_Compartment(float *Incorporations, float *Time, int nData, float param, 
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g);

int Two_Compartment(float *Incorporations, float *Time, int nData, float param1, float param2,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2);

float y_exmp[9];

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    int i;

    lbfgsfloatval_t fx = 0.0;

	if(0)
	{
		for (i = 0; i < n; i += 2) 
		{
			lbfgsfloatval_t t1 = 5.0 - x[i];

			lbfgsfloatval_t t2 = 10.0 * (x[i+1] - x[i] * x[i]);
        
			g[i+1] = 20.0 * t2;
        
			g[i] = -2.0 * (x[i] * g[i+1] + t1);
        
			fx += t1 * t1 + t2 * t2;
		}
		for(i =0; i < n; i++)
		{
			printf("g[%d] = %10.5f %10.5f\n", i, g[i], x[i]);
		}

		printf("FFFX = %10.5f\n", fx);
	}

	//A one-compartment decay model

	if(1)
	{
		 int j, k = 7;             //the number of time points 

		 float time[9];

		 char szLine[2048];

		 FILE *fp;

		/* y_exmp[0] =  -0.0002282408;  y_exmp[1] = 0.0875248723;   y_exmp[2] = 0.2133156562;  y_exmp[3] = 0.3816317807;
		 
		 y_exmp[4] = 0.6172371875;  y_exmp[5] = 0.8528870704;  y_exmp[6] = 0.9789444645;  y_exmp[7 ] = 0.9963421536;   y_exmp[8] = 0.9992965138;*/

		 time[0] = 0.00;  time[1] = 0.38; time[2] = 1.00; time[3] = 2.00; time[4] = 4.00; time[5] = 8.00; time[6] = 16.00;  time[7] = 24.00; time[8] = 32.00;
		 
		 if(0)
		 {
			 fp = fopen("inten.txt", "r");

			 if(NULL == fp)
			 {
				 printf("Cannot read the inputs. Exiting ...\n");

				 exit(1);
			 }

			 j = 0;

			 while(fgets(szLine, sizeof(szLine), fp) != NULL)
			 {
				 y_exmp[j] = atof(szLine);

				 //printf("Valuse y[%d]= %15.9f %s\n", j, y_exmp[j], szLine);

				 j++;
			 }

			 fclose(fp);
		 }

		 fx = 0;

		 if(1)     
		 {
			 if(n == 1)
			 //for (i = 0; i < n; i += 2)         //call for one-compartment model
			 {
				 lbfgsfloatval_t t1 = 0;

				 lbfgsfloatval_t t2 = 0; 
				 
				 One_Compartment(y_exmp, time, 9, x[0], &t1,  &t2);

				 //g[i+1] = 0.;

				 fx = t1; g[0] = t2;
			 }
			 else if(2 == n)     // two-paramer, two-compartment call
			 {
				  lbfgsfloatval_t t1 = 0, t2 = 0, t3 = 0.;

				  Two_Compartment(y_exmp, time, 9, x[0], x[1], 
					&t1, &t2, &t3);

				  fx = t1; g[0] = t2; g[1] = t3;
			 }

			  
		 }
		 else
		 {
			 for (i = 0; i < n; i += 2) 
			 {
				 lbfgsfloatval_t t1 = 0;

				 lbfgsfloatval_t t2 = 0; 

				 g[i] = 0.0;

				 for(j = 0; j < 9; j++)
				 {
					 t1 = 1.0 - exp(-x[i] * time[j]);

					 t2 = t2 + (y_exmp[j] - t1) * (y_exmp[j] - t1);

					 g[i] = g[i] - 2.0 * (y_exmp[j] - t1) * time[j] * exp(-x[i] * time[j]);

					 //g[i] = g[i] - 2.0* time[j] * exp(-x[i] * time[j]) / (y_exmp[j] - t1);x

					 //printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g[i], x[i],  exp(-x[i] * time[j]));
				 }

				 g[i+1] = 0.;

				 fx = fx + t2;

				 //printf("g[%d] = %10.5f %10.5f\n", i, g[i], x[i]);
			 }
		 }
	}

	return fx; 
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    /*printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n"); */  
    return 0;
}

//#define N   100

#define N   2       //determines the number of parameters to fit to the optimize the objective function

int main(int argc, char *argv[])
{
    int i, j, k, l, ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(N);
    lbfgs_parameter_t param;
	float Time[1024], Incorporations[1024];
	FILE *fp, *fp2;
	char szLine[1024], szTemp[1024];


	for(i = 0; i < 1024; i++)
	{
		Time[i] = Incorporations[i] = 0.0;
	}

	fp = fopen("Liver_Incorporations.txt", "r");

	fp2 = fopen("Rate_2Compart", "w");

	if(NULL == fp)
	{
		printf("Cannot read Guan_Liver_RIF.csv. Exiting ...\n");

		exit (1);
	}

	j = k = l = 0;

	while(fgets(szLine, sizeof(szLine), fp) != NULL)
	{
		j = 0; k = 0;

		for(i = 0; i < strlen(szLine); i++)
		{
			if(szLine[i] == ',')
			{
				szTemp[k] = '\0';

				Incorporations[j] = atof(szTemp);

				//printf("%10.5f, j = %d %s\n", Incorporations[j], j, szTemp);

				j++;

				szTemp[0] = '\0';

				k = 0;

			}
			else
			{
				szTemp[k] = szLine[i];

				k++;
			}
		}

		Incorporations[j] = atof(szTemp);

		//printf("%10.5f\n", Incorporations[j]);
		
		//printf("%s", szLine);

		for(i = 0; i < 9; i++)
			y_exmp[i] = Incorporations[i];

		/* Initialize the variables. */
		for (i = 0; i < N; i += 2) {
			x[i] = -1.2;

			x[i] = 0.1;

			if(N > 1)
			x[i+1] = 1.0;
		}

		/* Initialize the parameters for the L-BFGS optimization. */
		lbfgs_parameter_init(&param);
		/*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/

		/*
			Start the L-BFGS optimization; this will invoke the callback functions
			evaluate() and progress() when necessary.
			*/
		ret = lbfgs(N, x, &fx, evaluate, progress, NULL, &param);

		/* Report the result. */
		//printf("L-BFGS optimization terminated with status code = %d\n", ret);
		//printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);

		if(x[0] < x[1]) 
			fprintf(fp2, "%15.6f   %15.6f  %d\n", x[0], x[1], ret);
		else
			fprintf(fp2, "%15.6f   %15.6f  %d\n", x[1], x[0], ret);

		l++;
	}

	fclose(fp);

    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");

        return 1;
    }

	return 0;

	/* Initialize the variables. */
	for (i = 0; i < N; i += 2) {
		x[i] = -1.2;

		x[i] = 0.1;

		if(N > 1)
		x[i+1] = 1.0;
	}

	/* Initialize the parameters for the L-BFGS optimization. */
	lbfgs_parameter_init(&param);
	/*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/

	/*
		Start the L-BFGS optimization; this will invoke the callback functions
		evaluate() and progress() when necessary.
		*/
	ret = lbfgs(N, x, &fx, evaluate, progress, NULL, &param);

	/* Report the result. */
	printf("L-BFGS optimization terminated with status code = %d\n", ret);
	printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);


	lbfgs_free(x);

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
	(1 - exp(-param*t)) and experimental, Incorportations, data.
*/
int One_Compartment(float *Incorporations, float *Time, int nData, float param, 
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g)
{
	int j;

    lbfgsfloatval_t fx_current, g_current = 0.0;

	lbfgsfloatval_t t1 = 0;

	lbfgsfloatval_t t2 = 0; 

	for(j = 0; j < nData; j++)
	{
		t1 = 1.0 - exp(-param * Time[j]);

		t2 = t2 + (Incorporations[j] - t1) * (Incorporations[j] - t1);

		g_current = g_current - 2.0 * (Incorporations[j] - t1) * Time[j] * exp(-param * Time[j]);

		//printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g_current, param,  exp(-param * Time[j]));
	}

	fx_current = t2;

	*fx = fx_current; *g = g_current;

	//printf("g = %10.5f %10.5f\n", (*g), param);

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
int Two_Compartment(float *Incorporations, float *Time, int nData, float param1, float param2,
	lbfgsfloatval_t *fx, lbfgsfloatval_t *g1, lbfgsfloatval_t *g2)
{
	int j;

    lbfgsfloatval_t fx_current, g_current1 = 0.0, g_current2 = 0.0;

	lbfgsfloatval_t t1 = 0;

	lbfgsfloatval_t t2 = 0; 

	for(j = 0; j < nData; j++)
	{
		t1 = 1.0 - (param2 * exp(-Time[j]*param1)- param1 * exp(-Time[j]*param2))/
		                  (param2 - param1);

		t1 = Incorporations[j] - t1;

		t2 = t2 + t1 * t1;

		g_current1 = g_current1 - 2 * t1 *  ( (Time[j]*(param2 - param1) - 1) * exp(-param1 * Time[j]) + exp(-param2 * Time[j]) )*param2/
												((param2 - param1) * (param2 - param1));

		g_current2 = g_current2 - 2 * t1 * (exp(-param1 * Time[j]) - (1 + (param2 - param1) * Time[j])*exp(-param2 * Time[j]))*param1/
												((param2 - param1)*(param2 - param1));

		//printf("gi = %10.5f, xi = %15.9f  %10.5f\n", g_current1, param1,  exp(-param1 * Time[j]));
	}

	fx_current = t2;

	*fx = fx_current; *g1 = g_current1; *g2 = g_current2;

	//printf("g = %10.5f %10.5f\n", (*g1), param1);

	return 0;
}
