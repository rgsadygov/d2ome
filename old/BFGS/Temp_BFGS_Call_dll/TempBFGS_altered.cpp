#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace System;
using namespace LBFGS_Library_Call;

int WriteBFGSIsotopeIncorporationRates(ExperimentCollection ^allExperiments, 
	String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args)
{

	float *fTime, *fResponse;

	int i = 0, j = 0, k = 0, l = 0, nret;

	int ntimepoints = allExperiments->ExperimentsList->Count;

	fTime = (float *) calloc(ntimepoints, sizeof(float));

	fResponse = (float *) calloc(ntimepoints, sizeof(float));

	if(NULL == fTime || NULL == fResponse)
	{
		printf("Cannot Allocate Memory for fResponse\n");

		exit (1);
	}

	fTime[0] = 0.00;  
	
	fTime[1] = 4.00; 
	
	if(ntimepoints>2){fTime[2] = 10.00;} 
	
	if(ntimepoints>3){fTime[3] = 24.00;} 
	
	if(ntimepoints>4){fTime[4] = 48.00;} 
	
	if(ntimepoints>5){fTime[5] = 96.00;} 
	
	if(ntimepoints>6){fTime[6] = 168.00;}

	LBFGS ^ lbfgs = gcnew LBFGS(fTime, ntimepoints, 2, "Two_Compartment");

	lbfgs->InitializeTime();   //initialize the time points

	float Time[1024], Incorporations[1024];
	
	FILE *fp, *fp2;
	
	char szLine[1024], szTemp[1024];

	for(i = 0; i < 1024; i++)
	{
		Time[i] = Incorporations[i] = 0.0;
	}

	//fp = fopen("Liver_Incorporations.txt", "r");

	fp2 = fopen("Rate_2Compart.csv", "w");

	/*if(NULL == fp)
	{
		printf("Cannot read Guan_Liver_RIF.csv. Exiting ...\n");

		exit (1);
	}*/

	j = k = l = 0;

	while(fgets(szLine, sizeof(szLine), fp) != NULL)
	{
		j = 0; k = 0;

		for(i = 0; i < strlen(szLine); i++)
		{
			if(szLine[i] == ',')
			{
				szTemp[k] = '\0';

				fResponse[j] = atof(szTemp);

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

		fResponse[j] = atof(szTemp);

		nret = lbfgs->Optimize(fResponse);

		if(lbfgs->fParams[0] < lbfgs->fParams[1]) 
			fprintf(fp2, "%15.6f,   %15.6f,  %d\n", lbfgs->fParams[0], lbfgs->fParams[1], nret);
		else
			fprintf(fp2, "%15.6f,   %15.6f,  %d\n", lbfgs->fParams[1], lbfgs->fParams[0], nret);

	}

	lbfgs->Release_Memory();

	free(fResponse); free(fTime);

    return 0;
}
