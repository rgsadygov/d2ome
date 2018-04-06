// Temp_BFGS_Call_dll.cpp : main project file.

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace System;
using namespace LBFGS_Library_Call;
using namespace std;

int main(array<System::String ^> ^args)
{
   // Console::WriteLine(L"Hello World");
	float *fTime, *fResponse, tmax;

	//Added by JA 2016.10.05
	float rkd, rks, fx1;
	float kdfinal, ksfinal;
	int rrfinal;

	int i = 0, j = 0, k = 0, l = 0, nret, i1=0, ntimepoints=8;
	int imodel;
	int count=0;
	FILE *fp, *fp2;
	//String ^tUnits, ^templine;
	//char tUnits[1024], templine[1024];
	std::string templine, tUnits;
	fTime = (float *) calloc(ntimepoints, sizeof(float));
	fResponse = (float *) calloc(ntimepoints, sizeof(float));



	ifstream file1;//("timecourse.dat");
	file1.open("timecourse.dat");
	if(file1==NULL)
	{
		printf("Can not open time course file.\n");
		exit(-1);
	}

	//while(file1!=NULL)
	while(getline(file1,templine))
	{
		if(i1==0)
			{
				cout<<templine<<"\n";
				tUnits = templine;
				//getline(file1,tUnits);
				//in<<tUnits;
			}
		else
			{	//cout<<templine<<"\n";
				
				if("Day"==tUnits)
				{
					fTime[i1-1] = std::stof(templine);
				}
				else if("Hour"==tUnits)
				{
					fTime[i1-1] = std::stof(templine)/24;
				}
				cout<<fTime[i1-1]<<"\n";
				//cout<<templine<<"\n";
				//getline(file1,fTime[i1-1]);
				//in<<fTime[i1-1];
			}

		i1++;
	}
	file1.close();
	//fclose(fptime);
	ntimepoints = i1-1;
	cout<<ntimepoints<<"\n";




	//Liver
	//

	if(NULL == fTime || NULL == fResponse)
	{
		printf("Cannot Allocate Memory for fResponse\n");
		exit (1);
	}


	//Update Time normalization here
	/*
	tmax = fTime[ntimepoints-1];
	 for(i1=0;i1<ntimepoints;i1++)
	{
		fTime[i1] = fTime[i1]/fTime[ntimepoints-1];
	}
	*/

	imodel=2;

	//LBFGS ^ lbfgs = gcnew LBFGS(fTime, ntimepoints, imodel, "One_Compartment_exponential");

    LBFGS ^ lbfgs = gcnew LBFGS(fTime, ntimepoints, imodel, "Two_Compartment");

	//fp = fopen("Liver_Incorporations.csv", "r");
	//End Liver

	//fp = fopen("Brain_Incorporations_problem2.csv", "r");
	fp = fopen("Brain_Incorporations.csv", "r");
	//fp = fopen("APOB_HUMAN.Quant_NL.csv","r");

	//fp = fopen("Neomed_examples_Augustnew.csv","r");
		//fp = fopen("Neomed_examples_previous.csv","r");

/*
	//T2D3 data
	fTime = (float *) calloc(7, sizeof(float));
	fResponse = (float *) calloc(7, sizeof(float));
	if(NULL == fTime || NULL == fResponse)
	{
		printf("Cannot Allocate Memory for fResponse\n");
		exit (1);
	}

	fTime[0] = 0.00; fTime[1] = 4.00; fTime[2] = 10.00;
	fTime[3] = 24.00; fTime[4] = 48.00; fTime[5] = 96.00; fTime[6] = 168.00;

	
	for(i1=0;i1<7;i1++)
	{
		fTime[i1] = fTime[i1]/fTime[6];
	}
	
	LBFGS ^ lbfgs = gcnew LBFGS(fTime, 7, 1, "One_Compartment");
	fp = fopen("GI4502027_MPE_NetLabeling_TimeProfiles.csv","r");
	//End T2D3





	//fptime = fopen("timecourse.dat","r");

/*ifstream file1;// ("timecourse.dat");
  file1.open("timecourse.dat");
  if (file1.is_open())
  {
    while ( getline (file1,templine) )
    {
      cout << templine << '\n';
    }
    file1.close();
  }
  else cout << "Unable to open file"; 

 // return 0;

 */





	/*
	ofstream file2("timecourse.dat");
	file2<<tUnits<<"\n";
	*/


	lbfgs->InitializeTime();   //initialize the time points
	float Time[1024], Incorporations[1024];
	char szLine[1024], szTemp[1024];

	for(i = 0; i < 1024; i++)
	{
		Time[i] = Incorporations[i] = 0.0;
	}
	
	//srand(time(0));
	//srand( (unsigned)time( NULL ) );
	//rand();
	//srand(( unsigned )time( 0 ) * 1000 );
	

	//srand(( unsigned )time( 0 ) * 100 );

	//Update output filename here
	//fp2 = fopen("Neomed_examples_previous_no_norm_1c_rand1.csv", "w");
	//fp2 = fopen("APOB_HUMAN.Quant_NL_both_norm_3.csv", "w");
    fp2 = fopen("Rate_2Compart_Brain_no_norm_1_pt1_3.csv", "w");
	//fp2 = fopen("Rate_2Compart_Brain_no_norm_1_pt1_2.csv", "w");
	//fp2 = fopen("Rate_2Compart_Brain_notimenorm.csv", "w");
	//fp2 = fopen("Rate_1Compart_Liver.csv", "w");
	//fp2 = fopen("Rate_1Compart_expo.csv", "w");
	//fp2 = fopen("Rate_2Compart.csv", "w");

	if(NULL == fp)
	{
		printf("Cannot read Input file. Exiting ...\n");
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
		//printf("%10.5f\n",szTemp[6]);
		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f",fResponse[0],fResponse[1],fResponse[2],fResponse[3],fResponse[4],fResponse[5],fResponse[6]);
		
		
		
		//Update Response normalization here
		/*
		for(i1=0;i1<ntimepoints;i1++)
		{
			fResponse[i1] = fResponse[i1]/fResponse[ntimepoints-1];
		}
		*/

		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f",fResponse[0],fResponse[1],fResponse[2],fResponse[3],fResponse[4],fResponse[5],fResponse[6]);

		//srand(( unsigned )time( 0 ) * 100 );
		/*
		count=0;
		while(nret!=0 && count<=10)
		{
			nret = lbfgs->Optimize(fResponse,&rkd,&rks);
		    count++;
		}
		*/

		tmax = 1;

		/*
		nret = lbfgs->Optimize(fResponse,&rkd,&rks,&fx1);
				if(lbfgs->fParams[0] <= lbfgs->fParams[1]) 
				{
					kdfinal = lbfgs->fParams[0]/tmax;
					ksfinal = lbfgs->fParams[1]/tmax;
				}
				else
				{
					kdfinal = lbfgs->fParams[1]/tmax;
					ksfinal = lbfgs->fParams[0]/tmax;
				}
				*/


		
		rrfinal = 10;
		for(count=0;count<1;count++)
		{
			nret = lbfgs->Optimize(fResponse,&rkd,&rks,&fx1);
			//printf("%f",lbfgs->progress);
			
			if(Math::Abs(fx1)<Math::Abs(rrfinal))
			{
				rrfinal = fx1;
				if(lbfgs->fParams[0] < lbfgs->fParams[1]) 
				{
					kdfinal = lbfgs->fParams[0];///tmax;
					ksfinal = lbfgs->fParams[1];///tmax;
				}
				else
				{
					kdfinal = lbfgs->fParams[1];///tmax;
					ksfinal = lbfgs->fParams[0];///tmax;
				}
			}
		}
		
		

		//Update Readout here
		if(1==imodel)
		{
			//1 compartment
			//fprintf(fp2, "%15.6f,  %d\n", lbfgs->fParams[0]/tmax, nret);
			//fprintf(fp2, "%15.6f,  %d\n", lbfgs->fParams[0], nret);
			kdfinal = lbfgs->fParams[0]/tmax;
			fprintf(fp2, "%15.6f,  %d\n", kdfinal, nret);

		}
		else if(2==imodel)
		{
			//2 compartment
		/*
		if(lbfgs->fParams[0] < lbfgs->fParams[1]) 
		{
			//fprintf(fp2, "%15.6f,   %15.6f,  %d\n", lbfgs->fParams[0]/tmax, lbfgs->fParams[1]/tmax, nret);
			fprintf(fp2, "%15.6f,   %15.6f,  %d, %f, %f\n", lbfgs->fParams[0], lbfgs->fParams[1], nret, rkd, rks);
		}
		else
		{
			//fprintf(fp2, "%15.6f,   %15.6f,  %d\n", lbfgs->fParams[1]/tmax, lbfgs->fParams[0]/tmax, nret);
			fprintf(fp2, "%15.6f,   %15.6f,  %d, %f, %f\n", lbfgs->fParams[1], lbfgs->fParams[0], nret, rkd, rks);
		}
		*/
			fprintf(fp2, "%15.6f,   %15.6f,  %d, %f, %f\n", kdfinal, ksfinal, nret, rkd, rks);
		}
		else if(3==imodel)
		{
			//1 compartment exponential
			//fprintf(fp2, "%15.6f,  %d\n", (exp(lbfgs->fParams[0]))/tmax, nret);
			kdfinal = (exp(lbfgs->fParams[0]))/tmax;
			fprintf(fp2, "%15.6f,  %d\n", kdfinal, nret);
		}

	}

	lbfgs->Release_Memory();

	free(fResponse); free(fTime);

    return 0;
}
