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

	float *fTime, *fResponse;
	int i = 0, j = 0, k = 0, l = 0, nret, i1=0, ntimepoints;
	FILE *fp, *fp2;
	//String ^tUnits, ^templine;
	//char tUnits[1024], templine[1024];
	std::string templine, tUnits;

	ifstream file1("timecourse.dat");
	//file1.open("timecourse.dat");
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
				fTime[i1-1] = std::stod(templine);
				cout<<fTime[i1-1]<<"\n";
				//getline(file1,fTime[i1-1]);
				//in<<fTime[i1-1];
			}

		i1++;
	}
	file1.close();
	//fclose(fptime);
	ntimepoints = i1;





	//Liver
	//
	fTime = (float *) calloc(9, sizeof(float));
	fResponse = (float *) calloc(9, sizeof(float));
	if(NULL == fTime || NULL == fResponse)
	{
		printf("Cannot Allocate Memory for fResponse\n");

		exit (1);
	}

//	fTime[0] = 0.00;  fTime[1] = 0.38; fTime[2] = 1.00; fTime[3] = 2.00; fTime[4] = 4.00; fTime[5] = 8.00; 
//	fTime[6] = 16.00;  fTime[7] = 24.00; fTime[8] = 32.00;

	/*
		for(i1=0;i1<9;i1++)
	{
		fTime[i1] = fTime[i1]/fTime[8];
	}
		*/

	LBFGS ^ lbfgs = gcnew LBFGS(fTime, 9, 1, "One_Compartment");
	fp = fopen("Liver_Incorporations.csv", "r");
	//End Liver



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
	
	fp2 = fopen("Rate_1Compart_Liver_Test.csv", "w");
	//fp2 = fopen("Rate_1Compart_expo_Liver_notimenorm.csv", "w");
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
		
		for(i1=0;i1<7;i1++)
		{
			fResponse[i1] = fResponse[i1]/fResponse[8];
		}
		
		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f",fResponse[0],fResponse[1],fResponse[2],fResponse[3],fResponse[4],fResponse[5],fResponse[6]);

		nret = lbfgs->Optimize(fResponse);
		//fprintf(fp2, "%15.6f,  %d\n", (exp(lbfgs->fParams[0]))/32, nret);


		fprintf(fp2, "%15.6f,  %d\n", lbfgs->fParams[0]/32, nret);

/*		if(lbfgs->fParams[0] < lbfgs->fParams[1]) 
			fprintf(fp2, "%15.6f,   %15.6f,  %d\n", lbfgs->fParams[0], lbfgs->fParams[1], nret);
		else
			fprintf(fp2, "%15.6f,   %15.6f,  %d\n", lbfgs->fParams[1], lbfgs->fParams[0], nret);
*/
	}

	lbfgs->Release_Memory();

	free(fResponse); free(fTime);

    return 0;
}
