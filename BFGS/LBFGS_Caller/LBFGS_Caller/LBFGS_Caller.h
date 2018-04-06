// LBFGS_Caller.h

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "lbfgs.h"
#include <math.h>
#include <string.h>
#include <random>

using namespace System;

namespace LBFGS_Library_Call 
{

	
	public ref class LBFGS
	{
		//// TODO: Add your methods for this class here.
		public:

			int nData, nParams;

			float *fTime, *fParams;

			String ^ sMethod;

			LBFGS(float *fTime_exprnt, int nData_exprnt, int nParams_exprnt, String ^ sMethod_exprnt)
			{
				nData = nData_exprnt;

				sMethod = sMethod_exprnt;

				if(sMethod_exprnt->IndexOf("One_Compartment") > -1)
				{
					nParams = 1;
				}
				else if(sMethod_exprnt->IndexOf("Two_Compartment") > -1)
				{
					nParams = 2;
				}
				else if (sMethod_exprnt->IndexOf("Three_Compartment") > -1)
				{
					nParams = 3;
				}
				else 
				{
					printf("Unsupported Method for LBFGS ...\n");

					exit(1);
				}

				//nParams = nParams_exprnt;

				fTime = (float *)calloc(nData, sizeof(float));

				fParams = (float *)calloc(nData, sizeof(float));


				if(NULL == fTime || NULL == fParams)
				{
					printf("Cannot Allocate Memotry for Params\n");

					exit (1);
				}

				for(int i = 0; i < nData; i++)
				{
					fTime[i] = fTime_exprnt[i];
				}
			}

			int Optimize(float * fResponse,  float a1, float a2,  float * rkd, float *rks, float *fx1);

			void InitializeTime();

			void Release_Memory();

			//void RangedRandDemo( float range_min, float range_max, int n);
			//float RangedRandDemo( float range_min, float range_max, int n, float *rkd, float *rks );
	};
}