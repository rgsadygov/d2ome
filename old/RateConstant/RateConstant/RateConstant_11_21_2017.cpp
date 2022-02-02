// This is the main DLL file.

#include "stdafx.h"

#include "RateConstant.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

using namespace System;
using namespace System::IO;
using namespace System::Collections;
using namespace System::Collections::Generic;
using namespace RateConstant;
using namespace LBFGS_Library_Call;

using namespace std;


#define MASS_H2OH	 	 19.01784113
#define mProton			 1.00727642
#define mDeuterium		 2.014102
//#define mDeuterium     2.01355    //anohter source
#define mHydrogen	     1.007825035
#define M_PI             3.14159265358979323846

bool CheckForIonScores(int Number);


//double deltaC13 = 13.003354838 - 12.0;

double deltaC13 = 13.00281 - 11.99945;     //another source 13.00281

double dN14 = 14.003074, dN15 = 15.000109;

//double dO16 = 15.994915, dO17 = 16.999132, dO18 = 17.999161;

double dO16 = 15.99437, dO17 = 16.99858, dO18 = 17.99861;

//double dS32 = 31.972071, dS33 = 32.971459, dS34 = 33.967867, dS36 = 35.967081;

double dS32 = 31.97152, dS33 = 32.97091, dS34 = 33.96732, dS36 = 35.96653;    //Prot.Pros.


template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}


//local functions, used for houskeeping activities in the class

float QSelect_Float(List <float> ^lArray);

int compare_float (const void * a, const void * b);

unsigned Factorial(unsigned n);

void realft(array <double> ^data, const int isign);

void four1(array <double> ^data, const int n, const int isign);

void four1(array <double> ^data, const int isign);

double Log_Factorial(unsigned n);

float SingleNumberReader(char *szString, int nStart);

double Binomial(double prob, int nSuccess, int nTrials);

double AminoAcid_Sequence_Mass(String ^sSequence);


float ProteinRateConstant::ComputeRateConstant(String ^ sInputPeptide,
	array <float> ^TheoreticalIsotopes, array <double, 2> ^dExperimentIsotopes)
{
	//determine elemental compostion
	// and the number of exchangeable hydrogen atoms

	ElementalComposition(sInputPeptide);

	//uses approach where the rate constant is determine from
	// the decay of I0.
	//ComputeRateFromI0(TheoreticalIsotopes, dExperimentIsotopes);

	return 0.f;
}

/*
*
* Uses approach where I0's decay is used for rate constant determination
* I0_Natural is the expected value of the I0 isotope at the start of the
* labeling experiment, e.g., when there were no labeling except for the
* natural isotopes.
* I0_Asymptote is the value of the I0 isotope at the asymptote (long-term time) that 
* will be expected, given the body water enrichment and the number of exchangeable
* hydrogen atoms.
* nTimeCourse - the number of time-course data points 
*
*/
float ProteinRateConstant::ComputeRateFromI0(array <float> ^TheoreticalIsotopes, array <float> ^IonIntensities,  
						array <double> ^dLabelIncorporationTime, 
	int nTrueNumberTimePoints,float *fCorr, float *fRMRSS)
{
	float I0_AtAsymptote, I0_Natural, fBodyDueterium, I1_Natural;

	float rkd, rks, fx1,  fCorrelation;

	double fDegradationConstant;

	float fRootMeanRSSS = 0.0;

	int k, i, j, nTimeCourse, nRet;

	float fMaxCorrelation = 0.0, fDegradationConstant_MaxCorr=0.;

	float fDegradationConstant_Temp = 0.f, fRootMeanRSSS_MaxCorr=.0;

	nTimeCourse = nTrueNumberTimePoints;

	array <float> ^ I0_TimeCourse_Fit = gcnew array <float> (nTimeCourse);

	array <float> ^  fNaturalIsotopes = gcnew array <float> (nIsotopes);     //first 6 natural isotopes -computed using Poisson

	I0_AtAsymptote = fBodyDueterium = 0;

	k = fBodyWaterEnrichment->Length;
	
	//check if the experiments has isotopes - it is determined from the fTimeCourse;

	if(0 == k)
	{
		printf("BodyWaterEnrichment Syntax is not Accurate\n");

		exit (1);
	}

	//printf("k = %d\n", k);

	fBodyDueterium = fBodyWaterEnrichment[k - 1];

	if (0.00 == fBodyDueterium)
	{
		printf("O.0 Heavy Water Enrichment\n");

		printf("Will exit here...\n");

		exit(1);
	}

	I0_Natural = pow(1 - 0.011, nC) * pow(1 - 0.000115, nH) * pow(1 - 0.00366, nN);

	I0_Natural = I0_Natural * pow(1-0.00238, nO) * pow(1 - 0.0498, nS);

	/*
	*   Use the Poisson approximation for natural isotope distibutions calculations
	*/

    PoissonIsotopes(nIsotopes, fNaturalIsotopes, nH, nC, nN, nO, nS);

//	I1_Natural = lambda_HCN*exp(-lambda_HCN);

	I0_AtAsymptote = pow(1 - fBodyDueterium, nExHydrogen) * I0_Natural / pow(1 - 0.000115, (nH - nExHydrogen) );  
						//the last term accounts for the fact that there is now less hydrogens with "natural isotope",
						// nExHydrogen - are now exclusively exchangeable via heavy water hydrogens, with enrichment fBodyDueterium;

	float TimeCourseDates[500], TimeCourseI0Isotope[500];

	for(i = 0; i < nTrueNumberTimePoints; i++)
	{
		TimeCourseDates[i] = dLabelIncorporationTime[i];

		TimeCourseI0Isotope[i] = IonIntensities[i];

		//printf("IO_TIMECOURS[%d] %f, dates %f\n", i, TimeCourseI0Isotope[i], TimeCourseDates[i]/fScaleTime);
	}

	printf("nIsot= %d  I0_Natural = %5.3f  TheoAsympt. =  %5.3f ExpAsympt. = %5.3f\n", nIsotopes, I0_Natural, 
	I0_AtAsymptote, IonIntensities[nTrueNumberTimePoints - 1]);

	LBFGS ^ lbfgs = gcnew LBFGS(TimeCourseDates, nTrueNumberTimePoints, 1, "One_Compartment_exponential");

	lbfgs->InitializeTime();

	nRet = lbfgs->Optimize(TimeCourseI0Isotope, I0_AtAsymptote, (I0_Natural - I0_AtAsymptote), &rkd, &rks, &fx1);

	fDegradationConstant = exp(lbfgs->fParams[0]);     //take the exponent, as the paramer is in the exponent

	//printf("DF = %10.5f, fdegradation = %10.5f, fScale = %10.5f\n", lbfgs->fParams[0], fDegradationConstant*fScaleTime, fScaleTime);

	lbfgs->Release_Memory();

	for(i = 0; i < nTrueNumberTimePoints; i++)
	{
		I0_TimeCourse_Fit[i] = I0_AtAsymptote + (I0_Natural - I0_AtAsymptote) * exp( - fDegradationConstant * dLabelIncorporationTime[i]);

		printf("I0_Isotopes %f %f %10.5f\n", TimeCourseDates[i], TimeCourseDates[i]/fScaleTime, TimeCourseI0Isotope[i]);
	}

	if(false)
	for (i = 0; i < nTrueNumberTimePoints; i++)
	{
		if (i == nTrueNumberTimePoints - 1)
		{
			printf("%f\n", TimeCourseDates[i] / fScaleTime);
		}
		else
		printf("%f, ", TimeCourseDates[i] / fScaleTime);
	}

	for (i = 0; i < nTrueNumberTimePoints; i++)
	{
		if (i == nTrueNumberTimePoints - 1)
		{
			printf("%f\n", TimeCourseI0Isotope[i]);
		}
		else
		printf("%f, ",  TimeCourseI0Isotope[i]);
	}

	fCorrelation = PearsonCorrelation(IonIntensities, I0_TimeCourse_Fit, nTrueNumberTimePoints);

	fRootMeanRSSS = RMSS(IonIntensities, I0_TimeCourse_Fit, nTrueNumberTimePoints);

	/*
	*   if fCorrelation is less than 0.9, 
	*   drop one data point and see if it helps
	*   to improve the correlation and reduce the
	*   meanRSSS
	*   Max_Reduced_TimeCourse holds the time points for fit
	*   that produced maximum correlation.
	*/

	//for (i = 0; i < nTrueNumberTimePoints; i++)
		//printf(" %f,", IonIntensities[i]);

	//puts("");

	array <float> ^ Reduced_TimeCourse = gcnew array <float>(nTrueNumberTimePoints - 1);

	array <float> ^ Reduced_TimeCourseI0Isotope = gcnew array <float>(nTrueNumberTimePoints - 1);

	array <float> ^temp_I0_TimeCourse_Fit = gcnew array <float>(Reduced_TimeCourseI0Isotope->Length);

	array <float> ^ Max_Reduced_TimeCourse = gcnew array <float>(nTrueNumberTimePoints - 1);

	//if (fCorrelation > 0.0 && fCorrelation < 0.9)
	if (fCorrelation < 0.95)
	{
		for (i = 0; i < nTrueNumberTimePoints; i++)
		{
			k = 0;

			for (j = 0; j < nTrueNumberTimePoints; j++)
			{
				if (j != i)
				{
					//printf("k = %d, j = %d i = %d\n", k, j, i);

					Reduced_TimeCourse[k] = dLabelIncorporationTime[j];

					Reduced_TimeCourseI0Isotope[k] = IonIntensities[j];

					//printf("TimeC[%d] = %f I0 = %f\n", k, Reduced_TimeCourse[k], Reduced_TimeCourseI0Isotope[k]);

					k = k + 1;

				}

			}

			for (int ii = 0; ii < 100; ii++)
			{
				TimeCourseDates[ii] = 0.0;
			}

			for (int ii = 0; ii < Reduced_TimeCourse->Length; ii++)
			{
				TimeCourseDates[ii] = Reduced_TimeCourse[ii];
			}

			LBFGS ^ reduced_lbfgs = gcnew LBFGS(TimeCourseDates, Reduced_TimeCourse->Length, 1, "One_Compartment_exponential");

			reduced_lbfgs->InitializeTime();

			for (int ii = 0; ii < 100; ii++)
			{
				TimeCourseI0Isotope[ii] = 0.0;
			}

			for (int ii = 0; ii < Reduced_TimeCourseI0Isotope->Length; ii++)
			{
				TimeCourseI0Isotope[ii] = Reduced_TimeCourseI0Isotope[ii];
			}

			//for (int ii = 0; ii < Reduced_TimeCourseI0Isotope->Length; ii++)
				//printf("New values %f %f\n", Reduced_TimeCourse[ii], Reduced_TimeCourseI0Isotope[ii]);

			nRet = reduced_lbfgs->Optimize(TimeCourseI0Isotope, I0_AtAsymptote, (I0_Natural - I0_AtAsymptote), &rkd, &rks, &fx1);

			fDegradationConstant_Temp = exp(reduced_lbfgs->fParams[0]);     //take the exponent, as the paramer is in the exponent

																	   //printf("DF = %10.5f\n", lbfgs->fParams[0]); 

			reduced_lbfgs->Release_Memory();


			for (int ii = 0; ii < Reduced_TimeCourseI0Isotope->Length; ii++)
			{
				temp_I0_TimeCourse_Fit[ii] = I0_AtAsymptote + (I0_Natural - I0_AtAsymptote) * exp(-fDegradationConstant * Reduced_TimeCourse[ii]);

				//printf("%f %10.5f\n", TimeCourseDates[ii]/fScaleTime, TimeCourseI0Isotope[ii]);
			}

			float ftemp;

			ftemp = PearsonCorrelation(Reduced_TimeCourseI0Isotope, temp_I0_TimeCourse_Fit, Reduced_TimeCourseI0Isotope->Length);

			if (fMaxCorrelation < ftemp)
			{
				fMaxCorrelation = ftemp;

				Max_Reduced_TimeCourse = Reduced_TimeCourse;

				fDegradationConstant_MaxCorr = fDegradationConstant_Temp;

				fRootMeanRSSS_MaxCorr = RMSS(Reduced_TimeCourseI0Isotope, temp_I0_TimeCourse_Fit, Reduced_TimeCourseI0Isotope->Length);
				
			}

			delete reduced_lbfgs;

			//puts("");
		}
	}
	delete Reduced_TimeCourse, Reduced_TimeCourseI0Isotope, temp_I0_TimeCourse_Fit;

	int nTimePoinsUsed = dLabelIncorporationTime->Length;

	//printf("Fcorr %f MaxFcorr %f\n", fCorrelation, fMaxCorrelation);
	if (fMaxCorrelation > fCorrelation && fCorrelation > -1.0 )
	{
		//printf("Rate Constant Before: %f After = %f\n", fDegradationConstant* fScaleTime, fDegradationConstant_MaxCorr* fScaleTime);
		fCorrelation = fMaxCorrelation;

		fDegradationConstant = fDegradationConstant_MaxCorr;

		fRootMeanRSSS = fRootMeanRSSS_MaxCorr;

		nTimePoinsUsed = Max_Reduced_TimeCourse->Length;

	}

	//printf("RateConstant = %f, Correlation = %f, RMSS = %10.5f, nTimePoints Chosen = %d\n", fDegradationConstant * fScaleTime, fCorrelation, fRootMeanRSSS,
		//nTimePoinsUsed);

	(*fCorr) = fCorrelation; 

	(*fRMRSS) = fRootMeanRSSS;

	delete I0_TimeCourse_Fit, Max_Reduced_TimeCourse;

	return fDegradationConstant; 
}
/* 
*
*  Read the files in the current directory
*  find a file with quant.csv extension and 
*  compute Rate Constants for that protein
*  Keep the file name, just change the .quant.
*  to .Rate.
*
*/
void ProteinRateConstant::ReadFileFolder()
{
	int i,  iQuantFiles;

	cli::array<String ^,1> ^ QuantFiles, ^ FileList; 

	Directory::SetCurrentDirectory(workDirectory);

	FileList = Directory::GetFiles(workDirectory);	

	iQuantFiles = 0;

	//determine if the experiment contains replicates
	//
	printf("Starting Rate Constant Calculations\n");

	bReplicates = ProcessTimePoints(fTimeCourse);

	for(i=0; i < FileList->Length; i++)
	{
		if(FileList[i]->EndsWith(".Quant.csv"))
		{
			iQuantFiles++;
		}
	}

	QuantFiles = gcnew cli::array<String ^,1> (iQuantFiles);

	iQuantFiles = 0;

	for(i=0; i < FileList->Length; i++)
	{
		if(FileList[i]->EndsWith(".Quant.csv"))
		{
			QuantFiles[iQuantFiles] = FileList[i];

			iQuantFiles++;
		}
	}

	for(i=0; i < QuantFiles->Length; i++)
	{
		printf("FILE  %s\n", QuantFiles[i]);

		QuantFileReader(QuantFiles[i]);
	}
}

/*
*
*
*  A function to read and filter the conten of the Quant.csv files:
*  
*  It sends the data to ComputeRatesFromI0 function to compute the 
*  rate constants from the time course incorporation data
*
*/
void ProteinRateConstant::QuantFileReader(String ^ sQuantFile)
{
	int k, j, l, icomma, il, iP, itempComma, iPreviousComma, iU, nTrueIons;

	char szFile[2048], szRateFile[2048], szLine[18096];

	char szPeptide[2048], szUniquePepitde[1024];

	array <float> ^ fCurrentIsotopes = gcnew array <float> (6);

	array <float> ^fTheoreticalIsotopes = gcnew array <float> (6);

	array <double> ^aTimePointsForThisPeptde;

	array <float> ^fI0Intensities;

	RateConstResults ^cCurrentResult = gcnew RateConstResults();

	List <RateConstResults ^> ^AllResults = gcnew List <RateConstResults ^>;

	List <float> ^ fIonScores; 

	float fRateConstant, fCorr, fMeanRate, fSD, fRMRSS, fI0_IsotpeAccuracy;

	String ^ sRateConstFile;

	FILE *fp, *fpRateConst;

	szFile[0] = '\0';

	for(k = 0; k < sQuantFile->Length; k++)
	{
		szFile[k] = sQuantFile[k];
	}

	szFile[k] = '\0';

	j = sQuantFile->IndexOf(".Quant.csv");

	sRateConstFile = sQuantFile->Substring(0, j) + ".RateConst.csv";

	for(k = 0; k < sRateConstFile->Length; k++)
	{
		szRateFile[k] = sRateConstFile[k];
	}

	szRateFile[k] = '\0';

	fp = fopen(szFile, "r");

	if(NULL == fp)
	{
		printf("Cannot open the files %s for Reading Quant Values\n", szFile);

		exit(1);
	}

	fpRateConst = fopen(szRateFile, "w");

	if(NULL == fpRateConst)
	{
		printf("Cannot open the files %s for Reading Quant Values\n", szRateFile);

		exit(1);
	}

	fprintf(fpRateConst, "Peptides, PeptideUnique, RateConstants, Correlations, RootMeanRSS, AbsoluteIsotopeError\n");

	l = 0;

	//printf("FILE: %s\n", sQuantFile);

	szPeptide[0] = '\0';

	iPreviousComma = 0;

	while(fgets(szLine, sizeof(szLine), fp) != NULL)
	{

		// a new peptide and new time course isotopes

		SingleIsotopeCluster ^ OneIsotopeClusters = gcnew SingleIsotopeCluster();

		List <SingleIsotopeCluster ^> ^TimeCourseIsotopeCluster = gcnew List <SingleIsotopeCluster ^>;

		fIonScores = gcnew List <float>;

			//List <float> ^Intensity = gcnew List <float>

		if( 4 <= l)     //the first four lines in the quant.csv are not peptide information
		{
			iP = iU = 0;

			icomma = 0;   szUniquePepitde[0] = '\0';

			for(il = 0; il < strlen(szLine); il++)
			{
				//printf("%c ", szLine[il]);

				if(szLine[il] == ',')
				{
					icomma++;


					//printf("commad = %d %d\n", icomma, il);

					if(1 == icomma || CheckANumber(icomma))
					{
						itempComma = il;
					}
				}
				
				if(0 == icomma)    //peptide sequence comes before the first comma
				{
					szPeptide[iP] = szLine[il];

					iP++;
				}
				else if(1 == icomma && itempComma == il)
				{
					szPeptide[iP] = '\0';

					ElementalComposition(gcnew String(szPeptide) );     //assign the elemental composition of the peptide

					printf("Peptide: %s\n", szPeptide);
				}

				if (1 == icomma && szLine[il] != ',')       //this block assumes that in ",Yes ," once you hit the ",", until you hit the next
				{						// "," icomma will be 1. and you can copy the word ("Yes");
					szUniquePepitde[iU] = szLine[il];

					iU++;

					szUniquePepitde[iU] = '\0';
				}

				if (CheckANumber(icomma) && itempComma == il)
				{
					fCurrentIsotopes = NumberReader(szLine, il + 1);

					if (icomma < 16)
					{
						for (int ii = 0; ii < 6; ii++)
						{
							fTheoreticalIsotopes[ii] = fCurrentIsotopes[ii];
						}
					}

					delete OneIsotopeClusters;

					OneIsotopeClusters = gcnew SingleIsotopeCluster();

					if (icomma >= 16)
					{
						for (int ii = 0; ii < 6; ii++)
						{
							OneIsotopeClusters->fIsotopeCluster[ii] = fCurrentIsotopes[ii];
						}

						TimeCourseIsotopeCluster->Add(OneIsotopeClusters);

						if(fCurrentIsotopes[1] < 10 && fIonScores[fIonScores->Count - 1] > 10)
						{
							printf("This Does not Look Isotope Pattern\n");
						}
					}

					/*for(int ii=0; ii < 6; ii++)
					{
						printf("%10.1f ", fCurrentIsotopes[ii]);
					}

					puts("");*/
				}

				if (CheckForIonScores(icomma) && icomma != iPreviousComma)
				{
					fIonScores->Add(SingleNumberReader(szLine, il + 1));

					iPreviousComma = icomma;

					//printf("Ions Score = %f il = %d, icomma= %d\n", fIonScores[fIonScores->Count - 1], il, icomma);
				}
			}

			/*puts("WWWWWWWWWWWWWW");
			for(k = 0; k < TimeCourseIsotopeCluster->Count; k++)
			{
				for(j = 0; j < 6; j++)
				{
					printf("%10.1f ", TimeCourseIsotopeCluster[k]->fIsotopeCluster[j]);
				}

				puts("");
			}*/

			//printf("Lenght here = %d\n", TimeCourseIsotopeCluster->Count);

			// if there are replicates then process the TimeCourseIsotopes

			int iPeptide = 0;

			aTimePointsForThisPeptde = gcnew array <double>(TimeCourseIsotopeCluster->Count);

			fI0Intensities = gcnew array <float> (TimeCourseIsotopeCluster->Count);

			if (bReplicates)
			{
				iPeptide = AverageReplicates(TimeCourseIsotopeCluster, fIonScores, fTheoreticalIsotopes, 
					aTimePointsForThisPeptde, fI0Intensities, &nTrueIons);

				//for (int ii = 0; ii < nTrueIons; ii++)
					//printf("Inside the QuantRaad %f %f\n", aTimePointsForThisPeptde[ii], fI0Intensities[ii]);
			}
			else   //no replicates
			{
				iPeptide = NoReplicateExperiments(TimeCourseIsotopeCluster, fIonScores, fTheoreticalIsotopes,
					aTimePointsForThisPeptde, fI0Intensities, &nTrueIons);

				//printf("printf nTrueIons = %d  %d\n", nTrueIons, TimeCourseIsotopeCluster->Count);
			}

			if (-1 == iPeptide)
			{ 
				printf("No Analyzable Data for this Peptide\n");

				fRateConstant = -1;

				//continue;
			}
			else
			{
				//fRateConstant = ComputeRateFromI0(fTheoreticalIsotopes, TimeCourseIsotopeCluster, aTimePointsForThisPeptde,
					//&fCorr, &fRMRSS, &fI0_IsotpeAccuracy);

				//for (int ii = 0; ii < nTrueIons; ii++)
					//printf("Passing to COmpute[%d] = %f\n", ii, fI0Intensities[ii]);


				//printf("NUMBEROFIONS = %d\n", nTrueIons);

				fRateConstant = ComputeRateFromI0(fTheoreticalIsotopes, fI0Intensities, aTimePointsForThisPeptde, nTrueIons,
					&fCorr, &fRMRSS);
			}
			
			fI0_IsotpeAccuracy = ExperimentalTheoreticalIsotopeAccuracy(TimeCourseIsotopeCluster, fTheoreticalIsotopes);  

			delete OneIsotopeClusters, TimeCourseIsotopeCluster;


			//printf("Rate Constant = %f, Corr = %f\n", fRateConstant, fCorr);

			cCurrentResult = gcnew RateConstResults();

			cCurrentResult->sPeptide = gcnew String(szPeptide);

			cCurrentResult->sIsPeptideUnique = gcnew String(szUniquePepitde);

			cCurrentResult->fRateConst = fRateConstant;      //Dividing by ScaleTime, because the times were scaled.

			cCurrentResult->fCorrel = fCorr;

			cCurrentResult->fRMRSS = fRMRSS;

			cCurrentResult->fI0_IsotpeAccuracy = fI0_IsotpeAccuracy;

			AllResults->Add(cCurrentResult);

		}	

		l++;

		//delete OneIsotopeClusters, TimeCourseIsotopeCluster, cCurrentResult;

		delete cCurrentResult;

		delete fIonScores;
	}

	fclose(fp);

	

	List <float> ^ftempArray = gcnew List <float>;

	for(k = 0; k < AllResults->Count; k++)
	{
		if (AllResults[k]->fRateConst > -0.001) 
		{
			fprintf(fpRateConst, "%s, %s, %6.5f, %f, %f, %f\n", AllResults[k]->sPeptide, AllResults[k]->sIsPeptideUnique,
				AllResults[k]->fRateConst * fScaleTime, AllResults[k]->fCorrel, AllResults[k]->fRMRSS, AllResults[k]->fI0_IsotpeAccuracy);
		}
		
	}

	l = 0;

	fMeanRate = 0.f;

	for (k = 0; k < AllResults->Count; k++)
	{
		if (AllResults[k]->fRateConst > -0.001 && AllResults[k]->fCorrel >= fCorrThreshold &&
			AllResults[k]->fRMRSS < fRMSS_Threshold)
		{
			fMeanRate = fMeanRate + AllResults[k]->fRateConst * fScaleTime;

			ftempArray->Add(AllResults[k]->fRateConst * fScaleTime);

			l++;

		}

	}


	if(0 < l)
	{
		fMeanRate = fMeanRate/(float)l;
	}

	float fMedianRate  = 0.0;

	fMedianRate = QSelect_Float(ftempArray);

	fSD = 0.;

	for(k = 0; k < AllResults->Count; k++)
	{
		if(AllResults[k]->fRateConst > -0.001 && AllResults[k]->fCorrel >= fCorrThreshold &&
			AllResults[k]->fRMRSS < fRMSS_Threshold)
		{
			fSD = fSD + (fMeanRate - AllResults[k]->fRateConst * fScaleTime) * (fMeanRate - AllResults[k]->fRateConst * fScaleTime);
		}
	}

	if (l > 1)
	{
		fSD = fSD / (float)(l - 1);
	}
	
	//remove the outliers, the data points that are 3*sd away, and recalculated the Mean and SD again.

	double final_Mean, final_SD;

	j = 0;

	final_Mean = final_SD = 0.0;

	for(k = 0; k < AllResults->Count; k++)
	{
		if (AllResults[k]->fRateConst > -0.001 && AllResults[k]->fCorrel > fCorrThreshold &&
			AllResults[k]->fRMRSS < fRMSS_Threshold    &&
			Math::Abs(fMeanRate - AllResults[k]->fRateConst * fScaleTime) < 3.* sqrt(fSD))
		{
			final_Mean = final_Mean +  AllResults[k]->fRateConst * fScaleTime;

			j = j + 1;
		}

	}

	if (j <= 1)
	{
		j = l;

		final_Mean = fMeanRate;
	}
	else
	{
		final_Mean = final_Mean / (double)j;
	}
	
	final_SD = 0.0;
	
	if (j > 1)
	{
		ftempArray = gcnew List <float>;

		for (k = 0; k < AllResults->Count; k++)
		{
			if (AllResults[k]->fRateConst > -0.001 && AllResults[k]->fCorrel > fCorrThreshold && 
				AllResults[k]->fRMRSS < fRMSS_Threshold &&
				Math::Abs(fMeanRate - AllResults[k]->fRateConst * fScaleTime) < 3.* sqrt(fSD))
			{
				final_SD = final_SD + (final_Mean - AllResults[k]->fRateConst * fScaleTime) * (final_Mean - AllResults[k]->fRateConst * fScaleTime);

				ftempArray->Add(AllResults[k]->fRateConst * fScaleTime);

				//printf("%f %f  %f\n", AllResults[k]->fRateConst * fScaleTime, fMeanRate - AllResults[k]->fRateConst * fScaleTime, sqrt(fSD));
			}
		}

		final_SD = final_SD / (double)(j - 1);

		fMedianRate = QSelect_Float(ftempArray);
	}

	

	//printf("RawSD = %f, FinalSD %f, l = %d, j = %d\n", fSD, final_SD, l, j);



	if(1 < l)
	{	
		fprintf(fpRateConst, "\n\nMeanRateConst/CorrCutOff, %f, %f\n", final_Mean, fCorrThreshold);

		fprintf(fpRateConst, "\n\nMedianRateConst/RMSSCutOff, %f, %f\n", fMedianRate, fRMSS_Threshold);

		if (j > 1)
		{
			fprintf(fpRateConst, "StandDev/NumberPeptides, %f, %d\n", sqrt(final_SD), j);
		}
		else
		{
			fprintf(fpRateConst, "StandDev/NumberPeptides, %f, %d\n", sqrt(fSD), l);
		}
		
	}
	else
	{

	}


	delete ftempArray;

	fclose(fpRateConst);

	return;
}

/*
*  a simplifying check for numbers
*  It is used to determine the position of 
*  isotopes in the quant file.
*/

bool ProteinRateConstant::CheckANumber(int Number)
{
	bool bFound = false;

	if (5 == Number)
	{
		bFound = true;
	}
	else if ((Number - 1) % 16 == 0)
	{
		bFound = true;
	}

	return bFound;
}



/*
*  a simplifying check for numbers
*  It is used to determine the position of
*  IonScores in the Quant.csv file.
*/

bool CheckForIonScores(int Number)
{
	bool bFound = false;

	//if (Number == 12)
	if (Number == 13)
	{
		bFound = true;
	}

	if ((Number - 13) % 16 == 0)   //determines the remainder after devision by 12
	{
		bFound = true;
	}


	return bFound;
}
/*
* a small function to read the number between two commas
* from a line
*/
array <float> ^ ProteinRateConstant::NumberReader(char *szString, int nStart)
{
	int k, i, l;

	char szTemp[2048];

	array <float> ^ fIsotopes = gcnew array <float> (6);

	//printf("The number to Read %s\n", &szString[nStart]);

	i = l = 0;

	for(k = nStart; k < strlen(szString); k++)
	{
		if(szString[k] == ',')
		{
			szTemp[i] = '\0';

			//printf("TEMP = %s, l = %d\n", szTemp, l);

			if(i > 0)
			{
				fIsotopes[l] = atof(szTemp);
			}
			else
				fIsotopes[l] = 0.f;
			
			szTemp[0] = '\0';

			i = 0;
			
			if(5 == l)
			{
				break;     //overall there are 6 isotopes
			}

			l++;
		}
		else
		{
			szTemp[i] = szString[k];

			i++;
		}
	}

	return fIsotopes;
}

/*
* a small function to read the number between two commas
* from a line
*/
float SingleNumberReader(char *szString, int nStart)
{
	int k, i;

	char szTemp[2048];

	float ftemp = 0.0;

	i = 0;

	szTemp[0] = '\0';

	for (k = nStart; k < strlen(szString); k++)
	{
		if (szString[k] == ',')
		{
			szTemp[i] = '\0';

			break;
		}
		else
		{
			szTemp[i] = szString[k];

			i++;
		}
	}

	if (i > 0)
	{
		ftemp = atof(szTemp);
	}
	else
		ftemp = 0.0;
	

	return ftemp;
}
/*
*  this function reads a sequence of amino acids (a peptide)
*  assign the number of C, H, N, S, P atoms to class variables
*  nC, nH, nN, nO, nS, nP - which are public in class 
*/
int ProteinRateConstant::ElementalComposition(String ^sSequence)
{
	int iiC, iiN, iiH, iiO, iiS, iiP;

	int iLength, i, j; 
	
	char cTemp, cTemp2;

	float nHAA[256], ftemp; 

	int ElemenMatrix[256][5];

	char szSequence[1024];

	
	//assigns the szPeptide here

	szPeptide = sSequence;


	szSequence[0] = '\0';



	for(i = 0; i < sSequence->Length; i++)
	{
		szSequence[i] = sSequence[i];
	}

	szSequence[i] = '\0';


	for(i = 0; i < 256; i++)
	{
		for(j= 0; j < 5; j++)
		{
			ElemenMatrix[i][j] = 0;
		}
	}

	/* inittiate #1 is C, #2 is H, #3 is N, #4 is O, #5 is S*/
	ElemenMatrix[(int)'A'][0] = 3; ElemenMatrix['A'][1] = 5; 
	
	ElemenMatrix['A'][2] = 1; ElemenMatrix['A'][3] = 1;
	
	/*Glycine */
	ElemenMatrix['G'][0] = 2; ElemenMatrix['G'][1] = 3; 
	
	ElemenMatrix['G'][2] = 1; ElemenMatrix['G'][3] = 1;

	/* Serine */
	ElemenMatrix['S'][0] = 3; ElemenMatrix['S'][1] = 5; 
	
	ElemenMatrix['S'][2] = 1; ElemenMatrix['S'][3] = 2;

	/*Proline */
	ElemenMatrix['P'][0] = 5; ElemenMatrix['P'][1] = 7; 
	
	ElemenMatrix['P'][2] = 1; ElemenMatrix['P'][3] = 1;

	/* Valine */
	ElemenMatrix['V'][0] = 5; ElemenMatrix['V'][1] = 9; 
	
	ElemenMatrix['V'][2] = 1; ElemenMatrix['V'][3] = 1;


	/*Threonine */
	ElemenMatrix['T'][0] = 4; ElemenMatrix['T'][1] = 7; 
	
	ElemenMatrix['T'][2] = 1; ElemenMatrix['T'][3] = 2;

	/*Cystein */
	ElemenMatrix['C'][0] = 3; ElemenMatrix['C'][1] = 5; 
	
	ElemenMatrix['C'][2] = 1; ElemenMatrix['C'][3] = 1;  
	
	ElemenMatrix['C'][4] = 1;

	/* Leucine */
	ElemenMatrix['L'][0] = 6; ElemenMatrix['L'][1] = 11; 
	
	ElemenMatrix['L'][2] = 1; ElemenMatrix['L'][3] = 1;

	/* IsoLeucine */
	ElemenMatrix['I'][0] = 6; ElemenMatrix['I'][1] = 11; 
	
	ElemenMatrix['I'][2] = 1; ElemenMatrix['I'][3] = 1;

	/* Asparagine */
	ElemenMatrix['N'][0] = 4; ElemenMatrix['N'][1] = 6; 
	
	ElemenMatrix['N'][2] = 2; ElemenMatrix['N'][3] = 2;

	/* Aspartic Acid */

	ElemenMatrix['D'][0] = 4; ElemenMatrix['D'][1] = 5; 
	
	ElemenMatrix['D'][2] = 1; ElemenMatrix['D'][3] = 3;

	/* Glutamine */

	ElemenMatrix['Q'][0] = 5; ElemenMatrix['Q'][1] = 8; 
	
	ElemenMatrix['Q'][2] = 2; ElemenMatrix['Q'][3] = 2;

	/* Lysine */
	ElemenMatrix['K'][0] = 6; ElemenMatrix['K'][1] = 12; 
	
	ElemenMatrix['K'][2] = 2; ElemenMatrix['K'][3] = 1;


	/* Glutamic Acid */
	ElemenMatrix['E'][0] = 5; ElemenMatrix['E'][1] = 7; 
	
	ElemenMatrix['E'][2] = 1; ElemenMatrix['E'][3] = 3;

	/* Methinine */
	ElemenMatrix['M'][0] = 5; ElemenMatrix['M'][1] = 9; 
	
	ElemenMatrix['M'][2] = 1; ElemenMatrix['M'][3] = 1; 
	
	ElemenMatrix['M'][4] = 1;


	/* Histidine */
	ElemenMatrix['H'][0] = 6; ElemenMatrix['H'][1] = 7; 
	
	ElemenMatrix['H'][2] = 3; ElemenMatrix['H'][3] = 1; 

	/* Phenylalanine - mass is odd because odd number of N */
	ElemenMatrix['F'][0] = 9; ElemenMatrix['F'][1] = 9; 
	
	ElemenMatrix['F'][2] = 1; ElemenMatrix['F'][3] = 1; 

	/* Arginine - mass is even because even number of N */
	ElemenMatrix['R'][0] = 6; ElemenMatrix['R'][1] = 12; 
	
	ElemenMatrix['R'][2] = 4; ElemenMatrix['R'][3] = 1; 

	/* Tyrosine - mass is odd because odd number of N */
	ElemenMatrix['Y'][0] = 9; ElemenMatrix['Y'][1] = 9; 
	
	ElemenMatrix['Y'][2] = 1; ElemenMatrix['Y'][3] = 2;

	/* Tryptophan - mass is even because even number of N */
	ElemenMatrix['W'][0] = 11; ElemenMatrix['W'][1] = 10; 
	
	ElemenMatrix['W'][2] = 2; ElemenMatrix['W'][3] = 1; 

	//for modifications you will need to adjust from the original
	// make the upper and low case amino acid codes to have the same
	// number of H,C, N, O, S
	for(i = 65; i <= 87; i++)
	{
		for(j=0; j <= 4; j++)
		{
			ElemenMatrix[i+32][j] = ElemenMatrix[i][j];
		}
	}

	iiC = iiN = iiH = iiO = iiS = iiP = 0;

	//printf("%d %d %d %d %d %d\n", (int)'A', (int)'W', (int)'a', (int)'w', ElemenMatrix[65][0], ElemenMatrix[97][0]);

	iLength = strlen(szSequence);

	for(i = 0; i < iLength; i++)
	{
		j = (int) szSequence[i];


		if(i == iLength - 1)
		{
			if(szSequence[i] == '[')
			{
				continue;
			}
		}

		if(szSequence[i] == 'm')
		{

		}
		else if(szSequence[i] == 'c')
		{

		}
		else if(szSequence[i] == 'k')
		{

		}
		else if(j < 65 || j > 90)
		{
			printf("Cannot recognize the character %c in %s\n",
				szSequence[i], szSequence);

			puts("Unknown, or unexpected character... will exit");

			exit (1);
		}
		
		iiC += ElemenMatrix[szSequence[i]][0];

		iiH += ElemenMatrix[szSequence[i]][1];

		iiN += ElemenMatrix[szSequence[i]][2];

		iiO += ElemenMatrix[szSequence[i]][3];

		iiS += ElemenMatrix[szSequence[i]][4];

		//look for the modifications
		if(i < iLength - 1)
		{
			cTemp = szSequence[i];

			if(cTemp == 'S' || cTemp == 'T' ||
				cTemp == 'Y' || cTemp == 'C')
			{
				cTemp2 = szSequence[i+1];

				if(cTemp2 == '$' || cTemp2 == '#' ||
					cTemp2 == '*' || cTemp2 == '@')
					/* this site is phosphorylated
					* update the elemental composition
					* correspondingly, and increase the
					* index running over the sequence sites
					*/
				{
					iiO += 3;

					iiH += 1;

					iiP += 1;
					
					i++;
				}
			}
			else if(cTemp == 'M')
			{
				cTemp2 = szSequence[i+1];

				if(cTemp2 == '$' || cTemp2 == '#' ||
					cTemp2 == '*')
					/* this site is phosphorylated
					* update the elemental composition
					* correspondingly, and increase the
					* index running over the sequence sites
					*/
				{
					iiO += 1;
					
					i++;
				}
			}
			else if(cTemp == 'K')
			{
				cTemp2 = szSequence[i+1];

				if(cTemp2 == '$' || cTemp2 == '#' ||
					cTemp2 == '*')
					/* this site is phosphorylated
					* update the elemental composition
					* correspondingly, and increase the
					* index running over the sequence sites
					*/
				{
					iiO += 1;
					
					i++;
				}
			}
		}
	}

	nC = iiC; nH = iiH; nN = iiN; nO = iiO; nS = iiS; 

	nP = iiP;

	/* Correct for N- and C- termini*/

	nO = nO + 1; nH = nH + 3;

	// nubmer of exchangeable hydrogen atoms per AA

	for(int i = 0; i < 256; i++)
	{
		nHAA[i] = 0.0;
	}

	nHAA[(int)'A'] = nHAA[(int)'a']= 4.0f;  nHAA[(int)'C'] = nHAA[(int)'c'] = 1.62f; nHAA[(int)'D'] = nHAA[(int)'d']= 1.89f; nHAA[(int)'E'] = nHAA[(int)'e']= 3.95f; 

	nHAA[(int)'F'] = nHAA[(int)'f'] = 0.32f; nHAA[(int)'G'] = nHAA[(int)'g'] = 2.06f; nHAA[(int)'H'] = nHAA[(int)'h'] = 2.88f; nHAA[(int)'I'] = nHAA[(int)'i'] = 1.0f; 

	nHAA[(int)'L'] = nHAA[(int)'l'] = 0.6f; nHAA[(int)'K'] = nHAA[(int)'k'] = 0.54f; nHAA[(int)'M'] = nHAA[(int)'m'] = 1.12f; nHAA[(int)'N'] = nHAA[(int)'n'] = 1.89f; 

	nHAA[(int)'P'] = nHAA[(int)'p'] = 2.59f; nHAA[(int)'Q'] = nHAA[(int)'q'] = 3.95f; nHAA[(int)'R'] = nHAA[(int)'r'] = 3.43f; nHAA[(int)'S'] = nHAA[(int)'s'] = 2.61f; 

	nHAA[(int)'T'] = nHAA[(int)'t'] = 0.2f; nHAA[(int)'V'] = nHAA[(int)'v'] = 0.56f; nHAA[(int)'W'] = nHAA[(int)'w'] = 0.08f; nHAA[(int)'Y'] = nHAA[(int)'y'] = 0.42f;

	nExHydrogen = 0;

	ftemp = 0.f;

	for(i = 0; i < sSequence->Length; i++)
	{
		ftemp = ftemp + nHAA[sSequence[i]];
	}

	nExHydrogen = (int) ftemp;

	return 0;
}
/*
*
* Compute Resudial Mean Sum of Squares 
*
* Between two float arrays
*
*/

float ProteinRateConstant::RMSS(array <float> ^FirstArr, array <float> ^SecondArr, int nPoints)
{
	float fRMSS;

	int k;

	if (FirstArr->Length < nPoints || nPoints > SecondArr->Length)
	{
		printf("Pearson Correlation only Between Arrays of the same Size (larger than 1)\n");

		exit(1);
	}


	fRMSS = 0.0;

	for(k = 0; k < nPoints; k++)
	{
		fRMSS = fRMSS + (SecondArr[k]  - FirstArr[k]) * (SecondArr[k]  - FirstArr[k]);
	}

	fRMSS = fRMSS/(float)nPoints;

	return sqrt(fRMSS);
}

/*
*
* Compute a Pearson's Correlation 
*
* Between two float arrays
*
*/

float ProteinRateConstant::PearsonCorrelation(array <float> ^FirstArr, array <float> ^SecondArr, int nPoints)
{
	float fMean1, fMean2, ftemp, fsquare1, fsquare2;

	int k;



	fMean1 = fMean2 = ftemp = 0.0;

	for(k = 0; k <  nPoints; k++)
	{
		fMean1 = fMean1 + FirstArr[k];

		fMean2 = fMean2 + SecondArr[k];
	}

	fMean1 = fMean1/(float)nPoints;

	fMean2 = fMean2/(float)nPoints;

	fsquare1 = fsquare2 = 0.0;

	for(k = 0; k <  nPoints; k++)
	{
		fsquare1 = fsquare1 + (fMean1  - FirstArr[k]) * (fMean1  - FirstArr[k]);

		fsquare2 = fsquare2 + (fMean2  - SecondArr[k]) * (fMean2  - SecondArr[k]);

		ftemp = ftemp + (fMean1  - FirstArr[k]) * (fMean2  - SecondArr[k]);
	}

	if(fsquare1 < 0.00000001 || fsquare2 < 0.00000001)
		return -2.0f;

	ftemp = ftemp/sqrt(fsquare1 * fsquare2);

	return ftemp;
}


/*
* A recursive algorithm that returns a median
* of a list of float numbers
*/
float QSelect_Float(List <float> ^lArray)
{
	int i, nSize;

	float *fArray, fmedian;

	nSize = lArray->Count;

	fArray = (float *)malloc(nSize * sizeof(float));

	if(NULL == fArray)
	{
		printf("Not enough space for nArray\n");

		exit (1);
	}

	if(lArray->Count > 10000)
	{
		nSize = 10000;

		printf("Warning longer array Size is Sorting\n, Had to Cut the Size\n");
	}
	else
	{
		nSize = lArray->Count;
	}

	for(i=0; i < nSize; i++)
	{
		fArray[i] = lArray[i];
	}

	qsort(fArray, nSize, sizeof(float), compare_float);

	fmedian = fArray[nSize/2];

	free(fArray);

	return fmedian;
}
/*
*   The generic idiom for comparing two numerical values a and b 
*   for qsort looks as (a > b) - (a < b). 
*
*/
int compare_float (const void * a, const void * b)
{
	float fa = *(const float*) a;

	float fb = *(const float*) b;

	return (fa > fb) - (fa < fb);
}

/*
*   A factorial calculator
*/

unsigned Factorial(unsigned n)
{
    if (n == 1 || 0 == n)
        return 1;
    else
        return n * Factorial(n - 1);
}

/*
*   A Log_factorial calculator
*/

double Log_Factorial(unsigned n)
{
    if (n == 1 || 0 == n)
        return 0.;
    else
        return (log((double)n) + Log_Factorial(n - 1));
}
/*
*  Convolution of two real arrays of the same size
*  using FFT. The size should be a power of 2.
*
*  The arrays are arr1, arr2. They are of the type array <double> ^
*
*  The results will be return in the arr2
* 
*/

float TheoreticalIsotopeCalculator::Convolve(array <double> ^arr1, array <double> ^arr2)
{
	if ((arr1->Length != arr2->Length) || arr1->Length < 2 || arr1->Length & (arr1->Length-1) )
    {
		printf("Array Length must be power of 2 in four1");
	}

	int i;

	array <double> ^product = gcnew array <double> (arr1->Length);

	//fft of a real array - note that the 0th and 1 elements are the real valued
	// 1st and last elements

	realft(arr1, 1);

	realft(arr2, 1);

	// the 0th and 1 elements are the real valued
	// 1st and last elements
	product[0] = arr1[0] * arr2[0]; product[1] = arr1[1] * arr2[1];

	for(i = 2; i < arr1->Length; i = i + 2)
	{
		product[i] = arr1[i] * arr2[i] - arr1[i+1] * arr2[i+1];

		product[i+1] = arr1[i] * arr2[i+1] + arr1[i+1] * arr2[i];

		//printf("product[%d] = %f, data3[%d] = %f\n", i, product[i], i, data3[i]);
	}

	//inverse fft for the convolution
	realft(product, -1);

	//assign the arr2 for return values
	for(i = 0; i < arr1->Length; i++)
	{
		arr2[i] = product[i]  * 2./(double)product->Length;
	}

	delete(product);

	return 0.f;
}

/*
*
*   The method will self-convolve the arr1 nMultiplex times
*   the resutls will be stored in arr1
*   arr1 is of type array <double> ^
*   Since the convolution uses FFT, arr1 should be of length power of 2.
*
*/
float TheoreticalIsotopeCalculator::Self_Convolve(int nMultiplex, array <double> ^arr1)
{
	if (arr1->Length < 2 || arr1->Length & (arr1->Length-1) )
    {
		printf("Array Length must be power of 2 in four1");
	}

	int i, k;

	array <double> ^temp_array = gcnew array <double> (arr1->Length);

	array <double> ^temp_array2 = gcnew array <double> (arr1->Length);

	for(i = 0; i < arr1->Length; i++)
	{
		temp_array[i] = arr1[i];

		temp_array2[i] = arr1[i];
	}

	for(k = 1; k <= nMultiplex; k++)
	{
		for(i = 0; i < arr1->Length; i++)
		{
			temp_array[i] = arr1[i];
		}

		Convolve(temp_array, temp_array2);
	}

	for(i = 0 ; i < arr1->Length; i++)
	{
		arr1[i] = temp_array2[i];
	}

	delete(temp_array); delete(temp_array2);

	return 0.0;
}


/*
*  
*   Binomial probability - calculates the binomial
*   P, given the Ntrials, p, probability of success
*   NSucceses - the number of successes
*/
double Binomial(double prob, int nSuccess, int nTrials)
{
	double pResult;

	if (nSuccess > nTrials || nSuccess < 0 || nTrials == 0)
	{
		printf("Cannot have more succeses, %d, than trials, %d\n", nSuccess, nTrials);

		exit(1);
	}

	pResult = exp((double)nSuccess* log(prob) + (double)(nTrials - nSuccess)*log(1. - prob) );

	pResult = pResult * exp(Log_Factorial(nTrials) - Log_Factorial(nSuccess) - Log_Factorial(nTrials - nSuccess));

	return pResult;
}

/*
*
*  A method to check realft function for the sum of two Poissons with means of 4 and 7. it is just for test.
*
*/
float TheoreticalIsotopeCalculator::DoRealFFT()
{
	array <double> ^data = gcnew array <double> (32);

	array <double> ^data2 = gcnew array <double> (32);

	array <double> ^data3 = gcnew array <double> (32);

	array <double> ^product = gcnew array <double> (32);

	int i;

	for(i = 0; i < 32; i++)
	{
		data[i] = exp(-4.0 + (double)i*log(4.) - Log_Factorial((unsigned)i) );

		//data[i] = exp(-4.) * pow(4., (double)i)/(double)Factorial(i);

		data2[i] = exp(-7. + (double)i*log(7.) - Log_Factorial((unsigned)i) );


		data3[i] = exp(-11. + (double)i*log(11.) - Log_Factorial((unsigned)i) );

		//if(i > 12)
			//data[i] = 0.0;


		//printf("data[%d] = %f , data2 = %f %f\n", i, data[i], data2[i], data3[i]);
	}


	if(false)
	{


		realft(data, 1);

		realft(data2, 1);

		realft(data3, 1);

		printf("The Real Part\n");

		for(i=0; i < 32; i = i + 2)
			printf("data[%d] = %f  %f\n", i/2, data[i], data2[i]);

		printf("The Imaginary Part\n");

		for(i=1; i < 32; i = i + 2)
		  printf("data[%d] = %f  %f\n", (i-1)/2, data[i], data2[i]);

		puts("THE PRODUCT AND TRUE");

		product[0] = data[0] * data2[0]; 

		product[1] = data[1] * data2[1];

		for(i = 2; i < 32; i = i + 2)
		{
			product[i] = data[i] * data2[i] - data[i+1] * data2[i+1];

			product[i+1] = data[i] * data2[i+1] + data[i+1] * data2[i];

			printf("product[%d] = %f, data3[%d] = %f\n", i, product[i], i, data3[i]);
		}

		realft(product, -1);

		realft(data3, -1);

		for(i = 0; i < 32; i++)
		{
			printf("product[%d] = %f  Control[%d] = %f\n", i, product[i] * 2./32., i, data3[i] * 2./32.);
		}
	}

	printf("Convolution\n\n");

	Convolve(data, data2);

	for(i = 0; i < 32; i++)
	{
		printf("product[%d] = %f  Control[%d] = %f\n", i, data2[i], i, data3[i]);
	}


	delete (data); delete (data2); delete(data3); delete(product);
	return 0.;
}

/*
*
*
*  A program to compute the natural isotopes
*  using Poisson distribution
*  the numbers of atoms are given (self-explanotary),
*  the results are in fIsotopes, iIsotope is the number
*  of isotopolouges
*  Lambda_1, Lambda_2, Lambda_4 - expectation values in the Poisson
*  distributions, corresponding to 1 (H, C, N, O), 2 (O, S), and 4 (S) mass unit shifts
*
*
*/

float ProteinRateConstant::PoissonIsotopesWithMasses(int iIsotope, array <float> ^fIsotopes, int nH, int nC, int nN,
	int nO, int nS)
{
	double Lambda_1, Lambda_2, Lambda_4;

	double dtemp1, dtemp2, dtemp4, dtemp;

	double dLambda1_array[32], dLambda2_array[32], dLambda4_array[32];

	double dM1, dM2, dM4;

	double dM0 = 0;

	double probH1, probH2, probC12, probC13, probN14, probN15;

	double probO16, probO17, probO18, probS32, probS33, probS34, probS36;

	int i, i1, i2, i4;

	int iH, iC, iN, iO1, iO2, iS1, iS2;



	bool bHybrid = true;

	/*
	*   Use the Poisson approximation for natural isotope distibutions calculations
	*/

	probH1 = 0.99984426 / (0.99984426 + 0.00015574);      probH2 = 0.00015574 / (0.99984426 + 0.00015574);

	probC12 = 0.988922 / (0.988922 + 0.011078);           probC13 = 0.011078 / (0.988922 + 0.011078);

	probN14 = 0.996337 / (0.996337 + 0.003663);  probN15 = 0.003663 / (0.996337 + 0.003663);

	probO16 = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);

	probO17 = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	probO18 = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	probS32 = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);  

	probS33 = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS34 = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS36 = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;        //one mass shift isotopes of all elements; c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	Lambda_1 = Lambda_1 + nO *  probO17;

	Lambda_1 = Lambda_1 + nS * probS33;

	Lambda_2 = nO * probO18;

	Lambda_2 = Lambda_2 + nS * probS34;

	Lambda_4 = nS * probS36;

	dtemp1 = dtemp2 = dtemp4 = 0.0;

	for (i = 0; i < 32; i++)
	{
		dtemp = Log_Factorial(i);

		dLambda1_array[i] = exp(-Lambda_1) * pow(Lambda_1, i) / Factorial(i);

		dtemp1 = -Lambda_1 + i * Math::Log(Lambda_1) - dtemp;

		dLambda1_array[i] = exp(dtemp1);

		dLambda2_array[i] = exp(-Lambda_2) * pow(Lambda_2, i) / Factorial(i);

		dtemp2 = -Lambda_2 + i * Math::Log(Lambda_2) - dtemp;

		//dtemp2 = -Lambda_2 + i * Math::Log(Lambda_2) - dtemp;

		dLambda2_array[i] = exp(dtemp2);

		if (nS > 0 && i == 1 && false)
		{
			//dLambda2_array[i] = exp(dtemp2) - probS34 / sqrt(2.* exp(1.0) *M_PI);

			dLambda2_array[i] = nS * probS34 * pow(1 - probS34, nS - 1) * pow(probO16, nO);

			dLambda2_array[i] = dLambda2_array[i] + nO * pow(probO16, nO - 1) * probO18 * pow(1 - probS34, nS);
		}	
		else if (nS > 0 && i == 2 && false)
		{
			dLambda2_array[i] = nS * (nS - 1) * pow(probS34, 2)* pow(1 - probS34, nS - 2) * pow(probO16, nO)/2;

			dLambda2_array[i] = dLambda2_array[i] + nO * (nO - 1) *pow(probO18, 2) * pow(probO16, nO - 2) *pow(1 - probS34, nS)/2;

			dLambda2_array[i] = dLambda2_array[i] +
				nS * probS34 * pow(1 - probS34, nS - 1) * nO *  probO18 * pow(probO16, nO - 1);
		}

		dtemp4 = -Lambda_4 + i * Math::Log(Lambda_4) - dtemp;

		dLambda4_array[i] = exp(-Lambda_4) * pow(Lambda_4, i) / Factorial(i);

		dLambda4_array[i] = exp(dtemp4);

		dIsotopeMasses[i] = 0.0;
	}

	double ftemp_array[10000];

	for (i = 1; i < 100; i++)
		ftemp_array[i] = 0.0;

	double dSprob1 = 0.0, dOprob1 = 0.0, dOprob2;

	dtemp = dtemp1 = dtemp2 = dtemp4 = 0.0;

	//the only configuration of Sulfur providing input to 1st isotope

	if (nS > 0)
	{
		dSprob1 = Log_Factorial(nS) - Log_Factorial(nS - 1) + (nS - 1)*log(probS32) + 1. * log(probS33);

		dSprob1 = exp(dSprob1);
	}
	

	dOprob1 = Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO17);

	dOprob1 = exp(dOprob1);


	double dtempO, d_factorailnO;

	array <double, 2> ^ dOProbDist = gcnew array <double, 2>(nO + 1, nO + 1);

	d_factorailnO = Log_Factorial(nO);

	Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;

	for (i = 0; i < 32; i++)
	{
		dtemp = -Lambda_1 + i * Math::Log(Lambda_1) - Log_Factorial(i);

		dLambda1_array[i] = exp(dtemp);

		ftemp_array[i] = dIsotopeMasses[i] = 0.0;
	}

	for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
	{
		for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
		{
			dtempO = d_factorailnO - Log_Factorial(iO1) - Log_Factorial(iO2) - Log_Factorial(nO - iO1 - iO2);

			dtempO = dtempO + (nO - iO1 - iO2)* log(0.9976206) + iO1 * log(0.0003790) + iO2 * log(0.0020004);

			dtempO = exp(dtempO);

			dOProbDist[iO1, iO2] = dtempO;
		}
	}


	array <double> ^a_c, ^a_h, ^a_n;

	a_c = gcnew array <double>(nC + 1); a_h = gcnew array <double>(nH + 1); a_n = gcnew array <double>(nN + 1);

	for (i = 0; i <= nC; i++)
	{
		dtemp = -probC13*nC + i * Math::Log(probC13*nC) - Log_Factorial(i);

		a_c[i] = exp(dtemp);

		//a_c[i] = Binomial(probC13, i, nC);

	}

	for (i = 0; i <= nH; i++)
	{
		a_h[i] = Binomial(probH2, i, nH);

		dtemp = -probH2*nH + i * Math::Log(probH2*nH) - Log_Factorial(i);

		a_h[i] = exp(dtemp);
	}

	for (i = 0; i <= nN; i++)
	{
		a_n[i] = Binomial(probN15, i, nN);

		dtemp = -probN15*nN + i * Math::Log(probN15*nN) - Log_Factorial(i);

		a_n[i] = exp(dtemp);
	}

	double dThreeMasses[33];

	for (i = 0; i <= 32; i++)
		dThreeMasses[i] = 0;

	for (iC = 0; iC <= nC; iC++)
	{
		for (iH = 0; iH <= nH; iH++)
		{
			for (iN = 0; iN <= nN; iN++)
			{
				i = iC + iH + iN;

				if (i < 32)
				{
					dThreeMasses[i] = dThreeMasses[i] +
						(iC*deltaC13 + iH*(mDeuterium - mHydrogen) + iN*(dN15 - dN14)) *
						a_c[iC] * a_h[iH] * a_n[iN];

					ftemp_array[i] = ftemp_array[i] +
						a_c[iC] * a_h[iH] * a_n[iN];
				}
			}
		}
	}

	for (i = 0; i < 32; i++)
	{
		if (ftemp_array[i] > 0.00000001)
		{
			dThreeMasses[i] = dThreeMasses[i] / ftemp_array[i];

			//printf("THeMass[%d] = %10.5f\n", i, dThreeMasses[i]);
		}

		//dLambda1_array[i] = ftemp_array[i];

		ftemp_array[i] = 0.0;
	}

	dtemp = (deltaC13*probC13 + (mDeuterium - mHydrogen)*probH2 + (dN15 - dN14)*probN15) / (probC13 + probH2 + probN15);

	dtemp1 = (dO17 - dO16);

	dtemp2 = (dO18 - dO16);

	array <double> ^ a_o = gcnew array <double>(2 * nO + 1);

	array <double> ^ a_om = gcnew array <double>(2 * nO + 1);

	for (iO1 = 0; iO1 <= nO; iO1++)
	{
		for (iO2 = 0; iO2 <= nO; iO2++)
		{
			if (iO1 + 2 * iO2 <= nO)
			{
				a_o[iO1 + 2 * iO2] = a_o[iO1 + 2 * iO2] + dOProbDist[iO1, iO2];

				a_om[iO1 + 2 * iO2] = a_om[iO1 + 2 * iO2] +
					(iO1 * dtemp1 + iO2*dtemp2)* dOProbDist[iO1, iO2];
			}

		}

	}


	for (i = 0; i <= nO; i++)
	{
		if (a_o[i] > 0.00000001)
		{
			a_om[i] = a_om[i] / a_o[i];
		}
	}

	for (i1 = 0; i1 <= nC; i1++)
	{
		for (i2 = 0; i2 <= nO; i2++)
		{
			if (i1 + i2 < 32)
			{
				ftemp_array[i1 + i2] = ftemp_array[i1 + i2] +
					dLambda1_array[i1] * a_o[i2];

				dIsotopeMasses[i1 + i2] = dIsotopeMasses[i1 + i2] +
					(dThreeMasses[i1] + a_om[i2]) * dLambda1_array[i1] * a_o[i2];
			}
		}
	}

	for (i = 0; i < 32; i++)
	{
		if(ftemp_array[i] > 0.0000001)
		dIsotopeMasses[i] = dIsotopeMasses[i] / ftemp_array[i];
	}


	if (nS >= 1)
	{

		array <double, 3> ^a_sp = gcnew array <double, 3>(nS + 1, nS + 1, nS + 1);

		array <double>  ^a_s = gcnew array <double>(4 * nS + 1);

		array <double> ^a_sm = gcnew array <double>(4 * nS + 1);

		int iS1, iS2, iS4, iSum, iS;

		double d_factorial_nS, dtempS;

		double dTempProb[33], dTempIsotopeMass[33];

		d_factorial_nS = Log_Factorial(nS);

		for (iS1 = 0; iS1 <= nS && iS1 < 32; iS1++)
		{
			for (iS2 = 0; (iS2 + iS1) <= nS && (iS2 + iS1) < 32; iS2++)
			{
				for (iS4 = 0; (iS2 + iS1 + iS4) <= nS && (iS2 + iS1 + iS4) < 32; iS4++)
				{

					iSum = iS1 + 2 * iS2 + 4 * iS4;

					if (iSum < 32)
					{
						dtempS = d_factorial_nS - Log_Factorial(iS1) - Log_Factorial(iS2) -
							Log_Factorial(iS4) - Log_Factorial(nS - iS1 - iS2 - iS4);

						dtempS = dtempS + (nS - iS1 - iS2 - iS4) * log(0.9504074) + iS1 * log(0.0074869) +
							iS2 * log(0.0419599) + iS4 * log(0.0001458);

						dtempS = exp(dtempS);

						a_sp[iS1, iS2, iS4] = dtempS;

						a_s[iSum] = a_s[iSum] + dtempS;

						a_sm[iSum] = a_sm[iSum] +
							(iS1 * (dS33 - dS32) + iS2 * (dS34 - dS32) + iS4 * (dS36 - dS32)) * dtempS;
					}

				} //for (iS4 = 0; iS4 <= nS && iS4 < 32; iS4++)
			}
		}


		for (iS = 0; iS < 4 * nS; iS++)
		{
			if (a_s[iS] > 0.00000001)
			{
				a_sm[iS] = a_sm[iS] / a_s[iS];
			}
		}

		for (iS = 0; iS <= 4 * nS; iS++)
		{
			for (i = 0; i < 32; i++)
			{
				iSum = iS + i;

				if (iSum < 32)
				{
					dTempProb[iSum] = dTempProb[iSum] + a_s[iS] * ftemp_array[i];

					dTempIsotopeMass[iSum] = dTempIsotopeMass[iSum] + (dIsotopeMasses[i] + a_sm[iS]) * a_s[iS] * ftemp_array[i];
				}
			}
		}

		delete a_sp, a_s, a_sm;


		for (i = 0; i < 32; i++)
		{
			dIsotopeMasses[i] = dTempIsotopeMass[i];

			ftemp_array[i] = dTempProb[i];
		}


		dM0 = AminoAcid_Sequence_Mass(szPeptide);

		//dM0 = 2976.40658;

		for (int i = 0; i < iIsotope; i++)
		{
			fIsotopes[i] = ftemp_array[i];

			if (i <= 1)
			{
				dIsotopeMasses[i] = dM0 + dIsotopeMasses[i];

				//printf("The mass[%d] = %15.7f\n", i, dIsotopeMasses[i]);
			}
			else if (fIsotopes[i] > 0.00000001)
			{
				dIsotopeMasses[i] = dM0 + dIsotopeMasses[i] / fIsotopes[i];

				//printf("The mass[%d] = %15.7f\n", i, dIsotopeMasses[i]);
			}

		}


		return 0.f;


		if(false)
		for (i4 = 0; i4 <= nS && i4 < 32; i4++)
		{
			//dtemp4 = exp(-Lambda_4) * pow(Lambda_4, i4)/Factorial((unsigned)i4);

			//for(int i2 = 0; i2 <= nO && i2 <=2 ; i2++)      //only 5 isotopes, means can have only two of 18O
			for (i2 = 0; i2 <= nO && i2 < 32; i2++)
			{
				//dtemp2 = exp(-Lambda_2) * pow(Lambda_2, i2)/Factorial((unsigned)i2);

				dtemp = dLambda2_array[i2] * dLambda4_array[i4];

				//for(int i1 = 0; (i1 + 2 * i2 + 4 * i4) < iIsotope; i1++)
				//for (int i1 = 0; (i1 + 2 * i2 + 4 * i4) < 2*nO + 4*nS; i1++)
				for (i1 = 0; (i1 + 2 * i2 + 4 * i4) < 32; i1++)
				{
					if (i1 + 2 * i2 + 4 * i4 >= 2)
					{

						ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dLambda1_array[i1] * dtemp;


						//dIsotopeMasses[i1 + 2 * i2 + 4 * i4] = dIsotopeMasses[i1 + 2 * i2 + 4 * i4] +
						//dLambda1_array[i1] * dtemp * (i1 * dM1 + i2 * dM2 + i4 * dM4);

						if (i1 + 2 * i2 + 4 * i4 == 4)
							printf("Prob[%d, %d, %d] %f %f %f\n", i1, i2, i4, dLambda1_array[i1],
								dLambda2_array[i2], dLambda4_array[i4]);


						dIsotopeMasses[i1 + 2 * i2 + 4 * i4] = dIsotopeMasses[i1 + 2 * i2 + 4 * i4] +
							(i1*(deltaC13*probC13 + (mDeuterium - mHydrogen)*probH2 + (dN15 - dN14)*probN15 +
							(dS33 - dS32)*probS33 + (dO17 - dO16)*probO17) /
								(probC13 + probH2 + probN15 + probO17 + probS33) +
								i2*((dS34 - dS32)*probS34 + (dO18 - dO16)*probO18) / (probO18 + probS34) +
								i4 * (dS36 - dS32)) * dLambda1_array[i1] * dtemp;


					}
				}
			}
		}
	}
	else
	{

		if (bHybrid)
		{

			int iO1, iO2;

			double dtempO, d_factorailnO;

			array <double, 2> ^ dOProbDist = gcnew array <double, 2>(nO + 1, nO+1);

			d_factorailnO = Log_Factorial(nO);

			Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;

			for (i = 0; i < 32; i++)
			{
				dtemp = -Lambda_1 + i * Math::Log(Lambda_1) - Log_Factorial(i);

				dLambda1_array[i] = exp(dtemp);

				ftemp_array[i] = dIsotopeMasses[i] = 0.0;
			}

			for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
			{
				for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
				{
					dtempO = d_factorailnO - Log_Factorial(iO1) - Log_Factorial(iO2) - Log_Factorial(nO - iO1 - iO2);

					dtempO = dtempO + (nO - iO1 - iO2)* log(0.9976206) + iO1 * log(0.0003790) + iO2 * log(0.0020004);

					dtempO = exp(dtempO);

					dOProbDist[iO1, iO2] = dtempO;
				}
			}

			int iC, iH, iN;

			array <double> ^a_c, ^a_h, ^a_n;

			a_c = gcnew array <double>(nC+1); a_h = gcnew array <double>(nH+1); a_n = gcnew array <double>(nN+1);

			for (i = 0; i <= nC; i++)
			{
				dtemp = -probC13*nC + i * Math::Log(probC13*nC) - Log_Factorial(i);

				a_c[i] = exp(dtemp);

				//a_c[i] = Binomial(probC13, i, nC);

			}

			for (i = 0; i <= nH; i++)
			{
				a_h[i] = Binomial(probH2, i, nH);

				dtemp = -probH2*nH + i * Math::Log(probH2*nH) - Log_Factorial(i);

				a_h[i] = exp(dtemp);
			}

			for (i = 0; i <= nN; i++)
			{
				a_n[i] = Binomial(probN15, i, nN);

				dtemp = -probN15*nN + i * Math::Log(probN15*nN) - Log_Factorial(i);

				a_n[i] = exp(dtemp);
			}
	
			double dThreeMasses[33];

			for (i = 0; i <= 32; i++)
				dThreeMasses[i] = 0;

			for (iC = 0; iC <= nC; iC++)
			{
				for (iH = 0; iH <= nH; iH++)
				{
					for (iN = 0; iN <= nN; iN++)
					{
						i = iC + iH + iN;

						if (i < 32)
						{
							dThreeMasses[i] = dThreeMasses[i] +
								(iC*deltaC13 + iH*(mDeuterium - mHydrogen) + iN*(dN15 - dN14)) *
								a_c[iC] * a_h[iH] * a_n[iN];

							ftemp_array[i] = ftemp_array[i] +
								a_c[iC] * a_h[iH] * a_n[iN];
						}
					}
				}
			}

			for (i = 0; i < 32; i++)
			{
				if (ftemp_array[i] > 0.00000001)
				{
					dThreeMasses[i] = dThreeMasses[i] / ftemp_array[i];

					printf("THeMass[%d] = %10.5f\n", i, dThreeMasses[i]);
				}

				//dLambda1_array[i] = ftemp_array[i];

				ftemp_array[i] = 0.0;
			}


			if(false)
			for (iC = 0; iC <= nC; iC++)
			{
				for (iH = 0; iH <= nH; iH++)
				{
					for (iN = 0; iN <= nN; iN++)
					{
						for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
						{
							for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
							{
								if (iC + iH + iN + iO1 + 2 * iO2 < 32)
								{
									ftemp_array[iC + iH + iN + iO1 + 2 * iO2] = ftemp_array[iC + iH + iN + iO1 + 2 * iO2] +
										a_c[iC] * a_h[iH] * a_n[iN] *
										//Binomial(probC13, iC, nC) * Binomial(probH2, iH, nH) * Binomial(probN15, iN, nN) *
										dOProbDist[iO1, iO2];

									dIsotopeMasses[iC + iH + iN + iO1 + 2 * iO2] = dIsotopeMasses[iC + iH + iN + iO1 + 2 * iO2] +
										(iC*deltaC13 + iH*(mDeuterium - mHydrogen) + iN*(dN15 - dN14) +
											iO1 * (dO17 - dO16) + iO2*(dO18 - dO16) ) *
										a_c[iC] * a_h[iH] * a_n[iN] * 
										//Binomial(probC13, iC, nC) * Binomial(probH2, iH, nH) * Binomial(probN15, iN, nN) *
										dOProbDist[iO1, iO2];
								}
							}
						}
					}
				}
			}

			dtemp = (deltaC13*probC13 + (mDeuterium - mHydrogen)*probH2 + (dN15 - dN14)*probN15) / (probC13 + probH2 + probN15);

			dtemp1 = (dO17 - dO16);

			dtemp2 = (dO18 - dO16);

			array <double> ^ a_o = gcnew array <double>(2*nO + 1);

			array <double> ^ a_om = gcnew array <double>(2 * nO + 1);



			for (iO1 = 0; iO1 <= nO; iO1++)
			{
				for (iO2 = 0; iO2 <= nO; iO2++)
				{
					if (iO1 + 2 * iO2 <= nO)
					{
						a_o[iO1 + 2 * iO2]  = a_o[iO1 + 2 * iO2] + dOProbDist[iO1, iO2];

						a_om[iO1 + 2 * iO2] = a_om[iO1 + 2 * iO2] + 
							(iO1 * dtemp1 + iO2*dtemp2)* dOProbDist[iO1, iO2];
					}
					
				}

			}


			for (i = 0; i <= nO; i++)
			{
				if (a_o[i] > 0.00000001)
				{
					a_om[i] = a_om[i] / a_o[i];
				}
			}



			for (i1 = 0; i1 <= nC; i1++)
			{
				for (i2 = 0; i2 <= nO; i2++)
				{
					if (i1 + i2 < 32)
					{
						ftemp_array[i1 + i2] = ftemp_array[i1 + i2] +
							dLambda1_array[i1] * a_o[i2];

						dIsotopeMasses[i1 + i2] = dIsotopeMasses[i1 + i2] +
							(dThreeMasses[i1] + a_om[i2]) * dLambda1_array[i1] * a_o[i2];
					}
				}
			}



			if(false)
			for (i1 = 0; i1 < nC && i1 < 32; i1++)
			{
				for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
				{
					for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
					{
						if (i1 + iO1 + 2 * iO2 < 32)
						{
							ftemp_array[i1 + iO1 + 2 * iO2] = ftemp_array[i1 + iO1 + 2 * iO2] +
									dLambda1_array[i1] * dOProbDist[iO1, iO2];

							dIsotopeMasses[i1 + iO1 + 2 * iO2] = dIsotopeMasses[i1 + iO1 + 2 * iO2] +
								//(i1*dtemp +
								(dThreeMasses[i1] + 
										iO1 * dtemp1 + iO2*dtemp2 ) * dLambda1_array[i1] * dOProbDist[iO1, iO2];
						}
					}
				}
			}



			if(false) // this part reproduces accurate mass and isotope computation
			for (iC = 0; iC <= nC; iC++)
			{
				for (iH = 0; iH <= nH; iH++)
				{
					for (iN = 0; iN <= nH; iN++)
					{
						for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
						{
							for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
							{
								if (iC + iH + iN + iO1 + 2 * iO2 < 32)
								{
									ftemp_array[iC + iH + iN + iO1 + 2 * iO2] = ftemp_array[iC + iH + iN + iO1 + 2 * iO2] +
										Binomial(probC13, iC, nC) * Binomial(probH2, iH, nH) * Binomial(probN15, iN, nN) *
										dOProbDist[iO1, iO2];

									dIsotopeMasses[iC + iH + iN + iO1 + 2 * iO2] = dIsotopeMasses[iC + iH + iN + iO1 + 2 * iO2] +
										(iC*deltaC13 + iH*(mDeuterium - mHydrogen) + iN*(dN15 - dN14) +
											iO1 * (dO17 - dO16) + iO2*(dO18 - dO16) ) *
										Binomial(probC13, iC, nC) * Binomial(probH2, iH, nH) * Binomial(probN15, iN, nN) *
										dOProbDist[iO1, iO2];
								}
							}
						}
					}
				}
			}



			delete dOProbDist, a_c, a_h, a_n;
		}
		else
		{
			for (i2 = 0; i2 <= nO && i2 < 32; i2++)
			{
				for (i1 = 0; (i1 + 2 * i2) < 32; i1++)
				{
					if (i1 + 2 * i2 >= 2)
					{

						ftemp_array[i1 + 2 * i2] = ftemp_array[i1 + 2 * i2] + dLambda1_array[i1] * dLambda2_array[i2];

						dIsotopeMasses[i1 + 2 * i2] = dIsotopeMasses[i1 + 2 * i2] +
							(i1*(deltaC13*probC13 + (mDeuterium - mHydrogen)*probH2 + (dN15 - dN14)*probN15 +
							(dO17 - dO16)*probO17) /
								(probC13 + probH2 + probN15 + probO17) +
								i2*(dO18 - dO16))* dLambda1_array[i1] * dLambda2_array[i2];

					}
				}
			}
		}
		
	}

	printf("pC12 = %f pN14 = %f pH1 = %f pO16 = %f pS32 = %f pS36 = %f\n", probC12, probN14, probH1, probO16, probS32, probS36);

	printf("nC = %d nN = %d nH = %d nS = %d nO = %d\n", nC, nN, nH, nS, nO);

	

	//ftemp_array[0] = dLambda1_array[0] * dLambda2_array[0] * dLambda4_array[0];

	ftemp_array[0] = pow(probH1, nH) * pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

	ftemp_array[1] = Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO);

	ftemp_array[1] = ftemp_array[1] + Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

	ftemp_array[1] = ftemp_array[1] + Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO);

	ftemp_array[1] = ftemp_array[1] + 
		nO * probO17 *pow(probO16, nO-1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH);

	dIsotopeMasses[1] = deltaC13 * Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO);

	dIsotopeMasses[1] = dIsotopeMasses[1] +
						(mDeuterium - mHydrogen) * Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

	dIsotopeMasses[1] = dIsotopeMasses[1] +
		(dN15 - dN14) *  Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO);

	dIsotopeMasses[1] = dIsotopeMasses[1] +
		(dO17 - dO16) *  nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH);

	printf("MassOfOne = %15.7f %15.7f\n", dIsotopeMasses[1],
		
		dIsotopeMasses[1]/ (Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO) +
			Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO) +
			Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO) +
			nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH)));

	if (nS >= 1)
	{
		ftemp_array[0] = ftemp_array[0] * pow(probS32, nS);

		ftemp_array[0] = pow(probS32, nS) * pow(probH1, nH) * pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

		ftemp_array[1] = ftemp_array[1] * pow(probS32, nS);

		ftemp_array[1] = ftemp_array[1] +
			nS * probS33 * pow(probS32, nS - 1) * 
				pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) * pow(probO16, nO);


		dIsotopeMasses[1] = dIsotopeMasses[1] * pow(probS32, nS);

		dIsotopeMasses[1] = dIsotopeMasses[1] +
			(dS33- dS32)*nS * probS33 * pow(probS32, nS - 1) *
			pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) * pow(probO16, nO);

		dIsotopeMasses[1] = dIsotopeMasses[1] / (
			Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO)*pow(probS32, nS) +
			Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO)*pow(probS32, nS) +
			Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO)*pow(probS32, nS) +
			nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH)*pow(probS32, nS) +
			nS * probS33 * pow(probS32, nS - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) * pow(probO16, nO)  );
	}
	else //no Sulphur atoms
	{
		dIsotopeMasses[1] = dIsotopeMasses[1] / (
			Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO) +
			Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO) +
			Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO) +
			nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) );
	}

	printf("S = %10.7f, %10.7f, %10.7f, %10.7f\n ", probS32, probS33, probS34, probS36);

	//dM0 = AminoAcid_Sequence_Mass("CAHTNDEFKRMSGHQWYPAKLIVNM");

	//dM0 = 2976.40658;

	for (int i = 0; i < iIsotope; i++)
	{
		fIsotopes[i] = ftemp_array[i];

		if (i <= 1)
		{
			dIsotopeMasses[i] = dM0 + dIsotopeMasses[i];

			//printf("The mass[%d] = %15.7f\n", i, dIsotopeMasses[i]);
		}
		else if (fIsotopes[i] > 0.00000001)
		{
			dIsotopeMasses[i] = dM0 + dIsotopeMasses[i] / fIsotopes[i];

			//printf("The mass[%d] = %15.7f\n", i, dIsotopeMasses[i]);	
		}
			
	}

	return 0.f;
}


/*
*
*
*  A program to compute the natural isotopes
*  using Poisson distribution
*  the numbers of atoms are given (self-explanotary),
*  the results are in fIsotopes, iIsotope is the number
*  of isotopolouges
*  Lambda_1, Lambda_2, Lambda_4 - expectation values in the Poisson
*  distributions, corresponding to 1 (H, C, N, O), 2 (O, S), and 4 (S) mass unit shifts
*
*
*/

float ProteinRateConstant::PoissonIsotopes(int iIsotope, array <float> ^fIsotopes, int nH, int nC, int nN,
	int nO, int nS)
{
	double Lambda_1, Lambda_2, Lambda_4;

	double dtemp1, dtemp2, dtemp4, dtemp;

	double dLambda1_array[32], dLambda2_array[32], dLambda4_array[32];

	double dIsotopeMasses[32], dM1, dM2, dM4;

	double probH1, probH2, probC12, probC13, probN14, probN15;

	double probO16, probO17, probO18, probS32, probS33, probS34, probS36;

	int i, i1, i2, i4;

	int iH, iC, iN, iO1, iO2, iS1, iS2;
	/*
	*   Use the Poisson approximation for natural isotope distibutions calculations
	*/

	probH1 = 0.99984426 / (0.99984426 + 0.00015574);      probH2 = 0.00015574 / (0.99984426 + 0.00015574);

	probC12 = 0.988922 / (0.988922 + 0.011078);           probC13 = 0.011078 / (0.988922 + 0.011078);

	probN14 = 0.996337 / (0.996337 + 0.003663);  probN15 = 0.003663 / (0.996337 + 0.003663);

	probO16 = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);

	probO17 = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	probO18 = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	probS32 = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS33 = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS34 = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS36 = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;        //one mass shift isotopes of all elements; c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	Lambda_1 = Lambda_1 + nO *  probO17;

	Lambda_1 = Lambda_1 + nS * probS33;

	Lambda_2 = nO * probO18;

	Lambda_2 = Lambda_2 + nS * probS34;

	Lambda_4 = nS * probS36;

	dtemp1 = dtemp2 = dtemp4 = 0.0;

	for (i = 0; i < 32; i++)
	{
		dLambda1_array[i] = exp(-Lambda_1) * pow(Lambda_1, i) / Factorial(i);

		dLambda2_array[i] = exp(-Lambda_2) * pow(Lambda_2, i) / Factorial(i);

		dLambda4_array[i] = exp(-Lambda_4) * pow(Lambda_4, i) / Factorial(i);

		dIsotopeMasses[i] = 0.0;
	}

	float ftemp_array[10000];

	for (i = 1; i < 100; i++)
		ftemp_array[i] = 0.0;

	double dSprob1 = 0.0, dOprob1 = 0.0, dOprob2;



	//the only configuration of Sulfur providing input to 1st isotope

	if (nS > 0)
	{
		dSprob1 = Log_Factorial(nS) - Log_Factorial(nS - 1) + (nS - 1)*log(probS32) + 1. * log(probS33);

		dSprob1 = exp(dSprob1);
	}


	dOprob1 = Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO17);

	dOprob1 = exp(dOprob1);

	for (i4 = 0; i4 <= nS && i4 < 32; i4++)
	{
		//dtemp4 = exp(-Lambda_4) * pow(Lambda_4, i4)/Factorial((unsigned)i4);

		//for(int i2 = 0; i2 <= nO && i2 <=2 ; i2++)      //only 5 isotopes, means can have only two of 18O
		for (i2 = 0; i2 <= nO && i2 < 32; i2++)
		{
			//dtemp2 = exp(-Lambda_2) * pow(Lambda_2, i2)/Factorial((unsigned)i2);

			dtemp = dLambda2_array[i2] * dLambda4_array[i4];

			//for(int i1 = 0; (i1 + 2 * i2 + 4 * i4) < iIsotope; i1++)
			//for (int i1 = 0; (i1 + 2 * i2 + 4 * i4) < 2*nO + 4*nS; i1++)
			for (i1 = 0; (i1 + 2 * i2 + 4 * i4) < 32; i1++)
			{
				//dtemp1 = exp(-Lambda_1) * pow(Lambda_1, (int)i1)/Factorial(i1);

				//fIsotopes[i1 + 2 * i2 + 4 * i4] = fIsotopes[i1 + 2 * i2 + 4 * i4] + dtemp1 * dtemp2 * dtemp4;

				//ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dtemp1 * dtemp2 * dtemp4;

				if (i1 + 2 * i2 + 4 * i4 >= 2)
				{
					ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dLambda1_array[i1] * dtemp;

					dIsotopeMasses[i1 + 2 * i2 + 4 * i4] = dIsotopeMasses[i1 + 2 * i2 + 4 * i4] +
						dLambda1_array[i1] * dtemp * (i1 * dM1 + i2 * dM2 + i4 * dM4);
				}

			}
		}
	}


	//ftemp_array[0] = dLambda1_array[0] * dLambda2_array[0] * dLambda4_array[0];

	ftemp_array[0] = pow(0.99984426, nH) * pow(0.988922, nC) * pow(0.996337, nN);

	ftemp_array[0] = ftemp_array[0] * pow(0.9976206 / (0.9976206 + 0.0003790 + 0.0020004), nO) * pow(0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458), nS);



	for (int i = 0; i < iIsotope; i++)
	{
		fIsotopes[i] = ftemp_array[i];
	}

	return 0.f;
}

/*
*
*  Poisson Isotopes and Mean Mass Calculations
*
*/

float ProteinRateConstant::PoissonIsotopesMeanMasses(int iIsotope, array <float> ^fIsotopes, array <double> ^dMassesOfIsotopes, 
	          int nH, int nC, int nN, int nO, int nS)
{
	double Lambda_1, Lambda_2, Lambda_4;

	double dtemp1, dtemp2, dtemp4, dtemp;

	double dLambda1_array[32], dLambda2_array[32], dLambda4_array[32];

	double dIsotopeMasses[32], dM1, dM2, dM4;

	double probH1, probH2, probC12, probC13, probN14, probN15;

	double probO16, probO17, probO18, probS32, probS33, probS34, probS36;

	double M0;

	int i, i1, i2, i4;

	int iH, iC, iN, iO1, iO2, iS1, iS2;
	/*
	*   Use the Poisson approximation for natural isotope distibutions calculations
	*/

	probH1 = 0.99984426 / (0.99984426 + 0.00015574);      probH2 = 0.00015574 / (0.99984426 + 0.00015574);

	probC12 = 0.988922 / (0.988922 + 0.011078);           probC13 = 0.011078 / (0.988922 + 0.011078);

	probN14 = 0.996337 / (0.996337 + 0.003663);  probN15 = 0.003663 / (0.996337 + 0.003663);

	probO16 = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);

	probO17 = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	probO18 = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	probS32 = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS33 = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS34 = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS36 = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;        //one mass shift isotopes of all elements; c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	Lambda_1 = Lambda_1 + nO *  probO17;

	Lambda_1 = Lambda_1 + nS * probS33;

	Lambda_2 = nO * probO18;

	Lambda_2 = Lambda_2 + nS * probS34;

	Lambda_4 = nS * probS36;

	dtemp1 = dtemp2 = dtemp4 = 0.0;

	for (i = 0; i < 32; i++)
	{
		dtemp = -Lambda_1 + i * log(Lambda_1) - Log_Factorial(i);
		
		//dLambda1_array[i] = exp(-Lambda_1) * pow(Lambda_1, i) / Factorial(i);

		dLambda1_array[i] = exp(dtemp);

		dtemp = -Lambda_2 + i * log(Lambda_2) - Log_Factorial(i);

		//dLambda2_array[i] = exp(-Lambda_2) * pow(Lambda_2, i) / Factorial(i);

		dLambda2_array[i] = exp(dtemp);

		dtemp = -Lambda_4 + i * log(Lambda_4) - Log_Factorial(i);

		//dLambda4_array[i] = exp(-Lambda_4) * pow(Lambda_4, i) / Factorial(i);

		dLambda4_array[i] = exp(dtemp);

		dIsotopeMasses[i] = 0.0;
	}

	float ftemp_array[10000];

	for (i = 1; i < 100; i++)
		ftemp_array[i] = 0.0;

	double dOprob1 = 0.0, dOprob2;


	if (0 == nS)
	{
		dLambda4_array[0] = 1.0;
	}


	dOprob1 = Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO17);

	dOprob1 = exp(dOprob1);

	dM1 = nH * probH2 * (mDeuterium - mHydrogen) + nC * probC13 * deltaC13 + nN * probN15 * (dN15 - dN14);

	dM1 = dM1 + nO * probO17 * (dO17 - dO16) +  nS * probS33 * (dS33 - dS32);

	dM1 = dM1 / (nH * probH2 + nC * probC13 + nN * probN15 + nO * probO17 + nS * probS33);

	dM2 = nO * probO18 * (dO18 - dO16) + nS * probS34 * (dS34 - dS32);

	dM2 = dM2 / (nO * probO18 + nS * probS34);

	dM4 = dS36 - dS32;

	for (i4 = 0; i4 <= nS && i4 < 32; i4++)
	{
		//dtemp4 = exp(-Lambda_4) * pow(Lambda_4, i4)/Factorial((unsigned)i4);

		//for(int i2 = 0; i2 <= nO && i2 <=2 ; i2++)      //only 5 isotopes, means can have only two of 18O

		for (i2 = 0; i2 <= nO && i2 < 32; i2++)
		{
			//dtemp2 = exp(-Lambda_2) * pow(Lambda_2, i2)/Factorial((unsigned)i2);

			dtemp = dLambda2_array[i2] * dLambda4_array[i4];

			//for(int i1 = 0; (i1 + 2 * i2 + 4 * i4) < iIsotope; i1++)
			//for (int i1 = 0; (i1 + 2 * i2 + 4 * i4) < 2*nO + 4*nS; i1++)
			for (i1 = 0; (i1 + 2 * i2 + 4 * i4) < 32; i1++)
			{
				//dtemp1 = exp(-Lambda_1) * pow(Lambda_1, (int)i1)/Factorial(i1);

				//fIsotopes[i1 + 2 * i2 + 4 * i4] = fIsotopes[i1 + 2 * i2 + 4 * i4] + dtemp1 * dtemp2 * dtemp4;

				//ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dtemp1 * dtemp2 * dtemp4;

				if (i1 + 2 * i2 + 4 * i4 >= 1)
				{
					ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dLambda1_array[i1] * dtemp;

					dIsotopeMasses[i1 + 2 * i2 + 4 * i4] = dIsotopeMasses[i1 + 2 * i2 + 4 * i4] +
						dLambda1_array[i1] * dtemp * (i1 * dM1 + i2 * dM2 + i4 * dM4);
						//dLambda1_array[i1] * dtemp * (i1 * dM1 + i2 * dM2 + i4 * dM4);
				}

			}
		}
	}


	//ftemp_array[0] = dLambda1_array[0] * dLambda2_array[0] * dLambda4_array[0];

	ftemp_array[0] = pow(0.99984426, nH) * pow(0.988922, nC) * pow(0.996337, nN);

	ftemp_array[0] = ftemp_array[0] * pow(0.9976206 / (0.9976206 + 0.0003790 + 0.0020004), nO) * pow(0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458), nS);

	M0 = AminoAcid_Sequence_Mass(szPeptide);

	for (int i = 0; i < iIsotope; i++)
	{
		fIsotopes[i] = ftemp_array[i];

		dMassesOfIsotopes[i] = 0.0;

		if (ftemp_array[i] > 0.0000001)
		{
			//dMassesOfIsotopes[i] = dIsotopeMasses[i] / ftemp_array[i];

			dIsotopeMasses[i] = dIsotopeMasses[i] / ftemp_array[i];

			dIsotopeMasses[i] = M0 + dIsotopeMasses[i];

			dMassesOfIsotopes[i] = dIsotopeMasses[i];
			 
			//printf("%DISS[%d] = %10.7f\n", i, dIsotopeMasses[i]);
		}
		  
	}

	return 0.f;
}
/*
*
*  This method uses a hybrid - Poisson and exact distributions
*   to determine the isotope distributions and isotope masses.
*   the massses are  returned in dMassesOfIsotopes, Isotope intensities are 
*   in fIsotopes
*    nLimitIsotope - is the maximum number of isotopes to report
*
*/

float ProteinRateConstant::HybridIsotopes(int iIsotope, array <float> ^fIsotopes, array <double> ^dMassesOfIsotopes, 
											int nH, int nC, int nN, int nO, int nS)
{
	double Lambda_1, Lambda_2, Lambda_4;

	double dtemp1, dtemp2, dtemp4, dtemp;

	double dLambda1_array[32], dLambda2_array[32], dLambda4_array[32];

	double dM1, dM2, dM4;

	double dM0 = 0;

	double probH1, probH2, probC12, probC13, probN14, probN15;

	double probO16, probO17, probO18, probS32, probS33, probS34, probS36;

	int i, i1, i2, i4;

	int iH, iC, iN, iO1, iO2, iS1, iS2, nLimitIsotope;

	nLimitIsotope = 33;

	double d_CutOff = 0.0000000000001;

	double dSprob1 = 0.0, dOprob1 = 0.0, dOprob2;

	double ftemp_array[10000];

	double dtempO, d_factorailnO;

	array <double, 2> ^ dOProbDist = gcnew array <double, 2>(nO + 1, nO + 1);

	array <double> ^a_c, ^a_h, ^a_n;

	a_c = gcnew array <double>(nC + 1); a_h = gcnew array <double>(nH + 1); a_n = gcnew array <double>(nN + 1);

	array <double> ^ a_o = gcnew array <double>(2 * nO + 2);  array <double> ^ a_om = gcnew array <double>(2 * nO + 2);
	
	array <double> ^ dThreeMasses = gcnew array <double>(nLimitIsotope);

	/*
	*   Use the Poisson approximation for natural isotope distibutions calculations
	*/

	probH1 = 0.99984426 / (0.99984426 + 0.00015574);      probH2 = 0.00015574 / (0.99984426 + 0.00015574);

	probC12 = 0.988922 / (0.988922 + 0.011078);           probC13 = 0.011078 / (0.988922 + 0.011078);

	probN14 = 0.996337 / (0.996337 + 0.003663);  probN15 = 0.003663 / (0.996337 + 0.003663);

	probO16 = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);

	probO17 = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	probO18 = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	probS32 = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS33 = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS34 = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS36 = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;        //one mass shift isotopes of all elements; c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	Lambda_1 = Lambda_1 + nO *  probO17;

	Lambda_1 = Lambda_1 + nS * probS33;

	Lambda_2 = nO * probO18;

	Lambda_2 = Lambda_2 + nS * probS34;

	Lambda_4 = nS * probS36;

	dtemp1 = dtemp2 = dtemp4 = 0.0;

	for (i = 0; i < 32; i++)
	{
		dtemp = Log_Factorial(i);

		dtemp1 = -Lambda_1 + i * Math::Log(Lambda_1) - dtemp;

		dLambda1_array[i] = exp(dtemp1);

		dLambda2_array[i] = exp(-Lambda_2) * pow(Lambda_2, i) / Factorial(i);

		dtemp2 = -Lambda_2 + i * Math::Log(Lambda_2) - dtemp;

		dLambda2_array[i] = exp(dtemp2);

		dtemp4 = -Lambda_4 + i * Math::Log(Lambda_4) - dtemp;

		dLambda4_array[i] = exp(dtemp4);

		dIsotopeMasses[i] = 0.0;
	}

	for (i = 1; i < 100; i++)
		ftemp_array[i] = 0.0;



	dtemp = dtemp1 = dtemp2 = dtemp4 = 0.0;

	//the only configuration of Sulfur providing input to 1st isotope

	if (nS > 0)
	{
		dSprob1 = Log_Factorial(nS) - Log_Factorial(nS - 1) + (nS - 1)*log(probS32) + 1. * log(probS33);

		dSprob1 = exp(dSprob1);
	}


	dOprob1 = Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO17);

	dOprob1 = exp(dOprob1);

	d_factorailnO = Log_Factorial(nO);

	Lambda_1 = probH2 * nH + probC13 * nC + probN15 * nN;

	for (i = 0; i < 32; i++)
	{
		dtemp = -Lambda_1 + i * Math::Log(Lambda_1) - Log_Factorial(i);

		dLambda1_array[i] = exp(dtemp);

		ftemp_array[i] = dIsotopeMasses[i] = 0.0;
	}

	for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
	{
		for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
		{
			dtempO = d_factorailnO - Log_Factorial(iO1) - Log_Factorial(iO2) - Log_Factorial(nO - iO1 - iO2);

			dtempO = dtempO + (nO - iO1 - iO2)* log(0.9976206) + iO1 * log(0.0003790) + iO2 * log(0.0020004);

			dtempO = exp(dtempO);

			dOProbDist[iO1, iO2] = dtempO;
		}
	}

	for (i = 0; i <= nC; i++)
	{
		dtemp = -probC13*nC + i * Math::Log(probC13*nC) - Log_Factorial(i);

		a_c[i] = exp(dtemp);

		//a_c[i] = Binomial(probC13, i, nC);

	}

	for (i = 0; i <= nH; i++)
	{
		a_h[i] = Binomial(probH2, i, nH);

		dtemp = -probH2*nH + i * Math::Log(probH2*nH) - Log_Factorial(i);

		a_h[i] = exp(dtemp);
	}


	for (i = 0; i <= nN; i++)
	{
		dtemp = -probN15*nN + i * Math::Log(probN15*nN) - Log_Factorial(i);

		a_n[i] = exp(dtemp);
	}

	for (i = 0; i < nLimitIsotope; i++)
		dThreeMasses[i] = 0;


	for (iC = 0; iC <= nC; iC++)
	{
		for (iH = 0; iH <= nH; iH++)
		{
			for (iN = 0; iN <= nN; iN++)
			{
				i = iC + iH + iN;

				if (i < 32 && i < nLimitIsotope)
				{
					dThreeMasses[i] = dThreeMasses[i] +
						(iC*deltaC13 + iH*(mDeuterium - mHydrogen) + iN*(dN15 - dN14)) *
						a_c[iC] * a_h[iH] * a_n[iN];

					ftemp_array[i] = ftemp_array[i] +
						a_c[iC] * a_h[iH] * a_n[iN];
				}
			}
		}
	}
	


	for (i = 0; i < 32; i++)
	{
		if (ftemp_array[i] > d_CutOff)
		{
			dThreeMasses[i] = dThreeMasses[i] / ftemp_array[i];
		}
		//dLambda1_array[i] = ftemp_array[i];

		ftemp_array[i] = 0.0;
	}

	dtemp = (deltaC13*probC13 + (mDeuterium - mHydrogen)*probH2 + (dN15 - dN14)*probN15) / (probC13 + probH2 + probN15);

	dtemp1 = (dO17 - dO16);

	dtemp2 = (dO18 - dO16);

	for (iO1 = 0; iO1 <= nO; iO1++)
	{
		for (iO2 = 0; iO2 <= nO; iO2++)
		{
			if (iO1 + 2 * iO2 <= nO)
			{
				a_o[iO1 + 2 * iO2] = a_o[iO1 + 2 * iO2] + dOProbDist[iO1, iO2];

				a_om[iO1 + 2 * iO2] = a_om[iO1 + 2 * iO2] +
					(iO1 * dtemp1 + iO2*dtemp2)* dOProbDist[iO1, iO2];
			}
		}
	}

	for (i = 0; i <= nO; i++)
	{
		if (a_o[i] >  d_CutOff )
		{
			a_om[i] = a_om[i] / a_o[i];
		}
	}

	for (i1 = 0; i1 <= nC; i1++)
	{
		for (i2 = 0; i2 <= nO; i2++)
		{
			if (i1 + i2 < 32)
			{
				ftemp_array[i1 + i2] = ftemp_array[i1 + i2] +
					dLambda1_array[i1] * a_o[i2];

				dIsotopeMasses[i1 + i2] = dIsotopeMasses[i1 + i2] +
					(dThreeMasses[i1] + a_om[i2]) * dLambda1_array[i1] * a_o[i2];
			}
		}
	}
	

	if (nS >= 1)
	{
		array <double, 3> ^a_sp = gcnew array <double, 3>(nS + 1, nS + 1, nS + 1);

		array <double>  ^a_s = gcnew array <double>(4 * nS + 1);

		array <double> ^a_sm = gcnew array <double>(4 * nS + 1);

		array <double> ^ dTempProb = gcnew array <double>(nLimitIsotope); 
		
		array <double>  ^dTempIsotopeMass = gcnew array <double> (nLimitIsotope);

		for (i = 0; i < 32; i++)
		{
			if (ftemp_array[i] > d_CutOff )
				dIsotopeMasses[i] = dIsotopeMasses[i] / ftemp_array[i];
			else
			{
				ftemp_array[i] = 0.0;
			}
		}

		int iS1, iS2, iS4, iSum, iS;

		double d_factorial_nS, dtempS;

		for (i = 0; i < nLimitIsotope; i++)
		{
			dTempProb[i] = dTempIsotopeMass[i] = 0.0;
		}

		d_factorial_nS = Log_Factorial(nS);

		for (iS1 = 0; iS1 <= nS && iS1 < 32; iS1++)
		{
			for (iS2 = 0; (iS2 + iS1) <= nS && (iS2 + iS1) < 32; iS2++)
			{
				for (iS4 = 0; (iS2 + iS1 + iS4) <= nS && (iS2 + iS1 + iS4) < 32; iS4++)
				{

					iSum = iS1 + 2 * iS2 + 4 * iS4;

					if (iSum < 32)
					{
						dtempS = d_factorial_nS - Log_Factorial(iS1) - Log_Factorial(iS2) -
							Log_Factorial(iS4) - Log_Factorial(nS - iS1 - iS2 - iS4);

						dtempS = dtempS + (nS - iS1 - iS2 - iS4) * log(0.9504074) + iS1 * log(0.0074869) +
							iS2 * log(0.0419599) + iS4 * log(0.0001458);

						dtempS = exp(dtempS);

						a_sp[iS1, iS2, iS4] = dtempS;

						a_s[iSum] = a_s[iSum] + dtempS;

						a_sm[iSum] = a_sm[iSum] +
							(iS1 * (dS33 - dS32) + iS2 * (dS34 - dS32) + iS4 * (dS36 - dS32)) * dtempS;
					}

				} //for (iS4 = 0; iS4 <= nS && iS4 < 32; iS4++)
			}
		}


		for (iS = 0; iS <= 4 * nS; iS++)
		{
			if (a_s[iS] > d_CutOff)
			{
				a_sm[iS] = a_sm[iS] / a_s[iS];
			}
			else
			{
				a_s[iS] = 0.0;
			}

			//printf("Probs[%d] = %10.5f %10.5f, %15.10f\n", iS, a_s[iS], a_sm[iS], d_CutOff);
		}

		for (iS = 0; iS <= 4 * nS; iS++)
		{
			for (i = 0; i < 32; i++)
			{
				iSum = iS + i;

				if (iSum < 32)
				{
					dTempProb[iSum] = dTempProb[iSum] + a_s[iS] * ftemp_array[i];

					dTempIsotopeMass[iSum] = dTempIsotopeMass[iSum] + (dIsotopeMasses[i] + a_sm[iS]) * a_s[iS] * ftemp_array[i];
				}
			}
		}

		for (i = 0; i < 32; i++)
		{
			dIsotopeMasses[i] = dTempIsotopeMass[i];

			ftemp_array[i] = dTempProb[i];
		}

		delete a_sp, a_s, a_sm, dTempProb, dTempIsotopeMass;
	}

	//printf("pC12 = %f pN14 = %f pH1 = %f pO16 = %f pS32 = %f pS36 = %f\n", probC12, probN14, probH1, probO16, probS32, probS36);

	//printf("nC = %d nN = %d nH = %d nS = %d nO = %d\n", nC, nN, nH, nS, nO);


	//ftemp_array[0] = dLambda1_array[0] * dLambda2_array[0] * dLambda4_array[0];

	ftemp_array[0] = pow(probH1, nH) * pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

	ftemp_array[1] = Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO);

	ftemp_array[1] = ftemp_array[1] + Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

	ftemp_array[1] = ftemp_array[1] + Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO);

	ftemp_array[1] = ftemp_array[1] +
		nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH);

	dIsotopeMasses[1] = deltaC13 * Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO);

	dIsotopeMasses[1] = dIsotopeMasses[1] +
		(mDeuterium - mHydrogen) * Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

	dIsotopeMasses[1] = dIsotopeMasses[1] +
		(dN15 - dN14) *  Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO);

	dIsotopeMasses[1] = dIsotopeMasses[1] +
		(dO17 - dO16) *  nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH);

	if (nS >= 1)
	{
		ftemp_array[0] = ftemp_array[0] * pow(probS32, nS);

		ftemp_array[0] = pow(probS32, nS) * pow(probH1, nH) * pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO);

		ftemp_array[1] = ftemp_array[1] * pow(probS32, nS);

		ftemp_array[1] = ftemp_array[1] +
			nS * probS33 * pow(probS32, nS - 1) *
			pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) * pow(probO16, nO);


		dIsotopeMasses[1] = dIsotopeMasses[1] * pow(probS32, nS);

		dIsotopeMasses[1] = dIsotopeMasses[1] +
			(dS33 - dS32)*nS * probS33 * pow(probS32, nS - 1) *
			pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) * pow(probO16, nO);

		dIsotopeMasses[1] = dIsotopeMasses[1] / (
			Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO)*pow(probS32, nS) +
			Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO)*pow(probS32, nS) +
			Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO)*pow(probS32, nS) +
			nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH)*pow(probS32, nS) +
			nS * probS33 * pow(probS32, nS - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH) * pow(probO16, nO));
	}
	else //no Sulphur atoms
	{
		dIsotopeMasses[1] = dIsotopeMasses[1] / (
			Binomial(probC13, 1, nC)*pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO) +
			Binomial(probH2, 1, nH)*pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO) +
			Binomial(probN15, 1, nN)*pow(probC12, nC) * pow(probH1, nH) * pow(probO16, nO) +
			nO * probO17 *pow(probO16, nO - 1) * pow(probC12, nC) * pow(probN14, nN) * pow(probH1, nH));
	}


	dM0 = AminoAcid_Sequence_Mass(szPeptide);

	//dM0 = 2976.40658;

	for (int i = 0; i < iIsotope; i++)
	{
		fIsotopes[i] = ftemp_array[i];

		if (i <= 1)
		{
			dMassesOfIsotopes[i] = dM0 + dIsotopeMasses[i];

			dIsotopeMasses[i] = dM0 + dIsotopeMasses[i];

			//printf("The mass[%d] = %15.7f\n", i, dIsotopeMasses[i]);

		}
		else if (fIsotopes[i] > d_CutOff)
		{
			dMassesOfIsotopes[i] = dM0 + dIsotopeMasses[i] / fIsotopes[i];

			dIsotopeMasses[i] = dM0 + dIsotopeMasses[i] / fIsotopes[i];

			//printf("The mass[%d] = %15.7f\n", i, dIsotopeMasses[i]);
		}

	}


	delete dThreeMasses, a_c, a_h, a_n, a_o, dOProbDist;

	return 0.f;
}

/*
*
*
*  A program to compute the natural isotopes
*  using Poisson distribution
*  the numbers of atoms are given (self-explanotary),
*  the results are in fIsotopes, iIsotope is the number
*  of isotopolouges
*  Lambda_1, Lambda_2, Lambda_4 - expectation values in the Poisson
*  distributions, corresponding to 1 (H, C, N, O), 2 (O, S), and 4 (S) mass unit shifts
*
*
*/
float PoissonIsotopes_Fastest(int iIsotope, array <float> ^fIsotopes, int nH, int nC, int nN,
						int nO, int nS)
{
	double Lambda_1, Lambda_2, Lambda_4;

	double dtemp1, dtemp2, dtemp4, dtemp;

	double dLambda1_array[32], dLambda2_array[32], dLambda4_array[32];

	double dIsotopeMasses[32], dM1, dM2, dM4;

	double probH1, probH2, probC12, probC13, probN14, probN15;

	double probO16, probO17, probO18, probS32, probS33, probS34, probS36;

	int i, i1, i2, i4;

	int iH, iC, iN, iO1, iO2, iS1, iS2;
	/*
	*   Use the Poisson approximation for natural isotope distibutions calculations
	*/

	probH1 = 0.99984426 / (0.99984426 + 0.00015574);      probH2 = 0.00015574 / (0.99984426 + 0.00015574);

	probC12 = 0.988922 / (0.988922 + 0.011078);           probC13 = 0.011078 / (0.988922 + 0.011078);

	probN14 = 0.996337 / (0.996337 + 0.003663);  probN15 = 0.003663 / (0.996337 + 0.003663);

	probO16 = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);
	
	probO17 = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);
	
	probO18 = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	probS32 = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS33 = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	probS34 = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);
	
	probS36 = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	Lambda_1 = probH2 * nH + probC13 * nC +  probN15 * nN;        //one mass shift isotopes of all elements; c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	Lambda_1 = Lambda_1 + nO *  probO17;

	Lambda_1 = Lambda_1 + nS * probS33;

	Lambda_2 = nO * probO18;

	Lambda_2 = Lambda_2 + nS * probS34;

	Lambda_4 = nS * probS36;

	dtemp1 = dtemp2 =  dtemp4 = 0.0;

	for (i = 0; i < 32; i++)
	{
		dLambda1_array[i] = exp(-Lambda_1) * pow(Lambda_1, i) / Factorial(i);

		dLambda2_array[i] = exp(-Lambda_2) * pow(Lambda_2, i) / Factorial(i);

		dLambda4_array[i] = exp(-Lambda_4) * pow(Lambda_4, i) / Factorial(i);

		dIsotopeMasses[i] = 0.0;
	}

	float ftemp_array[10000];

	for (i = 1; i < 100; i++)
		ftemp_array[i] = 0.0; 

	// masses of isotopes:

	dM1 = deltaC13 * Binomial(probC13, 1, nC) * pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO) * pow(probS32, nS);

	dM1 = dM1 + (mDeuterium - mHydrogen) * Binomial(probH2, 1, nH) * pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO) * pow(probS32, nS);

	dM1 = dM1 + (dN15 - dN14) * Binomial(probN15, 1, nN) * pow(probC12, nC) * pow(probH1, nH) * pow(0.9976206, nO) * pow(probS32, nS);


	double dSprob1 = 0.0, dOprob1 = 0.0, dOprob2;
	
	//the only configuration of Sulfur providing input to 1st isotope

	if (nS >= 1)
	{
		dSprob1 = Log_Factorial(nS) - Log_Factorial(nS - 1) + (nS - 1)*log(probS32) + 1. * log(probS33);

		dSprob1 = exp(dSprob1);
	}
	else
	{
		dSprob1 = 0.0;
	}

	dOprob1 = Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO17);

	dOprob1 = exp(dOprob1);

	dM1 = dM1 + (dO17 - dO16) * dOprob1 * pow(probC12, nC) * pow(probH1, nH) * pow(probN14, nN) * pow(probS32, nS);

	dM1 = dM1 + (dS33 - dS32) * dSprob1  * pow(probC12, nC) * pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO);

	dM1 = dM1 / ( Binomial(probC13, 1, nC) * pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO) * pow(probS32, nS) +
			Binomial(probH2, 1, nH) * pow(probC12, nC) * pow(probN14, nN) * pow(probO16, nO) * pow(probS32, nS) +
			Binomial(probN15, 1, nN) * pow(probC12, nC) * pow(probH1, nH) *pow(probO16, nO) * pow(probS32, nS) +
			dOprob1 * pow(probC12, nC) * pow(probH1, nH) * pow(probN14, nN) * pow(probS32, nS) +
			dSprob1 * pow(probC12, nC) * pow(probH1, nH) * pow(probN14, nN) * pow(probO16, nO) );
	//dM1 = (dD)

	printf("dM1 = %15.7f\n", dM1);

	dM2 = (dO18 - dO16) * Binomial(0.0020004 / (0.9976206 + 0.0003790 + 0.0020004), 1, nO) * pow(1 - 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458), nS) +
		(dS34 - dS32) * Binomial(0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458), 1, nS) *
		pow(1 - 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004), nO);


	dM2 = (Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO18)) * (dO18 - dO16);

	dM2 = dM2 + (Log_Factorial(nS) - Log_Factorial(nS - 1) + (nS - 1)*log(probS32) + 1. * log(probS34) ) * (dS34 - dS32);


	dM2 = dM2 / ( (Log_Factorial(nO) - Log_Factorial(nO - 1) + (nO - 1)*log(probO16) + 1. * log(probO18)) * pow(probS32, nS) +
		(Log_Factorial(nS) - Log_Factorial(nS - 1) + (nS - 1)*log(probS32) + 1. * log(probS34)) * pow(probO16, nO)  );

	printf("dM2  = %10.7f\n", dM2);

	dM4 = dS36 - dS32;

	double dtemp_Prob = 0.0;

	dM2 = 0.0;

	for ( iH = 0; iH <= nH && iH <= 2; iH++)
	{
		for (iC = 0; iC <= nC && iC <= 2; iC++)
		{
			if (iC + iH > 2)
			{
				continue;
			}
			for (iN = 0; iN <= nN && iN <= 2; iN++)
			{
				if (iN + iC + iH > 2)
				{
					continue;
				}
				for (iO1 = 0; iO1 <= nO && iO1 <= 2; iO1++)
				{
					if (iO1 + iN + iC + iH > 2)
					{
						continue;
					}

					for (iO2 = 0; iO2 <= nO && iO2 <= 1; iO2++)
					{
						if (iO1 + iN + iC + iH + 2 * iO2 > 2)
						{
							continue;
						}
						int kO = iO2 + iO1;

						if (kO <= nO)
						{
							double dtempO;

							dtempO = Log_Factorial(nO) - Log_Factorial(iO1) - Log_Factorial(iO2) - Log_Factorial(nO - iO1 - iO2);

							dtempO = dtempO + (nO - iO1 - iO2)* log(0.9976206) + iO1 * log(0.0003790) + iO2 * log(0.0020004);

							dtempO = exp(dtempO);

							for (int iS1 = 0; iS1 <= nS && iS1 <= 2; iS1++)
							{
								for (int iS2 = 0; iS2 <= nS && iS2 <= 1; iS2++)
								{
									int kS;

									kS = iS1 + iS2;

									if (kS <= nS)
									{
										int iSum;

										iSum = iH + iC + iN + iO1 + iS1 + 2 * iO2 + 2 * iS2;

										double dtempS;

										if (iSum > 2)
										{
											continue;
										}
										else if (iSum == 2)
										{
											dtempS = Log_Factorial(nS) - Log_Factorial(iS1) - Log_Factorial(iS2) -
												 - Log_Factorial(nS - iS1 - iS2);

											dtempS = dtempS + (nS - iS1 - iS2) * log(0.9504074) + iS1 * log(0.0074869) +
												iS2 * log(0.0419599);

											dtempS = exp(dtempS);

											dtemp = Binomial(probH2, iH, nH) * Binomial(probC13, iC, nC) * Binomial(probN15, iN, nN);

											dtemp_Prob = dtemp_Prob + dtemp * dtempO * dtempS;

											dM2 = dM2 +
												(iH * (mDeuterium - mHydrogen) + iC * deltaC13 + iN * (dN15 - dN14) +
													iO1 * (dO17 - dO16) + iO2 * (dO18 - dO16) +
													iS1 * (dS33 - dS32) + iS2 * (dS34 - dS32)) * dtemp * dtempO * dtempS;
										}
									}

								}
							}
						}
					}
				}
			}
		}
	}
	
	//dIsotopeMasses[1] = 3514.66981015 + dM1;

	//dIsotopeMasses[2] = 3514.66981015 + dM2 / dtemp_Prob;

	printf("dM1 = %15.7f, dM2 = %15.7f dM4 = %15.7f, Mass = %20.7f\n", dM1, dM2/dtemp_Prob, dM4, (dM1 + 3514.66981015), (dM2 + 3514.66981015));

	dM2 = dM2 / dtemp_Prob; 

	for (i4 = 0; i4 <= nS && i4 < 32; i4++)
	{
		//dtemp4 = exp(-Lambda_4) * pow(Lambda_4, i4)/Factorial((unsigned)i4);

		//for(int i2 = 0; i2 <= nO && i2 <=2 ; i2++)      //only 5 isotopes, means can have only two of 18O
		for (i2 = 0; i2 <= nO && i2 < 32; i2++)
		{
			//dtemp2 = exp(-Lambda_2) * pow(Lambda_2, i2)/Factorial((unsigned)i2);

			dtemp = dLambda2_array[i2] * dLambda4_array[i4];

			//for(int i1 = 0; (i1 + 2 * i2 + 4 * i4) < iIsotope; i1++)
			//for (int i1 = 0; (i1 + 2 * i2 + 4 * i4) < 2*nO + 4*nS; i1++)
			for (i1 = 0; (i1 + 2 * i2 + 4 * i4) < 32; i1++)
			{
				//dtemp1 = exp(-Lambda_1) * pow(Lambda_1, (int)i1)/Factorial(i1);

				//fIsotopes[i1 + 2 * i2 + 4 * i4] = fIsotopes[i1 + 2 * i2 + 4 * i4] + dtemp1 * dtemp2 * dtemp4;

				//ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dtemp1 * dtemp2 * dtemp4;

				if (i1 + 2 * i2 + 4 * i4 >= 2)
				{
					ftemp_array[i1 + 2 * i2 + 4 * i4] = ftemp_array[i1 + 2 * i2 + 4 * i4] + dLambda1_array[i1] * dtemp;

					dIsotopeMasses[i1 + 2 * i2 + 4 * i4] = dIsotopeMasses[i1 + 2 * i2 + 4 * i4] +
						dLambda1_array[i1] * dtemp * (i1 * dM1 + i2 * dM2 + i4 * dM4);
				}
				
			}
		}
	}

	//adjust the masses of isotopes

	for (i = 0; i < 32; i++)
	{
		if (i <= 2)
		{
			printf("I[%d] = %15.7f\n", i, dIsotopeMasses[i]);
		}
		else if (ftemp_array[i] > 0.0000001)
		{
			dIsotopeMasses[i] = dIsotopeMasses[i] / ftemp_array[i];

			printf("I[%d] = %15.7f\n", i, (dIsotopeMasses[i] + 3514.66981015));
		}
			
	}

	
	//ftemp_array[0] = dLambda1_array[0] * dLambda2_array[0] * dLambda4_array[0];

	ftemp_array[0] = pow(0.99984426, nH) * pow(0.988922, nC) * pow(0.996337, nN);

	ftemp_array[0] = ftemp_array[0] * pow(0.9976206/(0.9976206 + 0.0003790 + 0.0020004), nO) * pow(0.9504074/(0.9504074 + 0.0074869 + 0.0419599 + 0.0001458), nS);

	printf("0th isotope  = %10.7f\n", ftemp_array[0]);

	AminoAcid_Sequence_Mass("ADEQRFTYIOPWERTYQIYTPOASDFGHLKCVNM");


	for (int i = 0; i < iIsotope; i++)
	{
		fIsotopes[i] = ftemp_array[i];
	}

	return 0.f;
}
/*
*   Compute Isotope Masses accurately
*   Uses full extention
*/
void TheoreticalIsotopeCalculator::IsotopeMasses(int nH, int nC, int nN, int nO, int nS)
{
	array <double> ^h_arr, ^ c_arr, ^ n_arr, ^ o_arr, ^s_arr, ^final_isotope;

	double delta_Deuterium, dM0;

	int iArraySize = 32;

	int i, j, iH, iC, iN, iO1, iO2, iS1, iS2, iS4, kO, kS, iSum;

	h_arr = gcnew array <double>(iArraySize);    c_arr = gcnew array <double>(iArraySize);

	n_arr = gcnew array <double>(iArraySize);    o_arr = gcnew array <double>(iArraySize);

	s_arr = gcnew array <double>(iArraySize);    final_isotope = gcnew array <double>(iArraySize);

	for (i = 0; i < iArraySize; i++)
	{
		h_arr[i] = c_arr[i] = n_arr[i] = o_arr[i] = s_arr[i] = final_isotope[i] = 0.0;
	}


	o_arr[0] = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);  o_arr[1] = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	o_arr[2] = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	s_arr[0] = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[1] = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	s_arr[2] = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[4] = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	for (j = 0; j < nH && j < iArraySize; j++)
	{
		h_arr[j] = Binomial(0.00015574 / (0.99984426 + 0.00015574), j, nH);
	}

	for (j = 0; j < nC && j < iArraySize; j++)
	{
		c_arr[j] = Binomial(0.011078 / (0.988922 + 0.011078), j, nC);
	}

	for (j = 0; j < nN && j < iArraySize; j++)
	{
		n_arr[j] = Binomial(0.003663 / (0.996337 + 0.003663), j, nN);
	}

	double dtemp = 0.0, dtempO, dtempS, d_factorailnO, d_factorial_nS;

	double dIsotopeMasses[32], dtemp_array[32];

	d_factorailnO = Log_Factorial(nO);

	d_factorial_nS = Log_Factorial(nS);  delta_Deuterium = mDeuterium - mHydrogen;

	for (i = 0; i < 32; i++)
	{
		dIsotopeMasses[i] = dtemp_array[i] = 0.0;
	}

	for (iH = 0; iH <= nH && iH < 32; iH++)
	{
		for (iC = 0; iC <= nC && iC < 32; iC++)
		{
			for (iN = 0; iN <= nN && iN < 32; iN++)
			{
				dtemp = h_arr[iH] * c_arr[iC] * n_arr[iN];

				for (iO1 = 0; iO1 <= nO && iO1 < 32; iO1++)
				{
					for (iO2 = 0; (iO2 + iO1) <= nO && iO2 < 32; iO2++)
					{
						dtempO = d_factorailnO - Log_Factorial(iO1) - Log_Factorial(iO2) - Log_Factorial(nO - iO1 - iO2);

						dtempO = dtempO + (nO - iO1 - iO2)* log(0.9976206) + iO1 * log(0.0003790) + iO2 * log(0.0020004);

						//dtempO = d_factorailnO - Log_Factorial(iO1) - Log_Factorial(iO2) - Log_Factorial(nO - kO);

						//dtempO = dtempO + (nO - kO)* log(0.9976206) + iO1 * log(0.0003790) + iO2 * log(0.0020004);

						dtempO = exp(dtempO);

						for (iS1 = 0; iS1 <= nS && iS1 < 32; iS1++)
						{
							for (iS2 = 0; (iS2 + iS1) <= nS && (iS2 + iS1) < 32; iS2++)
							{
								for (iS4 = 0; (iS2 + iS1 + iS4) <= nS && (iS2 + iS1 + iS4) < 32; iS4++)
								{

									iSum = iH + iC + iN + iO1 + iS1 + 2 * iO2 + 2 * iS2 + 4 * iS4;

									if (iSum < 32)
									{
										dtempS = d_factorial_nS - Log_Factorial(iS1) - Log_Factorial(iS2) -
											Log_Factorial(iS4) - Log_Factorial(nS - iS1 - iS2 - iS4);

										dtempS = dtempS + (nS - iS1 - iS2 - iS4) * log(0.9504074) + iS1 * log(0.0074869) +
												iS2 * log(0.0419599) + iS4 * log(0.0001458);

										dtempS = exp(dtempS);

										dtemp_array[iSum] = dtemp_array[iSum] + dtemp * dtempO * dtempS;

										dIsotopeMasses[iSum] = dIsotopeMasses[iSum] +
												(iH * delta_Deuterium + iC * deltaC13 + iN * (dN15 - dN14) +
													iO1 * (dO17 - dO16) + iO2 * (dO18 - dO16) +
													iS1 * (dS33 - dS32) + iS2 * (dS34 - dS32) + iS4 * (dS36 - dS32)) * dtemp * dtempO * dtempS;
									}

								}  //for (iS4 = 0; iS4 <= nS && iS4 < 32; iS4++)
							}
						}
					}
				}
			}
		}

	}

	dIsotopeMasses_FromCalculator = gcnew array <double>(64);

	dM0 = AminoAcid_Sequence_Mass(szPeptideInCalculator);

	for (i = 0; i < 32; i++)
	{
		if (dtemp_array[i] > 0.000001)
		{
			//printf("Accurate_IsotopeMasses[%d]_Corrected = %15.7f %15.7f\n", i, dtemp_array[i],
				//(dIsotopeMasses[i] / dtemp_array[i] + dM0));

			dIsotopeMasses[i] = dIsotopeMasses[i] / dtemp_array[i];

			dIsotopeMasses_FromCalculator[i] = dM0 + dIsotopeMasses[i];

			//printf("MASS = %f\n", dIsotopeMasses_FromCalculator[i]);
			
		}
	}


	return;
}

/*
   A method to compute isotope exactly using binomial and multinomial methods
*/

void TheoreticalIsotopeCalculator::Caller_StraightConvolution(array <float> ^fIsotopes, int nH, int nC, int nN, int nO, int nS)
{

	array <double> ^h_arr, ^ c_arr, ^ n_arr, ^ o_arr, ^s_arr, ^final_isotope;

	int iArraySize = fIsotopes->Length;

	int i, j;

	h_arr = gcnew array <double>(iArraySize);    c_arr = gcnew array <double>(iArraySize);

	n_arr = gcnew array <double>(iArraySize);    o_arr = gcnew array <double>(iArraySize);

	s_arr = gcnew array <double>(iArraySize);    final_isotope = gcnew array <double>(iArraySize);

	for (i = 0; i < iArraySize; i++)
	{
		h_arr[i] = c_arr[i] = n_arr[i] = o_arr[i] = s_arr[i] = final_isotope[i] = 0.0;
	}

	h_arr[0] = 0.99984426 / (0.99984426 + 0.00015574);  h_arr[1] = 0.00015574 / (0.99984426 + 0.00015574);

	c_arr[0] = 0.988922 / (0.988922 + 0.011078);       c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	n_arr[0] = 0.996337 / (0.996337 + 0.003663);     n_arr[1] = 0.003663 / (0.996337 + 0.003663);

	o_arr[0] = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);  o_arr[1] = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	o_arr[2] = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	s_arr[0] = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[1] = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	s_arr[2] = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[4] = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	//isotope masses
	IsotopeMasses(nH, nC, nN, nO, nS);
	
	

	for (j = 0; j < nH && j < iArraySize; j++)
	{
		h_arr[j] = Binomial(0.00015574 / (0.99984426 + 0.00015574), j, nH);
	}

	for (j = 0; j < nC && j < iArraySize; j++)
	{
		c_arr[j] = Binomial(0.011078 / (0.988922 + 0.011078), j, nC);
	}

	for (j = 0; j < nN && j < iArraySize; j++)
	{
		n_arr[j] = Binomial(0.003663 / (0.996337 + 0.003663), j, nN);
	}

	
	for (i = 0; i < iArraySize; i++)
	{
		final_isotope[i] = h_arr[i];
	}

	//convolve H & C

	StraightConvolution(c_arr, final_isotope);

	//convolve CH with N

	StraightConvolution(n_arr, final_isotope);
	
	//convolvle CHN with O

	for(i = 1; i <= nO; i++)
	{
		StraightConvolution(o_arr, final_isotope);
	}

	//convolvle CHNO with S

	for (i = 1; i <= nS; i++)
	{
		StraightConvolution(s_arr, final_isotope);
	}

	for (i = 0; i < fIsotopes->Length; i++)
	{
		fIsotopes[i] = final_isotope[i];
	}


	delete c_arr, n_arr,  o_arr, s_arr, final_isotope, h_arr;

	return;
}

/*

*/

void TheoreticalIsotopeCalculator::AssignPeptideSequence(String ^sPeptide)
{
	szPeptideInCalculator = sPeptide;
}

/*
	A method to do convolution of two arrays in a  straigforward way
	Used this to benchmark the FFT and Poisson method computations
	of Isotope distributios
*/

float TheoreticalIsotopeCalculator::StraightConvolution(array <double> ^fArray1, array <double> ^fArray2)
{
	int i, j, n;

	double dtemp;

	//array <double> ^dResult = gcnew array <double> (fArray2->Length);

	double dResult[32];

	for (n = 0; n < fArray2->Length; n++)
	{
		dtemp = 0.0;

		for (i = 0; i <= fArray1->Length && i <= n; i++)
		{
			dtemp = dtemp + fArray2[n - i] * fArray1[i];

			//printf("dtemp = %10.5f\n", dtemp);
		}

		dResult[n] = dtemp;
	}


	for (i = 0; i < fArray2->Length; i++)
	{
		fArray2[i] = dResult[i];

		//printf("New = %10.5f\n", fArray2[i]);
	}

	//delete dResult;

	return 0.;
}
/*
*
*  A program to compute the natural isotopes
*  using FFT
*  the numbers of atoms are given (self-explanotary),
*  the results are in fIsotopes, iIsotope is the number
*  of isotopolouges
*  The program uses Bionomial distributions for H, C and N
*  and the FFTs for Oxygen and Sulphur, and for the molecule
*
*  index - will hold the index of first atom  type (0 -H, 1-C, 2- N, 3 -O, 4 - S) that has non-zero number of atoms
*/
float TheoreticalIsotopeCalculator::FourierIsotopes(array <float> ^fIsotopes, int nH, int nC, int nN,
						int nO, int nS)
{

	if( 0 == nH &&  0 == nC &&  0 ==  nN &&
						0 == nO &&  0 == nS)
	{
		printf("There should be at least one atom for isotope calculation\n");

		exit (1);
	}

	int i, k, j, iIsotope;

	iIsotope = fIsotopes->Length;

	int iArraySize = 0;

	bool bFirst = false;

	for(i = 0; i < 32; i++)
	{
		if(pow(2., i) > iIsotope)
		{
			break;
		}
	}

	iArraySize = (int) pow(2., i);

	if(1 > iArraySize)
	{
		printf("Problem with ArraySize Assignment\n");

		exit (1);
	}

	if(iArraySize < 8)
	{
		iArraySize = 8;
	}

	/*
	*  assign the natural isotope probabilities   
	*/

	array <double> ^h_arr, ^ c_arr, ^ n_arr, ^ o_arr, ^s_arr, ^final_isotope;

	h_arr = gcnew array <double> (iArraySize);    c_arr = gcnew array <double> (iArraySize);

	n_arr = gcnew array <double> (iArraySize);    o_arr = gcnew array <double> (iArraySize);

	s_arr = gcnew array <double> (iArraySize);    final_isotope = gcnew array <double> (iArraySize);

	for(i = 0; i < iArraySize; i++)
	{
		h_arr[i] = c_arr[i] = n_arr[i] = o_arr[i] = s_arr[i] = final_isotope[i] = 0.0;
	}

	h_arr[0] = 0.99984426 / (0.99984426 + 0.00015574);  h_arr[1] = 0.00015574 / (0.99984426 + 0.00015574);

	c_arr[0] = 0.988922 / (0.988922 +  0.011078);       c_arr[1] = 0.011078 / (0.988922 +  0.011078);

	n_arr[0] = 0.996337 / (0.996337 + 0.003663);     n_arr[1] = 0.003663   / (0.996337 + 0.003663);

	o_arr[0] = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);  o_arr[1] = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	o_arr[2] = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	s_arr[0] = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[1] = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	s_arr[2] =  0.0419599/ (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[4] = 0.0001458  / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);


	//compute the Hydrogen Isotopes using Binomial Distribution.
	if (0 == nH || nC == 0 || 0 == nN || nO == 0)
	{
		printf("No Hydrogen Atoms in a biomolecule?  Exiting...\n");

		exit(1);
	}

	for (j = 0; j < nH && j < iArraySize; j++)
	{
		h_arr[j] = Binomial(0.00015574 / (0.99984426 + 0.00015574), j, nH);
	}

	//compute the Carbon Isotopes using Binomial Distribution.
	for (j = 0; j < nC && j < iArraySize; j++)
	{
		c_arr[j] = Binomial(0.011078 / (0.988922 + 0.011078), j, nC);
	}

	//compute the Nitrogen Isotopes using Binomial Distribution.
	for (j = 0; j < nN && j < iArraySize; j++)
	{
		n_arr[j] = Binomial(0.003663 / (0.996337 + 0.003663), j, nN);
	}

	
	//Oxygen Isotopes from self-convolution
	if(nO > 1)
	{
		Self_Convolve(nO - 1, o_arr);
	}


	//Sulfur Isotopes - from Self-convolution
	if(nS > 1)
	{
		if(nS > 1)
		{
			Self_Convolve(nS - 1, s_arr);
		}
	}

	for (i = 0; i < iArraySize; i++)
	{
		final_isotope[i] = h_arr[i];
	}

	//convolve H & C

	Convolve(c_arr, final_isotope);

	//convolve CH with N

	Convolve(n_arr, final_isotope);

	//convolvle CHN with O

	Convolve(o_arr, final_isotope);

	//convolvle CHNO with S
	if (nS >= 1)
	{
		Convolve(s_arr, final_isotope);
	}

	for (i = 0; i < fIsotopes->Length; i++)
	{
		fIsotopes[i] = final_isotope[i];
	}

	delete(h_arr); delete (n_arr); delete (c_arr); delete(o_arr); delete (s_arr);

	delete (final_isotope);

	return 0.f;
}

/*
*
*  A program to compute the natural isotopes
*  using FFT
*  the numbers of atoms are given (self-explanotary),
*  the results are in fIsotopes, iIsotope is the number
*  of isotopolouges
*  The program uses Bionomial distributions for H, C and N
*  and the FFTs for Oxygen and Sulphur, and for the molecule
*
*  index - will hold the index of first atom  type (0 -H, 1-C, 2- N, 3 -O, 4 - S) that has non-zero number of atoms

* if this function works, replace the FourierIsotopes with this function !!!!
*/
float TheoreticalIsotopeCalculator::FourierIsotopes_Temp(array <float> ^fIsotopes, int nH, int nC, int nN,
	int nO, int nS)
{

	if (0 == nH && 0 == nC && 0 == nN &&
		0 == nO && 0 == nS)
	{
		printf("There should be at least one atom for isotope calculation\n");

		exit(1);
	}

	int i, k, j, iIsotope;

	iIsotope = fIsotopes->Length;

	int iArraySize = 0;

	bool bFirst = false;

	double dtemp = 0.0;

	for (i = 0; i < 32; i++)
	{
		if (pow(2., i) > iIsotope)
		{
			break;
		}
	}

	iArraySize = (int)pow(2., i);

	if (1 > iArraySize)
	{
		printf("Problem with ArraySize Assignment\n");

		exit(1);
	}

	if (iArraySize < 8)
	{
		iArraySize = 8;
	}

	/*
	*  assign the natural isotope probabilities
	*/

	array <double> ^h_arr, ^ c_arr, ^ n_arr, ^ o_arr, ^s_arr, ^final_isotope;

	h_arr = gcnew array <double>(iArraySize);    c_arr = gcnew array <double>(iArraySize);

	n_arr = gcnew array <double>(iArraySize);    o_arr = gcnew array <double>(iArraySize);

	s_arr = gcnew array <double>(iArraySize);    final_isotope = gcnew array <double>(iArraySize);

	for (i = 0; i < iArraySize; i++)
	{
		h_arr[i] = c_arr[i] = n_arr[i] = o_arr[i] = s_arr[i] = final_isotope[i] = 0.0;
	}

	h_arr[0] = 0.99984426 / (0.99984426 + 0.00015574);  h_arr[1] = 0.00015574 / (0.99984426 + 0.00015574);

	c_arr[0] = 0.988922 / (0.988922 + 0.011078);       c_arr[1] = 0.011078 / (0.988922 + 0.011078);

	n_arr[0] = 0.996337 / (0.996337 + 0.003663);     n_arr[1] = 0.003663 / (0.996337 + 0.003663);

	o_arr[0] = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);  o_arr[1] = 0.0003790 / (0.9976206 + 0.0003790 + 0.0020004);

	o_arr[2] = 0.0020004 / (0.9976206 + 0.0003790 + 0.0020004);

	s_arr[0] = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[1] = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

	s_arr[2] = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);    s_arr[4] = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);


	//compute the Hydrogen Isotopes using Binomial Distribution.
	if (0 == nH || nC == 0 || 0 == nN || nO == 0)
	{
		printf("No Hydrogen Atoms in a biomolecule?  Exiting...\n");

		exit(1);
	}

	for (j = 0; j < nH && j < iArraySize; j++)
	{
		h_arr[j] = Binomial(0.00015574 / (0.99984426 + 0.00015574), j, nH);
	}

	//compute the Carbon Isotopes using Binomial Distribution.
	for (j = 0; j < nC && j < iArraySize; j++)
	{
		c_arr[j] = Binomial(0.011078 / (0.988922 + 0.011078), j, nC);
	}

	//compute the Nitrogen Isotopes using Binomial Distribution.
	for (j = 0; j < nN && j < iArraySize; j++)
	{
		n_arr[j] = Binomial(0.003663 / (0.996337 + 0.003663), j, nN);
	}


	//Oxygen Isotopes from self-convolution
	if (nO > 1)
	{
		Self_Convolve(nO - 1, o_arr);
	}


	//Sulfur Isotopes - from Self-convolution
	if (nS > 1)
	{
		if (nS > 1)
		{
			Self_Convolve(nS - 1, s_arr);
		}
	}

	for (i = 0; i < iArraySize; i++)
	{
		final_isotope[i] = h_arr[i];
	}

	//convolve H & C

	//Convolve(c_arr, final_isotope);


	array <double> ^product = gcnew array <double>(c_arr->Length);

	array <double> ^product2 = gcnew array <double>(c_arr->Length);

	//fft of a real array - note that the 0th and 1 elements are the real valued
	// 1st and last elements

	realft(c_arr, 1);

	realft(h_arr, 1);

	realft(n_arr, 1);

	realft(o_arr, 1);

	if(nS >= 1)
		realft(s_arr, 1);


	//C & H 
	// the 0th and 1 elements are the real valued
	// 1st and last elements
	product[0] = c_arr[0] * h_arr[0]; product[1] = c_arr[1] * h_arr[1];

	for (i = 2; i < c_arr->Length; i = i + 2)
	{
		product[i] = c_arr[i] * h_arr[i] - c_arr[i + 1] * h_arr[i + 1];

		product[i + 1] = c_arr[i] * h_arr[i + 1] + c_arr[i + 1] * h_arr[i];

		//printf("product[%d] = %f\n", i, product[i]);
	}

	product2[0] = product[0] * n_arr[0]; product2[1] = product[1] * n_arr[1];

	for (i = 2; i < c_arr->Length; i = i + 2)
	{
		product2[i] = product[i] * n_arr[i] - product[i + 1] * n_arr[i + 1];

		product2[i + 1] = product[i] * n_arr[i + 1] + product[i + 1] * n_arr[i];

		//printf("product2[%d] = %f %f\n", i, product2[i], n_arr[i]);
	}

	for (i = 0; i < c_arr->Length; i++)
	{
		product[i] = product2[i];
	}

	product2[0] = product[0] * o_arr[0]; product2[1] = product[1] * o_arr[1];

	for (i = 2; i < c_arr->Length; i = i + 2)
	{
		product2[i] = product[i] * o_arr[i] - product[i + 1] * o_arr[i + 1];

		product2[i + 1] = product[i] * o_arr[i + 1] + product[i + 1] * o_arr[i];

		//printf("product2[%d] = %f %f\n", i, product2[i], o_arr[i]);
	}


	//Include the S
	if (nS >= 1)
	{
		for (i = 0; i < c_arr->Length; i++)
		{
			product[i] = product2[i];
		}

		product2[0] = product[0] * s_arr[0]; product2[1] = product[1] * s_arr[1];

		for (i = 2; i < c_arr->Length; i = i + 2)
		{
			product2[i] = product[i] * s_arr[i] - product[i + 1] * s_arr[i + 1];

			product2[i + 1] = product[i] * s_arr[i + 1] + product[i + 1] * s_arr[i];

			//printf("product2[%d] = %f %f\n", i, product2[i], s_arr[i]);
		}
	}
	

	//inverse fft for the convolution
	realft(product2, -1);

	dtemp = 0.0;

	for (i = 0; i < c_arr->Length; i++)
	{
		dtemp = dtemp + product2[i];
	}

	for (i = 0; i < fIsotopes->Length && i < c_arr->Length; i++)
	{
		fIsotopes[i] = product2[i] / dtemp;
	}

	delete(h_arr); delete (n_arr); delete (c_arr); delete(o_arr); delete (s_arr);

	delete product, product2;

	delete (final_isotope);

	return 0.f;
}

/*
*  realft - fft of a real valued, one-dimensional array
*/
void realft(array <double> ^data, const int isign) {
	int i,i1,i2,i3,i4, n;
	
	//n=data.size();

	n = data->Length;

	double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp;
	
	double theta=3.141592653589793238/double(n>>1);
	
	if (isign == 1) {
		c2 = -0.5;
		four1(data,1);
	} else {
		c2=0.5;
		theta = -theta;
	}

	wtemp=sin(0.5*theta);
	
	wpr = -2.0*wtemp*wtemp;
	
	wpi=sin(theta);
	
	wr=1.0+wpr;
	
	wi=wpi;

	for (i=1;i<(n>>2);i++) {
		i2=1+(i1=i+i);
		i4=1+(i3=n-i1);
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r= -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[0] = (h1r=data[0])+data[1];
		data[1] = h1r-data[1];
	} else {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(data,-1);
	}
}

/*
*   fft routine from numerical recipes
*
*/
void four1(array <double> ^data, const int n, const int isign) 
{
	int nn,mmax,m,j,istep,i;

	double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

	double dtemp;

	if (n<2 || n&(n-1))
    {
		printf("n must be power of 2 in four1");
	}

	nn = n << 1;

	j = 1;
	
	for (i=1;i<nn;i+=2) 
	{
		if (j > i) 
		{
			//SWAP(data[j-1],data[i-1]);

			dtemp = data[j-1];

			data[j-1] = data[i-1];

			data[i-1] = dtemp;

			dtemp = data[j];

			data[j] = data[i];

			data[i] = dtemp;

			//SWAP(data[j],data[i]);
		}
	
		m=n;
		
		while (m >= 2 && j > m) 
		{
			j -= m;
			m >>= 1;
		}
		
		j += m;
	}
	
	mmax=2;
	
	while (nn > mmax) 
	{
		istep=mmax << 1;
	
		theta=isign*(6.28318530717959/mmax);
		
		wtemp=sin(0.5*theta);
		
		wpr = -2.0*wtemp*wtemp;
		
		wpi=sin(theta);
		
		wr=1.0;
		
		wi=0.0;
		
		for (m=1;m<mmax;m+=2) 
		{
			for (i=m;i<=nn;i+=istep) 
			{
				j=i+mmax;
				
				tempr=wr*data[j-1]-wi*data[j];
				
				tempi=wr*data[j]+wi*data[j-1];
				
				data[j-1]=data[i-1]-tempr;
				
				data[j]=data[i]-tempi;
				
				data[i-1] += tempr;
				
				data[i] += tempi;
			}
			
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			
			wi=wi*wpr+wtemp*wpi+wi;
		}
		
		mmax=istep;
	}
}
/*
*  overloaded call to the four1
*/
void four1(array <double> ^data, const int isign) 
{
	four1(data, data->Length/2, isign);
}

/*
*
*  The entry function to compute the theoretical isotopes - all routines start here.
*  A call to this function will trigger the isotope distributions
*
*/
float TheoreticalIsotopeCalculator::StartToComputeSequenceIsotopes(String ^strSequence, array <float> ^fIsotopes)
{
	int nHydrogen, nCarbon, nOxygen, nNitrogen, nSulfur, nPhosphor;

	nHydrogen = nCarbon = nOxygen = nNitrogen = nSulfur = nPhosphor = 0;

	bool bChoose = true;
	//compute the number of each atom type in the compound

	//printf("This is the Peptide: %s\n", strSequence);

	AtomicComposition(strSequence, &nHydrogen, &nCarbon, &nNitrogen, &nOxygen,  &nSulfur, &nPhosphor);

	szPeptideInCalculator = strSequence;

	if (!bChoose)
	{
		Caller_StraightConvolution(fIsotopes, nHydrogen, nCarbon, nNitrogen, nOxygen, nSulfur);

		//if(false)
		for (int i = 0; i < fIsotopes->Length; i++)
		{
			//printf("fIsotopes[%d] = %10.7f\n", i, fIsotopes[i]);

			fIsotopes[i] = 0.0;
		}
	}
	
	//Use FFT to compute the isotope distribution

	if (bChoose)
	{
		FourierIsotopes(fIsotopes, nHydrogen, nCarbon, nNitrogen, nOxygen, nSulfur);

		//FourierIsotopes_Temp(fIsotopes, nHydrogen, nCarbon, nNitrogen, nOxygen, nSulfur);
	}

	
	//printf("Sequence Mass = %20.8f\n", AminoAcid_Sequence_Mass(strSequence));

	return 0.;
}



/*
*
*  this function reads a sequence of amino acids (a peptide)
*  assign the number of C, H, N, S, P atoms to class variables
*  nC, nH, nN, nO, nS, nP - which are publicly available
*  this function does not have the Exchangeable hydrogens
*
*/
int TheoreticalIsotopeCalculator::AtomicComposition(String ^sSequence, int *nHydrogens, int *nCarbons, int *nNitrogens, 
	                         int *nOxygens, int *nSulfur, int *nPhosphor)
{
	int iiC, iiN, iiH, iiO, iiS, iiP;

	int iLength, i, j; 
	
	char cTemp, cTemp2;

	int ElemenMatrix[256][5];

	char szSequence[1024];

	szSequence[0] = '\0';

	for(i = 0; i < sSequence->Length; i++)
	{
		szSequence[i] = sSequence[i];
	}

	szSequence[i] = '\0';


	for(i = 0; i < 256; i++)
	{
		for(j= 0; j < 5; j++)
		{
			ElemenMatrix[i][j] = 0;
		}
	}

	/* inittiate #1 is C, #2 is H, #3 is N, #4 is O, #5 is S*/
	ElemenMatrix[(int)'A'][0] = 3; ElemenMatrix['A'][1] = 5; 
	
	ElemenMatrix['A'][2] = 1; ElemenMatrix['A'][3] = 1;
	
	/*Glycine */
	ElemenMatrix['G'][0] = 2; ElemenMatrix['G'][1] = 3; 
	
	ElemenMatrix['G'][2] = 1; ElemenMatrix['G'][3] = 1;

	/* Serine */
	ElemenMatrix['S'][0] = 3; ElemenMatrix['S'][1] = 5; 
	
	ElemenMatrix['S'][2] = 1; ElemenMatrix['S'][3] = 2;

	/*Proline */
	ElemenMatrix['P'][0] = 5; ElemenMatrix['P'][1] = 7; 
	
	ElemenMatrix['P'][2] = 1; ElemenMatrix['P'][3] = 1;

	/* Valine */
	ElemenMatrix['V'][0] = 5; ElemenMatrix['V'][1] = 9; 
	
	ElemenMatrix['V'][2] = 1; ElemenMatrix['V'][3] = 1;


	/*Threonine */
	ElemenMatrix['T'][0] = 4; ElemenMatrix['T'][1] = 7; 
	
	ElemenMatrix['T'][2] = 1; ElemenMatrix['T'][3] = 2;

	/*Cystein */
	ElemenMatrix['C'][0] = 3; ElemenMatrix['C'][1] = 5; 
	
	ElemenMatrix['C'][2] = 1; ElemenMatrix['C'][3] = 1;  
	
	ElemenMatrix['C'][4] = 1;

	/* Leucine */
	ElemenMatrix['L'][0] = 6; ElemenMatrix['L'][1] = 11; 
	
	ElemenMatrix['L'][2] = 1; ElemenMatrix['L'][3] = 1;

	/* IsoLeucine */
	ElemenMatrix['I'][0] = 6; ElemenMatrix['I'][1] = 11; 
	
	ElemenMatrix['I'][2] = 1; ElemenMatrix['I'][3] = 1;

	/* Asparagine */
	ElemenMatrix['N'][0] = 4; ElemenMatrix['N'][1] = 6; 
	
	ElemenMatrix['N'][2] = 2; ElemenMatrix['N'][3] = 2;

	/* Aspartic Acid */

	ElemenMatrix['D'][0] = 4; ElemenMatrix['D'][1] = 5; 
	
	ElemenMatrix['D'][2] = 1; ElemenMatrix['D'][3] = 3;

	/* Glutamine */

	ElemenMatrix['Q'][0] = 5; ElemenMatrix['Q'][1] = 8; 
	
	ElemenMatrix['Q'][2] = 2; ElemenMatrix['Q'][3] = 2;

	/* Lysine */
	ElemenMatrix['K'][0] = 6; ElemenMatrix['K'][1] = 12; 
	
	ElemenMatrix['K'][2] = 2; ElemenMatrix['K'][3] = 1;


	/* Glutamic Acid */
	ElemenMatrix['E'][0] = 5; ElemenMatrix['E'][1] = 7; 
	
	ElemenMatrix['E'][2] = 1; ElemenMatrix['E'][3] = 3;

	/* Methinine */
	ElemenMatrix['M'][0] = 5; ElemenMatrix['M'][1] = 9; 
	
	ElemenMatrix['M'][2] = 1; ElemenMatrix['M'][3] = 1; 
	
	ElemenMatrix['M'][4] = 1;


	/* Histidine */
	ElemenMatrix['H'][0] = 6; ElemenMatrix['H'][1] = 7; 
	
	ElemenMatrix['H'][2] = 3; ElemenMatrix['H'][3] = 1; 

	/* Phenylalanine - mass is odd because odd number of N */
	ElemenMatrix['F'][0] = 9; ElemenMatrix['F'][1] = 9; 
	
	ElemenMatrix['F'][2] = 1; ElemenMatrix['F'][3] = 1; 

	/* Arginine - mass is even because even number of N */
	ElemenMatrix['R'][0] = 6; ElemenMatrix['R'][1] = 12; 
	
	ElemenMatrix['R'][2] = 4; ElemenMatrix['R'][3] = 1; 

	/* Tyrosine - mass is odd because odd number of N */
	ElemenMatrix['Y'][0] = 9; ElemenMatrix['Y'][1] = 9; 
	
	ElemenMatrix['Y'][2] = 1; ElemenMatrix['Y'][3] = 2;

	/* Tryptophan - mass is even because even number of N */
	ElemenMatrix['W'][0] = 11; ElemenMatrix['W'][1] = 10; 
	
	ElemenMatrix['W'][2] = 2; ElemenMatrix['W'][3] = 1; 

	//for modifications you will need to adjust from the original
	// make the upper and low case amino acid codes to have the same
	// number of H,C, N, O, S
	for(i = 65; i <= 87; i++)
	{
		for(j=0; j <= 4; j++)
		{
			ElemenMatrix[i+32][j] = ElemenMatrix[i][j];
		}
	}

	iiC = iiN = iiH = iiO = iiS = iiP = 0;

	//printf("%d %d %d %d %d %d\n", (int)'A', (int)'W', (int)'a', (int)'w', ElemenMatrix[65][0], ElemenMatrix[97][0]);

	iLength = strlen(szSequence);

	for(i = 0; i < iLength; i++)
	{
		j = (int) szSequence[i];


		if(i == iLength - 1)
		{
			if(szSequence[i] == '[')
			{
				continue;
			}
		}

		if(szSequence[i] == 'm')
		{

		}
		else if(szSequence[i] == 'c')
		{

		}
		else if(szSequence[i] == 'k')
		{

		}
		else if(j < 65 || j > 90)
		{
			printf("Cannot recognize the character %c in %s\n",
				szSequence[i], szSequence);

			puts("Unknown, or unexpected character... will exit");

			exit (1);
		}
		
		iiC += ElemenMatrix[szSequence[i]][0];

		iiH += ElemenMatrix[szSequence[i]][1];

		iiN += ElemenMatrix[szSequence[i]][2];

		iiO += ElemenMatrix[szSequence[i]][3];

		iiS += ElemenMatrix[szSequence[i]][4];

		//look for the modifications
		if(i < iLength - 1)
		{
			cTemp = szSequence[i];

			if(cTemp == 'S' || cTemp == 'T' ||
				cTemp == 'Y' || cTemp == 'C')
			{
				cTemp2 = szSequence[i+1];

				if(cTemp2 == '$' || cTemp2 == '#' ||
					cTemp2 == '*' || cTemp2 == '@')
					/* this site is phosphorylated
					* update the elemental composition
					* correspondingly, and increase the
					* index running over the sequence sites
					*/
				{
					iiO += 3;

					iiH += 1;

					iiP += 1;
					
					i++;
				}
			}
			else if(cTemp == 'M')
			{
				cTemp2 = szSequence[i+1];

				if(cTemp2 == '$' || cTemp2 == '#' ||
					cTemp2 == '*')
					/* this site is phosphorylated
					* update the elemental composition
					* correspondingly, and increase the
					* index running over the sequence sites
					*/
				{
					iiO += 1;
					
					i++;
				}
			}
			else if(cTemp == 'K')
			{
				cTemp2 = szSequence[i+1];

				if(cTemp2 == '$' || cTemp2 == '#' ||
					cTemp2 == '*')
					/* this site is phosphorylated
					* update the elemental composition
					* correspondingly, and increase the
					* index running over the sequence sites
					*/
				{
					iiO += 1;
					
					i++;
				}
			}
		}
	}

	(*nCarbons) = iiC; 
	
	/* Correct for N- and C- termini*/
	(*nHydrogens) = iiH + 3; 
	
	(*nNitrogens) = iiN; 
	
	/* Correct for N- and C- termini*/
	(*nOxygens) = iiO + 1; 
	
	(*nSulfur) = iiS; 

	(*nPhosphor) = iiP;

	return 0;
}
/*
*  This method looks at the time course information to determine
*   if there are replicates, and assigns the unique time point
*   values to fUniqueTimePoints
*
*/
bool ProteinRateConstant::ProcessTimePoints(array <double> ^fTimePoints)
{
	int i = 0, nUniqTimePoints;

	bool bRepls = false;

	nUniqTimePoints = 1;

	//determine if there are any replicates
	for (i = 1; i < fTimePoints->Length; i++)
	{
		if (fabs(fTimePoints[i]/fScaleTime - fTimePoints[i - 1]/fScaleTime) <= 0.0001)
		{
			bRepls = true;
		}
		else
		{
			nUniqTimePoints++;
		}
	}

	if (bRepls)
	{
		fUniqueTimePoints = gcnew array <double> (nUniqTimePoints);

		aiReplicateStructure = gcnew array <int>(nUniqTimePoints);

		fUniqueTimePoints[0] = fTimePoints[0];

		aiReplicateStructure[0] = 1;

		nUniqTimePoints = 1;

		for (i = 1; i < fTimePoints->Length; i++)
		{
			if (fabs(fTimePoints[i] / fScaleTime - fTimePoints[i - 1] / fScaleTime) > 0.0001)
			{
				fUniqueTimePoints[nUniqTimePoints] = fTimePoints[i];

				nUniqTimePoints++;

				aiReplicateStructure[nUniqTimePoints - 1] = 1;
			}
			else
			{
				aiReplicateStructure[nUniqTimePoints-1] = aiReplicateStructure[nUniqTimePoints - 1] + 1;
			}
		}
	}
	else   //no replicates, single experiment per time point
	{
		fUniqueTimePoints = gcnew array <double>(fTimeCourse->Length);

		aiReplicateStructure = gcnew array <int>(fTimeCourse->Length);

		for (i = 0; i < fTimePoints->Length; i++)
		{
			fUniqueTimePoints[i] = fTimePoints[i];

			aiReplicateStructure[i] = 1;
		}
	}

	/*for (i = 0; i < fUniqueTimePoints->Length; i++)
	{
		printf("Uniq Time Poinst = t[%d] = %10.5f, #experiments %d\n", i, 
			fUniqueTimePoints[i] / fScaleTime, aiReplicateStructure[i]);
	}  */
	
	return bRepls;
}
/*
*
*  This method processes Isotopes in the case of replicates
*  It will use optimum data points - only values where the
*   data has actually been observed
*   It use the ProteinRateConstant class variable fUniqueTimePoints
*
*   Uses ProteinRateConstant class variables - fTimeCourse
*   It will use the information that was generated and stored by 
 *   ProcessTimePoints method, which identified the dUniqueTimes (unique time points), and
 *   aiReplicates (number of replicates at each time point);
 *
 *   First it will remove from the data the time points for which there is no experimental
 *   data!!
 *
 *  The function return 1, if after the filtering, 
 *   there is a one experimental data point for each time point
*/
int ProteinRateConstant::AnalyzeReplicates(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters, 
		List <float> ^fIonScores, array <float> ^fTheoreticalIsotopes, array <double> ^TimesForThisPeptide,
		array <float> ^fI0Ions)
{
	int i, j, k, l, m, nPassed;

	float ftemp, fI0, faverage;

	if (TimeCourseIsotopeClusters->Count != fIonScores->Count)
	{
		printf("RateConstat: Problem in AnalyzeReplicates with Scores and Clusters\n");

		exit(1);
	}

	array <bool> ^bFinal = gcnew array <bool>(fIonScores->Count);

	array <float> ^fTempI0 = gcnew array <float>(fIonScores->Count);

	fI0 = ftemp = 0.0;

	for (i = 0; i < fTheoreticalIsotopes->Length; i++)
	{
		ftemp = ftemp + fTheoreticalIsotopes[i];
	}

	if (ftemp < 0.001)
	{
		printf("Theoretical Isotopes aer zero??\n");

		exit(1);
	}

	fI0 = fTheoreticalIsotopes[0] / ftemp;


	//first filtering - cross out the data points that have 0 intensity, empty (for this peptide) experiments
	for (i = 0; i < fIonScores->Count; i++)
	{
		bFinal[i] = false;

		fTempI0[i] = 0.0;

		ftemp = 0.0;

		for (j = 0; j < TimeCourseIsotopeClusters[i]->fIsotopeCluster->Length; j++)
		{
			ftemp = ftemp + TimeCourseIsotopeClusters[i]->fIsotopeCluster[j];
		}
			
		if (ftemp < 1.0)     //first check
		{
			printf("No Signal %f\n", fTimeCourse[i]/fScaleTime);
		}
		else 
		{
			fTempI0[i] = TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp;

			//adjust for the differences in the body water enrichments at different time points

			if (fBodyWaterEnrichment[i] > 0.00001)
			{
				fTempI0[i] = fTempI0[i] * fBodyWaterEnrichment[fBodyWaterEnrichment->Length - 1] / fBodyWaterEnrichment[i];

				//fTempI0[i] = fTempI0[i] * pow(1 - fBodyWaterEnrichment[fBodyWaterEnrichment->Length - 1], nExHydrogen) /
					//pow(1 - fBodyWaterEnrichment[i], nExHydrogen);
			}

			//printf("Comparison with Original: Time = %10f I0 = %f  %f\n", fTimeCourse[i] / fScaleTime, fTempI0[i],
				//TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp);

			bFinal[i] = true;
		}
	}

	double dIsotopeWeihgted, d0sums;

	//second filtering: Averaged of Replicates is within 10% of each measurement, keep just the average for that time point 
	for (i = 0; i < fIonScores->Count; i++)
	{
		if (bFinal[i])
		{
			faverage = 0.0;


			dIsotopeWeihgted = d0sums = 0;

			k = 0;

			for (j = i; j < fIonScores->Count; j++)
			{
				if (Math::Abs(fTimeCourse[i] - fTimeCourse[j]) < 0.00001)
				{
					faverage = faverage + fTempI0[j];

					dIsotopeWeihgted = dIsotopeWeihgted + fTempI0[j] * TimeCourseIsotopeClusters[j]->fIsotopeCluster[0];

					k = k + 1;
				}
				else
				{
					break;
				}
			}

			if (k >= 2)    //there are replicates
			{

				faverage = faverage / double(k);

				k = 1;

				bool bAverageIsGoodApproximation = true;

				for (j = i; j < fIonScores->Count; j++)
				{
					if (Math::Abs(fTimeCourse[i] - fTimeCourse[j]) < 0.00001)
					{
						if (Math::Abs(faverage - fTempI0[j]) > 0.1*faverage)
						{
							bAverageIsGoodApproximation = false;

							break;
						}
					}
					else
					{
						break;
					}
					
				}

				if (bAverageIsGoodApproximation)   //replace the value by avearged value and keep only that value for the replicates
				{
					fTempI0[i] = faverage;

					fTempI0[i] = dIsotopeWeihgted / d0sums;

					for (j = i + 1; j < fIonScores->Count; j++)
					{
						if (Math::Abs(fTimeCourse[i] - fTimeCourse[j]) < 0.00001)
						{
							bFinal[j] = false;
						}
						else
						{
							break;
						}
					}
				}   //replace the value by avearged value ...
			}
		}
	}

	//third filtering: For the replicate points that are not filtering by averaging, remove the values that do not make the data monotonic 
	for (i = 0; i < fIonScores->Count; i++)
	{
		if (bFinal[i])
		{
			k = 0;

			for (j = i; j < fIonScores->Count; j++)
			{
				if (Math::Abs(fTimeCourse[i] - fTimeCourse[j]) < 0.00001)
				{
					k = k + 1;
				}
				else
				{
					break;
				}
			}

			// the data points from i to (j - 1) are replicates 
			if (k >= 2)    //there are replicates
			{

			}

		}
	}


	//final result

	//Check to see if each time point now has a single experimental
	// The total number of points could be equal or less than the number of uniquetime points
	//important thing is to have a single experimental value for time point

	bool bSingleDataAtEveryExperiment = true;
	
	bool bThisExperiment = false;

	int nTrueTimePointsForThisPeptide = 0;

	for (j = 0; j < fUniqueTimePoints->Length; j++)
	{
		k = 0;

		for (i = 0; i < fTimeCourse->Length; i++)
		{
			if (bFinal[i])
			{
				if (Math::Abs((fUniqueTimePoints[j] - fTimeCourse[i]) / fScaleTime) < 0.0001)
				{
					k = k + 1;
				}
				//printf("Time = %10f I0 = %f %f\n", fTimeCourse[i] / fScaleTime, fTempI0[i], fIonScores[i]);
			}

		}

		if (k > 1)
		{
			bSingleDataAtEveryExperiment = false;
		}


		if (k >= 1)
		{
			nTrueTimePointsForThisPeptide = nTrueTimePointsForThisPeptide + 1;
		}
	}


	//if only one data point per time send the data to BFGS
	// for good spectra it stops here.

	if (bSingleDataAtEveryExperiment)
	{
		TimesForThisPeptide = gcnew array <double>(nTrueTimePointsForThisPeptide);

		fI0Ions             = gcnew array <float>(nTrueTimePointsForThisPeptide);

		printf("This Peptide is ready for BFGS, nTrueNumberof Points = %d\n", nTrueTimePointsForThisPeptide);

		if (nTrueTimePointsForThisPeptide == fUniqueTimePoints->Length) //all time points are good for this peptide
		{
			TimesForThisPeptide = fUniqueTimePoints;
		}

		return 1;
	}
	

	double replicate[100][100], mSum[100][100];

	int nMonotone[100][100], ix[100][100], iy[100][100];


	for (i = 0; i < 100; i++)
	{
		for (j = 0; j < 100; j++)
		{
			replicate[i][j] = 0.0;

			nMonotone[i][j] = 0;

			mSum[i][j] = -1.0;
			
			ix[i][j] = 0;

			iy[i][j] = 0;
		}
	}

	//time points
	replicate[0][0] = fTempI0[0];

	k = 0;     //k will record the switch to new time point for each group of replicates

	for (i = 1; i < 100 && i < bFinal->Length; i++)
	{
		if (i >= 1 && bFinal[i])
		{
			if (Math::Abs(fTimeCourse[i - 1] - fTimeCourse[i]) < 0.00001)
			{
				//printf("This time is registered %f %f, k = %d, l = %d\n", fTimeCourse[i - 1]/fScaleTime, fTimeCourse[i]/fScaleTime, k, l);

				replicate[k][l] = fTempI0[i];

				l = l + 1;
			} 
			else
			{
				k = k + 1;

				replicate[k][0] = fTempI0[i];
				
				l = 1;

				//printf("A new time points %f, k = %d\n", fTimeCourse[i]/ fScaleTime, k);
			}
		}
	}

	for (i = 0; i < bFinal->Length; i++)
	{
		for (j = 0; j < bFinal->Length; j++)
		{
			if (replicate[i][j] > 0.0)
			{
				//printf("REplicates %f  ", replicate[i][j]);
			}
		}

	}

	//puts("");

	double fsum;

	int isum = 0;

	int i0, i1;


	//initiate the 0th row
	for (j = 0; j < bFinal->Length; j++)
	{
		mSum[5][j] = 0.0;
	}



	for (i = bFinal->Length - 1; i >= 0; i--)
	{
		for (m = 1; m < 10; m++)
		{
			i0 = i + m;    //step by step, check establish the best monotonic time path

			for (i1 = 0; i1 < bFinal->Length; i1++)
			{
				//for (j= 0; j < bFinal->Length; j++)  //for every value of time point replicates

				for (j = bFinal->Length - 1; j >= 0; j--)  //for every value of time point replicates
				{
					l = i1 + j;
					
					if (i0 < bFinal->Length && i1 < bFinal->Length)
					{
						if (replicate[i0][i1] > 0.0001 && replicate[i][j] > 0.0001)   //only from non-zero cells to non-zero cells
						{

							if (replicate[i0][i1] < replicate[i][j])   //decreasing I0
							{
								isum = nMonotone[i0][i1] + 1;

								//printf("replicate[%d][%d] = %f\n", i, i, replicate[i][j]);

								//printf("replicate[%d][%d] = %f\n", i0, i1, replicate[i0][i1]);
							}
							else                                       //increasing I0
							{
								//puts("CCC");

								isum = nMonotone[i0][i1] - 1;
							}

							fsum = mSum[i0][i1] + replicate[i0][i1] - replicate[i][j];

							if (fsum > mSum[i][j])
							{
								mSum[i][j] = fsum;

								ix[i][j] = i0;

								iy[i][j] = i1;
							}
						}
					}
				}

			}
		}
			

	}

	if(false)
	for (i = 0; i < bFinal->Length; i++)
	{
		for (j = 1; j < bFinal->Length; j++)
		{
			if(ix[i][j] != 0 && iy[i][j] != 0)
			printf(" ix[%d][%d] = %d. y[%d][%d] = %d\n", i, j, ix[i][j], i, j, iy[i][j]);
		}
	}

	for (i = 0; i < bFinal->Length; i++)
	{
		for (j = 0; j < 2; j++)
		{
			//printf("mSum[%d][%d] = %f\n", i, j, mSum[i][j]);
		}
	}

	return 0;
}

/*
	The function will remove the empty data points, 0, and average the rest of the replicats
	* the new data f0..., and new time points f... are saved.

*/

int ProteinRateConstant::AverageReplicates(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters,
	List <float> ^fIonScores, array <float> ^fTheoreticalIsotopes, array <double> ^TimesForThisPeptide,
	array <float> ^fI0Ions, int *nTrueIons)
{
	int i, j, k, l, m;

	float ftemp, fI0, faverage;

	if (TimeCourseIsotopeClusters->Count != fIonScores->Count)
	{
		printf("RateConstat: Problem in AverageReplicates with Scores and Clusters\n");

		exit(1);
	}

	array <bool> ^bFinal = gcnew array <bool>(fIonScores->Count);

	array <float> ^fTempI0 = gcnew array <float>(fIonScores->Count);

	fI0 = ftemp = 0.0;

	for (i = 0; i < fTheoreticalIsotopes->Length; i++)
	{
		ftemp = ftemp + fTheoreticalIsotopes[i];
	}

	if (ftemp < 0.001)
	{
		printf("Theoretical Isotopes aer zero??\n");

		exit(1);
	}

	fI0 = fTheoreticalIsotopes[0] / ftemp;


	//first filtering - remove  the data points that have 0 intensity, empty (for this peptide) experiments
	// and the time points that have ion score less than 1. (not identified)
	for (i = 0; i < fIonScores->Count; i++)
	{
		bFinal[i] = false;

		fTempI0[i] = 0.0;

		ftemp = 0.0;

		for (j = 0; j < TimeCourseIsotopeClusters[i]->fIsotopeCluster->Length; j++)
		{
			ftemp = ftemp + TimeCourseIsotopeClusters[i]->fIsotopeCluster[j];
		}

		//printf("II = %d\n", i);

		if (ftemp < 1.0)     //first check
		{
			printf("No Signal %f\n", fTimeCourse[i] / fScaleTime);
		}
		else
		{
			fTempI0[i] = TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp;

			//adjust for the differences in the body water enrichments at different time points

			if (fBodyWaterEnrichment[i] > 0.00001)
			{
				//fTempI0[i] = fTempI0[i] * (1 - nExHydrogen * fBodyWaterEnrichment[fBodyWaterEnrichment->Length - 1] / 
		        //									( 1 - fBodyWaterEnrichment[i]));

				fTempI0[i] = fTempI0[i] * pow(1 - fBodyWaterEnrichment[fBodyWaterEnrichment->Length - 1], nExHydrogen) /
					pow(1 - fBodyWaterEnrichment[i], nExHydrogen);
				;
			}

			//printf("Comparison with Original: Time = %10f I0 = %f  %f\n", fTimeCourse[i] / fScaleTime, fTempI0[i],
			//TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp);

			if (fIonScores[i] > 10)
			{
				bFinal[i] = true;

				//printf("Comparison with Original: Time = %10f I0 = %f  %f\n", fTimeCourse[i] / fScaleTime, fTempI0[i],
					//TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp);
			}
			//bFinal[i] = true;
		}
	}



	//second filtering: Averaged of Replicates is within 10% of each measurement, keep just the average for that time point 

	double dIsotopeWeihgted, d0sums;

	if(false)
	for (i = 0; i < fIonScores->Count; i++)
	{
		if (bFinal[i] && fTimeCourse[i] / fScaleTime > 0.0001)
		{
			faverage = 0.0;

			dIsotopeWeihgted = d0sums = 0.0;

			k = 0;

			for (j = i; j < fIonScores->Count; j++)
			{
				if (Math::Abs(fTimeCourse[i] - fTimeCourse[j]) < 0.00001)
				{
					faverage = faverage + fTempI0[j];

					dIsotopeWeihgted = dIsotopeWeihgted + fTempI0[j] * TimeCourseIsotopeClusters[j]->fIsotopeCluster[0];

					d0sums = d0sums + TimeCourseIsotopeClusters[j]->fIsotopeCluster[0];
						
					k = k + 1;
				}
				else
				{
					break;
				}
			}

			if (k >= 2)    //there are replicates, k >= 2 because the summation starts with i  !!
			{

				faverage = faverage / double(k);

				k = 1;

				bool bAverageIsGoodApproximation = true;

				for (j = i; j < fIonScores->Count; j++)
				{
					if (Math::Abs((fTimeCourse[i] - fTimeCourse[j]) / fScaleTime) < 0.001)
					{
						if (Math::Abs(faverage - fTempI0[j]) > 0.1*faverage)
						{
							bAverageIsGoodApproximation = false;

							break;
						}
					}
					else
					{
						break;
					}

				}

				if (bAverageIsGoodApproximation)   //replace the value by avearged value and keep only that value for the replicates
				{
					bool bNonzeroScore = false;

					for (j = i; j < fIonScores->Count; j++)
					{
						if (Math::Abs((fTimeCourse[i] - fTimeCourse[j]) / fScaleTime) < 0.001)
						{
							if (bNonzeroScore == false && fIonScores[j] > 1.)
							{
								fTempI0[j] = faverage;

								fTempI0[j] = dIsotopeWeihgted / d0sums;

								bNonzeroScore = true;
							}
							else 
							{
								bFinal[j] = false;
							}
							
						}
						else
						{
							break;
						}
					}
				}   //replace the value by avearged value ...
			}
		}
	}

	// third filtering :: O time point: if the replicates for the 0th time point are still not averaged, keep those that are closest to the
	// to theoretical value

	for (i = 0; i < fIonScores->Count; i++)
	{
		ftemp = 0.0;

		if (bFinal[i] && (fTimeCourse[i] / fScaleTime) < 0.0001 )   //0th time point
		{
			k = 0; l = 0;

			faverage = 0.0;

			for (j = i; j < fIonScores->Count; j++)
			{
				if (fIonScores[j] > 1 && (fTimeCourse[j]/ fScaleTime) < 0.0001)   //0th time point
				{
					if (Math::Abs(fTempI0[j] - fI0) < 0.5* fI0)    //fI0 is the theoretical Isotope's I0 
					{
						faverage = faverage + fTempI0[j];

						k = k + 1;

					}
				}

				ftemp = ftemp + fTempI0[j];

				l = l + 1;
			}


			if (faverage > 0.0001)   //replace the value by avearged value and keep only that value for the replicates
			{
				fTempI0[i] = faverage / (double)k;
			}
			else
			{

				//puts("MAYBE THIS");

				fTempI0[i] = ftemp / (double)l;
			}

			for (j = i + 1; j < fIonScores->Count; j++)
			{
				if (Math::Abs((fTimeCourse[j] - 0.) / fScaleTime) < 0.001)
				{
					bFinal[j] = false;
				}
				else
				{
					break;
				}
			}

			break;     // because it is a single time point

		}

	}


	//fourth filtering: For the replicate points that are not filtered by gentel averaging, remove the values that have 0 ion score
	//except the 0 time point
	
	for (i = 0; i < fIonScores->Count; i++)
	{
		if (bFinal[i] && fTimeCourse[i] > 0.0001)
		{
			k = 0;

			for (j = i + 1; j < fIonScores->Count; j++)
			{
				if (Math::Abs((fTimeCourse[i] - fTimeCourse[j]) / fScaleTime) < 0.001)
				{
					k++;
				}
				else
				{
					break;
				}
			}

			bool bNonZeroScore = false;
			if (k >= 1)      //replicates present, keep only the one that has a non-zero IonScore
			{
				for (j = i; j < fIonScores->Count; j++)
				{
					if (Math::Abs((fTimeCourse[i] - fTimeCourse[j]) / fScaleTime) < 0.001)
					{
						if (fIonScores[j] > 1.0)
						{
							bNonZeroScore = true;

							l = j;

							break;
						}
						else
						{
							break;
						}
					}
					
				}

				if (bNonZeroScore)
				{
					for (j = i; j < fIonScores->Count; j++)
					{
						if (Math::Abs((fTimeCourse[i] - fTimeCourse[j]) / fScaleTime) < 0.001)
						{
							if (j != l)                     //keep just one value for this replicates, the one that does not have 
							{								// 0 IonScore
								bFinal[j] = false;
							}
							else
							{
								break;
							}
						}
					}
				}
			}

		}
	}
	
	//fifth filtering: For the replicate points that are not filtered by gentel averaging +  0 ion scores,
	// do averaging based on the intensityweights

	for (i = 0; i < fIonScores->Count; i++)
	{
		if (bFinal[i]  && fIonScores[i] > 10.)
		{
			k = 0; 
			
			faverage = 0.0;

			d0sums = dIsotopeWeihgted = 0.0;

			for (j = i; j < fIonScores->Count; j++)
			{
				if (Math::Abs(fTimeCourse[i] - fTimeCourse[j]) / fScaleTime < 0.001 && bFinal[j] &&
					fIonScores[j] > 10.)
				{
					faverage = faverage + fTempI0[j];

					dIsotopeWeihgted = dIsotopeWeihgted + fTempI0[j] * TimeCourseIsotopeClusters[j]->fIsotopeCluster[0];
					
					d0sums = d0sums + TimeCourseIsotopeClusters[j]->fIsotopeCluster[0];

					k = k + 1;
				}
				else
				{
					//break;
				}
			}

			if (k > 0)   //this cluster is replicate, do averaging using isotope weights
			{

				fTempI0[i] = faverage / (double)k;

				fTempI0[i] = dIsotopeWeihgted / d0sums;

				//printf("THIS PEPTIDE and time = %f %f %f\n", fTimeCourse[i], faverage, fTempI0[i]);

				for (j = i + 1; j < fIonScores->Count; j++)
				{
					if (Math::Abs((fTimeCourse[i] - fTimeCourse[j]) / fScaleTime) < 0.001)
					{
						bFinal[j] = false;
					}
					else
					{
						break;
					}
				}
			}


		}
	}


	//final result

	// now get the number of time points, the time points array, and I0 array

	j = 0;

	for (i = 0; i < bFinal->Length; i++)
	{
		if (bFinal[i])
		{
			j++;
		}
	}

	j = 0;

	for (i = 0; i < bFinal->Length; i++)
	{
		if (bFinal[i] && fTempI0[i] > 0.00000001)
		{
			TimesForThisPeptide[j] = fTimeCourse[i];

			fI0Ions[j] = fTempI0[i];

			j++;
		}
	}
	
	(*nTrueIons) = j;

	if ((*nTrueIons) <= 2)
	{
		return -1;
	}

	/*for (i = 0; i < (*nTrueIons); i++)
	{
		printf("time = %f, I0 = %f, nTime = %d\n", TimesForThisPeptide[i] /fScaleTime, fI0Ions[i], (*nTrueIons));
	} */

	return 0;
}



/*
*   The function will work with single experiments - no replicates. It will return the value of
*   fI0Ions
*
*/

float ProteinRateConstant::ExperimentalTheoreticalIsotopeAccuracy(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters,
	array <float> ^fTheoreticalIsotopes)
{
	float fAccuracy = -1.0;

	float ftemp, ftemp2, ftemp3;

	int i;

	ftemp = ftemp2 = ftemp3 = 0.0;

	for (i = 0; i < fTheoreticalIsotopes->Length; i++)
	{
		ftemp = ftemp + fTheoreticalIsotopes[i];
	}

	if (ftemp < 0.001)
	{
		printf("Theoretical Isotopes are zero??\n");

		exit(1);
	}

	ftemp2 = ftemp3 = 0.0;

	for (i = 0; i < fTheoreticalIsotopes->Length; i++)
	{
		ftemp2 = ftemp2 + TimeCourseIsotopeClusters[0]->fIsotopeCluster[i]; // dExperimentIsotopes[k, i];
	}

	if (ftemp2 > 0.001)
	{
		ftemp3 = 0.0;

		for (i = 0; i < fTheoreticalIsotopes->Length; i++)
		{
			ftemp3 = ftemp3 +
				Math::Abs(TimeCourseIsotopeClusters[0]->fIsotopeCluster[i] / ftemp2 - fTheoreticalIsotopes[i]/ftemp); // dExperimentIsotopes[k, i];
		}

		fAccuracy = ftemp3;
	}
	else
	{
		fAccuracy = -100.0;
	}


	return fAccuracy;
}

/*
*   The function will work with single experiments - no replicates. It will return the value of  
*   fI0Ions
*
*/

int ProteinRateConstant::NoReplicateExperiments(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters,
	List <float> ^fIonScores, array <float> ^fTheoreticalIsotopes, array <double> ^TimesForThisPeptide,
	array <float> ^fI0Ions, int *nTrueIons)
{
	int i, j, k, l, m;

	float ftemp, fI0, faverage;

	if (TimeCourseIsotopeClusters->Count != fIonScores->Count)
	{
		printf("RateConstant: Problem in AverageReplicates with Scores and Clusters\n");

		exit(1);
	}

	array <bool> ^bFinal = gcnew array <bool>(fIonScores->Count);

	array <float> ^fTempI0 = gcnew array <float>(fIonScores->Count);

	fI0 = ftemp = 0.0;

	//for the 0times isotope cluster compute the match with 
	// with the theoretical isotope

	for (i = 0; i < fTheoreticalIsotopes->Length; i++)
	{
		ftemp = ftemp + fTheoreticalIsotopes[i];
	}

	if (ftemp < 0.001)
	{
		printf("Theoretical Isotopes are zero??\n");

		exit(1);
	}

	fI0 = fTheoreticalIsotopes[0] / ftemp;


	//only one filtering - remove  the data points that have 0 intensity, empty (for this peptide) experiments
	// and the time points that have ion score less than 1. (not identified)
	for (i = 0; i < fIonScores->Count; i++)
	{
		bFinal[i] = false;

		fTempI0[i] = 0.0;

		ftemp = 0.0;

		for (j = 0; j < TimeCourseIsotopeClusters[i]->fIsotopeCluster->Length; j++)
		{
			ftemp = ftemp + TimeCourseIsotopeClusters[i]->fIsotopeCluster[j];
		}

		//printf("II = %d\n", i);

		if (ftemp < 1.0)     //first check
		{
			printf("No Signal %f\n", fTimeCourse[i] / fScaleTime);
		}
		else
		{
			fTempI0[i] = TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp;

			//adjust for the differences in the body water enrichments at different time points

			if (fBodyWaterEnrichment[i] > 0.00001)
			{
				fTempI0[i] = fTempI0[i] * fBodyWaterEnrichment[fBodyWaterEnrichment->Length - 1] / fBodyWaterEnrichment[i];

				//fTempI0[i] = fTempI0[i] * pow(1 - fBodyWaterEnrichment[fBodyWaterEnrichment->Length - 1], nExHydrogen) /
				//pow(1 - fBodyWaterEnrichment[i], nExHydrogen);
			}

			//printf("Comparison with Original: Time = %10f I0 = %f  %f\n", fTimeCourse[i] / fScaleTime, fTempI0[i],
			//TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp);

			if (fIonScores[i] > 10)
			{
				bFinal[i] = true;

				//printf("Comparison with Original: Time = %10f I0 = %f  %f\n", fTimeCourse[i] / fScaleTime, fTempI0[i],
				//TimeCourseIsotopeClusters[i]->fIsotopeCluster[0] / ftemp);
			}
			//bFinal[i] = true;
		}
	}


	//final result

	// now get the number of time points, the time points array, and I0 array

	j = 0;

	for (i = 0; i < bFinal->Length; i++)
	{
		if (bFinal[i])
		{
			j++;
		}
	}

	j = 0;

	for (i = 0; i < bFinal->Length; i++)
	{
		if (bFinal[i] && fTempI0[i] > 0.00000001)
		{
			TimesForThisPeptide[j] = fTimeCourse[i];

			fI0Ions[j] = fTempI0[i];

			j++;
		}
	}

	//printf("J = %d\n", j);

	(*nTrueIons) = j;

	if ((*nTrueIons) <= 2)
	{
		return -1;
	}

	/*for (i = 0; i < (*nTrueIons); i++)
	{
	printf("time = %f, I0 = %f, nTime = %d\n", TimesForThisPeptide[i] /fScaleTime, fI0Ions[i], (*nTrueIons));
	} */

	return 0;
}
/*
    A method to compute the mass of a peptide sequence

*/

double AminoAcid_Sequence_Mass(String ^sSequence)
{
	int i;

	char c[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double aam[20] = { 57.02146372, 71.03711378, 87.03202840, 97.05276384, 99.06841390, 101.04767846, 103.00918451, 113.08406396, 113.08406396, 114.04292744, 115.02694302, 128.05857750, 128.09496300, 129.04259308, 131.04048463, 137.05891186, 147.06841390, 156.10111102, 163.06332852, 186.07931294 };
	//				  g            a            s            p            v            t             c             l             i             n             d             q			 k             e             m             h             f             r             y             w
	char aasymb[20] = { 'g', 'a', 's', 'p', 'v', 't', 'c', 'l', 'i', 'n', 'd', 'q', 'k', 'e', 'm', 'h', 'f', 'r', 'y', 'w' };

	double AA_Masses[128], dMass;

	for (i = 0; i < 128; i++)
	{
		AA_Masses[i] = 0.0;
	}

	dMass = MASS_H2OH;


	AA_Masses['g'] = AA_Masses['G'] = 57.02146372;   AA_Masses['a'] = AA_Masses['A'] = 71.03711378;  AA_Masses['s'] = AA_Masses['S'] = 87.03202840;

	AA_Masses['p'] = AA_Masses['P'] = 97.05276384;   AA_Masses['v'] = AA_Masses['V'] = 99.06841390;  AA_Masses['t'] = AA_Masses['T'] = 101.04767846;

	AA_Masses['c'] = AA_Masses['C'] = 103.00918451;   AA_Masses['i'] = AA_Masses['I'] = 113.08406396;  AA_Masses['l'] = AA_Masses['L'] = 113.08406396;

	AA_Masses['n'] = AA_Masses['N'] = 114.04292744;   AA_Masses['d'] = AA_Masses['D'] = 115.02694302;  AA_Masses['q'] = AA_Masses['Q'] = 128.05857750;

	AA_Masses['k'] = AA_Masses['K'] = 128.09496300;   AA_Masses['e'] = AA_Masses['E'] = 129.04259308;  AA_Masses['m'] = AA_Masses['M'] = 131.04048463;

	AA_Masses['h'] = AA_Masses['H'] = 137.05891186;   AA_Masses['f'] = AA_Masses['F'] = 147.06841390; AA_Masses['r'] = AA_Masses['R'] = 156.10111102;

	AA_Masses['y'] = AA_Masses['Y'] = 163.06332852;   AA_Masses['w'] = AA_Masses['W'] = 186.07931294;

	/*for (i = 0; i < 128; i++)
	{
		if (AA_Masses[i] > 0.001)
			printf("AA_Mass[%c] = %20.8f, i = %d\n", (char)i, AA_Masses[i], i);
	} */

	for (i = 0; i < sSequence->Length; i++)
	{
		dMass = dMass + AA_Masses[sSequence[i]];
	}

	return dMass;
}
