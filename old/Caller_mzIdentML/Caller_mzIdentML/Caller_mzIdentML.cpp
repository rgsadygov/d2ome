// Caller_mzIndentML.cpp : main project file.
//this project incorporates the mzIdentML class 
// which is a dll
/*
*  03.23.2017 to run the program 

*  calller_mzIdentML.exe files.txt 1 2 .
*/

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"
#include "ProteinsExperimentCollectorClasses.h"
#include <vector>
#include <list>
#include <sstream>
#include <iostream>
#include <fstream>



//using namespace std;
using namespace System;
using namespace System::Collections::Generic;
using namespace System::IO;
using namespace System::Globalization;
using namespace pepXML;                                   //for definitions in mzIdentML
using namespace ProteinCollector;                         //for definitions in this project
using namespace IsotopePeaks;                             //for computing isotopes
using namespace PeakDetectionIntegration;                 // peak integration
using namespace NNLS;                                     //for computing NNLS
using namespace RateConstant;                             // rate constant calculations, calls LBFGS


//using namespace System.Linq;
//using namespace System.Xml.Linq;


int Start_Process(array<System::String ^> ^args); // starts the Cleveland processing

void D2O_H2O(ExperimentCollection ^allResults,array<System::String ^> ^args);   //a pipeline to write data out for Cleveland Clinic

int ReadAbundances(array <String ^> ^smzML, ExperimentCollection ^AllExperiments); //reads isotopes from mzML files

void WriteTimeSeries(ExperimentCollection ^allResults,array<System::String ^> ^args);   //writes out the results in "time series" format

void WriteTimeSeriesPeptidesInProtein(ExperimentCollection ^allResults, 
	List  <String ^>  ^allAccessions, List  <bool>  ^CommonAccessions,array<System::String ^> ^args);

void ConsistentProteins(ExperimentCollection ^, int, float, int);

void ConsistentPeptides(ExperimentCollection ^allExperiments, float PeptideScoreThreshold, double dExpectationThreshold,
	int nExperimentThreshold); //PeptideScoreThreshold - threshold score for peptides, 
// nExperimentalThreshold - threshold number of experiments

void WriteQuantResults(ExperimentCollection ^allExperiments,array<System::String ^> ^args);

void WriteAProteinFile(ExperimentCollection ^allExperiments, 
	String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args);

void WriteAProteinFile1(ExperimentCollection ^allExperiments, 
	String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args);//Writing MPE info

void WriteAProteinFile2(ExperimentCollection ^allExperiments, 
	String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args);//Writing NNLS info

void WriteResults(ExperimentCollection ^allResults,array<System::String ^> ^args);

void WriteScansMasses(ExperimentCollection ^allResults);


int find_add_missing_time_points(ExperimentCollection ^allExperiments, int exp_no, int *charge,int *scan, float *retTime, 
	float *fstart, float *fend, double *dSeqMass, double *dSpecMass, double *dMassToCharge, String ^protein, String ^peptide,
	array<double, 2> ^dIsotopes, array<double, 1> ^w, List <int> ^m); //Find missing time points

bool check_matching_modlocations(List<int> ^m1, List<int> ^m2); //Checks each position and returns true/false based on match/mismatch

void WriteNNLSIsotopeDistribution(ExperimentCollection ^allExperiments, 
	String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args); //Writing NNLS

void WriteAProtein_O18(ExperimentCollection ^allExperiments, String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args);

void O18_Quantification(array<System::String ^> ^inputFiles);

int ReadmzID(array<System::String ^> ^args, ExperimentCollection ^allExperiments);

int Parse_InputTextFile(String ^strFile, array <String ^> ^ mzmlFiles, array <String ^> ^ mzidFiles, 
	array <float> ^fTimePoints, array <float> ^fBodyWater);

int NumberOfExperiments(String ^strFile);

int nTheoreticalIsotopes = 6;     //the default number of isotopes to calculate in Theoretical Isotope Generages

typedef struct _charArray
{
	char szChar[2048];
} charArray;

typedef std::vector <charArray> CharArray;

bool bWriteQuantCsv = true;

bool bScansMasses = true;

bool bO18 = false;

double mProton = 1.00727642, mHydrogen = 1.007825035, deltaC13 = 13.003354838 - 12.0;

double mDeuterium = 2.014102, dN14 = 14.003074, dN15 = 15.000109;

double dO16 = 15.994915, dO17 = 16.999132, dO18 = 17.999161;

double dS32 = 31.972071, dS33 = 32.971459, dS34 = 33.967867, dS36 = 35.967081;


//below is the list of the quant parameters and their default values
double ProteinScoreCut = 100, dExpectCut = 0.5, dIonScoreCut = 10.; // cut-off for Protein Scores
int nProteinConsistency = 4;   //number of experiments a protein needs to be observed in, for quantification
int nPeptideConsistency = 4;
float ElutionTimeWindow = 1.0;  // the time (min) 
double dMassAccuracy = 0.05, PeakWidth_M_over_Z = 0.04;
bool Compute_MPE = false;      // computes the MPE
float nHAA[256];               // an array that holds number of exchangeable H atoms for every AA
int Nparameter = 1;

unsigned short int mass_accuracy_unit = 0;    // mass accuracy unites, 0 - stands for ppm, 1 - stands for Da


/*Command line argument parameters */
int Isotope_Deconvolution = 0; //args[1]

std::string output_file_dir="";         //args[3]
/**/

int main(array<System::String ^> ^args)
{
	array<String ^> ^sFiles;

	PeptideHolder ^aPeptide = gcnew PeptideHolder();

	ProteinCollection ^currentExperiment = gcnew ProteinCollection();

	ExperimentCollection ^allExperiments = gcnew ExperimentCollection();

	allExperiments->ExperimentsList = gcnew List<ProteinCollection^>;

	ProteinSet ^aProtein;

	int i, k, l, j;

	i = args->Length;      //Fixing blank space in the output directory

	if (i < 1) {
		printf("please provide an input text file.\nYour command should look like:\n> d2ome.exe file.txt\n");
		return 2;
	}
	else if(i != 4){
		printf("The last 3 arguments are set to default values of 1, 1, and the output directory is the current location\n");
		array<System::String ^> ^newargs = {args[0],"1","1","."};
		args = newargs;
	}


	if(i > 4)
	{
		for(j=4; j < i; j++)
		{
			args[3] +=" "+args[j];
		}

	}

	Directory::SetCurrentDirectory(args[3]);

	Read_Params(args);

	//go to Cleveland Clinic Processing and exit
	if(bWriteQuantCsv)
	{
		//printf("\n\nPlease, cite this software as: \n\n");

		/*printf("Assessment of Cardiac Proteome Dynamics with Heavy Water: SlowerProtein \n");
		printf("Synthesis Rates in Interfibrillar than Subsarcolemmal Mitochondria, \n");
		printf("T. Kasumov, E.R. Dabkowski, K.C. Shekar, L. Li, R.F. Ribeiro Jr, \n");
		printf("K. Walsh, S.F. Previs, R.G. Sadygov, B. Willard, W.C. Stanley,\n");
		printf("Am J Physiol Heart Circ Physiol. 2013 May;304(9):H1201-14 \n\n");*/

		//System::Threading::Thread::Sleep(8000);

		Start_Process(args);

		/*printf("\n\nPlease, cite this software as: \n\n\n");

		printf("Assessment of Cardiac Proteome Dynamics with Heavy Water: SlowerProtein \n");
		printf("Synthesis Rates in Interfibrillar than Subsarcolemmal Mitochondria, \n");
		printf("T. Kasumov, E. R. Dabkowski, K. C. Shekar, L. Li, \n");
		printf("R. F. Ribeiro Jr, K. Walsh, S. F. Previs, R. G. Sadygov, \n");
		printf("B. Willard, W. C. Stanley, Am J Physiol Heart Circ Physiol. 2013 May;304(9):H1201-14 \n\n\n");*/

		return 0;
	}
	else if(bO18)
	{
		printf("O18 quantification\n");

		O18_Quantification(args);

		return 0;

	}



	if(args[0]->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase) == -1)
	{
		printf("Need an mzid file to run the program\n");

		exit (1);
	}

	if(args[0]->IndexOf("*.mzid", StringComparison::OrdinalIgnoreCase) != -1)
	{
		sFiles = System::IO::Directory::GetFiles("\\.");
	}
	else
	{
		//reads the mzML files from here
		sFiles = gcnew array<String ^> (1);

		sFiles[0] = args[0];
	}

	k = 0;

	for(l=0; l < sFiles->Length; l++)
	{
		if(sFiles[l]->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase) != -1)
		{
			k++;

			// process/read the mzid file
			mzIdentML ^mzI = gcnew mzIdentML(sFiles[l]);

			mzI->ReadProteins_mzIdent();

			//now create a storage for this experiment
			currentExperiment = gcnew ProteinCollection();

			currentExperiment->ProteinsList = gcnew List<ProteinSet^>;

			//count goes over every protein entry read in from mzid
			for(i=0; i < mzI->ProteinResultList->Count; i++)
			{
				//copy protein data
				aProtein = gcnew ProteinSet();

				aProtein->Peptides = gcnew List <PeptideHolder ^>;


				aProtein->accession          = mzI->ProteinResultList[i]->accession;
				aProtein->description        = mzI->ProteinResultList[i]->description;
				aProtein->ProteinScore       = mzI->ProteinResultList[i]->ProteinScore;
				aProtein->SeqCoverage        = mzI->ProteinResultList[i]->SeqCoverage;
				aProtein->nSpectralCount     = mzI->ProteinResultList[i]->nSpectralCount;
				aProtein->nDistinctSequences = mzI->ProteinResultList[i]->nDistinctSequences;
				aProtein->nSeqLength         = mzI->ProteinResultList[i]->nSeqLength;
				aProtein->ProteinMass        = mzI->ProteinResultList[i]->ProteinMass;

				for(int i0 = 0; i0 < mzI->ProteinResultList[i]->PeptideList->Count; i0++)
				{
					aPeptide = gcnew PeptideHolder();

					aPeptide->dIonscore = mzI->ProteinResultList[i]->PeptideList[i0]->dIonscore;
					aPeptide->bUniquePeptide = mzI->ProteinResultList[i]->PeptideList[i0]->bUniquePeptide;
					aPeptide->dRetTime  = mzI->ProteinResultList[i]->PeptideList[i0]->dRetTime;
					aPeptide->Peptide   = mzI->ProteinResultList[i]->PeptideList[i0]->Peptide;
					aPeptide->nScan     = mzI->ProteinResultList[i]->PeptideList[i0]->nScan;
					aPeptide->nCharge   = mzI->ProteinResultList[i]->PeptideList[i0]->nCharge;
					aPeptide->SeqMass  = mzI->ProteinResultList[i]->PeptideList[i0]->SeqMass;
					aPeptide->SpecMass = mzI->ProteinResultList[i]->PeptideList[i0]->SpecMass;
					aPeptide->dExpect = mzI->ProteinResultList[i]->PeptideList[i0]->dExpect;
					aPeptide->dIsotopes = gcnew array <double, 2>(2, 6);
					aPeptide->w = gcnew array <double, 1>(12);

					aPeptide->Protein = mzI->ProteinResultList[i]->PeptideList[i0]->Protein;

					//printf("Protein name : %s\n", aPeptide->Protein);

					aPeptide->dModMasses = gcnew List <double>;

					aPeptide->ModLocations = gcnew List <int>;

					aPeptide->dModMasses = mzI->ProteinResultList[i]->PeptideList[i0]->dModMasses;

					aPeptide->ModLocations = mzI->ProteinResultList[i]->PeptideList[i0]->ModLocations; 

					aProtein->Peptides->Add(aPeptide);

					delete(aPeptide);
				}

				currentExperiment->ProteinsList->Add(aProtein);

				delete (aProtein);
			}

			j = sFiles[l]->LastIndexOf("\\") + 1;

			if(j >= 1)
			{
				currentExperiment->sExperimentFile = sFiles[l]->Substring(j);
			}
			else
			{
				currentExperiment->sExperimentFile = sFiles[l];
			}

			allExperiments->ExperimentsList->Add(currentExperiment);

			delete(mzI);

			delete(currentExperiment->ProteinsList);
		} //if(sFiles[l]->IndexOf(".mzid") != -1)
	}

	printf("\nNumber of mzid files processed: %d\n", allExperiments->ExperimentsList->Count);

	if(allExperiments->ExperimentsList->Count == 0)
	{
		printf("Did not Read any Data to Process\n");

		exit (1);
	}

	if(bScansMasses)
	{
		printf("\nWrite Scans is ON\n");

		WriteScansMasses(allExperiments);
	}
	else
	{
		WriteResults(allExperiments,args);
	}

	delete(allExperiments);

	return 0;
}


bool check_matching_modlocations(List<int> ^m1, List<int> ^m2) //Checks each position and returns true/false based on match/mismatch
{
	int i;

	if(m1->Count != m2->Count)
		return false;
	else if(m1->Count == m2->Count)
	{
		for(i = 0; i < m1->Count; i++)
		{
			if(m1[i]!=m2[i])
				return false;
		}
	}

	return true;
}



int find_add_missing_time_points(ExperimentCollection ^allExperiments, int exp_no, int *Charge,int *nScan, float *fRetTime, float *fstart, float *fend, double *SeqMass, double *SpecMass, 
	double *dMassToCharge, String ^protein, String ^peptide,array<double, 2> ^dIsotopes,array<double, 1> ^w, List <int> ^m)
{
	int i1, j1, k1;
	double dtemp;

	//Go backward
	i1 = exp_no;
	while(i1 >= 0)
	{
		for(j1 = 0; j1 < allExperiments->ExperimentsList[i1]->ProteinsList->Count; j1++)
		{

			if(allExperiments->ExperimentsList[i1]->ProteinsList[j1]->ProteinScore > ProteinScoreCut)
			{

				for(k1 = 0; k1 < allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides->Count; k1++)
				{
					if(allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->bPeptidePassed)
					{

						if(peptide->Equals(allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->Peptide) &&
							protein->Equals(allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->Protein) &&         
							Math::Abs(*SpecMass - allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->SpecMass) <= 0.01 &&
							Math::Abs(*SeqMass - allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->SeqMass) <= 0.01 &&
							*Charge == allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->nCharge &&
							check_matching_modlocations(m,allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->ModLocations))
						{

							dtemp           = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->SeqMass;

							(*Charge)         =  allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->nCharge;

							//(*SeqMass)        = (double)(*Charge)*dtemp - (double)((*Charge) - 1)*mProton;
							(*SeqMass)           = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->SeqMass;

							(*SpecMass)       = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->SpecMass;

							(*nScan)           = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->nScan;

							(*dMassToCharge)   = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->SeqMass; 

							(*fRetTime)        = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->dRetTime; 

							(*fstart) = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->fFirstID_Elution;

							(*fend) =   allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->fLastID_Elution;

							//dIsotopes = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->dIsotopes;

							//w =  allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->w;
							m = allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides[k1]->ModLocations;

							printf("Found the missing time point Protein : %s, Peptide : %s, RetTime : %f, SpecMass : %f, SeqMass : %f\n",protein,peptide,*fRetTime/60.0, *SpecMass, *SeqMass);

							return 1;
							break;
						}


					}

				}

			}
			if(k1 < allExperiments->ExperimentsList[i1]->ProteinsList[j1]->Peptides->Count)
				break;


		}
		if(j1 < allExperiments->ExperimentsList[i1]->ProteinsList->Count)
			break;

		i1 = i1 - 1;

	}

	return 0;
}

/*
* a program to write results to an csv file
* generate a copy of the all experiments dataset
* that will hold the sorted/filtered data
*/
void WriteResults(ExperimentCollection ^allResults, array<System::String ^> ^args)
{
	FILE *fp;

	int i, j, k, m;

	bool bFound = false;

	List  <String ^>  ^allAccessions = gcnew List  <String ^>;

	List  <String ^>  ^allProteins = gcnew List  <String ^>;

	List  <bool>  ^bAccessions = gcnew List  <bool>;

	List  <int>  ^lLengths = gcnew List  <int>;

	char szTemp[10000];

	String ^sTemp;

	String ^filename4 = args[3];

	printf("Total Number of Experiments %d\n", allResults->ExperimentsList->Count);


	String ^tempFile;

	char szOutFile[2024];

	//file dir creation

	if(filename4 != "")

		tempFile = filename4 + "\\"+ "Proteins.csv";//change on 05/13/2016
	else

		tempFile = "Proteins.csv"; //change on 05/13/2016

	szOutFile[0] = '\0';

	for(k=0; k < tempFile->Length; k++)
	{
		szOutFile[k] = tempFile[k];
	}

	szOutFile[k] = '\0';


	fp = fopen(szOutFile, "w");

	if(NULL == fp)
	{
		printf("Cannot write/store in this folder ...\n");

		printf("Results were ready\n");

		exit (1);
	}

	//comprise the list of all proteins/accessions, allAccessions that have been
	//observed in anyone of the experiments
	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++)
		{
			bFound = false;

			for(k=0; k < allAccessions->Count; k++)
			{
				if(allAccessions[k]->Length != 
					allResults->ExperimentsList[i]->ProteinsList[j]->accession->Length)
				{
					continue;
				}
				else if(allResults->ExperimentsList[i]->ProteinsList[j]->accession->EndsWith(allAccessions[k]))
				{
					bFound = true;

					break;
				}
			}

			if(false == bFound &&
				allResults->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences >= 2 &&
				allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore > ProteinScoreCut)
			{
				allAccessions->Add(allResults->ExperimentsList[i]->ProteinsList[j]->accession);

				allProteins->Add(allResults->ExperimentsList[i]->ProteinsList[j]->description);

				lLengths->Add(allResults->ExperimentsList[i]->ProteinsList[j]->nSeqLength);

				bAccessions->Add(false);
			}
		}
	}

	printf("A Total of %d Proteins in %d Experiments\n", allAccessions->Count,
		allResults->ExperimentsList->Count);

	fprintf(fp, "Accession,Protein,Length");

	//write experiment designations
	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		j = allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".");

		if(j < 0)
		{
			j = allResults->ExperimentsList[i]->sExperimentFile->Length;
		}

		if(allResults->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, ",%s,\n", 
				allResults->ExperimentsList[i]->sExperimentFile->Substring(0, j));
		}
		else
		{
			fprintf(fp, ",%s,", allResults->ExperimentsList[i]->sExperimentFile->Substring(0, j));
		}
	}

	//write protein score and spectral count Headers
	//first one is shifted (the commas before ProteinScore, since two columns is needed for 
	// accessions, lengths;
	fprintf(fp, ",,,ProteinScore,SCN,");

	for(i=1; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		if(allResults->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, "ProteinScore, SCN\n");
		}
		else
		{
			fprintf(fp, "ProteinScore,SCN,");
		}
	}

	printf("Number of Proteins = %d\n", allProteins->Count);
	//Write the actual data

	for(k=0; k < allAccessions->Count; k++)
	{
		//remove the "," from the protein names;

		m = allProteins[k]->IndexOf(",");

		int l = 0;

		for(i=0; i < allProteins[k]->Length && i < 10000; i++)
		{
			if(allProteins[k][i] != ',')
			{
				szTemp[l] = allProteins[k][i];

				l++;
			}
		}

		szTemp[l] = '\0';

		//fprintf(fp, "%s, %-30.28s, %d,", allAccessions[k], sTemp, lLengths[k]);

		sTemp = gcnew String(szTemp);

		fprintf(fp, "%s, %s, %d,", allAccessions[k], sTemp, lLengths[k]);

		for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
		{
			bFound = false;

			for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++) // over proteins of a single experiment
			{
				if(allAccessions[k]->Equals(allResults->ExperimentsList[i]->ProteinsList[j]->accession))
				{

					if(  allResults->ExperimentsList->Count - 1 == i)
					{
						fprintf(fp,"%5.1f, %d\n", 
							allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore,
							allResults->ExperimentsList[i]->ProteinsList[j]->nSpectralCount);
					}
					else
					{
						fprintf(fp,"%5.1f, %d,", 
							allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore,
							allResults->ExperimentsList[i]->ProteinsList[j]->nSpectralCount);
					}

					bFound = true;

					break;
				}
			}

			if(false == bFound)
			{
				if( allResults->ExperimentsList->Count - 1 == i)
				{
					fprintf(fp, "0.0, 0\n");
				}
				else
				{
					fprintf(fp,"0.0, 0,");
				}
			}
		}
	}



	fclose(fp);

	delete(allProteins); 
	delete(allAccessions);
	delete(lLengths);

}



/*
* a program to write scans, corresponding
*  precursor masses, charge states, and scores
*/
void WriteScansMasses(ExperimentCollection ^allResults)
{
	int i, j, k, m, nCharge, l;

	int i_mod = 0, i_pep=0;

	double dtemp,  dSeqMass, dSpecMass;

	char szFile[2049];

	bool bMod = false;

	FILE *fp;

	i = k = j = 0;

	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		l = allResults->ExperimentsList[0]->sExperimentFile->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase);

		if(l < 1)
		{
			l = allResults->ExperimentsList[0]->sExperimentFile->Length;
		}

		szFile[0] = '\0';

		for(m=0; m < l; m++)
		{
			szFile[m] = allResults->ExperimentsList[0]->sExperimentFile[m];
		}

		szFile[m] = '\0';

		strcat(szFile, ".Scan.Charge.Score.Sequence");

		fp = fopen(szFile, "w");

		if(NULL == fp)
		{
			printf("Cannot write to %s  ...\n", szFile);

			exit (1);
		}

		fprintf(fp, "Scan\t\RetTime\tCharge\tScore\tExpectation\tSeqMass\tExpMass\Sequence\n");

		for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++) // over proteins of a single experiment
		{	
			for(int k=0; k < allResults->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; k++)
			{

				dtemp = allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass;
				nCharge =  allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge;

				dSeqMass  = (double)nCharge*dtemp - (double)(nCharge - 1)*mProton;

				dtemp = allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass;

				dSpecMass = (double)nCharge*dtemp - (double)(nCharge - 1)*mProton;

				fprintf(fp, "%6d %6.2f %2d %8.2f %e %10.5f %10.5f ", allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nScan,
					allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dRetTime,
					allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge, 
					allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIonscore,
					allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dExpect,
					dSeqMass, dSpecMass);
				/*allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass,
				allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass);*/

				if (allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations->Count == 0)
				{
					fprintf(fp, "%s\n", allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide); 
				}
				else
				{
					i_mod++;

					for(int ii = 0; ii < allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide->Length; ii++)
					{
						bMod = false;

						for(int q = 0; q < allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations->Count; q++)
						{
							if(ii == allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations[q] - 1)
							{
								fprintf(fp, "%c", tolower(allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide[ii]));

								bMod = true;

								break;
							}
						}

						if(false == bMod)
						{
							fprintf(fp, "%c", allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide[ii]);
						}
					}

					fprintf(fp, "\n");
				}

				i_pep++;

			}
		}

		fclose(fp);
	}

	printf("%d number of modified peptides out of %d peptides\n", i_mod, i_pep); 

	exit (1);
}

/*
* A method to write results for Heavy Water Labeling
* For every mzid file, we write a list proteins, their peptides,
*  peptide scores, masses, scans, and sequences.
*  calls WriteTimeSeries method
*/
void D2O_H2O(ExperimentCollection ^allResults,array<System::String ^> ^args)
{

	int i, j, k, l, nCharge,j1;

	double dSeqMass, dSpecMass,  dtemp;

	FILE *fp;

	String ^tempFile;

	String ^filename4 = args[3];

	char szOutFile[2024], szTempProt[2024];

	//sfprintf(fp, "Accession,Protein,Length");

	/*
	*    For every experiment, create a file source.Proteins.csv
	*    write into this file all proteins with scores greater than 100
	*    and for every proteins write all peptide sequences with scores
	*    higher than 20.
	*/
	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		j = allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf("\\")+1;
		j1 = allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".");

		if(j < 0)
		{
			j = allResults->ExperimentsList[i]->sExperimentFile->Length;
		}

		if(filename4 != "")
			tempFile = filename4 + "\\" + allResults->ExperimentsList[i]->sExperimentFile->Substring(j,j1-j) + ".Proteins.csv";  // previous it was Substring(0,j)
		else 
			tempFile = allResults->ExperimentsList[i]->sExperimentFile->Substring(j,j1-j) + ".Proteins.csv";
		//tempFile = ".Proteins.csv"; //Otherwise write into the current directory

		szOutFile[0] = '\0';

		for(k=0; k < tempFile->Length; k++)
		{
			szOutFile[k] = tempFile[k];
		}

		szOutFile[k] = '\0';

		fp = fopen(szOutFile, "w");

		if(NULL == fp)
		{
			printf("\nCannot write to %s\n", tempFile);

			printf("Exiting unfinished... \n");

			exit (1);
		}

		fprintf(fp, "Accession, Protein, ProtScore, ProtLength, Peptide, PepScore, Charge, Scan, RetTime,SeqMass, ExpMass\n");

		for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++) // over proteins of a single experiment
		{
			if(allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore > ProteinScoreCut)
			{
				for(k=0; k < allResults->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; k++)
				{
					if(0 == k)
					{
						szTempProt[0] = '\0';

						//if protein description contains a comma, remove it, as it confuses csv
						for(l=0; l < allResults->ExperimentsList[i]->ProteinsList[j]->description->Length; l++)
						{
							if(allResults->ExperimentsList[i]->ProteinsList[j]->description[l] == ',')
							{
								szTempProt[l] = ' ';
							}
							else
							{
								szTempProt[l] = allResults->ExperimentsList[i]->ProteinsList[j]->description[l];
							}
						}

						szTempProt[l] = '\0';

						fprintf(fp, "%s, %s, %6.1f, %8d,", allResults->ExperimentsList[i]->ProteinsList[j]->accession,
							szTempProt,
							allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore,
							allResults->ExperimentsList[i]->ProteinsList[j]->nSeqLength);
					}
					else
					{
						fprintf(fp, " , , , ,");
					}

					dtemp = allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass;
					nCharge =  allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge;

					dSeqMass  = (double)nCharge*dtemp - (double)(nCharge - 1)*mProton;

					dtemp = allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass;

					dSpecMass = (double)nCharge*dtemp - (double)(nCharge - 1)*mProton;

					fprintf(fp, "%s,%5.1f, %2d, %6d, %8.2f, %8.2f, %10.5f\n", 
						allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide,
						allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIonscore,
						allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge,
						allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nScan, 
						allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dRetTime/60.,
						dSeqMass, dSpecMass);
					/*allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass,
					allResults->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass);*/

				}
			}
		}


		fclose(fp);

	}

	//ConsistentProteins(allResults, 2, 100.0, allResults->ExperimentsList->Count -1, lProteins); 
	printf("calling ConsistentProteins\n");

	ConsistentProteins(allResults, 2, ProteinScoreCut, nProteinConsistency);

	printf("calling ConsistentPeptides\n");
	
	ConsistentPeptides(allResults, dIonScoreCut, dExpectCut, nPeptideConsistency);

	printf("calling WriteTimeSeries\n");
	//write the out the proteins that are in all experiments.
	WriteTimeSeries(allResults,args);

	printf("calling WriteResults\n");
	//write all proteins all times
	WriteResults(allResults,args);

	return;

}

/*
*  A workflow that starts Cleveland Clinic Processing
* smzID - holds names of mzid files
* mzML  - holds names of mzML files
*/
int Start_Process(array<System::String ^> ^args)
{

	// Stop the code if there is any csv file in the directory

	String^ CurrentFolder = Directory::GetCurrentDirectory();

	cli::array<String ^, 1> ^ FileList;
	FileList = Directory::GetFiles(CurrentFolder);
	for (int i = 0; i < FileList->Length; i++)
	{
		if (FileList[i]->EndsWith(".csv"))
		{
			printf("\n\n================\nPlease remove all the *.csv files from this directory and run again!\n================\n\n");
			exit(1);
		}
	}

	array <String ^> ^smzId, ^smzML;

	array <float> ^BodyWaterEnrich, ^ExperimentTime;

	PeptideHolder ^aPeptide = gcnew PeptideHolder();

	ProteinCollection ^currentExperiment = gcnew ProteinCollection();

	ExperimentCollection ^allExperiments = gcnew ExperimentCollection();

	allExperiments->ExperimentsList = gcnew List<ProteinCollection^>;

	ProteinSet ^aProtein;

	char szTemp[2024], szLine[2024], szTemp1[2024], szTemp2[2024];

	String ^strFile;           // name of the file holding information on mzML, mzID, time, and body water enrichment

	bool bSpectra = false;     //this varialbe controls reading and processing of mzML files if provided

	int i, k, l, j, m, num_of_experiments;

	if(args[0]->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase) == -1 && args[0]->IndexOf(".txt", StringComparison::OrdinalIgnoreCase) == -1)
	{
		printf("Need an mzid or txt file to run the program\n");

		exit (1);
	}

	if(args[0]->IndexOf("*.mzid", StringComparison::OrdinalIgnoreCase) != -1)              // caller_mzIdentML.exe *.mzid
	{
		smzId = System::IO::Directory::GetFiles("\.");
	}
	else if(args->Length >= 1 && args[0]->IndexOf(".txt") != -1)
	{
		Isotope_Deconvolution = int::Parse(args[1]); 

		printf("Isotope deconvolution option :%d\n",Isotope_Deconvolution);		
		
		bSpectra = true;

		strFile = args[0];      //asssumes that the argument is the files.txt (file name for holding time points, mzml, mzid and body water data
		

		num_of_experiments = NumberOfExperiments(args[0]);

		smzId = gcnew array<String ^> (num_of_experiments);

		smzML = gcnew array<String ^> (num_of_experiments);

		BodyWaterEnrich = gcnew array <float> (num_of_experiments);

		ExperimentTime = gcnew array <float> (num_of_experiments);      //default unit of time is a day!!!

		Parse_InputTextFile(args[0], smzML, smzId, ExperimentTime, BodyWaterEnrich);
	}
	else
	{
		smzId = gcnew array<String ^> (1);

		smzId[0] = args[0];
	}

	k = 0;

	for(l=0; l < smzId->Length; l++)
	{
		if(smzId[l]->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase) != -1)
		{
			k++;

			// process/read the mzid file
			mzIdentML ^mzI = gcnew mzIdentML(smzId[l]);

			mzI->ReadProteins_mzIdent_2();

			//now create a storage for this experiment
			currentExperiment = gcnew ProteinCollection();

			currentExperiment->ProteinsList = gcnew List<ProteinSet^>;

			//count goes over every protein entry read in from mzid
			for(i=0; i < mzI->ProteinResultList->Count; i++)
			{
				//copy protein data
				aProtein = gcnew ProteinSet();

				aProtein->Peptides = gcnew List <PeptideHolder ^>;


				aProtein->accession          = mzI->ProteinResultList[i]->accession;
				aProtein->description        = mzI->ProteinResultList[i]->description;
				aProtein->ProteinScore       = mzI->ProteinResultList[i]->ProteinScore;
				aProtein->SeqCoverage        = mzI->ProteinResultList[i]->SeqCoverage;
				aProtein->nSpectralCount     = mzI->ProteinResultList[i]->nSpectralCount;
				aProtein->nDistinctSequences = mzI->ProteinResultList[i]->nDistinctSequences;
				aProtein->nSeqLength         = mzI->ProteinResultList[i]->nSeqLength;
				aProtein->ProteinMass        = mzI->ProteinResultList[i]->ProteinMass;

				for(int i0 = 0; i0 < mzI->ProteinResultList[i]->PeptideList->Count; i0++)
				{

					/*if(mzI->ProteinResultList[i]->PeptideList[i0]->Peptide == "IGVELTGR" && k==2)
					{ 	*/
					aPeptide = gcnew PeptideHolder();

					aPeptide->dIonscore = mzI->ProteinResultList[i]->PeptideList[i0]->dIonscore;
					aPeptide->bUniquePeptide = mzI->ProteinResultList[i]->PeptideList[i0]->bUniquePeptide;
					aPeptide->dRetTime  = mzI->ProteinResultList[i]->PeptideList[i0]->dRetTime;
					aPeptide->Peptide   = mzI->ProteinResultList[i]->PeptideList[i0]->Peptide;
					aPeptide->nScan     = mzI->ProteinResultList[i]->PeptideList[i0]->nScan;
					aPeptide->nCharge   = mzI->ProteinResultList[i]->PeptideList[i0]->nCharge;
					aPeptide->SeqMass  = mzI->ProteinResultList[i]->PeptideList[i0]->SeqMass;
					aPeptide->SpecMass = mzI->ProteinResultList[i]->PeptideList[i0]->SpecMass;
					aPeptide->dExpect    = mzI->ProteinResultList[i]->PeptideList[i0]->dExpect;
					aPeptide->dIsotopes = gcnew array <double, 2>(2, 6);
					aPeptide->w = gcnew array <double, 1>(12);
					aPeptide->Protein = mzI->ProteinResultList[i]->accession;

					//printf("Protein name : %s\n", aPeptide->Protein);

					aPeptide->dModMasses = gcnew List <double>;

					aPeptide->ModLocations = gcnew List <int>;

					aPeptide->dModMasses = mzI->ProteinResultList[i]->PeptideList[i0]->dModMasses;

					aPeptide->ModLocations = mzI->ProteinResultList[i]->PeptideList[i0]->ModLocations; 


					aProtein->Peptides->Add(aPeptide);

					delete(aPeptide);

					/*}	*/



				}

				aProtein->ProteinIDnumber = -1;

				currentExperiment->ProteinsList->Add(aProtein);

				delete (aProtein);
			}

			j = smzId[l]->LastIndexOf("\\") + 1;

			if(j >= 1)
			{
				currentExperiment->sExperimentFile = smzId[l]->Substring(j);
			}
			else
			{
				currentExperiment->sExperimentFile = smzId[l];
			}

			currentExperiment->fExperimentTime = ExperimentTime[l];

			currentExperiment->fBWE = BodyWaterEnrich[l];


			allExperiments->ExperimentsList->Add(currentExperiment);

			j = smzML[l]->LastIndexOf("\\") + 1;

			j = -1;

			if(j >= 1)
			{
				currentExperiment->smzML = smzML[l]->Substring(j);
			}
			else
			{
				currentExperiment->smzML = smzML[l];
			}

			delete(mzI);

			delete(currentExperiment->ProteinsList);
		} //if(sFiles[l]->IndexOf(".mzid") != -1)




	}

	printf("\nNumber of mzid files processed: %d\n", allExperiments->ExperimentsList->Count);


	if(allExperiments->ExperimentsList->Count == 0)
	{
		printf("Did not Read any Data to Process\n");

		exit (1);
	}


	//determine the consistent peptides (e.g., peptides that have been identified in a certain number of 
	// experiemnts, and write some results - no quantification
	D2O_H2O(allExperiments,args);

	//read the intensity information from the mzML file, and
	// do the actual quantification
	if(bSpectra)
	{
		ReadAbundances(smzML, allExperiments);

		WriteQuantResults(allExperiments,args);   // RGS moved it here from ReadAbudnaces

		String^ path = Directory::GetCurrentDirectory();

		printf("Args %s   %s\n", args[3], path);

		//ProteinRates(strFile)

		ProteinRates(ExperimentTime, BodyWaterEnrich,Nparameter);
	}
	if (allExperiments) {
		delete(allExperiments);
	}
	
	delete(BodyWaterEnrich); delete(ExperimentTime);

	return 0;
}


/*
*   a method to read abundances and isotopes from mzML files
*/
int ReadAbundances(array <String ^> ^smzML, ExperimentCollection ^AllExperiments)
{
	int i, k, j, nCharge, q, count, Charge;

	int nScan;

	double dSeqMass, dSpecMass, dtemp, SeqMass, SpecMass, Prev_SpecMass;

	double dMassToCharge;

	float fRetTime, StartElution, EndElution;

	array <double, 1> ^w;

	String ^sTemp;
	String ^currentPeptide;
	bool peptideAssigned;
	int prev_peptide_count = 0;

	StreamWriter ^swOutput = gcnew StreamWriter("HeavyWater.log");



	for(i=0; i < AllExperiments->ExperimentsList->Count; i++) //over experiments
	{	
		printf("Processing %s\n", AllExperiments->ExperimentsList[i]->smzML);
		// initiate mzML file

		if(System::IO::FileInfo(AllExperiments->ExperimentsList[i]->smzML).Exists == false)
		{
			printf("Cannot read file: %s\n", AllExperiments->ExperimentsList[i]->smzML);

			printf("Exiting ..\n");

			exit (1);
		}

		DetectAndIntegratePeak ^DetectPeaks = gcnew  DetectAndIntegratePeak(AllExperiments->ExperimentsList[i]->smzML, 0.030, ElutionTimeWindow, dMassAccuracy, swOutput);
		
		DetectPeaks->BuildChromatogram(1);

		count = 0;
		for(j = 0; j < AllExperiments->ExperimentsList[i]->ProteinsList->Count; j++) // iterate over proteins of a single experiment
		{
			if(AllExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinScore > ProteinScoreCut)
			{
				for(k=0; k < AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; k++) // iterate over all peptides (including non-unique) of a protein
				{
	
					Prev_SpecMass = 0;

					//set this variables in any case, even if the are equal to 0.0
					AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes = gcnew array <double, 2>(2, 6);
					w = gcnew array <double, 1>(12);


					if(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptidePassed /*&& AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide=="ALVDTLK"*/) 
					{					
						dtemp           = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass;

						nCharge         =  AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge;

						dSeqMass        = (double)nCharge*dtemp - (double)(nCharge - 1)*mProton;

						dSpecMass       = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass;

						nScan           = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nScan;

						dMassToCharge   = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass; 

						fRetTime        = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dRetTime; 


						Prev_SpecMass = dSpecMass; 

						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes = gcnew array <double, 2> (2, 6);
						w = gcnew array <double, 1> (12);

						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->w = gcnew array <double, 1> (12);


						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptideTimePoint = true;
						//If at least one peptide of the protein has the time point available, we set the protein time point availability to true
						AllExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinTimePoint = AllExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinTimePoint | AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptideTimePoint; 

						if(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations->Count > 0)
						{
							for(q=0; q < AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations->Count; q++)
							{
								if(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations[q] == 
									AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide->Length + 1)
								{
									//printf("MAss = %10.5f\n", AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dModMasses[q]);

									dMassToCharge = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass - 
										AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dModMasses[q]/(double)nCharge;
								}
							}
						}

						//printf("Scan: %d, Peptide %s\n", nScan, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide);
						//fprintf(fpLog, 
							//"Scan: %d, Peptide %s\n", nScan, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide);

						swOutput->Write("Scan: {0}, Peptide {1} \n", nScan, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide);

		

						StartElution = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->fFirstID_Elution;
						EndElution =   AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->fLastID_Elution;

						count = count+1;

						DetectPeaks->IntegratePeakTimeMz(fRetTime, nScan -1, dMassToCharge, nCharge, &StartElution, &EndElution,
						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes,w);

						//printf("\nAfter DetectPeaks: Ret time = %f, Spec Mass = %f, Scan = %d, Start Elution = %f, End Elution = %f\n",fRetTime,dMassToCharge,nScan,StartElution,EndElution);
						//						DetectPeaks->IntegratePeakTimeMz(fRetTime, nScan -1, dMassToCharge, nCharge, &StartElution, &EndElution,
						//	AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes);


						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->fStartElution = StartElution;

						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->fEndElution = EndElution;

						peptideAssigned = true;

						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->w = w;

						//fprintf(ftemp,"%f,%f,%f,%f,%f,%f\n",w[0],w[1],w[2],w[3],w[4],w[5]);

						//fprintf(ftemp,"%f,%f,%f,%f,%f,%f\n",w[6],w[7],w[8],w[9],w[10],w[11]);


						/*printf("MS Scan = %d MS2_Scan = %d %10.5f %10.5f\n", nScan,
						AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nScan,
						mzml->Spectrum[0, i0], mzml->Spectrum[1, i0]);*/

					}
					// it needs additional look most of the time, 09/08/2017
					else if(!AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptidePassed && i-1 >= 0 /*&& AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide=="ALVDTLK"*/) 
					{
						//checks in the peptide list if this peptide has elution profile before or after this current missed eluted peptide
						String ^protein = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Protein;
						String ^peptide = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide;
						SpecMass = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass;
						Charge = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge;
						SeqMass = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass;
						List <int> ^m = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations;

						peptideAssigned = false;
						
						for(int k1 = 0; k1 < AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; k1++)
						{
							if(peptide->Equals(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->Peptide) && 
								Charge == AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->nCharge &&
								Math::Abs(SpecMass - AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->SpecMass) <= 0.01 &&
								Math::Abs(SeqMass - AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->SeqMass) <= 0.01 &&
								check_matching_modlocations(m, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->ModLocations) &&
								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->bPeptidePassed)
							{
								peptideAssigned = true;
								
								break;
							}
							else if(peptide->Equals(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->Peptide) && 
								Charge == AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->nCharge &&
								check_matching_modlocations(m, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->ModLocations) &&
								Math::Abs(SpecMass - AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k1]->SpecMass) >= 0.01) //If there is SpecMass difference 
							{

								peptideAssigned = true;
								
								break;

							}
						}

						if(!peptideAssigned) 
						{
							printf("Assigning time points in !peptideAssigned: Protein :%s, Peptide : %s, RetTime : %f, SpecMass :%f, SeqMass :%f, Score = %f\n",
								protein, peptide, fRetTime/60.0, SpecMass, SeqMass, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIonscore);
						}

						if(!peptideAssigned)   // remove this block
						{
							//if(i-1 >= 0) //checking if the time point is > 0
							//{
							//AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptidePassed = true;


							//printf("Assigning time points: Protein :%s, Peptide : %s, RetTime : %f, SpecMass :%f, SeqMass :%f\n",protein, peptide, fRetTime/60.0, SpecMass, SeqMass);

							//Calling the function that will assign return the previous time points

							array<double, 2> ^dIsotopes = gcnew array <double, 2> (2,6);
							array<double, 1> ^w = gcnew array <double, 1> (12);


							if(find_add_missing_time_points(AllExperiments, i-1, &Charge, &nScan, &fRetTime, &StartElution, &EndElution, &SeqMass, &SpecMass, &dMassToCharge, protein, peptide, dIsotopes, w, m))
							{

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptidePassed = true;
								
								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nCharge = Charge;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SpecMass = SpecMass;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->nScan = nScan;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass = SeqMass; 

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dRetTime = fRetTime; 


								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations = m;

								Prev_SpecMass = SpecMass;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes = gcnew array <double, 2> (2, 6);
								w = gcnew array <double, 1> (12);

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->w = gcnew array <double, 1> (12);


								/*AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes = dIsotopes;								

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->w = w;*/
								//addjust the mononoisotopic mass to that of the original molecular in the cases of stable-isotope
								// labeling. E.g., for O18 labeling, take out 4 Da or 2 Da dependending on the incorporation level
								// Modification location for O18 peptide comes out as 1 + peptide.length
								if(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations->Count > 0)
								{
									for(q=0; q < AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations->Count; q++)
									{
										if(AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->ModLocations[q] == 
											AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide->Length + 1)
										{
											//printf("MAss = %10.5f\n", AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dModMasses[q]);

											dMassToCharge = AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->SeqMass - 
												AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dModMasses[q]/(double)Charge;
										}
									}
								}

								//printf("Scan: %d, Peptide %s\n", nScan, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide);
								//fprintf(fpLog, 
								//"Scan: %d, Peptide %s\n", nScan, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide);
								swOutput->Write(
									"Scan: {0}, Peptide {1}\n", nScan, AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->Peptide);

								DetectPeaks->IntegratePeakTimeMz(fRetTime, nScan -1, dMassToCharge, Charge, &StartElution, &EndElution,
									AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->dIsotopes,w);


								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->fStartElution = StartElution;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->fEndElution = EndElution;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->w = w;

								AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptideTimePoint = true;
								//If at least one peptide of the protein has the time point available, we set the protein time point availability to true
								AllExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinTimePoint = AllExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinTimePoint | AllExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[k]->bPeptideTimePoint; 

		
								peptideAssigned = true;
							}
							

						}//if (!peptideAssigned)
					}
				}
			}
		}

		DetectPeaks->Close_mzML_file();

		delete(DetectPeaks);

		printf("Finished Processing %s\n", AllExperiments->ExperimentsList[i]->smzML);

		swOutput->Flush();
	}

	//fclose(fpLog);
	
	swOutput->Close();

	//WriteQuantResults(AllExperiments,args);  RGS 03/22/2017

	return 0;
}

/*
*  A method to write out all informations
*  about proteins in a "time series" format.
*  It will write out only those proteins that
*  have been observed in all experiments.
*/
void WriteTimeSeries(ExperimentCollection ^allResults,array<System::String ^> ^args)
{
	int i, j, k, l, iExp, m = 0;

	//From here below find allAccessions, then find those accessions that are
	// common to all experiments, and print their peptides


	List  <String ^>  ^allAccessions = gcnew List  <String ^>;

	List  <bool>  ^CommonAccessions = gcnew List  <bool>;

	bool bFound = false;

	String ^filename4 = args[3];

	String ^tempFile;

	char szTempFile[2046], c;	

	m = 0; szTempFile[0] = '\0';

	if(filename4 != "")
		tempFile = filename4 + "\\" + "Time.Series.AllProteins.csv";
	else
		tempFile = "Time.Series.AllProteins.csv";

	szTempFile[0] = '\0';

	for(k=0; k < tempFile->Length; k++)
	{
		szTempFile[k] = tempFile[k];
	}

	szTempFile[k] = '\0';


	//comprise the list of all proteins/accessions, allAccessions that have been
	//observed in anyone of the experiments
	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++)
		{
			bFound = false;

			for(k=0; k < allAccessions->Count; k++)
			{
				if(allAccessions[k]->Length != 
					allResults->ExperimentsList[i]->ProteinsList[j]->accession->Length)
				{
					continue;
				}
				else if(allResults->ExperimentsList[i]->ProteinsList[j]->accession->EndsWith(allAccessions[k]))
				{
					bFound = true;

					break;
				}
			}

			if(false == bFound &&
				allResults->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences >= 2 &&
				allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore > ProteinScoreCut)
			{
				allAccessions->Add(allResults->ExperimentsList[i]->ProteinsList[j]->accession);

				CommonAccessions->Add(false);
			}
		}

	}

	//index proteins that pass the protein cut-off score and are observed in all experiments
	iExp = 0;                            

	for(k=0; k < allAccessions->Count; k++)  //over all proteins
	{
		iExp = 0;

		for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
		{
			for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++)
			{
				if(allAccessions[k]->Length == 
					allResults->ExperimentsList[i]->ProteinsList[j]->accession->Length &&
					allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore > ProteinScoreCut)
				{
					if(allResults->ExperimentsList[i]->ProteinsList[j]->accession->EndsWith(allAccessions[k]))
					{
						iExp++;

						break;
					}
				}
			}
		}

		if(allResults->ExperimentsList->Count == iExp) // if a given protein is observed in all experiments, its common = true; 
			//if(iExp == nProteinConsistency)
		{
			CommonAccessions[k] = true;
		}
	}

	bFound = false;

	for(i=0; i < CommonAccessions->Count; i++)
	{
		if(CommonAccessions[i])
		{
			bFound = true;

			break;
		}
	}

	//write common proteins

	FILE *fp1 = fopen(szTempFile, "w");
	if(fp1==NULL)
	{
		printf("cannot open out put directory\n");
		exit(1);
	}
	fprintf(fp1, " , , , ,");

	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		if( allResults->ExperimentsList->Count - 1  == i)
		{
			fprintf(fp1, "%s, , , \n", allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
				allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );
		}
		else
		{
			fprintf(fp1, "%s, , , ,", 
				allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
				allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );
		}
	}

	fprintf(fp1, "Number, Accession, Description, Length,"); 

	for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
	{
		if( allResults->ExperimentsList->Count - 1  == i)
		{
			fprintf(fp1, "Score, SpecCount, UniquePeptides, Coverage\n");
		}
		else
		{
			fprintf(fp1, "Score, SpecCount, UniquePeptides, Coverage,");
		}

	}

	l = 0;

	//writes most of the Time.Series.Proteins.csv
	for(k=0; k < allAccessions->Count; k++)  //over all proteins
	{
		if(true == CommonAccessions[k])
		{
			for(i=0; i < allResults->ExperimentsList->Count; i++) //over experiments
			{
				for(j = 0; j < allResults->ExperimentsList[i]->ProteinsList->Count; j++)
				{
					bFound = false;

					if(allAccessions[k]->Length == 
						allResults->ExperimentsList[i]->ProteinsList[j]->accession->Length)
					{
						if(allResults->ExperimentsList[i]->ProteinsList[j]->accession->EndsWith(allAccessions[k]))
						{
							if(0 == i)
							{
								fprintf(fp1, "%d, %s,", (l+1), allResults->ExperimentsList[i]->ProteinsList[j]->accession);

								for(int ii=0; ii < allResults->ExperimentsList[i]->ProteinsList[j]->description->Length; ii++)
								{
									if(allResults->ExperimentsList[i]->ProteinsList[j]->description[ii] != ',')
									{
										fprintf(fp1, "%c", allResults->ExperimentsList[i]->ProteinsList[j]->description[ii]);
									}

									if(ii == allResults->ExperimentsList[i]->ProteinsList[j]->description->Length - 1)
									{
										fprintf(fp1, ",");
									}
								}


								fprintf(fp1, "%d,",  
									allResults->ExperimentsList[i]->ProteinsList[j]->nSeqLength);
							}

							if( allResults->ExperimentsList->Count - 1  == i)
							{
								fprintf(fp1, "%5.1f, %d, %d, %5.1f\n", allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore,
									allResults->ExperimentsList[i]->ProteinsList[j]->nSpectralCount, 
									allResults->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences,
									allResults->ExperimentsList[i]->ProteinsList[j]->SeqCoverage);
							}
							else
							{
								fprintf(fp1, "%5.1f, %d, %d, %4.1f, ", allResults->ExperimentsList[i]->ProteinsList[j]->ProteinScore,
									allResults->ExperimentsList[i]->ProteinsList[j]->nSpectralCount, 
									allResults->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences,
									allResults->ExperimentsList[i]->ProteinsList[j]->SeqCoverage);
							}

							bFound = true;

							break;
						}
					}
				}
			}

			l++;
		}
	}

	fclose(fp1);

	WriteTimeSeriesPeptidesInProtein(allResults, allAccessions, CommonAccessions,args);

	delete(CommonAccessions);

	delete(allAccessions);

	return;
}

/*
*  A method to write out separate information
*  about each protein and its peptides in a "time series" format.
*  It will write out only those proteins that
*  have been observed in all experiments.
*/
void WriteTimeSeriesPeptidesInProtein(ExperimentCollection ^allResults, 
	List  <String ^>  ^allAccessions, List  <bool>  ^CommonAccessions,array<System::String ^> ^args)
{
	int i, j, k, l, m, t, n, nCharge, iExp = 0, imax, s, q;

	double dSeqMass, dSpecMass,  dtemp, dError;

	bool bFound = false;

	List  <String ^>  ^ProteinPeptides;

	List <int> ^PeptideCharges;

	List  <int>  ^PeptideIndices;

	List <double> ^PeptideMass;

	List  <bool>  ^PeptideInAllExperiments;

	SingleList ^aList;

	List <SingleList ^> ^AllLists;

	//List <ModLocs> ^PeptideModLocs; //modification locations for all peptides;
	GeneralListCollector ^ ListCollector = gcnew GeneralListCollector();

	FILE *fp1 = stdout;

	//writes most of the Time.Series.Proteins.csv

	for(k=0; k < allAccessions->Count; k++)  //over all proteins
	{
		if(true == CommonAccessions[k])
		{
			ProteinPeptides = gcnew List  <String ^>;

			PeptideCharges = gcnew List <int>;

			PeptideIndices = gcnew List  <int>;

			PeptideMass = gcnew List <double>;

			PeptideInAllExperiments = gcnew List  <bool>;


			AllLists = gcnew List<SingleList^>;

			for(j = 0; j < allResults->ExperimentsList[0]->ProteinsList->Count; j++)   // if a protein is in all experiments, it is in the first experiment
			{
				if(allAccessions[k]->Length == 
					allResults->ExperimentsList[0]->ProteinsList[j]->accession->Length &&
					allResults->ExperimentsList[0]->ProteinsList[j]->accession->EndsWith(allAccessions[k]))
				{
					for(m=0; m < allResults->ExperimentsList[0]->ProteinsList[j]->Peptides->Count; m++) // over all peptides
					{
						aList = gcnew SingleList(); 

						aList->IntegerList = gcnew List <int>;

						aList->DoubleList = gcnew List <double>;

						aList->StringList = gcnew List <String ^>;


						bFound = false; dtemp = 0;

						for(t=0; t < ProteinPeptides->Count; t++)
						{
							if(allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->Peptide->Length == ProteinPeptides[t]->Length &&
								allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->Peptide->EndsWith(ProteinPeptides[t]) &&
								PeptideCharges[t] == allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->nCharge &&
								PeptideMass[t] == allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->SeqMass)	
							{
								if(allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->dIonscore >
									allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[PeptideIndices[t]]->dIonscore)
								{
									PeptideIndices[t] = m;
								}

								bFound = true;

								break;
							}
						}

						//collecting all peptides, such that one peptide per sequence (i.e. if two spectra matched the same peptide
						// store the sequence only once.
						if(false == bFound)
						{
							ProteinPeptides->Add(allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->Peptide);

							PeptideCharges->Add(allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->nCharge);

							PeptideMass->Add(allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->SeqMass);

							PeptideIndices->Add(m);

							aList->IntegerList = allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->ModLocations;

							aList->DoubleList = allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->dModMasses;

							aList->StringList->Add(allResults->ExperimentsList[0]->ProteinsList[j]->Peptides[m]->Peptide);

							AllLists->Add(aList);

							delete(aList);

							PeptideInAllExperiments->Add(false);
						}
					}

					//search and identify which peptides of the protein accessions[k] is present in 
					// all experiments, starts with 1, because the 0th experiment has already been analyzed

					for(m = 0; m < ProteinPeptides->Count; m++)
					{
						s = 0; //will not the number of experiments in which the m-the peptide is present
						// if that equals the total number of experiments, then that peptide is present in all experiments.

						for(i=1; i < allResults->ExperimentsList->Count; i++)
						{
							bFound = false;

							for(t = 0; t < allResults->ExperimentsList[i]->ProteinsList->Count; t++)   // if a protein is in all experiments, it is in the first experiment
							{
								if(allAccessions[k]->Length == 
									allResults->ExperimentsList[i]->ProteinsList[t]->accession->Length &&
									allResults->ExperimentsList[i]->ProteinsList[t]->accession->EndsWith(allAccessions[k]))
								{
									for(n=0; n < allResults->ExperimentsList[i]->ProteinsList[t]->Peptides->Count; n++)
									{
										if(allResults->ExperimentsList[i]->ProteinsList[t]->Peptides[n]->Peptide->Length == ProteinPeptides[m]->Length &&
											allResults->ExperimentsList[i]->ProteinsList[t]->Peptides[n]->Peptide->EndsWith(ProteinPeptides[m]) &&
											allResults->ExperimentsList[i]->ProteinsList[t]->Peptides[n]->nCharge == PeptideCharges[m] &&
											allResults->ExperimentsList[i]->ProteinsList[t]->Peptides[n]->SeqMass == PeptideMass[m])
										{
											//check to make sure that the mods are equal as well

											bool bSameMods = true;

											for(q=0; q < AllLists[m]->IntegerList->Count; q++)
											{
												if(AllLists[m]->IntegerList[q] != 
													allResults->ExperimentsList[i]->ProteinsList[t]->Peptides[n]->ModLocations[q])
												{
													bSameMods = false;

													break;
												}
											}


											if(bSameMods)
											{
												s++;

												bFound = true;

												break;
											}
										}
									}

								}
							}
						}

						// (s + 1) because, s is the number of experiments in which the peptide is observed, not counting
						// the first experiment.
						if(allResults->ExperimentsList->Count == (s + 1) )
						{
							PeptideInAllExperiments[m] = true;
						}
					}

					char szTempFile[2094], c; 

					s = 0;

					String ^filename4 = args[3];

					if(filename4 != "")
					{

						for(s =  0; s < filename4->Length; s++)
							szTempFile[s] = filename4[s];

						szTempFile[s++] = '\\';
					}


					for(l=0; l < allAccessions[k]->Length; l++)
					{
						c = allAccessions[k][l];

						if(c != '|' && c != ',' && c != '\n' && c != '-' &&
							c != '\t' && c != ' ')
						{
							szTempFile[s] = allAccessions[k][l];
							s++;
						}
					}

					szTempFile[s] = '\0';

					strcat(szTempFile, ".csv");

					FILE *fp2 = fopen(szTempFile, "w");

					if(NULL == fp2)
					{
						printf("Cannot write out %s\n", szTempFile);

						puts("Exiting ...");

						exit (1);
					}

					for(i=0; i < allResults->ExperimentsList->Count; i++)
					{
						if(0 == i)
						{
							fprintf(fp2, ",%s, , , , , , ,", allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
								allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );

							if(allResults->ExperimentsList->Count - 1 == i)
							{
								fprintf(fp2, "\n");
							}
						}
						else if(allResults->ExperimentsList->Count - 1 == i)
						{
							/*fprintf(fp2, "%s, , , , , , \n", allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
							allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );*/
							fprintf(fp2, "%s, , , , \n", allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
								allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );
						}
						else
						{
							/*fprintf(fp2, "%s, , , , , , ,", allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
							allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );*/
							fprintf(fp2, "%s, , , , ,", allResults->ExperimentsList[i]->sExperimentFile->Substring(0,
								allResults->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) );
						}
					}

					for(i=0; i < allResults->ExperimentsList->Count; i++)
					{
						if(0 == i)
						{
							//fprintf(fp2, "Peptide, SpecMass, Charge, SeqMass, IonScore, Expectn, Error(ppm), Scan,");
							fprintf(fp2, "Peptide, Charge, SeqMass, SpecMass, IonScore, Expectn, Error(ppm), Scan,");

							if(allResults->ExperimentsList->Count - 1 == i) 
							{
								fprintf(fp2, "\n");
							}
						}
						else if(allResults->ExperimentsList->Count - 1 == i)
						{
							//fprintf(fp2, "SpecMass, Charge, SeqMass, IonScore, Expectn, Error(ppm), Scan\n");
							fprintf(fp2, "SpecMass, IonScore, Expectn, Error(ppm), Scan\n");
						}
						else
						{
							//fprintf(fp2, "SpecMass, Charge, SeqMass, IonScore, Expectn, Error(ppm), Scan,");
							fprintf(fp2, "SpecMass, IonScore, Expectn, Error(ppm), Scan,");
						}
					}

					for(m=0; m < ProteinPeptides->Count; m++)
					{
						if(PeptideInAllExperiments[m])
						{

							if(AllLists[m]->IntegerList->Count == 0)
							{
								fprintf(fp2, "%s,", ProteinPeptides[m]);
							}
							else
							{
								for(int ii=0; ii < ProteinPeptides[m]->Length; ii++)
								{
									bool bModified = false;

									for(q=0; q < AllLists[m]->IntegerList->Count; q++)
									{
										if(AllLists[m]->IntegerList[q] == ii + 1)   //the mod position is counted starting from 1.
										{
											bModified = true;

											break;
										}
									}

									if(bModified)
									{
										fprintf(fp2, "%c", tolower(ProteinPeptides[m][ii]));
									}
									else
									{
										fprintf(fp2, "%c", ProteinPeptides[m][ii]);
									}
								}
								fprintf(fp2, ",");
							}

							for(l=0; l < allResults->ExperimentsList->Count; l++)
							{
								for(n=0; n < allResults->ExperimentsList[l]->ProteinsList->Count; n++)
								{
									//if(allResults->ExperimentsList[l]->ProteinsList[n]->accession->EndsWith(allAccessions[k]))
									if(allResults->ExperimentsList[l]->ProteinsList[n]->accession->Equals(allAccessions[k]))
									{
										imax = 0;

										dtemp = 0.0;

										for(t=0; t < allResults->ExperimentsList[l]->ProteinsList[n]->Peptides->Count; t++)
										{
											if(allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[t]->Peptide->Length == ProteinPeptides[m]->Length &&
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[t]->Peptide->EndsWith(ProteinPeptides[m]) &&
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[t]->nCharge == PeptideCharges[m] &&
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[t]->SeqMass == PeptideMass[m])
											{
												if(allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[t]->dIonscore > dtemp)
												{
													dtemp = allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[t]->dIonscore;

													imax = t;
												}
											}
										}

										nCharge = allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->nCharge;

										dSeqMass = allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->SeqMass;

										dSpecMass = allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->SpecMass;

										dSeqMass = (dSeqMass - mProton)*(double)nCharge + mProton;

										dtemp = (dSpecMass - mProton)*(double)nCharge + mProton;

										dError = (dSeqMass - dtemp)/dSeqMass*pow(10., 6);

										//set the bQuant to true as this peptide has been observed in ALL experiments

										allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->bQuant = true;

										if(0 == l)
										{
											fprintf(fp2, "%d, %10.5f,  %10.5f, %10.5f, %e, %4.1f, %d,",
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->nCharge, 
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->SeqMass,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->SpecMass, 
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->dIonscore,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->dExpect,
												dError,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->nScan);

											if(allResults->ExperimentsList->Count -1 == l)
											{
												fprintf(fp2, "\n");
											}
										}
										else if(allResults->ExperimentsList->Count -1 == l)
										{
											fprintf(fp2, "%10.5f,  %10.5f, %e, %4.1f, %d\n", 
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->SpecMass,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->dIonscore,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->dExpect,
												dError,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->nScan);
										}
										else
										{
											fprintf(fp2, "%10.5f, %10.5f, %e, %4.1f, %d,", 
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->SpecMass, 
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->dIonscore,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->dExpect,
												dError,
												allResults->ExperimentsList[l]->ProteinsList[n]->Peptides[imax]->nScan);
										}
									}
								}
							}


						}

					}

					fclose(fp2);

					break;
				}

			}

		}
	}

	delete (PeptideCharges);
	delete (ProteinPeptides);
	delete (PeptideIndices);
	delete (PeptideInAllExperiments);
	delete (PeptideMass);

	/*int c;
	c=string(c,12);
	void( );
	if (c==69)
	printf(" the under lined text is",&c);*/

	return;
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
*  A method to identify which proteins have been
*  observed in a certain percentage of all experiments
*  nDistinctPeptides - is the number of distinct peptides that a protein needs to have
*  Protein score is the cut off score for proteins
*  nExperimentThreshold is the number of Experiments the protein observed in
*  The Proteins that are found are stored in allAccessions;
*  ProteinIDnumber will be equal to the number at which that protein passed the threshold in the 
*   list of proteins.
*/

void ConsistentProteins(ExperimentCollection ^allExperiments, int nDistinctPeptides,
	float ProteinScoreThreshold, int nExperimentThreshold)
{
	ProteinList ^currentProtein;

	List <ProteinList ^> ^lProteins = gcnew List <ProteinList ^>;

	int i, j, k;

	bool bNew = false;

	//comprise the list of all proteins/accessions, allAccessions that have been
	//observed in anyone of the experiments and passes the thresholds for experiments 
	for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
	{
		for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
		{
			if( allExperiments->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences >= nDistinctPeptides &&
				allExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinScore >= ProteinScoreThreshold )
			{

				bNew = false;    //shouldn't this be true and the rest reversed?

				for(k=0; k < lProteins->Count; k++)
				{
					if(allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(lProteins[k]->Accession))
					{
						bNew = true;

						lProteins[k]->ProteinInExperiments->Add(i);		

						lProteins[k]->ExperimentTime->Add(allExperiments->ExperimentsList[i]->fExperimentTime);

						break;
					}
				}

				// Protein filtering - score and the number of distinct peptides
				if(false == bNew)
				{
					currentProtein = gcnew ProteinList;

					currentProtein->Accession    = gcnew String(allExperiments->ExperimentsList[i]->ProteinsList[j]->accession);

					currentProtein->Protein      = gcnew String(allExperiments->ExperimentsList[i]->ProteinsList[j]->description);

					currentProtein->bProteinPassed = false;

					currentProtein->ProteinInExperiments = gcnew List <int>;

					currentProtein->ExperimentTime = gcnew List <float>;

					currentProtein->ProteinInExperiments->Add(i);

					currentProtein->ExperimentTime->Add(allExperiments->ExperimentsList[i]->fExperimentTime);

					lProteins->Add(currentProtein);

					allExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinIDnumber = lProteins->Count - 1;  //the count starts from 0.

				}
			} // if( allExperiments->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences >= nDistinctPeptides &&
		}
	}
	/*
	*  check to see if the proteins
	*/

	//Proteins that passed the threshold in at least one experiment;
	// the block of the code below determines the number of
	// protein observation in non-replicate experiments (counts
	// proteins one per experimental time point).

	j = 0;

	for(k = 0; k < lProteins->Count; k++)
	{
		j = 1;     //starts with 1, since the first element has to be a new non-replicata

		for (i = 1; i < lProteins[k]->ProteinInExperiments->Count; i++)
		{
			if (lProteins[k]->ExperimentTime[i] != lProteins[k]->ExperimentTime[i - 1])
			{
				j = j + 1;
			}

		}

		if(j >= nExperimentThreshold)
		{
			lProteins[k]->bProteinPassed = true;

			//printf("%s MMHHHnot passed j = %d\n", lProteins[k]->Accession, j);
		}
		else
		{
			lProteins[k]->bProteinPassed = false;

			//printf("%s SSHHHnot pass j = %d\n", lProteins[k]->Accession, j);
		}
	}

	for(k = 0; k < lProteins->Count; k++)
	{
		if(lProteins[k]->bProteinPassed)
		{
			for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
			{
				for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
				{
					if(lProteins[k]->Accession->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->accession) )
					{
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinScore >= ProteinScoreThreshold &&
							allExperiments->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences >= nDistinctPeptides )
						{
							allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed = true;

							//printf("Protein %s passed\n", allExperiments->ExperimentsList[i]->ProteinsList[j]->accession);
						}

						break;
					}
				}

			}
		}
	}


	//printf("ListLength = %d\n", lProteins->Count);

	delete currentProtein, lProteins;

	printf("Finished with Consistent Proteins\n");

	return;
}

/*
*  a method to filter peptides and 
*  PeptideScoreThreshold is the threshold score - currently Mascot ion score
*  nExperimentThreshold number of experiments for a peptide to be observed in for
*  subsequent processing
*  Need to create a list of new allExperiments - where each peptide will be represented in
*  a protein by a single entry - instead of what is happening now
*/


void ConsistentPeptides(ExperimentCollection ^allExperiments, float PeptideScoreThreshold, 
	double dExpectationScoreThreshold, int nExperimentThreshold)
{
	PeptideList ^currentPeptide;

	List <PeptideList ^> ^lPeptides = gcnew List <PeptideList ^>;

	int i, j, l, k, q;

	int i1, j1, k1, l1, imax2;


	bool bDuplicate = false, bFound = false, bFoundProtein = false;

	bool bPeptide_ThisExperiment = false;

	float fFirstID_Elution, fLastID_Elution;     //first and last times that a peptide has been identified via MS/MS


	//comprise the list of all proteins/accessions, allAccessions that have been
	//observed in anyone of the experiments and passes the thresholds for experiments 
	for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
	{
		for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
		{
			if( allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed)
			{
				for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
				{
					//check to see if the l-th peptide is already in the list
					// if already in the list then do nothing, just add quantify true;
					// if not in the list then add it
					if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore >= PeptideScoreThreshold &&
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect <= dExpectationScoreThreshold)
					{
						bDuplicate = false;

						for(k=0; k < lPeptides->Count; k++)
						{
							if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide->Equals(lPeptides[k]->Peptide->Peptide)  &&
								allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge == lPeptides[k]->Peptide->nCharge &&
								allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass == lPeptides[k]->Peptide->SeqMass)
							{
								bDuplicate = true;

								for(q=0; q < lPeptides[k]->Peptide->ModLocations->Count; q++)
								{
									if(lPeptides[k]->Peptide->ModLocations[q] != 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
									{
										bDuplicate = false;

										//printf("Not Duplicates %s %s %d %d\n", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide,
											//lPeptides[k]->Peptide->Peptide,
											//lPeptides[k]->Peptide->ModLocations[q], allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q]);

										break;
									}
								}


								//if Duplicate true, replace the peptide if its score is higher
								if (bDuplicate)
								{
									bool bPeptideInThisExperiment = false;

									for (q = 0; q < lPeptides[k]->PeptideInExperiments->Count; q++)
									{
										//already holding this peptide in this experiment
										if (lPeptides[k]->PeptideInExperiments[q] == i)
										{
											bPeptideInThisExperiment = true;

											break;
										}
									}

									if (!bPeptideInThisExperiment)
									{
										lPeptides[k]->PeptideInExperiments->Add(i);

										lPeptides[k]->ExperimentTime->Add(allExperiments->ExperimentsList[i]->fExperimentTime);
									}

									break;
								}

							}

						} //for(k=0; k < lPeptides->Count; k++)


						// add new peptides
						if(false == bDuplicate)
						{
							currentPeptide = gcnew PeptideList;

							currentPeptide->Peptide = gcnew PeptideHolder();

							currentPeptide->Peptide->ModLocations = gcnew List <int>;

							currentPeptide->Peptide->dModMasses = gcnew List <double >;

							currentPeptide->Peptide = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l];

							currentPeptide->Peptide->ModLocations = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations;

							currentPeptide->Peptide->dModMasses = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dModMasses;

							currentPeptide->Peptide->dRetTime = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime;

							currentPeptide->PeptideInExperiments = gcnew List <int>;

							currentPeptide->ExperimentTime = gcnew List <float>;

							currentPeptide->PeptideInExperiments->Add(i);

							currentPeptide->ExperimentTime->Add(allExperiments->ExperimentsList[i]->fExperimentTime);

							lPeptides->Add(currentPeptide);

							//printf("CUrrent Length %d %s\n", lPeptides->Count, lPeptides[[lPeptides->Count-1]->Peptide->Peptide);
						}
					} //if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore >= PeptideScoreThreshold

				}

			} // if( allExperiments->ExperimentsList[i]->ProteinsList[j]->nDistinctSequences >= nDistinctPeptides &&
		}
	}

	//for a peptide to pass in overall, it has to be observed in a certain number of experiments
	//only after that set its pass variable to true

	// check for the replicate experiments - if a peptide observed
	// in several EXPERIMENTS at single time of label incorporation,
	// count the peptide once only.

	j = 1; 

	for(k=0; k < lPeptides->Count; k++)
	{
		j = 1;

		for (i = 1; i < lPeptides[k]->ExperimentTime->Count; i++)
		{
			if (lPeptides[k]->ExperimentTime[i] != lPeptides[k]->ExperimentTime[i - 1])
			{
				j = j + 1;
			}

		}


		//if(lPeptides[k]->PeptideInExperiments->Count >= nExperimentThreshold)

		if (j >= nExperimentThreshold)
		{
			lPeptides[k]->bPeptidePassed = true;

			//printf("%s PASSED , j = %d\n", lPeptides[k]->Peptide->Peptide, j);
		}
		else 
		{
			//printf("Did not pass %s\n", lPeptides[k]->Peptide->Peptide);

			lPeptides[k]->bPeptidePassed = false;

		}

	}

	// find peptides that have passed and designate them bPassedPeptide
	// if the sequence is not unique in the peptide, mark as passed the 
	// peptide with the maximum score
	int imax = 0; double dmaxScore = 0.0;

	ProteinSet ^aProtein; //To save previously passed Protein info
	PeptideHolder ^aPeptide; // To save previously passed Peptide info;
	int exp_no;
	String ^prot_name;
	imax2 = 0;
	List <int > ^count_prot = gcnew List<int >;
	for(k=0; k < lPeptides->Count; k++)
	{
		if(lPeptides[k]->bPeptidePassed /*&& lPeptides[k]->Peptide->Peptide == "YNALDLTNNGK"*/)
		{

			prot_name = lPeptides[k]->Peptide->Protein;

			for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
			{

			
				exp_no = i;

				bPeptide_ThisExperiment = false;    //this boolean variable tracks if this current peptide (which passed the experimental threshold on number of experiments) has been seen in
				                                    // in this experiment at all (with any score). If not, it will be added with the score 0, and m/z values from previous experiment.

				for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)  //over all proteins of an experiment
				{

					/*if(lPeptides[k]->Peptide->Protein->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->accession))
					{
						printf("COMMON PROTEIN %s\n", lPeptides[k]->Peptide->Protein);
					}*/

					bFoundProtein = false;

					dmaxScore = 0.0;

					fFirstID_Elution = 1000000.0;

					fLastID_Elution = -1.0;


					for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
					{
						bFound = false;

						//prot_name = allExperiments->ExperimentsList[i]->ProteinsList[j]->accession;

						//prot_name = lPeptides[k]->Peptide->Protein;


						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide->Equals(lPeptides[k]->Peptide->Peptide)  &&
							allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge == lPeptides[k]->Peptide->nCharge &&
							/*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass == lPeptides[k]->Peptide->SeqMass &&
							allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore >= PeptideScoreThreshold && 
							allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect <= dExpectationScoreThreshold)*/
							allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass == lPeptides[k]->Peptide->SeqMass)   //12.21.2016 - no need for other filters
						{	                                                                                                               // as lPeptides  lPeptides[k]->bPeptidePassed - 



							bFound = true;

							//find if there is a difference in modifications
							for(q=0; q < lPeptides[k]->Peptide->ModLocations->Count; q++)
							{

								/*printf("lPeptides[k]->Peptide->ModLocations[%d] = %d, %d\n", q,
								lPeptides[k]->Peptide->ModLocations[q], 
								allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q] );*/

								if(lPeptides[k]->Peptide->ModLocations[q] != 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
								{ 
									bFound = false;

									//imax2 = l;

									break;
								}
							}

							if(bFound)  //determines the maximum score for peptide, if it has been identified several times in
								// a single experiment
							{
								if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore > dmaxScore)
								{

									if(dmaxScore > 0.001)
									{
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->bPeptidePassed = false;   //check if this is needed!!!
									}

									imax = l;

									dmaxScore = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore;
								}

								bFoundProtein = true;  //it is enough to set it to true once per protein to
								//mean that the corresponding peptide has been found in this peptide


								if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime < fFirstID_Elution)
								{
									fFirstID_Elution = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime;
								}

								if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime > fLastID_Elution)
								{
									fLastID_Elution = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime;
								}

								bPeptide_ThisExperiment = true;
							}

						}

						//designate the peptide (if it had duplicates) with the maximum score
						// as the one that passes

						//printf("imax = %d\n", imax);
						if(bFound)
						{
							/*printf("Copies %s %s\n", lPeptides[k]->Peptide->Peptide, 
							allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide);

							printf("%s  Last Detection Elution 10.5%f %10.5f %10.5f\n", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide,
							fFirstID_Elution/60, fLastID_Elution/60, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime/60);*/
						}
					} // for(l = 0; 

					//once found the peptide inside the protein protein, exit
					if(bFoundProtein)
					{
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->bPeptidePassed = true;
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->fFirstID_Elution = fFirstID_Elution;
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->fLastID_Elution = fLastID_Elution;

						/*printf("experime = %d, Protein = %s, peptide %s  %d Score %10.5f  %10.5f %10.5f %d\n", i, allExperiments->ExperimentsList[i]->ProteinsList[j]->accession, 
							        allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->Peptide, 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->nScan, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->dIonscore,
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->fFirstID_Elution/60, 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->fLastID_Elution/60,
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[imax]->nCharge); */
						
						count_prot->Add(i);
						
						break;
					}

				} //for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)


				if(!bPeptide_ThisExperiment && lPeptides[k]->bPeptidePassed)
				{
					//printf("Peptide %s, protein %s not found in experi %d\n", lPeptides[k]->Peptide->Peptide, lPeptides[k]->Peptide->Protein, i);

					for(int ll=0; ll < allExperiments->ExperimentsList[i]->ProteinsList->Count; ll++)
					{
						if(lPeptides[k]->Peptide->Protein->Equals(allExperiments->ExperimentsList[i]->ProteinsList[ll]->accession))
						{
							/*printf("Adding Peptide %s to protein %s, exp = %d\n", lPeptides[k]->Peptide->Peptide, 
								allExperiments->ExperimentsList[i]->ProteinsList[ll]->accession, i);*/


							PeptideHolder ^ cPeptide = gcnew PeptideHolder();

							cPeptide->ModLocations = gcnew List <int>;

							cPeptide->dModMasses = gcnew List <double >;

							cPeptide->ModLocations = lPeptides[k]->Peptide->ModLocations;

							cPeptide->dModMasses = lPeptides[k]->Peptide->dModMasses;

							cPeptide->bUniquePeptide = lPeptides[k]->Peptide->bUniquePeptide;

							cPeptide->dIonscore = 0.0;
							
							cPeptide->dRetTime  = 0.0;
							
							cPeptide->Peptide   = lPeptides[k]->Peptide->Peptide;

							cPeptide->Protein   = lPeptides[k]->Peptide->Protein;
							
							cPeptide->nScan     =  lPeptides[k]->Peptide->nScan;
							
							cPeptide->nCharge   = lPeptides[k]->Peptide->nCharge;
							
							cPeptide->SeqMass  =  lPeptides[k]->Peptide->SeqMass;
							
							cPeptide->SpecMass = lPeptides[k]->Peptide->SpecMass;

							cPeptide->dIsotopes = gcnew array <double, 2>(2, 6);

							cPeptide->w = gcnew array <double, 1>(12);
							

							cPeptide->dExpect = 0.0;

							cPeptide->dRetTime = lPeptides[k]->Peptide->dRetTime;
							
							cPeptide->fFirstID_Elution = lPeptides[k]->Peptide->dRetTime;

							cPeptide->fLastID_Elution  = lPeptides[k]->Peptide->dRetTime;


							//printf("Retime is %10.5f\n", lPeptides[k]->Peptide->dRetTime);

							cPeptide->bPeptidePassed  = true;                   //this is the key for processing it later for Abundance

							cPeptide->Protein = lPeptides[k]->Peptide->Protein;

							allExperiments->ExperimentsList[i]->ProteinsList[ll]->Peptides->Add(cPeptide);


							/*printf("This Peptide  = %s , from lpep %s, the dRettime = %10.5f\n", lPeptides[k]->Peptide->Peptide,
								lPeptides[k]->Peptide->Peptide, lPeptides[k]->Peptide->dRetTime/60.);


							for(int kk = 0; kk < allExperiments->ExperimentsList[i]->ProteinsList[ll]->Peptides->Count; kk++)
							{
								 printf("Protein in %d exp %s   %s %10.5f\n", i, allExperiments->ExperimentsList[i]->ProteinsList[ll]->accession, 
									    allExperiments->ExperimentsList[i]->ProteinsList[ll]->Peptides[kk]->Peptide,
										allExperiments->ExperimentsList[i]->ProteinsList[ll]->Peptides[kk]->dRetTime/60);
							}*/

							delete cPeptide;

							break;
							//break;   // could be useful to uncheck this
						}
					}
					
				}
			}// for(i = 0; i <...)						
		}

	}


	//debug
	if (false)
	{
		for (i = 0; i < allExperiments->ExperimentsList->Count; i++)
		{
			for (j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
			{
				prot_name = allExperiments->ExperimentsList[i]->ProteinsList[j]->accession;

				for (l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
				{

					printf("Exp no: %d, Scan: %d, RetTime : %f Start Elution :%f End Elution :%f\n", i,
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan,
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dRetTime / 60.0,
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fFirstID_Elution / 60.0,
						allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fLastID_Elution / 60.0);

				}

			}

		}

		for (k = 0; k < lPeptides->Count; k++)
		{
			if (lPeptides[k]->bPeptidePassed)
			{
				printf("PASSES: %s\n", lPeptides[k]->Peptide->Peptide);
			}
		}
	}


	delete (currentPeptide);

	delete lPeptides;

	printf("Finished with Consistent Peptides\n");

	return;
}

/*
*
*   The function call the WriteProteins function to write Quant.csv file
*    per every protein.
*
*/

void WriteQuantResults(ExperimentCollection ^allExperiments,array<System::String ^> ^args)
{
	List <PeptideHolder ^> ^QuantPeptides = gcnew List <PeptideHolder ^>;

	int i, j, k, l, q, m, m1, i1;

	bool bNewPeptide = false, bThisPeptide= false;

	PeptideHolder ^currentPeptide = gcnew PeptideHolder;

	List <String ^> ^PassedAccessions = gcnew List <String ^>;

	for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
	{
		for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
		{
			if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed)
			{
				bool bNewAccession = true;

				for(k=0; k < PassedAccessions->Count; k++)
				{
					if(PassedAccessions[k]->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->accession) )
					{
						bNewAccession = false;

						break;
					}
				}

				if(bNewAccession)
				{
					//among the accession there could be ones that did not have enough peptides in them
					// even though the scores were consistently high, so remove those accessions 
					// passed peptide means it has passed the threshold in several (specified number of) experiments
					for(l=0; l <  allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
					{
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bPeptidePassed)
						{
							PassedAccessions->Add(allExperiments->ExperimentsList[i]->ProteinsList[j]->accession);

							break;
						}
					}
				} //if(bNewAccession)
			}
		}
	} //for(i=0; i < allExperiments->ExperimentsList->Count; i++) 

	//new loop, for every protein that that passed the threshold determine all
	// peptides (that passed the threshold) and print them

	for(k = 0; k < PassedAccessions->Count; k++)
	{
		List <PeptideHolder ^> ^ProteinPeptides = gcnew List <PeptideHolder ^>;

		for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
		{
			for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
			{
				if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
					allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(PassedAccessions[k]) )
				{
					for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
					{ 
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bPeptidePassed)
						{
							bool bNewPeptide = true;

							bool bDuplicate = true;

							bDuplicate = false;

							for (m = 0; m < ProteinPeptides->Count; m++)
							{
								//this if statement misses the cases with the same number of modificaitons, but on different position
								if (ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide) &&
									ProteinPeptides[m]->nCharge == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge &&
									ProteinPeptides[m]->SeqMass == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass &&
									ProteinPeptides[m]->ModLocations->Count == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations->Count)
								{

									bDuplicate = true;

									bNewPeptide = false;

									for (q = 0; q < ProteinPeptides[m]->ModLocations->Count; q++)
									{
										if (ProteinPeptides[m]->ModLocations[q] ==
											allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
										{
											//bNewPeptide = false;

											bDuplicate = true;
										}

										if (ProteinPeptides[m]->ModLocations[q] !=
											allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
										{
											//printf("Location = %d Location %d\n", ProteinPeptides[m]->ModLocations[q],
												//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q]);

											//bNewPeptide = true;

											bDuplicate = false;

											break;
										}

									} //for (q = 0; q < ProteinPeptides[m]->ModLocations->Count; q++)

								}  //if (ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide)

								if (bDuplicate)  //found the duplicate in the peptide list break out.
								{
									bNewPeptide = false;

									break;
								}
	
							}  //for (m = 0; m < ProteinPeptides->Count; m++)

							if(bNewPeptide)
							{
								ProteinPeptides->Add(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]);

								ProteinPeptides[ProteinPeptides->Count - 1]->bUniquePeptide = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bUniquePeptide;
							}
						}
					}
				}
			}
		}

		//now have the list of peptides that are used for this protein:
		// send them for printing into a file

		if(bWriteQuantCsv)
		{
			
			if(Isotope_Deconvolution==0)
				WriteAProteinFile(allExperiments, PassedAccessions[k], ProteinPeptides,args); //No MPE or NNLS
			if(Isotope_Deconvolution==1)
			{
				WriteAProteinFile1(allExperiments, PassedAccessions[k], ProteinPeptides,args);//Writing MPE info
				//WriteAProteinRate(allExperiments, PassedAccessions[k], ProteinPeptides,args,fprate,fTime,tmax);
			}
			if(Isotope_Deconvolution==2)
			{
				WriteAProteinFile2(allExperiments, PassedAccessions[k], ProteinPeptides,args);//Writing NNLS info
			}
				
			//			if(Isotope_Deconvolution==3)
			//			WriteAProteinFile3(allExperiments, PassedAccessions[k], ProteinPeptides,args);//Writing MPE info along with BFGS calculated degradation rate
			//	if(Isotope_Deconvolution==4)
			//	WriteAProteinFile4(allExperiments, PassedAccessions[k], ProteinPeptides,args);//Writing NNLS info along with BFGS calculated degradation rate

		}
		else if(bO18)
		{
			WriteAProtein_O18(allExperiments, PassedAccessions[k], ProteinPeptides,args);
		}

		delete ProteinPeptides; 
	} //for(k = 0; k < PassedAccessions->Count; k++)

	return;
}

/*
*  Writes time course results for each protein separately.
*
*   nExchblH is the number of exchangeable hydrogens in a peptide - obtained from the appropriate number in AAs
*   
*  nTimePoints - number of time points 
*/

void WriteAProteinFile(ExperimentCollection ^allExperiments, String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args)
{
	int i, j, l, q, m, nCharge;

	bool bNewPeptide = false, bThisPeptide= false;

	double dSeqMass, dSpecMass, dtemp, dError;

	double dTheorIsotopes[10];

	char szTempPeptide[2046];

	float fExchblH;

	String  ^sTemp;

	PeptideHolder ^aPeptide;

	FILE *fp;

	String ^filename4 = args[3];

	char szTempFile[2046], c;	

	m = 0; szTempFile[0] = '\0';

	if(filename4 != "")
	{
		for(m=0; m<filename4->Length; m++)
			szTempFile[m] = filename4[m];

		szTempFile[m++] = '\\';
	}

	for(l=0; l < sProtein->Length; l++)
	{
		c = sProtein[l];

		if(c != '|' && c != ',' && c != '\n' && c != '-' &&
			c != '\t' && c != ' ')
		{
			szTempFile[m] = sProtein[l];

			m++;
		}
	}

	szTempFile[m] = '\0';

	strcat(szTempFile, ".Quant.csv");

	fp = fopen(szTempFile, "w");

	if(NULL == fp)
	{
		printf("Cannot write to %s\n", szTempFile);

		exit (1);
	}

	TheoreticalIsotopeCalculator ^TheoIsotopePeaks = gcnew TheoreticalIsotopeCalculator();

	array <float> ^fFourierIsotopes = gcnew array <float> (6);


	//print the first line of the csv file
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{
		sTemp = allExperiments->ExperimentsList[i]->sExperimentFile->Substring(0,
			allExperiments->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) ;

		if(0 == i)
		{
			//fprintf(fp, ", , , , , , , , ,%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, ", , , , , , , , , ,%s, , , , , , , , , , , , ,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			//fprintf(fp, "%s, , , , , , , , , , \n", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , \n", sTemp);
		}
		else
		{
			//fprintf(fp, "%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , , ,", sTemp);
		}
	}

	sTemp = gcnew String("SpecMass, IonScore, Expectn, Error(ppm), Scan, I0, I1, I2, I3, I4, I5, Start Elution (min), End Elution (min)");

	//print the file header
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{

		if(allExperiments->ExperimentsList->Count - 1 == i && 0 == i) // in the case when there is a single experiment
		{
			//fprintf(fp, "\nPeptide, Charge, SeqMass,  M0, M1, M2, M3, M4, M5, %s\n", sTemp);

			fprintf(fp, "\nPeptide, Exchangeable Hydrogens, Charge, SeqMass,  M0, M1, M2, M3, M4, M5, %s\n", sTemp);
		}
		else if(0 == i)
		{
			//fprintf(fp, "Peptide, Charge, SeqMass, M0, M1, M2, M3, M4, M5, %s,", sTemp);

			fprintf(fp, "Peptide, Exchangeable Hydrogens, Charge, SeqMass, M0, M1, M2, M3, M4, M5, %s,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, "%s\n", sTemp);
		}
		else
		{
			fprintf(fp, "%s,", sTemp);
		}
	}

	//print results for every peptide in the protein

	for(m = 0; m < ProteinPeptides->Count; m++)
	{
		if(ProteinPeptides[m]->ModLocations->Count == 0)
		{
			fprintf(fp, "%s,", ProteinPeptides[m]->Peptide);
		}
		else
		{
			for(int ii=0; ii < ProteinPeptides[m]->Peptide->Length; ii++)
			{
				bool bModified = false, bCterm = false;   //for O18 Mascot writes the modification outside of peptide length

				for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
				{
					if(ProteinPeptides[m]->ModLocations[q] == ii + 1)   //the mod position is counted starting from 1.
					{
						bModified = true;

						break;
					}
					else if (ProteinPeptides[m]->ModLocations[q] == ProteinPeptides[m]->Peptide->Length + 1 )
					{
						bCterm = true;

						break;
					}
				}

				if(bModified)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else if(bCterm == true && ii == ProteinPeptides[m]->Peptide->Length - 1)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else
				{
					fprintf(fp, "%c", ProteinPeptides[m]->Peptide[ii]);
				}
			}
			fprintf(fp, ",");
		}

		//compute the theoretical isotope

		//Isotopes ^isotope = gcnew Isotopes();
		

		szTempPeptide[0] = '\0';

		for(q=0; q < ProteinPeptides[m]->Peptide->Length; q++)
		{
			szTempPeptide[q] = ProteinPeptides[m]->Peptide[q];
		}

		szTempPeptide[q] = '\0';

		//printf("Seqence = %s\n", szTempPeptide);

		printf("First Six Theoretical Isotopes for %s\n", ProteinPeptides[m]->Peptide);

		//isotope->ComputeIsotopePeaks(szTempPeptide, dTheorIsotopes);

		TheoIsotopePeaks ->StartToComputeSequenceIsotopes(ProteinPeptides[m]->Peptide, fFourierIsotopes);

		for(int ii=0; ii < 10 && ii < fFourierIsotopes->Length; ii++)
		{
			dTheorIsotopes[ii] = 100.*fFourierIsotopes[ii]; 
		}

		//delete isotope;

		aPeptide = gcnew PeptideHolder;

		for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
		{
			bool bFoundPeptide = false;

			for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)    //over proteins in experiment
			{
				if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
					allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(sProtein) )
				{
					for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)    //over peptides of a protein
					{ 
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bPeptidePassed &&
							ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide) &&
							ProteinPeptides[m]->nCharge == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge &&
							ProteinPeptides[m]->SeqMass == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass &&
							ProteinPeptides[m]->ModLocations->Count == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations->Count)
						{
							bool bThisPeptide = true;

							for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
							{
								if(ProteinPeptides[m]->ModLocations[q] != 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
								{ 
									bThisPeptide = false;

									break;
								}
							}

							if(bThisPeptide)
							{
								bFoundPeptide = true;

								nCharge = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge;

								dSeqMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass;

								dSpecMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass;

								dSeqMass = (dSeqMass - mProton)*(double)nCharge + mProton;

								dtemp = (dSpecMass - mProton)*(double)nCharge + mProton;

								dError = (dSeqMass - dtemp)/dSeqMass*pow(10., 6);		

								fExchblH = 0.0;

								for(int iH = 0; iH < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide->Length; iH++)
								{
									fExchblH = fExchblH + nHAA[allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide[iH]];
								}

								if(0 == i)
								{
									fprintf(fp, "%d, %d, %10.5f,",
										(int)fExchblH, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass);

									//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fTheoreticIsotopes =
										//gcnew array <float> (6);
									//print theoretical isotopes											
									for(q=0; q < 6; q++)
									{
										//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fTheoreticIsotopes[q] =
											//(float)dTheorIsotopes[q];

										fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
									}

									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,",
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									for(q = 0; q < 6; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);

									if(allExperiments->ExperimentsList->Count -1 == i) // in the case when there is only one experiment
									{
										fprintf(fp, "\n");
									}
								}
								else if(allExperiments->ExperimentsList->Count -1 == i)
								{
									fprintf(fp, "%10.5f,  %10.5f, %e, %4.1f, %d, ", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2 - 1; q++)
									for(q = 0; q < 5; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}



									fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]);

									fprintf(fp,"%f, %f\n",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);

									//fprintf(fp, "\n");
								}
								else
								{
									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									for(q = 0; q < 6; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);

								}
							}

						}
					} //for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
				}
			}

			//this peptide was not observed in all experiments;
			if(false == bFoundPeptide)
			{
				if(0 == i)
				{

					fExchblH = 0;

					for(int iH = 0; iH < ProteinPeptides[m]->Peptide->Length; iH++)
					{
						fExchblH = fExchblH + nHAA[ProteinPeptides[m]->Peptide[iH]];
					}

					fprintf(fp, "%d, %d, %10.5f,   ",
						(int)fExchblH, ProteinPeptides[m]->nCharge, 
						ProteinPeptides[m]->SeqMass);

					//print theoretical isotopes											
					for(q=0; q < 6; q++)
					{
						fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
					}

					fprintf(fp, ", , , , , , , , , , , , ,");
				}
				else if(allExperiments->ExperimentsList->Count -1 == i)
				{
					fprintf(fp, ", , , , , , , , , , ,  \n");
				}
				else
				{
					fprintf(fp, ", , , , , , , , , , , , , ");
				}
			}

		} // for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	} //for(m = 0; m < ProteinPeptides->Count; m++)

	fclose(fp);

	delete(TheoIsotopePeaks);  delete(fFourierIsotopes);

	return;

}

/*
*
*   Writing MPE informatio
*   Currently this function is used to print out resutls, 03/28/2017
*   Writes the Quant.csv files for each protein.
*
*
*/
void WriteAProteinFile1(ExperimentCollection ^allExperiments, String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args)
{
	int i, j, l, q, m, nCharge, nret, FoundPeptide, i1, jtemp, ltemp;

	bool bNewPeptide = false, bThisPeptide= false;

	double mass_sum, intensity_sum, frac_intensity, initial_label;

	double factor[6], mass_sum2;

	double w[12], intensity_sum2;

	char szTempProt[2024];


	String ^szTempMean;

	double dSeqMass, dSpecMass, dtemp, dError;

	double dTheorIsotopes[10];

	float fExchblH;

	String  ^sTemp;

	String ^sTemp2;

	PeptideHolder ^aPeptide;

	String ^filename4 = args[3]; //command line output 

	char szTempFile[2046], c;

	m = 0; szTempFile[0] = '\0';

	int ntimepoints = allExperiments->ExperimentsList->Count;

	float *MPE_vector = (float*)calloc(ntimepoints,sizeof(float));

	//float *meanMPE_vector = (float*)calloc(ntimepoints,sizeof(float));

	//float *medianMPE_vector = (float*)calloc(ntimepoints,sizeof(float));

	FILE *fp, *fp2;

	std::string templine, tUnits;
	//	array<float> ^MPE_vector = gcnew array <float> (ntimepoints);

	float tmax;

	if(NULL == MPE_vector)
	{
		printf("Cannot Allocate Memory for Response\n");
		exit (1);
	}


	TheoreticalIsotopeCalculator ^TheoIsotopePeaks = gcnew TheoreticalIsotopeCalculator();

	array <float> ^fFourierIsotopes = gcnew array <float> (6);

	if(filename4 != "")
	{
		for(m=0; m<filename4->Length; m++)
			szTempFile[m] = filename4[m];

		szTempFile[m++] = '\\';
	}

	for(l=0; l < sProtein->Length; l++)
	{
		c = sProtein[l];

		if(c != '|' && c != ',' && c != '\n' && c != '-' &&
			c != '\t' && c != ' ')
		{
			szTempFile[m] = sProtein[l];

			m++;
		}
	}

	szTempFile[m] = '\0';

	strcat(szTempFile, ".Quant.csv");

	fp = fopen(szTempFile, "w");

	if(NULL == fp)
	{
		printf("Cannot write to %s\n", szTempFile);

		exit (1);
	}

	//
	for(j = 0; j < allExperiments->ExperimentsList[0]->ProteinsList->Count; j++)    //over proteins over the experiment at the zeroth time point
	{
		if(allExperiments->ExperimentsList[0]->ProteinsList[j]->bProteinPassed &&
			allExperiments->ExperimentsList[0]->ProteinsList[j]->accession->Equals(sProtein) )
		{

			szTempProt[0] = '\0';

			//if protein description contains a comma, remove it, as it confuses csv
			for(l=0; l < allExperiments->ExperimentsList[0]->ProteinsList[j]->description->Length; l++)
			{
				if(allExperiments->ExperimentsList[0]->ProteinsList[j]->description[l] == ',')
				{
					szTempProt[l] = ' ';
				}
				else
				{
					szTempProt[l] = allExperiments->ExperimentsList[0]->ProteinsList[j]->description[l];
				}
			}

			szTempProt[l] = '\0';


			//print the first and second line of the Quant.csv file
			fprintf(fp,"%s\n",allExperiments->ExperimentsList[0]->ProteinsList[j]->accession);
			fprintf(fp,"%s\n",szTempProt);
		}
	}
	//

	//print the third line of the Quant.csv file
	for (i = 0; i < allExperiments->ExperimentsList->Count; i++)
	{
		sTemp = allExperiments->ExperimentsList[i]->sExperimentFile->Substring(0,
			allExperiments->ExperimentsList[i]->sExperimentFile->LastIndexOf("."));

		if (0 == i)
		{
			//fprintf(fp, ", , , , , , , , ,%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, ", , , , , , , , , , , ,%s, , , , , , , , , , , , , , , ,", sTemp);
		}
		else if (allExperiments->ExperimentsList->Count - 1 == i)
		{
			//fprintf(fp, "%s, , , , , , , , , , \n", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , , , , , \n", sTemp);
		}
		else
		{
			//fprintf(fp, "%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , , , , ,", sTemp);
		}
	}

	sTemp = gcnew String("SpecMass, IonScore, Expectn, Error(ppm), Scan, I0, I1, I2, I3, I4, I5, Start Elution (min), End Elution (min), I0 Peak Width, Total Labeling, Net Labeling");
	sTemp2 = gcnew String("SpecMass, IonScore, Expectn, Error(ppm), Scan, I0, I1, I2, I3, I4, I5, Start Elution (min), End Elution (min), I0 Peak Width, Total Labeling, Net Labeling");//, Degradation Rate (BFGS), Synthesis Rate (BFGS)");
	//print the file header
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{

		if(allExperiments->ExperimentsList->Count - 1 == i && 0 == i) // in the case when there is a single experiment
		{
			//fprintf(fp, "\nPeptide, Charge, SeqMass,  M0, M1, M2, M3, M4, M5, %s\n", sTemp);

			//fprintf(fp, "\nPeptide, Exchangeable Hydrogens, Charge, SeqMass,  M0, M1, M2, M3, M4, M5, Total Labeling, %s\n", sTemp2);

			fprintf(fp, "\nPeptide, UniqueToProtein, Exchangeable Hydrogens, Charge,  SeqMass,  M0, M1, M2, M3, M4, M5, Total Labeling, %s\n", sTemp2);
		}
		else if(0 == i)
		{
			//fprintf(fp, "Peptide, Exchangeable Hydrogens, Charge, SeqMass, M0, M1, M2, M3, M4, M5, Total Labeling, %s,", sTemp);
			fprintf(fp, "Peptide, UniqueToProtein, Exchangeable Hydrogens, Charge, SeqMass, M0, M1, M2, M3, M4, M5, Total Labeling, %s,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, "%s\n", sTemp2);
		}
		else
		{
			fprintf(fp, "%s,", sTemp);
		}
	}


	for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
	{
		for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)    //over proteins an experiment
		{
			if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
				allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(sProtein) )
			{
				allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE = 0;
				allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE = 0;
				allExperiments->ExperimentsList[i]->ProteinsList[j]->NetLabeling = gcnew List<float>;
				allExperiments->ExperimentsList[i]->ProteinsList[j]->PeptideIndex = gcnew List<int>;
			}
		}
	}

	jtemp = -1;

	for(m = 0; m < ProteinPeptides->Count; m++)
	{
		FoundPeptide = 0;

		for(i=0;i<ntimepoints;i++)
		{
			MPE_vector[i] = 0.00;
		}

		if(ProteinPeptides[m]->ModLocations->Count == 0)
		{
			fprintf(fp, "%s,", ProteinPeptides[m]->Peptide);
		}
		else
		{
			for (int ii = 0; ii < ProteinPeptides[m]->Peptide->Length; ii++)
			{
				bool bModified = false, bCterm = false;   //for O18 Mascot writes the modification outside of peptide length

				for (q = 0; q < ProteinPeptides[m]->ModLocations->Count; q++)
				{
					if (ProteinPeptides[m]->ModLocations[q] == ii + 1)   //the mod position is counted starting from 1.
					{
						bModified = true;

						break;
					}
					else if (ProteinPeptides[m]->ModLocations[q] == ProteinPeptides[m]->Peptide->Length + 1)
					{
						bCterm = true;

						break;
					}
				}

				if (bModified)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else if (bCterm == true && ii == ProteinPeptides[m]->Peptide->Length - 1)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else
				{
					fprintf(fp, "%c", ProteinPeptides[m]->Peptide[ii]);
				}

			}
			fprintf(fp, ",");
		}  //if (ProteinPeptides[m]->ModLocations->Count == 0), else

		printf("WRP1:: First Six Theoretical Isotopes for %s\n", ProteinPeptides[m]->Peptide);

		TheoIsotopePeaks->StartToComputeSequenceIsotopes(ProteinPeptides[m]->Peptide, fFourierIsotopes);

		for(int ii=0; ii < 10 && ii < fFourierIsotopes->Length; ii++)
		{
			dTheorIsotopes[ii] = 100.*fFourierIsotopes[ii]; 

			if(ii < 6)
			{
				printf("   %10.5f\n", dTheorIsotopes[ii]);
			}
		}

		String ^prot_name;

		aPeptide = gcnew PeptideHolder;

		for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
		{
			bool bFoundPeptide = false;


			for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)    //over proteins an experiment
			{
				prot_name = allExperiments->ExperimentsList[i]->ProteinsList[j]->accession;

				if((allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
					allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(sProtein)))
				{

					for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)    //over peptides of a protein
					{ 
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bPeptidePassed &&
							ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide) &&
							ProteinPeptides[m]->nCharge == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge &&
							ProteinPeptides[m]->SeqMass == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass &&
							ProteinPeptides[m]->ModLocations->Count == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations->Count)
						{
							bool bThisPeptide = true;

							for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
							{
								if(ProteinPeptides[m]->ModLocations[q] != 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
								{ 
									bThisPeptide = false;

									break;
								}
							}

							if(bThisPeptide)
							{
								bFoundPeptide = true;

								FoundPeptide++;
								
								nCharge = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge;

								dSeqMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass;

								dSpecMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass;

								dSeqMass = (dSeqMass - mProton)*(double)nCharge + mProton;

								dtemp = (dSpecMass - mProton)*(double)nCharge + mProton;

								dError = (dSeqMass - dtemp)/dSeqMass*pow(10., 6);		

								fExchblH = 0.0;

								for(int iH = 0; iH < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide->Length; iH++)
								{
									fExchblH = fExchblH + nHAA[allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide[iH]];
								}
								
								if(0 == i)
								{
									if (allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bUniquePeptide)
									{
										fprintf(fp, "Yes, %d, %d, %10.5f,",
											(int)fExchblH, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge,
											allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass);
									}
									else
									{
										fprintf(fp, "No, %d, %d, %10.5f,",
											(int)fExchblH, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge,
											allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass);
									}
									//fprintf(fp, "%d, %d, %10.5f,",
										//(int)fExchblH, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge, 
										//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass);

									//print theoretical isotopes											
									
									
									//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fTheoreticIsotopes = gcnew array <float> (6);

									for(q=0; q < 6; q++)
									{
										//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fTheoreticIsotopes[q] =
											//(float)dTheorIsotopes[q];

										fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
									}

									mass_sum = 0;
									for(q=1; q < 6; q++)
									{
										mass_sum += q*dTheorIsotopes[q];
									}

									fprintf(fp, "%5.3f,",mass_sum);

									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,",
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)

									intensity_sum = 0;
									for(q = 0; q < 6; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);

										intensity_sum += allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}


									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);

									frac_intensity = 0;
									intensity_sum2 = 0;

									/*
									for(q = 1; q < 6; q++)
									{																				
									frac_intensity += q * (allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum);
									}
									*/

									/*
									mass_sum2 = 0;
									for(q=0; q < 6; q++)
									{
									mass_sum2 += dTheorIsotopes[q];
									}
									*/



									if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0)
									{
										factor[0]=factor[1]=factor[2]=factor[3]=factor[4]=factor[5]=1;
									}
									else
									{
										factor[2]=factor[3]=factor[4]=factor[5]=0;
										factor[0]=factor[1]=1;
									}

									for(q=0;q<6;q++)
									{
										intensity_sum2 += factor[q]*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0)
									{
										factor[0]=0;
										factor[1]=1;
										factor[2]=2;
										factor[3]=3;
										factor[4]=4;
										factor[5]=5;
									}
									else
									{
										factor[2]=factor[3]=factor[4]=factor[5]=0;
										factor[0]=0;
										factor[1]=1;
									}

									if(intensity_sum2==0)
									{
										frac_intensity = 0;
									}
									else
									{
										for(q=0;q<6;q++)
										{
											frac_intensity += factor[q] * (allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum2);
										}
									}
									fprintf(fp,"%f,",2*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->w[0]);
									fprintf(fp,"%f,",frac_intensity);


									initial_label = frac_intensity; //keeps the total intesity at 0;

									fprintf(fp,"%f,",frac_intensity-initial_label);


									MPE_vector[i] = frac_intensity-initial_label;

									//if(!(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 2]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 3]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 4]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]==0))
									//{
									//if(NULL!=MPE_vector[i])
									//{
									allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE = allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE +  MPE_vector[i];
									allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE = allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE +  1;
									allExperiments->ExperimentsList[i]->ProteinsList[j]->NetLabeling->Add(MPE_vector[i]);
									//allExperiments->ExperimentsList[i]->ProteinsList[j]->TimePointsAvailable->Add(i);
									allExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinIndex = j;
									allExperiments->ExperimentsList[i]->ProteinsList[j]->PeptideIndex->Add(l);

									jtemp = j;
									ltemp = l;

									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->NetLabeling = MPE_vector[i];

									//}
									//}

									if(allExperiments->ExperimentsList->Count -1 == i) // in the case when there is only one experiment
									{
										fprintf(fp, "\n");
									}
								}
								else if(allExperiments->ExperimentsList->Count -1 == i)
								{
									fprintf(fp, "%10.5f,  %10.5f, %e, %4.1f, %d, ", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2 - 1; q++)
									intensity_sum = 0;
									for(q = 0; q < 5; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);

										intensity_sum += allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]);

									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);

									intensity_sum += allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5];

									frac_intensity = 0;
									intensity_sum2 = 0;

									/*
									for(q = 1; q < 5; q++)
									{																				
									frac_intensity += q*(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum);
									}
									frac_intensity += 5*(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]/intensity_sum);
									*/

									/*
									mass_sum2 = 0;
									for(q=0; q < 6; q++)
									{
									mass_sum2 += dTheorIsotopes[q];
									}
									*/



									if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0)
									{
										factor[0]=factor[1]=factor[2]=factor[3]=factor[4]=factor[5]=1;
									}
									else
									{
										factor[2]=factor[3]=factor[4]=factor[5]=0;
										factor[0]=factor[1]=1;
									}

									for(q=0;q<6;q++)
									{
										intensity_sum2 += factor[q]*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0)
									{
										factor[0]=0;
										factor[1]=1;
										factor[2]=2;
										factor[3]=3;
										factor[4]=4;
										factor[5]=5;
									}
									else
									{
										factor[2]=factor[3]=factor[4]=factor[5]=0;
										factor[0]=0;
										factor[1]=1;

									}


									if(intensity_sum2==0)
									{
										frac_intensity = 0;
									}
									else
									{
										for(q=0;q<6;q++)
										{
											frac_intensity += factor[q] * (allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum2);
										}  
									}

									fprintf(fp,"%f,",2*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->w[0]);

									fprintf(fp,"%f,",frac_intensity);
									fprintf(fp,"%f,",frac_intensity-initial_label);

									MPE_vector[i] = frac_intensity-initial_label;


									//if(!(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 2]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 3]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 4]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]==0))
									//{
									//if(NULL!=MPE_vector[i])
									//{
									allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE = allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE +  MPE_vector[i];
									allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE = allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE +  1;
									allExperiments->ExperimentsList[i]->ProteinsList[j]->NetLabeling->Add(MPE_vector[i]);
									allExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinIndex = jtemp;
									allExperiments->ExperimentsList[i]->ProteinsList[j]->PeptideIndex->Add(ltemp);

									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->NetLabeling = MPE_vector[i];
									//}
									//}

									//fprintf(fp, "\n");
								}
								else
								{
									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									intensity_sum = 0;
									for(q = 0; q < 6; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);

										intensity_sum += allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}


									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);


									frac_intensity = 0;
									intensity_sum2 = 0;
									/*
									for(q = 1; q < 6; q++)
									{																				
									frac_intensity += q*(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum);
									}
									*/

									/*
									mass_sum2 = 0;
									for(q=0; q < 6; q++)
									{
									mass_sum2 += dTheorIsotopes[q];
									}
									*/

									if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0)
									{
										factor[0]=factor[1]=factor[2]=factor[3]=factor[4]=factor[5]=1;
									}
									else
									{
										factor[2]=factor[3]=factor[4]=factor[5]=0;
										factor[0]=factor[1]=1;
									}

									for(q=0;q<6;q++)
									{
										intensity_sum2 += factor[q]*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0)
									{
										factor[0]=0;
										factor[1]=1;
										factor[2]=2;
										factor[3]=3;
										factor[4]=4;
										factor[5]=5;
									}
									else
									{
										factor[2]=factor[3]=factor[4]=factor[5]=0;
										factor[0]=0;
										factor[1]=1;
									}

									if(intensity_sum2==0)
									{
										frac_intensity = 0;
									}
									else
									{
										for(q=0;q<6;q++)
										{
											frac_intensity += factor[q] * (allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum2);
										}

									}

									fprintf(fp,"%f,",2*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->w[0]);

									fprintf(fp,"%f,",frac_intensity);
									fprintf(fp,"%f,",frac_intensity-initial_label);

									MPE_vector[i] = frac_intensity-initial_label;

									//if(!(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 1]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 2]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 3]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 4]==0 && allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]==0))
									//{
									//if(NULL!=MPE_vector[i])
									//{
									allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE = allExperiments->ExperimentsList[i]->ProteinsList[j]->sumMPE +  MPE_vector[i];
									allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE = allExperiments->ExperimentsList[i]->ProteinsList[j]->countMPE +  1;
									allExperiments->ExperimentsList[i]->ProteinsList[j]->NetLabeling->Add(MPE_vector[i]);
									allExperiments->ExperimentsList[i]->ProteinsList[j]->ProteinIndex = jtemp;
									allExperiments->ExperimentsList[i]->ProteinsList[j]->PeptideIndex->Add(ltemp);

									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->NetLabeling = MPE_vector[i];
									//}
									//}
								}
							}

						}
					} //for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)

				} //if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
				//allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(sProtein) )

			}//for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)

			//this peptide was not observed in all experiments;
			if(false == bFoundPeptide)
			{
				if(0 == i)
				{

					fExchblH = 0;

					for(int iH = 0; iH < ProteinPeptides[m]->Peptide->Length; iH++)
					{
						fExchblH = fExchblH + nHAA[ProteinPeptides[m]->Peptide[iH]];
					}
					
					//if (allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bUniquePeptide)
					if(ProteinPeptides[m]->bUniquePeptide)
					{
						fprintf(fp, "Yes, %d, %d, %10.5f,   ",
							(int)fExchblH, ProteinPeptides[m]->nCharge,
							ProteinPeptides[m]->SeqMass);
					}
					else
					{
						fprintf(fp, "No, %d, %d, %10.5f,   ",
							(int)fExchblH, ProteinPeptides[m]->nCharge,
							ProteinPeptides[m]->SeqMass);
					}

					//fprintf(fp, "%d, %d, %10.5f,   ",
						//(int)fExchblH, ProteinPeptides[m]->nCharge, 
						//ProteinPeptides[m]->SeqMass);

					//print theoretical isotopes	
					//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fTheoreticIsotopes = gcnew array <float> (6);

					for(q=0; q < 6; q++)
					{
						//allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fTheoreticIsotopes[q] =
							//				(float)dTheorIsotopes[q];

						fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
					}

					mass_sum = 0;
					for(q=1; q < 6; q++)
					{
						mass_sum += q*dTheorIsotopes[q];
					}

					fprintf(fp, "%5.3f,",mass_sum);


					fprintf(fp, " , , , , , , , , , , , , , , , , ");
				}
				else
				{
					//Exceptional case : where peptide info has been added only for this time

					fprintf(fp, ", , , , , , , , , , , , , , , ,");
				}
			}

		} // for(i=0; i < allExperiments->ExperimentsList->Count; i++)


		fprintf(fp, "\n");

	} //for(m = 0; m < ProteinPeptides->Count; m++)

	if(MPE_vector)
	{
		delete[] MPE_vector;
	}


	fclose(fp);

	delete(TheoIsotopePeaks); delete(fFourierIsotopes);

	return;

}


/*Writing NNLS information*/
void WriteAProteinFile2(ExperimentCollection ^allExperiments, String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args)
{
	int i, j, l, q, m, nCharge;

	bool bNewPeptide = false, bThisPeptide= false;

	double dSeqMass, dSpecMass, dtemp, dError;

	double dTheorIsotopes[10];

	char szTempPeptide[2046];

	char szTempProt[2024];
	//

	float fExchblH;

	String  ^sTemp;

	PeptideHolder ^aPeptide;

	FILE *fp;

	String ^filename4 = args[3];

	char szTempFile[2046], c;

	m = 0; szTempFile[0] = '\0';

	if(filename4 != "")
	{
		for(m=0; m<filename4->Length; m++)
			szTempFile[m] = filename4[m];

		szTempFile[m++] = '\\';
	}

	for(l=0; l < sProtein->Length; l++)
	{
		c = sProtein[l];

		if(c != '|' && c != ',' && c != '\n' && c != '-' &&
			c != '\t' && c != ' ')
		{
			szTempFile[m] = sProtein[l];

			m++;
		}
	}

	szTempFile[m] = '\0';

	strcat(szTempFile, ".Quant.csv");

	fp = fopen(szTempFile, "w");

	if(NULL == fp)
	{
		printf("Cannot write to %s\n", szTempFile);

		exit (1);
	}

	TheoreticalIsotopeCalculator ^TheoIsotopePeaks = gcnew TheoreticalIsotopeCalculator();

	array <float> ^fFourierIsotopes = gcnew array <float> (6);

	//
	for(j = 0; j < allExperiments->ExperimentsList[0]->ProteinsList->Count; j++)    //over proteins over the experiment at the zeroth time point
	{
		if(allExperiments->ExperimentsList[0]->ProteinsList[j]->bProteinPassed &&
			allExperiments->ExperimentsList[0]->ProteinsList[j]->accession->Equals(sProtein) )
		{

			szTempProt[0] = '\0';

			//if protein description contains a comma, remove it, as it confuses csv
			for(l=0; l < allExperiments->ExperimentsList[0]->ProteinsList[j]->description->Length; l++)
			{
				if(allExperiments->ExperimentsList[0]->ProteinsList[j]->description[l] == ',')
				{
					szTempProt[l] = ' ';
				}
				else
				{
					szTempProt[l] = allExperiments->ExperimentsList[0]->ProteinsList[j]->description[l];
				}
			}

			szTempProt[l] = '\0';


			//print the first and second line of the Quant.csv file
			fprintf(fp,"%s\n",allExperiments->ExperimentsList[0]->ProteinsList[j]->accession);
			fprintf(fp,"%s\n",szTempProt);
		}
	}
	//


	//print the third line of the csv file
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{
		sTemp = allExperiments->ExperimentsList[i]->sExperimentFile->Substring(0,
			allExperiments->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) ;

		if(0 == i)
		{
			//fprintf(fp, ", , , , , , , , ,%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, ", , , , , , , , , ,%s, , , , , , , , , , , , , ,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			//fprintf(fp, "%s, , , , , , , , , , \n", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , ,\n", sTemp);
		}
		else
		{
			//fprintf(fp, "%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , , ,", sTemp);
		}
	}

	sTemp = gcnew String("SpecMass, IonScore, Expectn, Error(ppm), Scan, I0, I1, I2, I3, I4, I5, Start Elution (min), End Elution (min), I0 Peak Width");

	//print the file header
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{

		if(allExperiments->ExperimentsList->Count - 1 == i && 0 == i) // in the case when there is a single experiment
		{
			//fprintf(fp, "\nPeptide, Charge, SeqMass,  M0, M1, M2, M3, M4, M5, %s\n", sTemp);

			fprintf(fp, "\nPeptide, Exchangeable Hydrogens, Charge, SeqMass,  M0, M1, M2, M3, M4, M5, %s\n", sTemp);
		}
		else if(0 == i)
		{
			//fprintf(fp, "Peptide, Charge, SeqMass, M0, M1, M2, M3, M4, M5, %s,", sTemp);

			fprintf(fp, "Peptide, Exchangeable Hydrogens, Charge, SeqMass, M0, M1, M2, M3, M4, M5, %s,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, "%s\n", sTemp);
		}
		else
		{
			fprintf(fp, "%s,", sTemp);
		}
	}

	//print results for every peptide in the protein

	for(m = 0; m < ProteinPeptides->Count; m++)
	{
		if(ProteinPeptides[m]->ModLocations->Count == 0)
		{
			fprintf(fp, "%s,", ProteinPeptides[m]->Peptide);
		}
		else
		{
			for(int ii=0; ii < ProteinPeptides[m]->Peptide->Length; ii++)
			{
				bool bModified = false, bCterm = false;   //for O18 Mascot writes the modification outside of peptide length

				for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
				{
					if(ProteinPeptides[m]->ModLocations[q] == ii + 1)   //the mod position is counted starting from 1.
					{
						bModified = true;

						break;
					}
					else if (ProteinPeptides[m]->ModLocations[q] == ProteinPeptides[m]->Peptide->Length + 1 )
					{
						bCterm = true;

						break;
					}
				}

				if(bModified)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else if(bCterm == true && ii == ProteinPeptides[m]->Peptide->Length - 1)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else
				{
					fprintf(fp, "%c", ProteinPeptides[m]->Peptide[ii]);
				}
			}
			fprintf(fp, ",");
		}

		//compute the theoretical isotope

		//Isotopes ^isotope = gcnew Isotopes();

		
		szTempPeptide[0] = '\0';

		for(q=0; q < ProteinPeptides[m]->Peptide->Length; q++)
		{
			szTempPeptide[q] = ProteinPeptides[m]->Peptide[q];
		}

		szTempPeptide[q] = '\0';

		TheoIsotopePeaks ->StartToComputeSequenceIsotopes(ProteinPeptides[m]->Peptide, fFourierIsotopes);

		for(int ii=0; ii < 10 && ii < fFourierIsotopes->Length; ii++)
		{
			dTheorIsotopes[ii] = 100.*fFourierIsotopes[ii]; 
		}
		

		//printf("Seqence = %s\n", szTempPeptide);

		printf("WRP2: First Six Theoretical Isotopes for %s\n", ProteinPeptides[m]->Peptide);

		//isotope->ComputeIsotopePeaks(szTempPeptide, dTheorIsotopes);

		//delete isotope;

		aPeptide = gcnew PeptideHolder;

		for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
		{
			bool bFoundPeptide = false;

			for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)    //over proteins an experiment
			{
				if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
					allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(sProtein) )
				{
					for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)    //over peptides of a protein
					{ 
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bPeptidePassed &&
							ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide) &&
							ProteinPeptides[m]->nCharge == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge &&
							ProteinPeptides[m]->SeqMass == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass &&
							ProteinPeptides[m]->ModLocations->Count == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations->Count)
						{
							bool bThisPeptide = true;

							for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
							{
								if(ProteinPeptides[m]->ModLocations[q] != 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
								{ 
									bThisPeptide = false;

									break;
								}
							}

							if(bThisPeptide)
							{
								bFoundPeptide = true;

								nCharge = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge;

								dSeqMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass;

								dSpecMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass;

								dSeqMass = (dSeqMass - mProton)*(double)nCharge + mProton;

								dtemp = (dSpecMass - mProton)*(double)nCharge + mProton;

								dError = (dSeqMass - dtemp)/dSeqMass*pow(10., 6);		

								fExchblH = 0.0;

								for(int iH = 0; iH < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide->Length; iH++)
								{
									fExchblH = fExchblH + nHAA[allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide[iH]];
								}

								if(0 == i)
								{
									fprintf(fp, "%d, %d, %10.5f,",
										(int)fExchblH, allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass);

									//print theoretical isotopes											
									for(q=0; q < 6; q++)
									{
										fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
									}

									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,",
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									for(q = 0; q < 6; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);
									fprintf(fp,"%f,",2*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->w[0]);

									if(allExperiments->ExperimentsList->Count -1 == i) // in the case when there is only one experiment
									{
										fprintf(fp, "\n");
									}
								}
								else if(allExperiments->ExperimentsList->Count -1 == i)
								{
									fprintf(fp, "%10.5f,  %10.5f, %e, %4.1f, %d, ", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2 - 1; q++)
									for(q = 0; q < 5; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}



									fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]);

									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);
									fprintf(fp,"%f\n",2*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->w[0]);

									//fprintf(fp, "\n");
								}
								else
								{
									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									//for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									for(q = 0; q < 6; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp,"%f, %f,",allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fStartElution,allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->fEndElution);
									fprintf(fp,"%f,",2*allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->w[0]);

								}
							}

						}
					} //for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
				}
			}

			//this peptide was not observed in all experiments;
			if(false == bFoundPeptide)
			{
				if(0 == i)
				{

					fExchblH = 0;

					for(int iH = 0; iH < ProteinPeptides[m]->Peptide->Length; iH++)
					{
						fExchblH = fExchblH + nHAA[ProteinPeptides[m]->Peptide[iH]];
					}

					fprintf(fp, "%d, %d, %10.5f,   ",
						(int)fExchblH, ProteinPeptides[m]->nCharge, 
						ProteinPeptides[m]->SeqMass);

					//print theoretical isotopes											
					for(q=0; q < 6; q++)
					{
						fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
					}

					fprintf(fp, ", , , , , , , , , , , , ,");
				}
				else if(allExperiments->ExperimentsList->Count -1 == i)
				{
					fprintf(fp, ", , , , , , , , , , , \n");
				}
				else
				{
					fprintf(fp, ", , , , , , , , , , , , , ");
				}
			}

		} // for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	} //for(m = 0; m < ProteinPeptides->Count; m++)

	fclose(fp);

	delete(TheoIsotopePeaks);  delete(fFourierIsotopes);

	WriteNNLSIsotopeDistribution(allExperiments, sProtein, ProteinPeptides, args);

	//printf("FROM WTIE2\n");

	return;

}



//Writes the isotope distribution of a protein into a csv file
void WriteNNLSIsotopeDistribution(ExperimentCollection ^allExperiments, String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args)
{
	const int m1=6;
	const int n1=6;
	int count, m, ret, row_counter, q, nCharge, FoundPeptide, nret;
	double *wp;
	double *zzp;
	int *indexp;
	FILE *fp;//, *fptime;
	String  ^sTemp;
	int j = 0;
	int i =  0;
	int k = 0;
	int l = 0;
	double sumx = 0;

	int q1 = 0;
	int i1 = 0;

	array<double,2> ^a = gcnew array <double,2> (m1,n1);
	double *a_array = (double *)malloc(n1*m1*sizeof(double));
	double *b = (double*)calloc(m1,sizeof(double));
	double *x = (double*)calloc(n1,sizeof(double));

	double temp;
	double temp2;
	double temp3;
	/*for (count=0; count<m1; count++)
	{
	a[count] = (double*)malloc(n1*sizeof(double));
	}*/
	wp = NULL;
	zzp = NULL;
	indexp= NULL; 

	char szTempFile[2046], szIsotopeFile[2046], c;

	char szTempPeptide[2046];

	double dSeqMass, dSpecMass, dtemp, dError;

	double dTheorIsotopes[10];

	double b0[10];

	PeptideHolder ^aPeptide;

	double mass_sum, intensity_sum, frac_intensity, initial_label, intensity_sum0;


	int ntimepoints = allExperiments->ExperimentsList->Count;

	float *MPE_vector = (float*)calloc(ntimepoints,sizeof(float));


	float rkd, rks, fx1;

	std::string templine, tUnits;

	if(NULL == MPE_vector)
	{
		printf("Cannot Allocate Memory for Response\n");

		exit (1);
	}


	TheoreticalIsotopeCalculator ^TheoIsotopePeaks = gcnew TheoreticalIsotopeCalculator();

	array <float> ^fFourierIsotopes = gcnew array <float> (6);

	i =0; m = 0; szTempFile[0] = '\0'; szIsotopeFile[0] ='\0'; 

	String ^filename4 = args[3];		

	if(filename4 != "")
	{
		for(m=0; m<filename4->Length; m++)
		{
			szTempFile[m] = filename4[m];
			szIsotopeFile[i++] = filename4[m];
		}

		szTempFile[m++] = '\\';
		szIsotopeFile[i++] ='\\';
	}

	for(l=0; l < sProtein->Length; l++)
	{
		c = sProtein[l];

		if(c != '|' && c != ',' && c != '\n' && c != '-' &&
			c != '\t' && c != ' ')
		{
			szTempFile[m] = sProtein[l];
			szIsotopeFile[i] = sProtein[l];
			m++;
			i++;
		}
	}

	szTempFile[m] = '\0';
	szIsotopeFile[i] = '\0';

	strcat(szTempFile, ".Quant.csv");
	strcat(szIsotopeFile, ".Kinetics.csv");

	fp = fopen(szIsotopeFile,"w");
	if(fp == NULL)
	{

		printf("Can not open isotope distribution file.\n");
		exit(-1);
	}


	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{
		sTemp = allExperiments->ExperimentsList[i]->sExperimentFile->Substring(0,
			allExperiments->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) ;

		if(0 == i)
		{
			//fprintf(fp, ", , , , , , , , ,%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, ", , , , , , , , ,%s, , , , , , , , , , , ,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			//fprintf(fp, "%s, , , , , , , , , , \n", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , , \n", sTemp);
		}
		else
		{
			//fprintf(fp, "%s, , , , , , , , , , ,", sTemp);

			fprintf(fp, "%s, , , , , , , , , , , ,", sTemp);
		}
	}

	//
	// Print file headers
	fprintf(fp,"%s,","Peptide,Charge,SeqMass,M0,M1,M2,M3,M4,M5");

	for(i=0; i<allExperiments->ExperimentsList->Count; i++)
	{
		if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp,"%s,","I0,I1,I2,I3,I4,I5,X0,X1,X2,X3,X4,X5,Degradation Rate (BFGS),Synthesis Rate (BFGS)");
		}
		else
		{
			fprintf(fp,"%s,","I0,I1,I2,I3,I4,I5,X0,X1,X2,X3,X4,X5");
		}
	}
	fprintf(fp,"%s","\n");



	//**********************************************************************************************
	//print results for every peptide in the protein

	for(m = 0; m < ProteinPeptides->Count; m++)
	{


		FoundPeptide = 0;
		for(i=0;i<ntimepoints;i++)
		{
			MPE_vector[i] = 0.00;
		}
		//

		if(ProteinPeptides[m]->ModLocations->Count == 0)
		{
			fprintf(fp, "%s,", ProteinPeptides[m]->Peptide);
		}
		else
		{
			for(int ii=0; ii < ProteinPeptides[m]->Peptide->Length; ii++)
			{
				bool bModified = false, bCterm = false;   //for O18 Mascot writes the modification outside of peptide length

				for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
				{
					if(ProteinPeptides[m]->ModLocations[q] == ii + 1)   //the mod position is counted starting from 1.
					{
						bModified = true;

						break;
					}
					else if (ProteinPeptides[m]->ModLocations[q] == ProteinPeptides[m]->Peptide->Length + 1 )
					{
						bCterm = true;

						break;
					}
				}

				if(bModified)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else if(bCterm == true && ii == ProteinPeptides[m]->Peptide->Length - 1)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else
				{
					fprintf(fp, "%c", ProteinPeptides[m]->Peptide[ii]);
				}
			}
			fprintf(fp, ",");
		}

		//compute the theoretical isotope

		Isotopes ^isotope = gcnew Isotopes();

		szTempPeptide[0] = '\0';

		for(q=0; q < ProteinPeptides[m]->Peptide->Length; q++)
		{
			szTempPeptide[q] = ProteinPeptides[m]->Peptide[q];
		}

		szTempPeptide[q] = '\0';

		//printf("Seqence = %s\n", szTempPeptide);

		printf("NNLS:First Six Theoretical Isotopes for %s\n", ProteinPeptides[m]->Peptide);

		TheoIsotopePeaks ->StartToComputeSequenceIsotopes(ProteinPeptides[m]->Peptide, fFourierIsotopes);

		for(int ii=0; ii < 10 && ii < fFourierIsotopes->Length; ii++)
		{
			dTheorIsotopes[ii] = 100.*fFourierIsotopes[ii]; 
		}

		//isotope->ComputeIsotopePeaks(szTempPeptide, dTheorIsotopes);
		//delete isotope;

		aPeptide = gcnew PeptideHolder;
		for(i1=0; i1 < allExperiments->ExperimentsList->Count; i1++) //over experiments
		{
			bool bFoundPeptide = false;

			for(j = 0; j < allExperiments->ExperimentsList[i1]->ProteinsList->Count; j++)    //over proteins an experiment
			{
				if(allExperiments->ExperimentsList[i1]->ProteinsList[j]->bProteinPassed &&
					allExperiments->ExperimentsList[i1]->ProteinsList[j]->accession->Equals(sProtein) )
				{
					for(l = 0; l < allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides->Count; l++)    //over peptides of a protein
					{ 
						if(allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->bPeptidePassed &&
							ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->Peptide) &&
							ProteinPeptides[m]->nCharge == allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->nCharge &&
							ProteinPeptides[m]->SeqMass == allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->SeqMass &&
							ProteinPeptides[m]->ModLocations->Count == allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->ModLocations->Count)
						{
							bool bThisPeptide = true;
							NNLS_Class ^call_NNLS = gcnew NNLS_Class();

							for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
							{
								if(ProteinPeptides[m]->ModLocations[q] != 
									allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->ModLocations[q])
								{ 
									bThisPeptide = false;

									break;
								}
							}

							if(bThisPeptide)
							{
								bFoundPeptide = true;
								FoundPeptide++;
								nCharge = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->nCharge;

								dSeqMass = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->SeqMass;

								dSpecMass = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->SpecMass;

								dSeqMass = (dSeqMass - mProton)*(double)nCharge + mProton;

								dtemp = (dSpecMass - mProton)*(double)nCharge + mProton;

								dError = (dSeqMass - dtemp)/dSeqMass*pow(10., 6);		


								if(0 == i1)
								{
									fprintf(fp, "%d, %10.5f,",
										allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->nCharge, 
										allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->SeqMass);

									intensity_sum = 0;
									//print theoretical isotopes											
									for(q=0; q < 6; q++)
									{
										fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
										intensity_sum += allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									//a = gcnew array <double, 2> (2, 6);
									//b = gcnew array <double> (6);

									for(q=0; q < 6; q++)
									{
										a[q,0] = dTheorIsotopes[q]/100; 
										b[q] = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum;
										b0[q] = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum;
									}

									//Assemble a
									a[0,1] = 0;
									a[1,1] = a[0,0];
									a[2,1] = a[1,0];
									a[3,1] = a[2,0];	
									a[4,1] = a[3,0];
									a[5,1] = a[4,0];

									a[0,2] = 0;
									a[1,2] = 0;
									a[2,2] = a[0,0];
									a[3,2] = a[1,0];
									a[4,2] = a[2,0];
									a[5,2] = a[3,0];

									a[0,3] = 0;
									a[1,3] = 0;
									a[2,3] = 0;
									a[3,3] = a[0,0];
									a[4,3] = a[1,0];
									a[5,3] = a[2,0];

									a[0,4] = 0;
									a[1,4] = 0;
									a[2,4] = 0;
									a[3,4] = 0;
									a[4,4] = a[0,0];
									a[5,4] = a[1,0];

									a[0,5] = 0;
									a[1,5] = 0;
									a[2,5] = 0;
									a[3,5] = 0;
									a[4,5] = 0;
									a[5,5] = a[0,0];					

									for(q=0; q < m1; q++)
									{
										for(q1=0; q1 < n1; q1++)
										{
											a_array[q+q1*m1] = a[q,q1];
										}
									}

									ret = call_NNLS->NonNegativeLeastSquares_2(a_array, m1, n1, b, x, NULL);
									sumx = 0;
									for(q=0; q<n1; q++)
									{
										sumx +=x[q];
										fprintf(fp,"%e,",allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									//if(sumx!=0 && !Double::IsInfinity(x[0])) 
									if(sumx!=0)
									{
										for(q=0; q<n1; q++)
										{

											temp3 = x[q]/sumx;
											fprintf(fp,"%10.5f,",temp3);
										}										
									}
									else
									{
										for(q=0; q<n1; q++)
										{
											temp3 = x[q]; fprintf(fp,"%10.5f,",temp3);
										}
									}

									if(ret != 0)
									{
										printf("NNLS does not calculate the isotope distribution correctly.\n");
									}

									frac_intensity = 0;
									if(sumx!=0)
									{
										for(q = 0; q < n1; q++)
										{																				
											frac_intensity += q*x[q]/sumx;
										}
									}

									MPE_vector[i1] = frac_intensity;

									if(allExperiments->ExperimentsList->Count -1 == i1) // in the case when there is only one experiment
									{
										fprintf(fp, "\n");
									}
								}
								else if(allExperiments->ExperimentsList->Count -1 == i1)
								{

									intensity_sum = 0;
									for(q = 0; q < 5; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
										intensity_sum += allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									fprintf(fp, "%e,", allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5]);
									intensity_sum += allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 5];



									for(q=0; q < 6; q++)
									{
										a[q,0] = b0[q];
									}
									for(q=0; q < 6; q++)
									{
										b[q] = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum;
									}

									//Assemble a
									a[0,1] = 0;
									a[1,1] = a[0,0];
									a[2,1] = a[1,0];
									a[3,1] = a[2,0];	
									a[4,1] = a[3,0];
									a[5,1] = a[4,0];

									a[0,2] = 0;
									a[1,2] = 0;
									a[2,2] = a[0,0];
									a[3,2] = a[1,0];
									a[4,2] = a[2,0];
									a[5,2] = a[3,0];

									a[0,3] = 0;
									a[1,3] = 0;
									a[2,3] = 0;
									a[3,3] = a[0,0];
									a[4,3] = a[1,0];
									a[5,3] = a[2,0];

									a[0,4] = 0;
									a[1,4] = 0;
									a[2,4] = 0;
									a[3,4] = 0;
									a[4,4] = a[0,0];
									a[5,4] = a[1,0];

									a[0,5] = 0;
									a[1,5] = 0;
									a[2,5] = 0;
									a[3,5] = 0;
									a[4,5] = 0;
									a[5,5] = a[0,0];					

									for(q=0; q<m1; q++)
									{
										for(q1=0; q1<n1; q1++)
										{
											a_array[q+q1*m1] = a[q,q1];
										}
									}

									ret = call_NNLS->NonNegativeLeastSquares_2(a_array, m1, n1, b, x, NULL);
									sumx = 0;
									for(q=0; q<n1; q++)
									{
										sumx +=x[q];
									}
									//if(sumx!=0 && !Double::IsInfinity(x[0])) 
									if(sumx!=0)
									{
										for(q=0; q<n1; q++)
										{
											temp3 = x[q]/sumx;
											fprintf(fp,"%10.5f,",temp3);
										}										
									}
									else
									{
										for(q=0; q<n1; q++)
										{
											temp3 = x[q]; fprintf(fp,"%10.5f,",temp3);
										}
									}
									//fprintf(fp, "\n");

									if(ret != 0)
									{
										printf("NNLS does not calculate the isotope distribution correctly.\n");
									}
									//fprintf(fp,"%s","\n");


									frac_intensity = 0;
									if(sumx!=0)
									{
										for(q = 0; q < n1; q++)
										{																				
											frac_intensity += q*x[q]/sumx;
										}

									}
									MPE_vector[i1] = frac_intensity;

								}
								else
								{
									intensity_sum = 0;
									for(q = 0; q < 6; q++)
									{
										intensity_sum += allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q];
									}

									for(q=0; q < 6; q++)
									{
										a[q,0] = b0[q];
									}

									for(q=0; q < 6; q++)
									{
										b[q] = allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]/intensity_sum;
									}
									//Assemble a
									a[0,1] = 0;
									a[1,1] = a[0,0];
									a[2,1] = a[1,0];
									a[3,1] = a[2,0];	
									a[4,1] = a[3,0];
									a[5,1] = a[4,0];

									a[0,2] = 0;
									a[1,2] = 0;
									a[2,2] = a[0,0];
									a[3,2] = a[1,0];
									a[4,2] = a[2,0];
									a[5,2] = a[3,0];

									a[0,3] = 0;
									a[1,3] = 0;
									a[2,3] = 0;
									a[3,3] = a[0,0];
									a[4,3] = a[1,0];
									a[5,3] = a[2,0];

									a[0,4] = 0;
									a[1,4] = 0;
									a[2,4] = 0;
									a[3,4] = 0;
									a[4,4] = a[0,0];
									a[5,4] = a[1,0];

									a[0,5] = 0;
									a[1,5] = 0;
									a[2,5] = 0;
									a[3,5] = 0;
									a[4,5] = 0;
									a[5,5] = a[0,0];

									for(q=0; q < m1; q++)
									{
										for(q1=0; q1 < n1; q1++)
										{
											a_array[q+q1*m1] = a[q,q1];
										}
									}
									ret = call_NNLS->NonNegativeLeastSquares_2(a_array, m1, n1, b, x, NULL);
									sumx = 0;
									for(q=0; q<n1; q++)
									{
										sumx +=x[q];
										fprintf(fp,"%e,",allExperiments->ExperimentsList[i1]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}
									//if(sumx!=0 && !Double::IsInfinity(x[0])) 
									if(sumx!=0)
									{
										for(q=0; q<n1; q++)
										{

											temp3 = x[q]/sumx;
											fprintf(fp,"%10.5f,",temp3);
										}										
									}
									else
									{
										for(q=0; q<n1; q++)
										{
											temp3 = x[q]; fprintf(fp,"%10.5f,",temp3);
										}
									}

									if(ret != 0)
									{
										printf("NNLS does not calculate the isotope distribution correctly.\n");
									}

									frac_intensity = 0;
									if(sumx!=0)
									{
										for(q = 0; q < n1; q++)
										{																				
											frac_intensity += q*x[q]/sumx;
										}

									}
									MPE_vector[i1] = frac_intensity;
								}
							}
							//fprintf(fp,"%s","\n");
							delete(call_NNLS);

						}
					}

				}
			}

			//this peptide was not observed in all experiments;
			if(false == bFoundPeptide)
			{
				if(0 == i1)
				{
					fprintf(fp, "%d, %10.5f,   ",
						ProteinPeptides[m]->nCharge, 
						ProteinPeptides[m]->SeqMass);

					//print theoretical isotopes											
					for(q=0; q < 6; q++)
					{
						fprintf(fp, "%5.3f,",  dTheorIsotopes[q]);
					}

					fprintf(fp, ", , , , , , , , , , , ,");
				}
				else
				{
					fprintf(fp, ", , , , , , , , , , , ,");
				}
			}
		}//for i1

		if(FoundPeptide>1)
		{
			//Passing normalized MPE_vector as input to the BFGS dll
			for(q=0;q<ntimepoints;q++){
				MPE_vector[q] = MPE_vector[q]/MPE_vector[ntimepoints-1];
			}



			//nret = lbfgs->Optimize(MPE_vector,&rkd,&rks,&fx1);
		}
		else if(FoundPeptide<=1)
		{
			fprintf(fp, "\n");
		}

	}//for m


	//**********************************************************************************************
	//write isotope distribution into the csv file

	//Deletion

	if(b)
		delete[] b;
	if(x)
		delete[] x;
	if(a_array)
		delete[] a_array;

	if(MPE_vector)
	{
		delete[] MPE_vector;
	}

	//free(fTime);

	fclose(fp);

	delete(TheoIsotopePeaks); delete(fFourierIsotopes); 


	return;


}


void O18_Quantification(array<System::String ^> ^inputFiles)
{
	ExperimentCollection ^allProteinData = gcnew ExperimentCollection;

	Read_Params(inputFiles);

	printf("The current params setting: %10.1f %10.1f\n", ProteinScoreCut, dIonScoreCut);

	ReadmzID(inputFiles, allProteinData);

	ConsistentProteins(allProteinData, 1, ProteinScoreCut, nProteinConsistency);

	ConsistentPeptides(allProteinData, dIonScoreCut, dExpectCut, nPeptideConsistency);

	array <String ^> ^szML;

	ReadAbundances(szML, allProteinData);

	WriteQuantResults(allProteinData,inputFiles);

	delete (allProteinData);

	return;
}

void ReadFileNames(array<System::String ^> ^inputs, List <String ^> ^asmzML, List <String ^> ^asmzID,array<System::String ^> ^args)
{
	int i;

	char szLine[10046], szFileName[10046], szTemp1[2046], szTemp2[2046];

	String ^sTemp;

	FILE *fp;

	for(i=0; i < inputs[0]->Length; i++)
	{
		szFileName[i] = inputs[0][i];
	}

	szFileName[i] = '\0';

	fp = fopen(szFileName, "r");

	if(NULL == fp)
	{
		printf("Cannot read %s\n", szFileName);

		exit (1);
	}


	while(NULL != fgets(szLine, sizeof(szLine), fp) )
	{
		if(2 != sscanf(szLine, "%s %s", szTemp1, szTemp2) )
		{
			printf("Unexpected file format: %s in %s\n", szLine, szFileName);

			exit (1);
		}

		sTemp = gcnew String(szTemp1);

		asmzML->Add(sTemp);

		sTemp = gcnew String(szTemp2);

		asmzID->Add(sTemp);
	}

	fclose(fp);

	return;
}

/*
*
* A program to read mzID files,
* to combine them in a single ExperimentCollection collection List
* mzML filenames are used to add to the collection
* mzML files are actually read somewhere else
* at the end there is a ready to use file allExperiments which has all necessary
* protein and peptide identification information
*
*/
int ReadmzID(array<System::String ^> ^args, ExperimentCollection ^allExperiments)
{
	array<String ^> ^smzId, ^smzML;

	charArray szFile;

	PeptideHolder ^aPeptide = gcnew PeptideHolder();

	ProteinCollection ^currentExperiment = gcnew ProteinCollection();

	//ExperimentCollection ^allExperiments = gcnew ExperimentCollection();

	allExperiments->ExperimentsList = gcnew List<ProteinCollection^>;

	ProteinSet ^aProtein;

	char szTemp[2024], szLine[2024], szTemp1[2024], szTemp2[2024];

	bool bSpectra = false;     //this varialbe controls reading and processing of mzML files if provided

	int i, k, l, j, m;

	if(args[0]->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase) == -1 && args[0]->IndexOf(".txt", StringComparison::OrdinalIgnoreCase) == -1)
	{
		printf("Need an mzid or txt file to run the program\n");

		exit (1);
	}

	if(args[0]->IndexOf("*.mzid", StringComparison::OrdinalIgnoreCase) != -1)              // caller_mzIdentML.exe *.mzid
	{
		smzId = System::IO::Directory::GetFiles("\.");
	}
	else if(args->Length == 1 && args[0]->IndexOf(".txt") != -1)
	{
		bSpectra = true;
		//read file names in pairs from the XXXX.txt file

		szTemp[0] = '\0';

		for(i=0; i < args[0]->Length; i++)
		{
			szTemp[i] = args[0][i];
		}

		szTemp[i] = '\0';

		FILE *fp = fopen(szTemp, "r");

		if(fp == NULL)
		{
			printf("Cannot read %s\n", szTemp); 

			exit (1);
		}

		k = 0;

		while(NULL != fgets(szLine, sizeof(szLine), fp) )
		{
			k++;      //determine number of files/experiments
		}

		smzId = gcnew array<String ^> (k);

		smzML = gcnew array<String ^> (k);

		rewind(fp); k = 0;

		while(NULL != fgets(szLine, sizeof(szLine), fp) )
		{
			sscanf(szLine, "%s %s", szTemp1, szTemp2);

			szFile.szChar[0] = '\0';

			for(m = 0; m < strlen(szTemp1); m++)
			{
				szFile.szChar[m] = szTemp1[m];
			}

			szFile.szChar[m] = '\0';

			smzML[k] = gcnew String(szTemp1);

			smzId[k] = gcnew String(szTemp2);

			k++;      
		}

		fclose(fp);
	}
	else
	{
		smzId = gcnew array<String ^> (1);

		smzId[0] = args[0];
	}

	k = 0;

	for(l=0; l < smzId->Length; l++)
	{
		if(smzId[l]->IndexOf(".mzid", StringComparison::OrdinalIgnoreCase) != -1)
		{
			k++;

			// process/read the mzid file
			mzIdentML ^mzI = gcnew mzIdentML(smzId[l]);

			mzI->ReadProteins_mzIdent_2();

			//now create a storage for this experiment
			currentExperiment = gcnew ProteinCollection();

			currentExperiment->ProteinsList = gcnew List<ProteinSet^>;

			//count goes over every protein entry read in from mzid
			for(i=0; i < mzI->ProteinResultList->Count; i++)
			{
				//copy protein data
				aProtein = gcnew ProteinSet();

				aProtein->Peptides = gcnew List <PeptideHolder ^>;


				aProtein->accession          = mzI->ProteinResultList[i]->accession;
				aProtein->description        = mzI->ProteinResultList[i]->description;
				aProtein->ProteinScore       = mzI->ProteinResultList[i]->ProteinScore;
				aProtein->SeqCoverage        = mzI->ProteinResultList[i]->SeqCoverage;
				aProtein->nSpectralCount     = mzI->ProteinResultList[i]->nSpectralCount;
				aProtein->nDistinctSequences = mzI->ProteinResultList[i]->nDistinctSequences;
				aProtein->nSeqLength         = mzI->ProteinResultList[i]->nSeqLength;
				aProtein->ProteinMass        = mzI->ProteinResultList[i]->ProteinMass;

				for(int i0 = 0; i0 < mzI->ProteinResultList[i]->PeptideList->Count; i0++)
				{
					aPeptide = gcnew PeptideHolder();

					aPeptide->dIonscore = mzI->ProteinResultList[i]->PeptideList[i0]->dIonscore;
					aPeptide->bUniquePeptide = mzI->ProteinResultList[i]->PeptideList[i0]->bUniquePeptide;
					aPeptide->dRetTime  = mzI->ProteinResultList[i]->PeptideList[i0]->dRetTime;
					aPeptide->Peptide   = mzI->ProteinResultList[i]->PeptideList[i0]->Peptide;
					aPeptide->nScan     = mzI->ProteinResultList[i]->PeptideList[i0]->nScan;
					aPeptide->nCharge   = mzI->ProteinResultList[i]->PeptideList[i0]->nCharge;
					aPeptide->SeqMass  = mzI->ProteinResultList[i]->PeptideList[i0]->SeqMass;
					aPeptide->SpecMass = mzI->ProteinResultList[i]->PeptideList[i0]->SpecMass;
					aPeptide->dExpect    = mzI->ProteinResultList[i]->PeptideList[i0]->dExpect;
					aPeptide->dIsotopes = gcnew array <double, 2>(2, 6);
					aPeptide->w = gcnew array <double, 1>(12);
					

					aPeptide->dModMasses = gcnew List <double>;

					aPeptide->ModLocations = gcnew List <int>;

					aPeptide->dModMasses = mzI->ProteinResultList[i]->PeptideList[i0]->dModMasses;

					aPeptide->ModLocations = mzI->ProteinResultList[i]->PeptideList[i0]->ModLocations; 

					/*printf("Peptide %s Expct = %10.5f Score = %10.5f\n", aPeptide->Peptide,  aPeptide->dExpect,
					aPeptide->dIonscore);*/

					aProtein->Peptides->Add(aPeptide);

					delete(aPeptide);
				}

				currentExperiment->ProteinsList->Add(aProtein);

				delete (aProtein);
			}

			j = smzId[l]->LastIndexOf("\\") + 1;

			if(j >= 1)
			{
				currentExperiment->sExperimentFile = smzId[l]->Substring(j);
			}
			else
			{
				currentExperiment->sExperimentFile = smzId[l];
			}

			allExperiments->ExperimentsList->Add(currentExperiment);

			j = smzML[l]->LastIndexOf("\\") + 1;

			if(j >= 1 && false)
			{
				currentExperiment->smzML = smzML[l]->Substring(j);
			}
			else
			{
				currentExperiment->smzML = smzML[l];
			}

			delete(mzI);

			delete(currentExperiment->ProteinsList);
		} //if(sFiles[l]->IndexOf(".mzid") != -1)
	}

	printf("\nNumber of mzid files processed: %d\n", allExperiments->ExperimentsList->Count);

	if(allExperiments->ExperimentsList->Count == 0)
	{
		printf("Did not Read any Data to Process\n");

		exit (1);
	}


	//read the intensity information from the mzML file
	if(bSpectra)
	{
		//ReadAbundances(smzML, allExperiments);
	}

	delete aPeptide, aProtein, currentExperiment;

	return 0; 
}


/*
*
* reads parameters from quant.state file
*
*/
void Read_Params(array<System::String ^> ^args)
{
	FILE *fp;

	char szLine[10046], szTemp1[256], szTemp2[256], szTemp3[256], szTemp4[256];

	char *szFile = "Quant.state";

	char *word1, *word2, *word3;

	fp = fopen(szFile, "r");

	if(NULL == fp)
	{
		printf("No quant.state file, will use default settings for filtering\n");

		printf("Protein Score cut-off: %10.5f\n", ProteinScoreCut);

		printf("Peptide Score cut-off: %10.5f\n", dIonScoreCut);

		printf("Protein Consistency cut-off = (Number.Of.Experiments - 1)\n");

		return;
	}

	while(NULL != fgets(szLine, sizeof(szLine), fp) )
	{


		szLine[strcspn(szLine, "\n")] = '\0';	// remove newline character(s)

		//word1 = strtok(szLine, " =");			// split characters are ' ' and '='

		if(NULL == szLine || '#' == szLine[0])
		{
			continue;
		}


		if(strstr(szLine, "mass_accuracy") != NULL)
		{
			word1 = strtok(szLine, " =");			// split characters are ' ' and '='

			word2 = strtok(NULL, " =");

			word3 = strtok(NULL, " =");

			//printf("Word %s %s %s\n", word1, word2, word3);

			dMassAccuracy = atof(word2);

			if(0 == strcmp(word3, "ppm") || 0 == strcmp(word3, "PPM") )
			{
				mass_accuracy_unit = 0;

				printf("Mass Accuracy in ppm\n");

			}
			else if(0 == strcmp(word3, "Da") || 0 == strcmp(word3, "DA"))
			{
				mass_accuracy_unit = 1;

				//printf("Mass Accuracy in Da\n");
			}
			else
			{
				printf("The program cannot recognize the mass accuracy units: Please, use ppm or Da for mass accuracy units\n");

				printf("Exiting now\n");
				exit (1);
			}

		}
		else if(strstr(szLine, "protein_score") != NULL)
		{
			if(4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read Protein Score: Error in quant.state file\n");

				exit (1);
			}

			ProteinScoreCut = atof(szTemp3);
		}
		else if(strstr(szLine, "peptide_score") != NULL)
		{
			if(4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read Protein Score: Error in quant.state file\n");

				exit (1);
			}

			dIonScoreCut = atof(szTemp3);
		} 
		else if(strstr(szLine, "peptide_expectation") != NULL)
		{
			if(4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read Peptide Expectation: Error in quant.state file\n");

				exit (1);
			}

			dExpectCut = atof(szTemp3);
		}
		else if(strstr(szLine, "elutiontimewindow") != NULL)
		{
			if(4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read Protein Score: Error in quant.state file\n");

				exit (1);
			}

			ElutionTimeWindow = atof(szTemp3);
			//ElutionTimeWindow = float::Parse(args[2]); //New addition; add from the command line
			printf("Elution time window : %f\n",ElutionTimeWindow);
		}
		else if(strstr(szLine, "protein_consistency") != NULL)
		{
			if(4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read Protein Score: Error in quant.state file\n");

				exit (1);
			}

			nProteinConsistency = atoi(szTemp3);
		}
		else if(strstr(szLine, "peptide_consistency") != NULL)
		{
			if(4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read Protein Score: Error in quant.state file\n");

				exit (1);
			}

			nPeptideConsistency = atoi(szTemp3);
		}
		else if (strstr(szLine, "NParam_RateConst_Fit") != NULL)
		{
			if (4 != sscanf(szLine, "%s %s %s %s", szTemp1, szTemp2, szTemp3, szTemp4) )
			{
				printf("Cannot Read NParam_RateConst_Fit: Error in quant.state file\n");

				exit(1);
			}

			if(atoi(szTemp3) != 1 && atoi(szTemp3)!=2 && atoi(szTemp3)!=3 && atoi(szTemp3) != -1)
			{
				printf("Error in quant.state file. Please set NParam_RateConst_Fit to 1 or 2 or 3\n");
				exit(1);
			}
			Nparameter = atoi(szTemp3);
		}
	}

	fclose(fp);

	printf("Params: MassAccuracy = %10.3f PeptideConsistency = %d, ProteinConsistency = %d, PeptideScore = %5.2f PeptideExpectation = %e\n#Params in rate constant: %d\n", 
		dMassAccuracy, nPeptideConsistency, nProteinConsistency, dIonScoreCut, dExpectCut, Nparameter);

	//initialize the number of exchangeable H atoms for AA

	for(int i = 0; i < 256; i++)
	{
		nHAA[i] = 0.0;
	}

	nHAA[(int)'A'] = nHAA[(int)'a']= 4.0f;  nHAA[(int)'C'] = nHAA[(int)'c'] = 1.62f; nHAA[(int)'D'] = nHAA[(int)'d']= 1.89f; nHAA[(int)'E'] = nHAA[(int)'e']= 3.95f; 

	nHAA[(int)'F'] = nHAA[(int)'f'] = 0.32f; nHAA[(int)'G'] = nHAA[(int)'g'] = 2.06f; nHAA[(int)'H'] = nHAA[(int)'h'] = 2.88f; nHAA[(int)'I'] = nHAA[(int)'i'] = 1.0f; 

	nHAA[(int)'L'] = nHAA[(int)'l'] = 0.6f; nHAA[(int)'K'] = nHAA[(int)'k'] = 0.54f; nHAA[(int)'M'] = nHAA[(int)'m'] = 1.12f; nHAA[(int)'N'] = nHAA[(int)'n'] = 1.89f; 

	nHAA[(int)'P'] = nHAA[(int)'p'] = 2.59f; nHAA[(int)'Q'] = nHAA[(int)'q'] = 3.95f; nHAA[(int)'R'] = nHAA[(int)'r'] = 3.43f; nHAA[(int)'S'] = nHAA[(int)'s'] = 2.61f; 

	nHAA[(int)'T'] = nHAA[(int)'t'] = 0.2f; nHAA[(int)'V'] = nHAA[(int)'v'] = 0.56f; nHAA[(int)'W'] = nHAA[(int)'w'] = 0.08f; nHAA[(int)'Y'] = nHAA[(int)'y'] = 0.42f;


	return;
}

/*
*  will write results from O18 labeling experiments
*/

void WriteAProtein_O18(ExperimentCollection ^allExperiments, String ^sProtein, List <PeptideHolder ^> ^ProteinPeptides,array<System::String ^> ^args)
{
	int i, j, l, q, m, nCharge;

	bool bNewPeptide = false, bThisPeptide= false;

	double dSeqMass, dSpecMass, dtemp, dError, dRatio;

	double dTheorIsotopes[10];

	char szTempPeptide[2046];

	String ^sTemp;

	PeptideHolder ^aPeptide;

	FILE *fp;

	String ^filename4 = args[3];

	char szTempFile[2046], c;

	m = 0; szTempFile[0] = '\0';


	if(filename4 != "")
	{
		for(m=0; m<filename4->Length; m++)
			szTempFile[m] = filename4[m];

		szTempFile[m++] = '\\';
	}


	for(l=0; l < sProtein->Length; l++)
	{
		c = sProtein[l];

		if(c != '|' && c != ',' && c != '\n' && c != '-' &&
			c != '\t' && c != ' ')
		{
			szTempFile[m] = sProtein[l];

			m++;
		}
	}

	szTempFile[m] = '\0';

	strcat(szTempFile, ".Quant.csv");

	fp = fopen(szTempFile, "w");

	if(NULL == fp)
	{
		printf("Cannot write to %s\n", szTempFile);

		exit (1);
	}

	//print the first line of the csv file
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{
		sTemp = allExperiments->ExperimentsList[i]->sExperimentFile->Substring(0,
			allExperiments->ExperimentsList[i]->sExperimentFile->LastIndexOf(".")) ;

		if(0 == i)
		{
			fprintf(fp, ",%s, , , , , , , , , , , , , ,", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, "%s, , , , , , , , , , ,\n", sTemp);
		}
		else
		{
			fprintf(fp, "%s, , , , , , , , , , , ,", sTemp);
		}
	}

	sTemp = gcnew String("SpecMass, IonScore, Expectn, Error(ppm), Scan, M0, M1, M2, M3, M4, M5, Ratio");

	//print the file header
	for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	{
		if(allExperiments->ExperimentsList->Count - 1 == i && 0 == i) // in the case when there is a single experiment
		{
			fprintf(fp, "\nPeptide, Charge, SeqMass, %s\n", sTemp);
		}
		else if(allExperiments->ExperimentsList->Count - 1 == i)
		{
			fprintf(fp, "%s\n", sTemp);
		}
		else if(0 == i)
		{
			fprintf(fp, "Peptide, Charge, SeqMass, %s,", sTemp);
		}
		else
		{
			fprintf(fp, "%s,", sTemp);
		}
	}

	//print results for every peptide in the protein

	for(m = 0; m < ProteinPeptides->Count; m++)
	{
		if(ProteinPeptides[m]->ModLocations->Count == 0)
		{
			fprintf(fp, "%s,", ProteinPeptides[m]->Peptide);
		}
		else
		{
			//print peptide sequence, AA residues
			for(int ii=0; ii < ProteinPeptides[m]->Peptide->Length; ii++)
			{
				bool bModified = false, bCterm = false;   //for O18 Mascot writes the modification outside of peptide length

				for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
				{
					if(ProteinPeptides[m]->ModLocations[q] == ii + 1)   //the mod position is counted starting from 1.
					{
						bModified = true;

						break;
					}
					else if (ProteinPeptides[m]->ModLocations[q] == ProteinPeptides[m]->Peptide->Length + 1 )
					{
						bCterm = true;

						break;
					}
				}

				if(bModified)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else if(bCterm == true && ii == ProteinPeptides[m]->Peptide->Length - 1)
				{
					fprintf(fp, "%c", tolower(ProteinPeptides[m]->Peptide[ii]));
				}
				else
				{
					fprintf(fp, "%c", ProteinPeptides[m]->Peptide[ii]);
				}
			}
			fprintf(fp, ",");
		}

		Isotopes ^isotope = gcnew Isotopes();

		szTempPeptide[0] = '\0';

		for(q=0; q < ProteinPeptides[m]->Peptide->Length; q++)
		{
			szTempPeptide[q] = ProteinPeptides[m]->Peptide[q];
		}

		szTempPeptide[q] = '\0';

		//printf("Seqence = %s\n", szTempPeptide);

		printf("First Five Theoretical Isotopes for %s\n", ProteinPeptides[m]->Peptide);

		isotope->ComputeIsotopePeaks(szTempPeptide, dTheorIsotopes);

		delete isotope;

		aPeptide = gcnew PeptideHolder;

		for(i=0; i < allExperiments->ExperimentsList->Count; i++) //over experiments
		{
			bool bFoundPeptide = false;

			for(j = 0; j < allExperiments->ExperimentsList[i]->ProteinsList->Count; j++)
			{
				if(allExperiments->ExperimentsList[i]->ProteinsList[j]->bProteinPassed &&
					allExperiments->ExperimentsList[i]->ProteinsList[j]->accession->Equals(sProtein) )
				{
					for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
					{ 
						if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->bPeptidePassed &&
							ProteinPeptides[m]->Peptide->Equals(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->Peptide) &&
							ProteinPeptides[m]->nCharge == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge &&
							ProteinPeptides[m]->SeqMass == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass &&
							ProteinPeptides[m]->ModLocations->Count == allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations->Count)
						{
							bool bThisPeptide = true;

							for(q=0; q < ProteinPeptides[m]->ModLocations->Count; q++)
							{
								if(ProteinPeptides[m]->ModLocations[q] != 
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->ModLocations[q])
								{ 
									bThisPeptide = false;

									break;
								}
							}

							if(bThisPeptide)
							{
								bFoundPeptide = true;

								nCharge = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge;

								dSeqMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass;

								dSpecMass = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass;

								dSeqMass = (dSeqMass - mProton)*(double)nCharge + mProton;

								dtemp = (dSpecMass - mProton)*(double)nCharge + mProton;

								dError = (dSeqMass - dtemp)/dSeqMass*pow(10., 6);		

								if(allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0] < 1.0)
								{
									allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0] = 1.0;
								}

								//compute the ratio:

								double d0, d2, d4, d2_temp;

								d0 = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 0];

								d2 = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 2];

								d4 = allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, 4];

								d2_temp = d2 - d0*dTheorIsotopes[2]/dTheorIsotopes[0];

								if(d2_temp < 0.0)
								{
									d2_temp = 0.0;
								}

								dRatio = d0 / (d4 - d0*dTheorIsotopes[4]/dTheorIsotopes[0] + d2_temp * ( 1. - dTheorIsotopes[2]/dTheorIsotopes[0]) );

								if(dRatio < 0)
								{
									dRatio = 0.0;
								}


								if(0 == i && allExperiments->ExperimentsList->Count -1 == i) // in the case when there is only one experiment
								{
									fprintf(fp, "%d, %10.5f,  %10.5f, %10.5f, %e, %4.1f, %d,",
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									{
										fprintf(fp, "%e,", 
											allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp, "%5.1f\n", dRatio);
								}
								else if(0 == i)
								{
									fprintf(fp, "%d, %10.5f,  %10.5f, %10.5f, %e, %4.1f, %d,",
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nCharge, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SeqMass,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									{
										fprintf(fp, "%e,", 
											allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp, "%5.1f,", dRatio);
								}
								else if(allExperiments->ExperimentsList->Count -1 == i)
								{
									fprintf(fp, "%10.5f,  %10.5f, %e, %4.1f, %d, ", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp, "%5.1f\n", dRatio);
								}
								else
								{
									fprintf(fp, "%10.5f, %10.5f, %e, %4.1f, %d,", 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->SpecMass, 
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIonscore,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dExpect,
										dError,
										allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->nScan);

									for(q = 0; q < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes->Length/2; q++)
									{
										fprintf(fp, "%e,", allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides[l]->dIsotopes[1, q]);
									}

									fprintf(fp, "%5.3f,", dRatio);
								}
							}

						}
					} //for(l = 0; l < allExperiments->ExperimentsList[i]->ProteinsList[j]->Peptides->Count; l++)
				}
			}

			//this peptide was not observed in all experiments;
			if(false == bFoundPeptide)
			{
				if(0 == i)
				{
					fprintf(fp, "%d, %10.5f, , , , , , , , , , , , ,",
						ProteinPeptides[m]->nCharge, 
						ProteinPeptides[m]->SeqMass);
				}
				else if(allExperiments->ExperimentsList->Count -1 == i)
				{
					fprintf(fp, ", , , , , , , , , , ,\n");
				}
				else
				{
					fprintf(fp, ", , , , , , , , , , , , ");
				}
			}

		} // for(i=0; i < allExperiments->ExperimentsList->Count; i++)
	} //for(m = 0; m < ProteinPeptides->Count; m++)

	fclose(fp);

	return;

}

/*
*  This function calls the RateConstant routine
*  To call RateConstant, the function needs the time course information (time array),
*  and the BodyWater enrichment of Heavy Water. This information is read from the files.txt
*  Also needed is the working directory information.
*  strFile - is the file containing mzml, mzid, time and body water enrichment information
*  DeuEnrichment is a List that will hold the body water enrichment
*
*/
//void  ProteinRates(String ^strFile)
void ProteinRates(array <float> ^fTimeCourse,  array <float> ^fBodyWaterEnrichment, int Nparameter)
{
	FILE *fpFilesTxt;

	List <float> ^Times = gcnew List <float>;

	//List <float> ^DeuEnrichment = gcnew List <float>;

	float ftemp;

	int i;

	//array <float> ^ fTimeCourse, ^fBodyWaterEnrichment;

	/*char szFile[3024];

	for(i = 0; i < strFile->Length; i++)
	{
		szFile[i] = strFile[i];
	}

	szFile[i] = '\0';

	fpFilesTxt = fopen(szFile, "r");*/

	String^ CurrentFolder = Directory::GetCurrentDirectory();

	//printf("Current directory is %s\n", CurrentFolder);

	//if(NULL == fpFilesTxt)
	//{
	//	printf("Cannot find in the %s the \"files.txt\" file\n", CurrentFolder);

	//	exit (1);
	//}

	//ftemp = 0.0;

	//while(fgets(szLine, sizeof(szLine), fpFilesTxt) != NULL)
	//{
	//	if(szLine[0] != '#')
	//	{
	//		if (4 != sscanf(szLine, "%s %s %s %s", &szTemp1, &szTemp2, &szTemp3, &szTemp4) )
	//		{
	//			printf("A wrong format for files.txt\n");
	//		}
	//		else
	//		{
	//			ftemp = atof(szTemp1);

	//			Times->Add(ftemp);

	//			ftemp = atof(szTemp4);

	//			DeuEnrichment->Add(ftemp);
	//		}
	//	}
	//}

	//fclose(fpFilesTxt);

	//nTimePoints = Times->Count;

	//fBodyWaterEnrichment = gcnew array <float> (nTimePoints);

	//fTimeCourse = gcnew array <float> (nTimePoints);

	//for(i = 0; i < nTimePoints; i++)
	//{
	//	fTimeCourse[i] =  Times[i];

	//	fBodyWaterEnrichment[i] = DeuEnrichment[i];
	//}

	/*
	*   call the protein rate constant calculation
	*/

	ProteinRateConstant ^ProtRate;


	ProtRate = gcnew ProteinRateConstant(CurrentFolder, fTimeCourse, fBodyWaterEnrichment,Nparameter);

	ProtRate->ReadFileFolder();

	 
	//delete(fTimeCourse); delete(Times);

	//delete(DeuEnrichment); delete(fBodyWaterEnrichment);

}

/*
*
*    A function to process the input text file - most of the time
*    it is the files.txt
*    A sample line from this file is:
*
*    14.0 Z:\RawFiles\UCLA\Mice\C57\CTRL\Cyto\mzML\C57_CTRL_Heart_Cyto_d14_A03.mzML recent_search\C57_CTRL_Heart_Cyto_d14_A03.mzid   0.05
*
*    the first number is the time point of incorporation experiment, the string is the mzml file, the third string is the mzid file,
*    the fourth number is the enrichment of D2O in the water.
*
*    # - is the escape/comment line
* 
*
*
*/
int Parse_InputTextFile(String ^strFile, array <String ^> ^ mzmlFiles, array <String ^> ^ mzidFiles, array <float> ^fTimePoints, array <float> ^fBodyWaterD2O)
{
	char szFile[2034], szLine[4050], szTemp[2030];

	int i, k, j, nExperiments, len;

	String ^strTemp;

	FILE *fpInpTxt;

	StringComparison ^ComparisonType = gcnew StringComparison ();   

	StringComparer ^Comparea;

	//String ^ mzmlNoCase = gcnew String(".mzml"); // (StringComparer.InvariantCultureIgnoreCase);

	szFile[0] = '\0';

	for(i = 0; i < strFile->Length; i++)
	{
		szFile[i] = strFile[i];
	}

	szFile[i] = '\0';

	fpInpTxt = fopen(szFile, "r");

	if(NULL == fpInpTxt)
	{
		printf("Cannot open the files %s for reading\n", szFile);

		exit (1);
	}

	nExperiments = 0;

	k = 0;

	while(fgets(szLine, sizeof(szLine), fpInpTxt) != NULL)
	{
		//printf("szLine = %s\n", szLine);

		if(strlen(szLine) > 1)       //only the lines that are longer than 1 character
		{
			if(szLine[0] != '#')     //ignore the comment line
			{
				strTemp = gcnew String(szLine);

				if(strTemp->LastIndexOf(".mzid", StringComparison::OrdinalIgnoreCase) > 1 && strTemp->LastIndexOf(".mzml", StringComparison::OrdinalIgnoreCase) > 1)
				{
					//printf("SzLine = %s\n", szLine);

					k = 0;  szTemp[0] = '\0';

					for(j = 0; j < strlen(szLine); j++)
					{
						if(isspace(szLine[j]))
						{
							szTemp[k] = '\0';

							fTimePoints[nExperiments] = atof(szTemp);    //number of time points at the start of the line

							break;     //break when you reach the first blank character
						}
						else
						{
							szTemp[k] = szLine[j];

							k++;
						}
					}

					//find a new non-space character:

					for(k = j; k < strlen(szLine); k++)
					{
						if(!isspace(szLine[k]))
						{
							break;
						}
					}

					len = strTemp->LastIndexOf(".mzml", StringComparison::OrdinalIgnoreCase) - k + 5;

					//printf("%s  %d %f\n", strTemp->Substring(k, len), strTemp->IndexOf(".mzml"), atof(szTemp));

					mzmlFiles[nExperiments] = strTemp->Substring(k, len);
					//nFilePath =  PathFileExists(mzmlFiles[nExperiments]);

					//get the mzid FIles

					len = strTemp->LastIndexOf(".mzml", StringComparison::OrdinalIgnoreCase);

					k = 0;

					for(k = len + 6; k < strlen(szLine); k++)
					{
						if(!isspace(szLine[k]))
						{
							break;
						}
					}

					len = strTemp->LastIndexOf(".mzid", StringComparison::OrdinalIgnoreCase) - k + 5;

					//len = strTemp->LastIndexOf(".mzid") - k + 5;

					mzidFiles[nExperiments] = strTemp->Substring(k, len);

					//parse the body water enrichment

					len = strTemp->LastIndexOf(".mzid", StringComparison::OrdinalIgnoreCase) + 5;

					k = 0; szTemp[0] = '\0';
					
					for(j = len + 1; j < strlen(szLine); j++)
					{
						if(szLine[j] == '#')
						{
							break;
						}
						//else if(!iss(szLine[j]))
						else if(isdigit(szLine[j]) || szLine[j] == '.')
						{
							szTemp[k] = szLine[j];

							k++;
						}
					}

					szTemp[k + 1] = '\0';

					fBodyWaterD2O[nExperiments] = atof(szTemp);

					//printf("nExperiments  = %d  %s %s %f %f\n", nExperiments, mzmlFiles[nExperiments], mzidFiles[nExperiments], 
						//fTimePoints[nExperiments], fBodyWaterD2O[nExperiments]);

					nExperiments++;
				}
			}
		}
	}

	fclose(fpInpTxt);

	//check if the files can be read
	for (i = 0; i < mzmlFiles->Length; i++)
	{
		FILE *fpCheck, *fpmzID;

		char szTempFileName[4058], szTempmzID[5048];

		szTempFileName[0] = '\0';

		szTempmzID[0] = '\0';
		

		for (k = 0; k < mzmlFiles[i]->Length; k++)
		{
			szTempFileName[k] = mzmlFiles[i][k];
		}

		szTempFileName[k] = '\0';

		for (k = 0; k < mzidFiles[i]->Length; k++)
		{
			szTempmzID[k] = mzidFiles[i][k];
		}
		
		szTempmzID[k] = '\0';

		fpCheck = fopen(szTempFileName, "r");

		fpmzID = fopen(szTempmzID, "r");

		if (NULL == fpCheck)
		{
			printf("Cannot Access %s %s\n", mzmlFiles[i], szTempFileName);

			printf("stopping program prematurely...\n");

			exit(1); 
		}
		else if (NULL == fpmzID)
		{
			printf("Cannot Access %s %s\n", mzidFiles[i], szTempmzID);

			printf("Stopping program prematurely, Parse_Input_Text_File ...\n");

			exit(1);
		}

		fclose(fpCheck); fclose(fpmzID);
	}


	if (nExperiments == mzmlFiles->Length)
	{
		printf("Successfully Parsed the %s file - %d experiments\n", strFile, nExperiments);
	}
	else
	{
		printf("There was a problem with parsing %s\n", strFile);
	}

	return nExperiments;
}



/*
*  the function returns the number of experiments
*  in the input file
*
*/
int NumberOfExperiments(String ^strFile)
{
	char szFile[2034], szLine[4050];

	int i, nExperiments;

	String ^strTemp;

	FILE *fpInpTxt;

	szFile[0] = '\0';

	for(i = 0; i < strFile->Length; i++)
	{
		szFile[i] = strFile[i];
	}

	szFile[i] = '\0';

	fpInpTxt = fopen(szFile, "r");

	if(NULL == fpInpTxt)
	{
		printf("Cannot open the files %s for reading\n", szFile);

		exit (1);
	}

	nExperiments = 0;

	//determine the number of experiments:
	while(fgets(szLine, sizeof(szLine), fpInpTxt) != NULL)
	{
		//printf("szLine %s %d\n", szLine, strlen(szLine));

		if(strlen(szLine) > 1)       //only the lines that are longer than 1 character
		{
			if(szLine[0] != '#')     //ignore the comment line
			{
				strTemp = gcnew String(szLine);

				if(strTemp->LastIndexOf(".mzid", StringComparison::OrdinalIgnoreCase) > 1 && strTemp->LastIndexOf(".mzml", StringComparison::OrdinalIgnoreCase) > 1)
				{
					nExperiments++;
				}
			}
		}
	}

	fclose(fpInpTxt);

	return nExperiments;
}