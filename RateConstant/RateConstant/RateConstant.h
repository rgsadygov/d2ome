// RateConstant.h
/*
*
*  A group of methods that will incorporation rate of a peptide
*  and return the decay rate constant, and correlation between the
*  the fit and the experimental data.
*
*  The basic workflow of the algorithm is as follows:
*  1) The class is initiated, and the path to the directory holding
*  .qaunt.csv files and the time points at which the experimental data were collected
*   are passed to it.
*
*   2) Algorithm reads quant.csv files line by line. For each peptide,
*   it reads the sequence information (sets the number of nC, nH, .., nS
*   nExHydrogen) and the isotopes at every time point. 
*
*   The isotope information is passed to ComputeRateFromI0. This method converts the isotope
*   information to necessary I0_Isotopes form - which are the monoisotope at each time point.
*   
*   3) The I0_Isotopes, along with the initial and asymptotic values of I0 are passed to
*   LBFGS algorithm lbfgs->Optimize(TimeCourseI0Isotope, I0_AtAsymptote, (I0_Natural - I0_AtAsymptote), &rkd, &rks, &fx1)
*   to determine the rate constant.
*   
*   4) PearsonCorrelation computes the correlation between the fit (theoretical) and experimental I0_Isotopes.
*
*   5) The results (peptides, rate constants, correlations) are written to "Protein.RateConst.csv" files. 
*
*
*
*  nExHydrogen is the number of exchangeable hydrogen atoms. It is
*  calculated in this algorithm in the function, ElementalComposition(String ^sSequence)
*  
*  nIsotopes - the number of isotopes that are used to
*  compute the relative value of I0. It is set to 5, 
*  but could possibly be made adaptable to a specific 
*  experiment.
*  fTimeScale - scale the time so that there are no time values larger than 10. Otherwise, LBFGS has difficulities to to converge.
*
*  Generic array:
* 
*  array< MyClass^ >^ local = gcnew array< MyClass^ >(ARRAY_SIZE);
*
*   bReplicates - if replicates are present set to true
*
*/

#include "stdio.h"
#include "stdlib.h"
#include <map>
#include <vector>

using namespace System;
using namespace System::IO;
using namespace System::Collections;
using namespace System::Collections::Generic;


namespace RateConstant {



	//this class will hold the isotopes for each 
	// time point in experiment
	public ref class SingleIsotopeCluster
	{
		public:
			array <float> ^ fIsotopeCluster;

			SingleIsotopeCluster()
			{
				fIsotopeCluster = gcnew array <float> (6);
			}
	};

	public ref class RateConstResults
	{
		public:
			float fRateConst, fCorrel, fRMRSS, fI0_IsotpeAccuracy;

			String ^sPeptide, ^sIsPeptideUnique;
	};

	public ref class ProteinRateConstant
	{
		// TODO: Add your methods for this class here.

		public:

			int nC, nH, nN, nO, nS, nP;

			int fNparameter;

			int nExHydrogen, nIsotopes;

			array <float>  ^fBodyWaterEnrichment;

			array <double> ^fTimeCourse, ^ fUniqueTimePoints; //if the Replicates is true, determine the unique/actual experimental time points;

			double fScaleTime;     //scale the time so that there are not time values larger than 10. Otherwise, LBFGS has difficulities to to converge. 

			float fRateConstant, fCorrThreshold, fRMSS_Threshold;

			String ^workDirectory;   //need to set this

			bool bReplicates;       //set to true if replicates are present

			array <int> ^aiReplicateStructure;  //will store how may replicates at each time point 

			array <double> ^dIsotopeMasses;

			String ^szPeptide;

			ProteinRateConstant(String ^sDirectory, array <float> ^ExperimentTimeCourse, 
				array <float> ^ BodyWaterDeuterium, int Nparameter)
			{
				fNparameter = Nparameter;

				workDirectory = sDirectory;

				fTimeCourse = gcnew array <double> (ExperimentTimeCourse->Length);

				for (int i = 0; i < ExperimentTimeCourse->Length; i++)
				{
					fTimeCourse[i] =  (double) ExperimentTimeCourse[i];
				}

				fBodyWaterEnrichment = BodyWaterDeuterium;

				bReplicates = false;

				nC = nH = nN = nO = nS = nP = 0;
				
				nExHydrogen = 0;

				nIsotopes = 6; // may need to make it truely parametrized

				fCorrThreshold = 0.90f; // may need to make it truely parametrized

				fRMSS_Threshold = 0.05;

				fScaleTime = 0.0;

				for(int i = 0; i < fTimeCourse->Length; i++)
				{
					if(fTimeCourse[i] > fScaleTime)
					{
						fScaleTime = fTimeCourse[i];
					}
				}


				if(fScaleTime == 0)
				{
					printf("No Time Course. All Time points are 0\n");

					exit (1);
				}

				fScaleTime = 13. / fScaleTime;

				//fScaleTime = 1.0;

				for(int i = 0; i < fTimeCourse->Length; i++)
				{
					fTimeCourse[i] = fScaleTime * fTimeCourse[i];
				}

				dIsotopeMasses = gcnew array <double>(64);
			}

			int ElementalComposition(String ^sSequence);

			float ComputeRateConstant(String ^ sInputPeptide, array <float> ^TheoreticalIsotopes, 
				 array <double, 2> ^dExperimentIsotopes);

                        //std::map<std::string, std::vector<float>> ComputeRateFromI0(array <float> ^TheoreticalIsotopes, array <float> ^IonIntensities,
                        //      array <double> ^adLabelIncorporationTime, int nNumberOfTimesforThePeptie, float *fCorr,
                        //      float *fRMRSS);

            float ComputeRateFromI0(array <float> ^TheoreticalIsotopes, array <float> ^IonIntensities,
                                array <double> ^adLabelIncorporationTime, int nNumberOfTimesforThePeptie, float *fCorr,
                                float *fRMRSS);

			void ReadFileFolder();

			void QuantFileReader(String ^ sQuantFile);

			bool CheckANumber(int Number);

			array <float> ^ NumberReader(char *szString, int nStart);

			float PearsonCorrelation(array <float> ^FirstArr, array <float> ^SecondArr, int nPoints);

			float RMSS(array <float> ^FirstArr, array <float> ^SecondArr, int nPoints);

			float PoissonIsotopesWithMasses(int iIsotope, array <float> ^fIsotopes, int nH, int nC, int nN,
						int nO, int nS);

			float PoissonIsotopes(int iIsotope, array <float> ^fIsotopes, int nH, int nC, int nN,
				int nO, int nS);

			float PoissonIsotopesMeanMasses(int iIsotope, array <float> ^fIsotopes, array <double> ^dMassesOfIsotopes, 
				int nH, int nC, int nN, int nO, int nS);

			float HybridIsotopes(int iIsotope, array <float> ^fIsotopes, array <double> ^MassesofIsotopes,
				int nH, int nC, int nN, int nO, int nS);

			bool ProcessTimePoints(array <double> ^fTimePoints);   //will determine if there are replicates

			int AnalyzeReplicates(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters, 
				List <float> ^fIonScores, array <float> ^fTheoreticalIsotopes, array <double> ^,
				array <float> ^);


			int AverageReplicates(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters,
				List <float> ^fIonScores, array <float> ^fTheoreticalIsotopes, array <double> ^,
				array <float> ^, int *);

			int NoReplicateExperiments(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters,
				List <float> ^fIonScores, array <float> ^fTheoreticalIsotopes, array <double> ^TimesForThisPeptide,
				array <float> ^fI0Ions, int *nTrueIons);

			float ExperimentalTheoreticalIsotopeAccuracy(List <SingleIsotopeCluster ^> ^TimeCourseIsotopeClusters,
				array <float> ^fTheoreticalIsotopes);

	};

	public ref class TheoreticalIsotopeCalculator
	{	
		public:


			array <double> ^dIsotopeMasses_FromCalculator;

			String ^szPeptideInCalculator;

			//int nH, nC, nN, nO, nS, nP;
			
			TheoreticalIsotopeCalculator()
			{
				//nH = nC = nN = nO = nS = nP = 0;

				dIsotopeMasses_FromCalculator = gcnew array <double>(64);
			}

			//float PoissonIsotopes(int iIsotope, array <float> ^fIsotopes, int nH, int nC, int nN,
				//		int nO, int nS);

			int AtomicComposition(String ^sSequence, int *nHydrogens, int *nCarbons, int *nNitrogens, 
	                         int *nOxygens, int *nSulfur, int *nPhosphor);

			float Convolve(array <double> ^arr1, array <double> ^arr2);

			float Self_Convolve(int nMultiplex, array <double> ^arr1);

			float FourierIsotopes(array <float> ^fIsotopes, int nH, int nC, int nN,
						int nO, int nS);

			float StartToComputeSequenceIsotopes(String ^sSequence, array <float> ^fIsotopes);

			float StraightConvolution(array <double> ^fArray1, array <double> ^fArray22);

			void Caller_StraightConvolution(array <float> ^fIsotopes, int nH, int nC, int nN, int nO, int nS);

			float DoRealFFT();

			void IsotopeMasses(int nH, int nC, int nN, int nO, int nS);

			float TheoreticalIsotopeCalculator::FourierIsotopes_Temp(array <float> ^fIsotopes, int nH, int nC, int nN,
				int nO, int nS);

			void AssignPeptideSequence(String ^sPeptide);
	};
}
