#pragma once

#include <vector>

using namespace System;
using namespace System::Collections::Generic;

//using namespace std;


//function definitions
void ReadFileNames(array<System::String ^> ^inputs, List <String ^> ^asmzML, List <String ^> ^asmzID);

void ProteinRates(array <float> ^, array <float> ^, int);

int compare_float (const void * a, const void * b);



float QSelect_Float(List <float> ^lArray);

void Read_Params(array<System::String ^> ^args);

namespace ProteinCollector
{
	// a peptide entry
	public ref class PeptideHolder
	{
		public:
		Byte nCharge;
		Byte nDuplicity;
		double FirstScore;
		double deltaScore;
		double SeqMass;
		double SpecMass;
		String^ Peptide;
		String^ Protein;
		float SecondScore;
		int RankSecondScore;
		long int nScan;
		float dRetTime;
		float fStartElution, fEndElution;
		float fFirstID_Elution, fLastID_Elution;     //first and last times that a peptide has been identified via MS/MS

		short int nRank;

		bool bUniquePeptide;           //it is true for a peptide that is unique to a protein, false if peptide sequence is shared between proteins
		
		String ^sProtId;     //example "DBSeq_1_FIBB_HUMAN"
		String ^sPepID;     // "peptide_4_1"
		String ^sEvidence;   // <PeptideEvidence id="PE_2_4_LV106_HUMAN_0_47_51

		array <double, 2> ^dIsotopes;

		array <float> ^fTheoreticIsotopes;

		//Added by JA 2016.09.16
		//Miscellaneous items for passing (Peak half Width, Left Width, Right Width, Retention Time Peak)
		array <double, 1> ^w;
		//

		//Added by JA 2016.12.07
		float NetLabeling;

		//Added by JA 2016.12.07
		bool bPeptideTimePoint;


		bool bQuant, bPeptidePassed;        // if set true, will read peptide abundance from MS1

		double dIonscore, dHomscore, dIdenscore, dExpect;

		List <int> ^ ModLocations;  //peptide modification locations
		List <double> ^dModMasses;  //modification masses, <Modification location="4" residues="M" monoisotopicMassDelta="15.994915">
	};
/*
*  a reference class that will hold protein results of a database
*  search. Every instance of this class will hold results from a s
*  single experiment. A list defined on this class will hold results
*  from multiple experiments
*
*/
	public ref class ProteinSet   //holds information about a single protein
		{
		public:
			double ProteinScore, ProteinMass;
			double SeqCoverage;    //sequence coverage
			int nSeqLength, nDistinctSequences;
			String^ accession, ^description;       //protein description from Mascot
			int nSpectralCount;

			bool bProteinPassed;

			int ProteinIDnumber;

			//Added by JA 2016.12.07
			bool bProteinTimePoint;

			List <PeptideHolder ^> ^Peptides;

			List <float> ^NetLabeling;

			int ProteinIndex;

			List <int> ^PeptideIndex;

			//Added by JA 2016.10.14
		double sumMPE;
		int countMPE;
		};

	public ref class ProteinCollection   //holder for all proteins from a single experiment
	{
	public:
		System::Collections::Generic::List<ProteinSet^> ^ ProteinsList;

		String ^sExperimentFile, ^smzML;

		float fExperimentTime;    //the time that data was collected, eg. 0 days, 0.3 days, 1 day, 4 days, ....60 days

		float fBWE;               // Body water enrichment at the time that the sample was collected.
	};

	public ref class ExperimentCollection   //holder for all proteins from mulitple experiments
	{
	public:
		System::Collections::Generic::List<ProteinCollection^> ^ ExperimentsList;
	};

	public ref class SingleList   //holder for a single IntegerList
	{
	public:
		
		List <int> ^ IntegerList;

		List <double> ^DoubleList;

		List <String ^> ^StringList;
	};


	public ref class GeneralListCollector
	{
	public:
		List <SingleList ^> ^aList;

	};

	public ref class ProteinList
	{
	public:
		String ^Accession, ^Protein;

		bool bProteinPassed;

		List <int> ^ ProteinInExperiments; 

		List <float> ^ ExperimentTime;        //this is the incorporation time, like 0, 3, 5, ... days

		//Added by JA 2016.10.14
		double sumMPE;
		int countMPE;
	};


	public ref class PeptideList
	{
	public:

		PeptideHolder ^Peptide;

		bool bPeptidePassed;

		List <int> ^ PeptideInExperiments;

		List <float> ^ ExperimentTime;        //this is the incorporation time, like 0, 3, 5, ... days
	};
} //name psace ProteinCollector

typedef struct _OneSet
{
	int Value;
} OneSet;

typedef std::vector <OneSet> ListCollector;
