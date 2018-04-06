// pepXML.h

#pragma once

using namespace System;
using namespace System::Collections::Generic;

namespace pepXML {

	public ref class PeptideEntry
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

		short int nRank;

		bool bUniquePeptide;
		
		String ^sProtId;     //example "DBSeq_1_FIBB_HUMAN"
		String ^sPepID;     // "peptide_4_1"
		String ^sEvidence;   // <PeptideEvidence id="PE_2_4_LV106_HUMAN_0_47_51

		double dIonscore, dHomscore, dIdenscore, dExpect;

		List <int> ^ ModLocations;  //peptide modification locations
		List <double> ^dModMasses;  //modification masses,  <Modification location="4" residues="M" monoisotopicMassDelta="15.994915">
	};

	//temporary peptide storage, stores peptides reported at the start of mzid file
   //this information is not copied to the MzIdent caller
	public ref class PeptideAndID
	{
	public:
		String ^sPeptide;     // <peptideSequence>DILMK</peptideSequence>
		String ^PeptideID;   // <Peptide id="peptide_4_2"
		List <int> ^ ModLocations;  //peptide modification locations
		List <double> ^dModMasses;  //modification masses, <Modification location="4" residues="M" monoisotopicMassDelta="15.994915">
	};
	/*
	*  temporal protein list to hold protein description, and length 
	*  from mzid file. It is necessary, because in mzid protein length
	*  and description appear much before the protein score information
	*  and in a different area.
	*/
	public ref class ProteinDescrip
	{
	public:
		int nSeqLength;
		String^ accession;
		String ^sDescription;       //protein description from Mascot
		String ^sId;

		//System::Collections::Generic::List<PeptideEntry^> ^ PeptideList;
};
	/*
	*  entry to store protein summary information at the end of the 
	*  mzid file
	*/
	public ref class ProteinEntry
	{
public:
	double ProteinScore;
	double SeqCoverage;    //sequence coverage
	double ProteinMass;
	int nSeqLength, nDistinctSequences;
	String^ accession, ^description;       //protein description from Mascot
	int nSpectralCount;

	System::Collections::Generic::List<PeptideEntry^> ^ PeptideList;

	List <String ^> ^PeptideEvidence;
};

	public ref class pepX 
	{
	public:

			System::Collections::Generic::List<PeptideEntry^> ^ ResultList;
			array<double, 2> ^ Spectrum;
			array<double, 2> ^ aaMod;

			String ^ sFilepepxml;

			pepX(String ^ sInputFile)
			{
				sFilepepxml = sInputFile;
			}

			int ReadPepXML(double minScore, double minDeltaScore, int minSecondRank);

			int ReadMascotPepXML(double minScore, double minDeltaScore, int minSecondRank);
	};

	public ref class mzIdentML 
	{
	public:

			System::Collections::Generic::List<PeptideEntry^> ^ ResultList;

			System::Collections::Generic::List<ProteinEntry^> ^ ProteinResultList;

			array<double, 2> ^ Spectrum;
			array<double, 2> ^ aaMod;

			String ^ smzIdentML;

			mzIdentML(String ^ sInputFile)
			{
				smzIdentML = sInputFile;
			}

			void ReadProteins_mzIdent();

			void ReadProteins_mzIdent_2();

			void mzIdent();
	};
}


