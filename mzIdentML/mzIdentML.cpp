// pepXML.cpp : main project file.
// This is the main DLL file.

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "mzIdentML.h"
#include "math.h"
#using <System.Xml.dll> 

using namespace System;
using namespace System::IO;
using namespace System::Xml;
using namespace System::Collections::Generic;
using namespace pepXML;

int ReadScanNumber(String ^sScanLine);

String^ ReadNode(XmlTextReader^ reader,String^ nodeName) 
{
	String^ nodeValue;

	while (reader->MoveToNextAttribute()) 
	{
	  
	 //if (reader->Name==nodeName) //modified by Mahbubur Rahman date : 10/02/2015
		if(System::String::Compare(reader->Name,nodeName,true) == 0) //modified by Mahbubur Rahman date : 10/02/2015
		{ 
			nodeValue=reader->Value;

			break;
		}
		
	}
	return nodeValue; 
}

//*************

// a method returns the lenght of the spectrum
// (scan number) nScan

 int pepX::ReadPepXML(double minScore, double minDeltaScore, int minSecondRank)

{

float fragment_ion_tol;
float peptide_mass_tol;
Byte hit_rank;
String ^ nvalue;
String ^ base_name;
String ^ spectrum;
String ^ peptide;
String ^ peptide_prev_aa;
String ^ peptide_next_aa;
String ^ protein;
String ^ local_path;
String ^ symbol;
int num_tot_proteins, num_matched_ions, nScan, nHit_Rank; 
int tot_num_ions, num_tol_term, num_missed_cleavages, is_rejected;
double calc_neutral_pep_mass, massdiff, FirstScore, deltaScore;
double SpecMass, SeqMass;
float SecondScore;
Byte nCharge, nDuplicity;
int RankSecondScore;
Boolean cterMod=false, bFirst = false;





ResultList = gcnew System::Collections::Generic::List<PeptideEntry^> ;
PeptideEntry ^ Results = gcnew PeptideEntry;
XmlTextReader^ reader = nullptr;
aaMod = gcnew array<double, 2> (2,10);
	try
	{
      reader = gcnew XmlTextReader(sFilepepxml);
      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
		 {
         case XmlNodeType::Element:

			if("search_summary" == reader->Name )
			{
				nvalue=ReadNode(reader,"base_name");
				base_name=nvalue;
 			}                     

			if("search_database" == reader->Name )
			{
				nvalue=ReadNode(reader,"local_path");
				local_path=nvalue;
 			}                     


/*			if("aminoacid_modification" == reader->Name )
			{
							nvalue=ReadNode(reader,"variable");
							varMod =nvalue;
							reader->MoveToAttribute(0);
							while (reader->MoveToNextAttribute())  
							{		
								if (reader->Name == "aminoacid" && varMod == "Y")
								{ aminoacid	=	reader->Value;
									aaMod
								}
								if (reader->Name == "massdiff")					massdiff	=	double::Parse(reader->Value);
								if (reader->Name == "mass")						mMass		=	double::Parse(reader->Value);
								if (reader->Name == "symbol")					varSymbol	=	reader->Value;
								

							}
			}
*/

			if (reader->Name == "terminal_modification")	
			{
				nvalue=ReadNode(reader,"symbol"); 
				symbol=nvalue; 
			}			

			if("parameter" == reader->Name ) 
			{

				nvalue=ReadNode(reader,"name");

				if (nvalue=="peptide_mass_tol")
				{
					nvalue = reader->Value;  
					peptide_mass_tol=float::Parse(nvalue); 

				}

				if (nvalue=="fragment_ion_tol")
				{
					nvalue = reader->Value;  
					fragment_ion_tol=float::Parse(nvalue); 

				}

			}
		
			if("spectrum_query" == reader->Name ) 
			{
							reader->MoveToAttribute(0);
							while (reader->MoveToNextAttribute())  
							{		
								if (reader->Name == "spectrum")					spectrum = reader->Value;
								if (reader->Name == "start_scan")				nScan = int::Parse(reader->Value);
								if (reader->Name == "precursor_neutral_mass")	SpecMass = double::Parse(reader->Value);
								if (reader->Name == "assumed_charge")			nCharge = Byte::Parse(reader->Value);

							}
			}

 
			if("search_hit" == reader->Name)  
			{
				cterMod = false;

				reader->MoveToAttribute(0);

				//printf("%s\n", reader->Name);

				if (reader->Name == "hit_rank")
				{
					hit_rank = Byte::Parse(reader->Value);

					nHit_Rank = int::Parse(reader->Value);
				}

				//if(1 == nHit_Rank)
				{
					while (reader->MoveToNextAttribute())  
					{	
						if (reader->Name == "peptide")					peptide = reader->Value;
						if (reader->Name == "peptide_prev_aa")			peptide_prev_aa = reader->Value;
						if (reader->Name == "peptide_next_aa")			peptide_next_aa = reader->Value;
						if (reader->Name == "protein")					protein = reader->Value;
						if (reader->Name == "num_tot_proteins")			num_tot_proteins = int::Parse(reader->Value);
						if (reader->Name == "num_matched_ions")			num_matched_ions = int::Parse(reader->Value);
						if (reader->Name == "tot_num_ions")				tot_num_ions = int::Parse(reader->Value);
						if (reader->Name == "calc_neutral_pep_mass")	calc_neutral_pep_mass = double::Parse(reader->Value);
						if (reader->Name == "massdiff")					massdiff= double::Parse(reader->Value);
						if (reader->Name == "num_tol_term")				num_tol_term = int::Parse(reader->Value);
						if (reader->Name == "num_missed_cleavages")		num_missed_cleavages = int::Parse(reader->Value);
						if (reader->Name == "is_rejected")				is_rejected = int::Parse(reader->Value);
						if (reader->Name == "num_tot_proteins")			nDuplicity = Byte::Parse(reader->Value);
						if (reader->Name == "calc_neutral_pep_mass")	SeqMass = double::Parse(reader->Value);
						if (reader->Name == "num_tot_proteins")			nDuplicity = Byte::Parse(reader->Value);

					}
				}
			}

			//printf("nHi %d %s\n", nHit_Rank, peptide);
			
			if (reader->Name == "modification_info" && nHit_Rank == 1)	
			{
				while(reader->MoveToNextAttribute())
				{	
					if(reader->Name->IndexOf("mod_cterm_mass") != -1)
					{
						cterMod = true;

						break;
					}
					else
					{
						cterMod = false;
					}
				}

				
				/*if(reader->Name->IndexOf("mod_cterm_mass") != -1)
				{
					cterMod = true;
				}
				else
					cterMod = false;*/

				//nvalue=ReadNode(reader,"mod_cterm_mass"); 
				
				/*if (double::Parse(nvalue) > 0.1)
				{
					cterMod = true; 

					printf("MODDEDDD PEPE = %d %s %s\n", nHit_Rank, peptide, nvalue);
				}
				else
				{
					cterMod = false;
				}*/
			}		

			if("search_score" == reader->Name && nHit_Rank == 1) 	
			{
						nvalue = ReadNode(reader,"name"); 
						if ( nvalue=="xcorr")
						{
							FirstScore = double::Parse(reader->Value);
						}			
						if ( nvalue=="deltacn")
						{
							deltaScore  = double::Parse(reader->Value);
						}			
						if ( nvalue=="spscore")
						{
							SecondScore  = float::Parse(reader->Value); 
						}			
						if ( nvalue=="sprank")
						{
							RankSecondScore  = int::Parse(reader->Value); 
						}			
			}

			break;

		 case XmlNodeType::EndElement:

			 if (reader->Name=="search_hit" && nHit_Rank == 1)    
			 {
				 Results = gcnew PeptideEntry;
				 if (cterMod) 
				 {
					peptide += symbol;
				 }
				 cterMod=false;
				 peptide_prev_aa += "."; 
				 peptide += "."; 
				 peptide = peptide_prev_aa + peptide;
				 peptide += peptide_next_aa;
				 Results->Peptide = peptide; 
				 Results->Protein = protein;
				 Results->FirstScore = FirstScore;
				 Results->deltaScore = deltaScore;
				 Results->nCharge = nCharge;
				 Results->RankSecondScore = RankSecondScore;
				 Results->SecondScore = SecondScore;
				 Results->SeqMass = SeqMass;
				 Results->SpecMass = SpecMass;
				 Results->nDuplicity = nDuplicity;
				 Results->nCharge = nCharge;
				 Results->nScan = nScan;
				 
// Filters.... ????
					if (FirstScore>minScore)
					{
						ResultList->Add(Results);
					}
				
			 }
            break; 
		 } 

} 
	   while (reader->Read());    

return 0;

   }
   finally
   { 

      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilepepxml);
      if ( reader != nullptr )
            reader->Close();
   }
	return 0;
 }




 //a method that reads result from MascotPepXML file

 //*************

// a method returns the lenght of the spectrum
// (scan number) nScan

 int pepX::ReadMascotPepXML(double minScore, double minDeltaScore, int minSecondRank)

{

float fragment_ion_tol;
float peptide_mass_tol;
Byte hit_rank;
String ^ nvalue;
String ^ base_name;
String ^ spectrum;
String ^ peptide;
String ^ peptide_prev_aa;
String ^ peptide_next_aa;
String ^ protein;
String ^ local_path;
String ^ symbol;
int num_tot_proteins, num_matched_ions, nScan, nHit_Rank; 
int tot_num_ions, num_tol_term, num_missed_cleavages, is_rejected;
double calc_neutral_pep_mass, massdiff, FirstScore, deltaScore;
double SpecMass, SeqMass;
float SecondScore;
Byte nCharge, nDuplicity;
double dExpection;
Boolean cterMod=false, bFirst = false;

ResultList = gcnew System::Collections::Generic::List<PeptideEntry^> ;
PeptideEntry ^ Results = gcnew PeptideEntry;
XmlTextReader^ reader = nullptr;
aaMod = gcnew array<double, 2> (2,10);
	try
	{
      reader = gcnew XmlTextReader(sFilepepxml);
      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
		 {
         case XmlNodeType::Element:

			if("search_summary" == reader->Name )
			{
				nvalue=ReadNode(reader,"base_name");
				base_name=nvalue;
 			}                     

			if("search_database" == reader->Name )
			{
				nvalue=ReadNode(reader,"local_path");
				local_path=nvalue;
 			}                     


/*			if("aminoacid_modification" == reader->Name )
			{
							nvalue=ReadNode(reader,"variable");
							varMod =nvalue;
							reader->MoveToAttribute(0);
							while (reader->MoveToNextAttribute())  
							{		
								if (reader->Name == "aminoacid" && varMod == "Y")
								{ aminoacid	=	reader->Value;
									aaMod
								}
								if (reader->Name == "massdiff")					massdiff	=	double::Parse(reader->Value);
								if (reader->Name == "mass")						mMass		=	double::Parse(reader->Value);
								if (reader->Name == "symbol")					varSymbol	=	reader->Value;
								

							}
			}
*/

			if (reader->Name == "terminal_modification")	
			{
				nvalue=ReadNode(reader,"symbol"); 
				symbol=nvalue; 
			}			

			if("parameter" == reader->Name ) 
			{

				nvalue=ReadNode(reader,"name");

				if (nvalue=="peptide_mass_tol")
				{
					nvalue = reader->Value;  
					peptide_mass_tol=float::Parse(nvalue); 

				}

				if (nvalue=="fragment_ion_tol")
				{
					nvalue = reader->Value;  
					fragment_ion_tol=float::Parse(nvalue); 

				}

			}
		
			if("spectrum_query" == reader->Name ) 
			{
							reader->MoveToAttribute(0);
							while (reader->MoveToNextAttribute())  
							{		
								if (reader->Name == "spectrum")					spectrum = reader->Value;
								if (reader->Name == "start_scan")				nScan = int::Parse(reader->Value);
								if (reader->Name == "precursor_neutral_mass")	SpecMass = double::Parse(reader->Value);
								if (reader->Name == "assumed_charge")			nCharge = Byte::Parse(reader->Value);

							}
			}

 
			if("search_hit" == reader->Name)  
			{
				cterMod = false;

				reader->MoveToAttribute(0);

				//printf("%s\n", reader->Name);

				if (reader->Name == "hit_rank")
				{
					hit_rank = Byte::Parse(reader->Value);

					nHit_Rank = int::Parse(reader->Value);
				}

				//if(1 == nHit_Rank)
				{
					while (reader->MoveToNextAttribute())  
					{	
						if (reader->Name == "peptide")					peptide = reader->Value;
						if (reader->Name == "peptide_prev_aa")			peptide_prev_aa = reader->Value;
						if (reader->Name == "peptide_next_aa")			peptide_next_aa = reader->Value;
						if (reader->Name == "protein")					protein = reader->Value;
						if (reader->Name == "num_tot_proteins")			num_tot_proteins = int::Parse(reader->Value);
						if (reader->Name == "num_matched_ions")			num_matched_ions = int::Parse(reader->Value);
						if (reader->Name == "tot_num_ions")				tot_num_ions = int::Parse(reader->Value);
						if (reader->Name == "calc_neutral_pep_mass")	calc_neutral_pep_mass = double::Parse(reader->Value);
						if (reader->Name == "massdiff")					massdiff= double::Parse(reader->Value);
						if (reader->Name == "num_tol_term")				num_tol_term = int::Parse(reader->Value);
						if (reader->Name == "num_missed_cleavages")		num_missed_cleavages = int::Parse(reader->Value);
						if (reader->Name == "is_rejected")				is_rejected = int::Parse(reader->Value);
						if (reader->Name == "num_tot_proteins")			nDuplicity = Byte::Parse(reader->Value);
						if (reader->Name == "calc_neutral_pep_mass")	SeqMass = double::Parse(reader->Value);
						if (reader->Name == "num_tot_proteins")			nDuplicity = Byte::Parse(reader->Value);

					}
				}
			}

			//printf("nHi %d %s\n", nHit_Rank, peptide);
			
			if (reader->Name == "modification_info" && nHit_Rank == 1)	
			{
				while(reader->MoveToNextAttribute())
				{	
					if(reader->Name->IndexOf("mod_cterm_mass") != -1)
					{
						cterMod = true;

						break;
					}
					else
					{
						cterMod = false;
					}
				}

				
				/*if(reader->Name->IndexOf("mod_cterm_mass") != -1)
				{
					cterMod = true;
				}
				else
					cterMod = false;*/

				//nvalue=ReadNode(reader,"mod_cterm_mass"); 
				
				/*if (double::Parse(nvalue) > 0.1)
				{
					cterMod = true; 

					printf("MODDEDDD PEPE = %d %s %s\n", nHit_Rank, peptide, nvalue);
				}
				else
				{
					cterMod = false;
				}*/
			}		

			if("search_score" == reader->Name && nHit_Rank == 1) 	
			{
						nvalue = ReadNode(reader,"name");

						if ( nvalue=="ionscore")
						{
							FirstScore = double::Parse(reader->Value);
						}			
						if ( nvalue=="identityscore")
						{
							deltaScore  = double::Parse(reader->Value);
						}			
						if ( nvalue=="homologyscore")
						{
							SecondScore  = float::Parse(reader->Value); 
						}			
						if ( nvalue=="expect")
						{
							dExpection  = double::Parse(reader->Value); 
						}			
			}

			break;

		 case XmlNodeType::EndElement:

			 if (reader->Name=="search_hit" && nHit_Rank == 1)    
			 {
				 Results = gcnew PeptideEntry;
				 if (cterMod) 
				 {
					peptide += symbol;
				 }
				 cterMod=false;
				 peptide_prev_aa += "."; 
				 peptide += "."; 
				 peptide = peptide_prev_aa + peptide;
				 peptide += peptide_next_aa;
				 Results->Peptide = peptide; 
				 Results->Protein = protein;
				 Results->dIonscore = FirstScore;
				 Results->dIdenscore = deltaScore;
				 Results->nCharge = nCharge;
				 Results->dExpect = dExpection;
				 Results->dHomscore = SecondScore;
				 Results->SeqMass = SeqMass;
				 Results->SpecMass = SpecMass;
				 Results->nDuplicity = nDuplicity;
				 Results->nCharge = nCharge;
				 Results->nScan = nScan;
				 
// Filters.... ????
					if (FirstScore>minScore)
					{
						ResultList->Add(Results);
					}
				
			 }
            break; 
		 } 

} 
	   while (reader->Read());    

return 0;

   }
   finally
   { 

      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilepepxml);
      if ( reader != nullptr )
            reader->Close();
   }
	return 0;
 }


/*
* A method to read mzIdentML files
*
* 
*/

void mzIdentML::mzIdent()
{

	float fragment_ion_tol;
	float peptide_mass_tol;
	Byte hit_rank;
	String ^ nvalue;
	String ^ base_name;
	String ^ spectrum;
	String ^ peptide;
	String ^ peptide_prev_aa;
	String ^ peptide_next_aa;
	String ^ protein;
	String ^ local_path;
	String ^ symbol;
	int num_tot_proteins, num_matched_ions, nScan, nHit_Rank; 
	int tot_num_ions, num_tol_term, num_missed_cleavages, is_rejected;
	double calc_neutral_pep_mass, massdiff, FirstScore, deltaScore;
	double SpecMass, SeqMass;
	float SecondScore;
	Byte nCharge, nDuplicity;
	double dExpection;
	Boolean cterMod=false, bFirst = false, bProteinAmb = false;

	int RankSecondScore;

	double minScore = 0.0;

	ResultList = gcnew System::Collections::Generic::List<PeptideEntry^> ;
	PeptideEntry ^ Results = gcnew PeptideEntry;
	XmlTextReader^ reader = nullptr;
	aaMod = gcnew array<double, 2> (2,10);
	try
	{
      reader = gcnew XmlTextReader(smzIdentML);
      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
		 {
         case XmlNodeType::Element:

			if("search_summary" == reader->Name )
			{
				nvalue=ReadNode(reader,"base_name");
				base_name=nvalue;
 			}                     

			if("search_database" == reader->Name )
			{
				nvalue=ReadNode(reader,"local_path");
				local_path=nvalue;
 			}                     

			if("ProteinAmbiguityGroup" == reader->Name)
			{
				ReadNode(reader,"id");

				bProteinAmb = true;

				/*while (reader->MoveToNextAttribute())  
				{
					printf("%s\n", reader->Value);

					if (reader->Name == "value")	
					{
						FirstScore = double::Parse(reader->Value);

						break;
					}
				}*/
			}

			if("cvParam" == reader->Name && bProteinAmb ) 
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1001171")
				{
					reader->MoveToAttribute(0);

					while (reader->MoveToNextAttribute())  
					{
						if (reader->Name == "value")	
						{
							FirstScore = double::Parse(reader->Value);

							printf("Print1 %s %10.5f\n", reader->Value, FirstScore);

							break;
						}
					}
				}
			}

			break;

		 case XmlNodeType::EndElement:

			 if (reader->Name=="search_hit" && nHit_Rank == 1)    
			 {
				 Results = gcnew PeptideEntry;
				 if (cterMod) 
				 {
					peptide += symbol;
				 }
				 cterMod=false;
				 peptide_prev_aa += "."; 
				 peptide += "."; 
				 peptide = peptide_prev_aa + peptide;
				 peptide += peptide_next_aa;
				 Results->Peptide = peptide; 
				 Results->Protein = protein;
				 Results->FirstScore = FirstScore;
				 Results->deltaScore = deltaScore;
				 Results->nCharge = nCharge;
				 Results->RankSecondScore = RankSecondScore;
				 Results->SecondScore = SecondScore;
				 Results->SeqMass = SeqMass;
				 Results->SpecMass = SpecMass;
				 Results->nDuplicity = nDuplicity;
				 Results->nCharge = nCharge;
				 Results->nScan = nScan;
				 
// Filters.... ????
					if (FirstScore>minScore)
					{
						ResultList->Add(Results);
					}
				
			 }

			 if("ProteinAmbiguityGroup" == reader->Name)
			 {
				bProteinAmb = false;
				
				//printf("Falsigying protein\n");
			 }

            break; 
		 } 

} 
	   while (reader->Read());    

return ;

   }
   finally
   { 

      Console::WriteLine( "\nProcessing of the file {0} complete.\n", smzIdentML);
      if ( reader != nullptr )
            reader->Close();
   }
}

/*
* A method to read mzIdentML files
* read proteins
* num_tot_proteins - total number of proteins
* TempProtCollection  temporaroly holds  the protein description,
*     length and other information, which is stored at the start of mzid
*
* n1Scan is the scan number read in the Start of the SpectrumIdentificationResult, e.g.,
* SpectrumIdentificationResult id="SIR_10632" spectrumID="mzMLid=controllerType=0 controllerNumber=1 scan=13434" spectraData_ref="SD_1"
* 
* n2Scan is the scan number read in at the End of the SpectrumIdentificationResult block from the line:
cvParam accession="MS:1000796" name="spectrum title"  cvRef="PSI-MS" value="controllerType=0 controllerNumber=1 scan=13434"
* the should be equal.
*/

void mzIdentML::ReadProteins_mzIdent_2()
{
	String ^ nvalue, ^sPeptideEvidence;

	//String ^ protein, ^peptide, ^sTemp, ^sTempProtID;
	String  ^protein, ^sTemp, ^sTempProtID;
	
	String ^rank1_peptide, ^rank1_sTemPortID;

	List <String ^> ^ListEvidence;

	int num_tot_proteins = 0, iRank = 0, iCharge, iTemp, n1Scan, n2Scan;

	int nSpectralCount= 0, nDistinctSeqs = 0, nModLocation = 0;
	
	double dProteinScore, dSeqCoverage = 0;

	double dSeqMass, dExpMass, dRetTime, dModMass= 0.0;

	Boolean bProteinAmb = false, bRankOne = false, bSpectumId = false, bSpecResult = false;

	bool bDBSequence = false, bPeptideEvidence = false;


	Boolean sportAddCheck = false; //Check if the db sequence can be added


	PeptideAndID ^sTempPeptideID = gcnew PeptideAndID();

	List <PeptideAndID ^> ^lTempPeptides = gcnew List <PeptideAndID ^>;

	ProteinDescrip ^tempProtDes = gcnew ProteinDescrip();

	PeptideEntry ^tempPeptide = gcnew PeptideEntry();

	List <PeptideEntry ^> ^ListPeptides = gcnew List <PeptideEntry ^>;

	System::Collections::Generic::List<ProteinDescrip^> ^ TempProtCollection = 
		gcnew System::Collections::Generic::List<ProteinDescrip^>;

	ProteinEntry ^currProtein;

	ProteinResultList = gcnew System::Collections::Generic::List<ProteinEntry^>;
	
	XmlTextReader^ reader = nullptr;
	
	try
	{
      reader = gcnew XmlTextReader(smzIdentML);

      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
		 {
         case XmlNodeType::Element:  

			//section to read protein description, accession, start of mzid

			if("DBSequence" == reader->Name)
			{
				bDBSequence = true;

				tempProtDes = gcnew ProteinDescrip();

				//tempProtDes->PeptideList = gcnew List<PeptideEntry ^>;

				nvalue = ReadNode(reader, "id");

				if(String::IsNullOrEmpty(nvalue))
				{
					printf("FOUND IT\n");
				}

				tempProtDes->sId = nvalue;

				nvalue = ReadNode(reader, "length");

				tempProtDes->nSeqLength = int::Parse(nvalue);

				nvalue = ReadNode(reader,"accession");

				tempProtDes->accession = nvalue;
			}

			if(bDBSequence)
			{
				if("cvParam" == reader->Name)
				{
					nvalue = ReadNode(reader,"accession");

					if("MS:1001088" == nvalue)   //for protein description
					{
						tempProtDes->sDescription = ReadNode(reader, "value");
						//store the protein infor
						//TempProtCollection->Add(tempProtDes);
					}
				}
			}
			//section to read protein description, accession

			//read and save peptides
			if("Peptide" == reader->Name)
			{
				nvalue = ReadNode(reader, "id");

				sTempPeptideID = gcnew PeptideAndID();

				sTempPeptideID->ModLocations = gcnew List <int>;

				sTempPeptideID->dModMasses  = gcnew List <double>;
				
				sTempPeptideID->PeptideID = nvalue;

				//printf("Temp %s\n", sTempPeptideID->PeptideID);
			}
			else if("peptideSequence" == reader->Name || "PeptideSequence" == reader->Name)
			{
				
				nvalue = reader->ReadString();

				sTempPeptideID->sPeptide = nvalue;

				//printf("PeptideSequence and nvalue %s\n",nvalue);
				//lTempPeptides->Add(sTempPeptideID);    //MAY NEED TO REINSTATE THIS
			}
			
			else if("Modification" == reader->Name)
			{
				
				//printf("Modification\n");

				nvalue = ReadNode(reader, "location");

				nModLocation = int::Parse(nvalue);

				//printf("nModLocation %d\n",nModLocation);

				sTempPeptideID->ModLocations->Add(nModLocation);

				nvalue = ReadNode(reader, "monoisotopicMassDelta");

				dModMass = float::Parse(nvalue);

				//printf("monoisotopicMassDelta %f\n",dModMass);

				sTempPeptideID->dModMasses->Add(dModMass);

				//printf("ModLoc = %d %s\n", int::Parse(nvalue), sTempPeptideID->sPeptide);

			}

			//end of peptide reading
			
			//if("PeptideEvidence" == reader->Name) //Modified by Mahbubur Rahman to capture peptide evidence and sportid; date: 10/06/2015 
			//{

			//	nvalue = ReadNode(reader, "id");

			//	sPeptideEvidence = nvalue;

			//	//nvalue = ReadNode(reader, "DBSequence_Ref");

			//	////tempPeptide->sProtId = nvalue;
			//	//sTempProtID = nvalue;

			//	//tempPeptide->sEvidence = sPeptideEvidence;
			//	//ListPeptides->Add(tempPeptide); //Adding peptide

			//	//ListPeptides[ListPeptides->Count-1]->sProtId = sTempProtID; //adding db sequence


			//}



			//section to read peptide entries;
            if("SpectrumIdentificationResult" == reader->Name)
			{
				bSpecResult = true;

				nvalue = ReadNode(reader, "SpectrumID");

				n1Scan = ReadScanNumber(nvalue);

				//printf("Scan Number = %d\n", n1Scan);
			}

			if("SpectrumIdentificationItem" == reader->Name)
			{
				bSpectumId = true;
				
				tempPeptide = gcnew PeptideEntry();

				nvalue = ReadNode(reader, "calculatedMassToCharge");

				dSeqMass = float::Parse(nvalue);

				nvalue = ReadNode(reader, "chargeState");

				iCharge = int::Parse(nvalue);

				nvalue = ReadNode(reader, "experimentalMassToCharge");
				
				dExpMass = float::Parse(nvalue);

				nvalue = ReadNode(reader, "Peptide_ref");

				tempPeptide->sPepID = nvalue;

				nvalue = ReadNode(reader, "rank");

				iRank = int::Parse(nvalue);

				if(1 == iRank)
				{
					tempPeptide->nScan   = n1Scan;

					tempPeptide->SeqMass = dSeqMass;

					tempPeptide->nCharge = iCharge;

					tempPeptide->SpecMass = dExpMass;

					tempPeptide->nRank = iRank;

					bRankOne = true;


					//////////////modification 10/09/2015////////////////
					
					
				    sPeptideEvidence = gcnew String("NONE");

					sTempProtID = gcnew String("NONE");      //some times the is no evidence reported for a peptide in mzid    change : 10/06/2015
				}
				else
				{
					bRankOne = false;
				}
			}

			if(bSpectumId && "cvParam" == reader->Name && bRankOne)
			{
				nvalue = ReadNode(reader,"accession");
				
				if("MS:1001171" == nvalue)   //scores, Mascot
				{
					tempPeptide->dIonscore = float::Parse(ReadNode(reader, "value"));
				}
				else if("MS:1001172" == nvalue)
				{
					tempPeptide->dExpect = double::Parse(ReadNode(reader, "value"));
				}
				else if ("MS:1001175" == nvalue)    // shared peptide
				{
					tempPeptide->bUniquePeptide = false;
				}
				else if ("MS:1001363" == nvalue)    // unique peptide
				{
					tempPeptide->bUniquePeptide = true;
				}

				//printf("MASCOT SCORES %10.5f\n", dIonScore);

				bPeptideEvidence = false;   //read about bPeptideEvidence
			}

			
			if(bSpecResult && "cvParam" == reader->Name)
			{
				nvalue = ReadNode(reader,"accession");

				if("MS:1001114" == nvalue || "MS:1000894" == nvalue || "MS:1000016" == nvalue)   //retention time <cvParam accession="MS:1001114" name="retention time(s)"  cvRef="PSI-MS" value="655.774" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" 
				{
					dRetTime  = float::Parse(ReadNode(reader, "value"));
				}
				//else if("MS:1000797" == nvalue)   //cvParam accession="MS:1000797" name="peak list scans"  cvRef="PSI-MS" value="1194" 
				else if("MS:1000796" == nvalue)     //name="spectrum title"  cvRef="PSI-MS" value="15apr2604.21098.21098.5 File:&quot;15apr2604.raw&quot;, 
                {                                    //NativeID:&quot;controllerType=0 controllerNumber=1 scan=21098&quot;" 
					//nvalue = ReadNode(reader, "value");

					String^ stemp;

					stemp = ReadNode(reader, "value");

					n2Scan = ReadScanNumber(stemp);
				}
			}
			//read the peptide evidence information;

			// <PeptideEvidence id="PE_2_4_LV106_HUMAN_0_47_51" start="47" end="51" pre="K" post="D" missedCleavages="0" 
			// added bPeptideEvidence, so that to read the evidence only once. It happens that in some proteins, the same
			// sequence is encountered multiple times, and then its peptide evidence is multiplexed as well, e.g.,
			//  <PeptideEvidence id="PE_481_1_SPA_STAA8_0_62_69" start="62" end="69" pre="R" post="D" missedCleavages="0" isDecoy="false" DBSequence_Ref="DBSeq_1_SPA_STAA8" />
			//  <PeptideEvidence id="PE_481_1_SPA_STAA8_0_123_130" start="123" end="130" pre="R" post="D" missedCleavages="0" isDecoy="false" DBSequence_Ref="DBSeq_1_SPA_STAA8" />
						
			
			if(bRankOne && "PeptideEvidence" == reader->Name && false == bPeptideEvidence)   //in Mascot 2.5.1 outputs do not reach here date:02.18.2016
			{
				nvalue = ReadNode(reader, "id");

				sPeptideEvidence = nvalue;

				nvalue = ReadNode(reader, "DBSequence_Ref");

				//tempPeptide->sProtId = nvalue;
				sTempProtID = nvalue;

				bPeptideEvidence = true;
			}

			if(bRankOne && "PeptideEvidenceRef" == reader->Name && false == bPeptideEvidence) //Mahbubur Rahman
			{

				nvalue = ReadNode(reader, "peptideEvidence_ref");

				//printf("peptideEvidence_ref %s\n", nvalue);

				sPeptideEvidence = nvalue;

				rank1_peptide = nvalue;

				//nvalue = ReadNode(reader, "DBSequence_Ref");

				//tempPeptide->sProtId = nvalue;
				//sTempProtID = nvalue;

				bPeptideEvidence = true;
				sportAddCheck = true;
			}
			//end of section for peptide entries

			//read protein information and scores
			if("ProteinAmbiguityGroup" == reader->Name)
			{

				nvalue = ReadNode(reader, "id");

				nSpectralCount = 0;

				bProteinAmb = true;

				ListEvidence = gcnew List <String ^>;
			}

			//read the protein ambiguity node
			if(bProteinAmb)
			{	
				if("cvParam" == reader->Name)
				{
					nvalue=ReadNode(reader,"accession");
				
					if (nvalue=="MS:1001171") // protein score
					{
						dProteinScore = double::Parse(ReadNode(reader, "value"));
					}
					else if (nvalue == "MS:1001093")    //sequence coverage
					{
						dSeqCoverage = double::Parse(ReadNode(reader, "value"));
					}
					else if (nvalue == "MS:1001097")    //number of distinct sequences
					{
						nDistinctSeqs = int::Parse(ReadNode(reader, "value"));
					}
					
				} // end of if("cvparams" == reader->name ...)
				else if("PeptideHypothesis" == reader->Name ) 
				{
					nvalue = ReadNode(reader, "PeptideEvidence_Ref");

					if(rank1_peptide == nvalue)
					{
						sTempProtID = rank1_sTemPortID; //Adding port id 10-12-2015

						/*if(sportAddCheck==true) //10/12/2015
						{ */ 
				            ListPeptides[ListPeptides->Count-1]->sProtId = sTempProtID;  
							//sportAddCheck = false;
						//}

					}

					ListEvidence->Add(nvalue);

					nSpectralCount++;
				}
				else if("ProteinDetectionHypothesis" == reader->Name)
				{
					protein = ReadNode(reader,"id");
					
					rank1_sTemPortID = ReadNode(reader, "dBSequence_ref"); //keeping Db_sequence
						//sTempProtID = nvalue;
						//sportAddCheck = false;
					
					//printf("rank1_sTemp = %s\n", rank1_sTemPortID);	      
				}
			}// if(bProteinAmb)

			break;


		 case XmlNodeType::EndElement:


			 if("ProteinDetectionHypothesis" == reader->Name)
			 {
				//date:02.18.2016 Commented out RGS. Because some protein ambiguity groups
				 //contain several proteins. In that case, the option only reads
				 //the last protein of several contained in ProteinAmbiguity
				if(nSpectralCount > 0)                           
				{
					currProtein = gcnew ProteinEntry();

					currProtein->nSpectralCount = nSpectralCount;

					currProtein->accession = protein;

					currProtein->ProteinScore = dProteinScore;

					currProtein->SeqCoverage = dSeqCoverage;

					currProtein->nDistinctSequences = nDistinctSeqs;

					currProtein->PeptideEvidence = ListEvidence;

					ProteinResultList->Add(currProtein);

                    /*printf("Protein %s SC %d ProtScore %10.2f Cov = %10.1f %d\n", 
							protein, nSpectralCount, dProteinScore, dSeqCoverage, num_tot_proteins);*/

					num_tot_proteins++;
				}
			 }
			 else if("ProteinAmbiguityGroup" == reader->Name)
			 {
				 //ReadNode(reader, "id");

				 //printf("Protein AmbigID %s\n", reader->Name);

				bProteinAmb = false;
				

				/*if(nSpectralCount > 0)      //date:02.18.2016 Commented out RGS. Because some protein ambiguity groups
				{									//contain several proteins. In that case, the option only reads
					currProtein = gcnew ProteinEntry();                //the last protein of several contained in ProteinAmbiguity

					currProtein->nSpectralCount = nSpectralCount;

					currProtein->accession = protein;

					currProtein->ProteinScore = dProteinScore;

					currProtein->SeqCoverage = dSeqCoverage;

					currProtein->nDistinctSequences = nDistinctSeqs;

					currProtein->PeptideEvidence = ListEvidence;

					ProteinResultList->Add(currProtein);

                    printf("Protein %s SC %d ProtScore %10.2f Cov = %10.1f %d\n", 
							protein, nSpectralCount, dProteinScore, dSeqCoverage, num_tot_proteins);

					num_tot_proteins++;
				}*/
			 }
			 else if("SpectrumIdentificationResult" == reader->Name)
			 {
				 bSpecResult = false;
                //update the ret time, as it is per spectrum Identification Result
                if(ListPeptides->Count > 0) //added on 05/16/2016
				{
					//printf("Scan Before = %d\n", ListPeptides[ListPeptides->Count-1]->nScan);


					if(n2Scan == ListPeptides[ListPeptides->Count-1]->nScan)
					{

						ListPeptides[ListPeptides->Count-1]->dRetTime = dRetTime;

						ListPeptides[ListPeptides->Count-1]->sProtId = sTempProtID;  
					}
				
					//ListPeptides[ListPeptides->Count-1]->nScan   = nScan + 1;

					//ListPeptides[ListPeptides->Count-1]->nScan   = n1Scan;

					/*if(ListPeptides[ListPeptides->Count-1]->sPepID->Equals("peptide_10630_1") )
					{
						printf("Scan = %d RetTime = %10.5f Count# = %d\n", 
							ListPeptides[ListPeptides->Count-1]->nScan, ListPeptides[ListPeptides->Count-1]->dRetTime/60.0, 
							ListPeptides->Count);
					}*/
				}
		
				//delete (tempPeptide);
			 }
			 else if("SpectrumIdentificationItem" == reader->Name && bRankOne)
			 {
				tempPeptide->sEvidence = sPeptideEvidence;

				ListPeptides->Add(tempPeptide);

				/*printf("Peptid %s %d %d %s  Count# = %d\n", tempPeptide->Peptide, tempPeptide->nScan, tempPeptide->nCharge,
					tempPeptide->sPepID, ListPeptides->Count);*/

				delete (tempPeptide);

				bRankOne = false;

				bSpectumId = false;
			 }
			 else if("DBSequence" == reader->Name)
			 {
				 bDBSequence = false;

				 //store the protein infor
				TempProtCollection->Add(tempProtDes);
			 }
			 else if("Peptide" == reader->Name)    //finished reading the peptide information 
			 {                                    // add it to the temporary list
				 //printf("End of peptide\n");
				 lTempPeptides->Add(sTempPeptideID);
			 }

            break; 
		 } 

} 
	   while (reader->Read());    

	   if(ProteinResultList->Count != TempProtCollection->Count)
	   {
			printf("Error in the mzid file, %s, different number of proteins %d %d\n", smzIdentML,
				ProteinResultList->Count, TempProtCollection->Count);

			puts("Exitting ...\n");

			exit (1);
	   }

	   printf("Number of Proteins in %s is %d\n", smzIdentML, num_tot_proteins);

	   //combine results from temporal protein holder and protein depository
	   // some results are in protein holder and some are in protein depository
	   // that is why it is needed to combine them

	   printf("Combining Protein Results\n");

	   for(int i=0; i < TempProtCollection->Count; i++)
	   {
			if(ProteinResultList[i]->accession->IndexOf(TempProtCollection[i]->accession) == -1)
			{
				printf("Problem for %s %s\n", ProteinResultList[i]->accession,
						TempProtCollection[i]->accession);

				exit(1);
			}
			else
			{
				ProteinResultList[i]->accession = TempProtCollection[i]->accession;

				ProteinResultList[i]->description = TempProtCollection[i]->sDescription;

				ProteinResultList[i]->nSeqLength = TempProtCollection[i]->nSeqLength;
			}
	   }

	   printf("Finished Combining Protein Results\n");

	   delete (TempProtCollection);

return ;

   }
   finally
   { 

	   printf("Combining Peptide Results\n");

      for(int i=0; i < lTempPeptides->Count; i++)
	  {
		  for(int j=0; j < ListPeptides->Count; j++)
		  {
			if(lTempPeptides[i]->PeptideID->Equals(ListPeptides[j]->sPepID) )
			{
				ListPeptides[j]->Peptide = lTempPeptides[i]->sPeptide;

				ListPeptides[j]->dModMasses = gcnew List <double>;

				ListPeptides[j]->dModMasses = lTempPeptides[i]->dModMasses;

				ListPeptides[j]->ModLocations = gcnew List <int>;

				ListPeptides[j]->ModLocations = lTempPeptides[i]->ModLocations;

				/*printf("Pepitde = %s  %d  %d\n", ListPeptides[j]->Peptide, ListPeptides[j]->nCharge, ListPeptides[j]->nScan);

				printf("Evidences %s %s\n", lTempPeptides[i]->PeptideID, ListPeptides[j]->sPepID);*/

				break;
			}
		  }
	  }

	  printf("Finished Combining Peptide Results\n");
	   

	  delete(lTempPeptides);


	/*  for(int j=0; j < ListPeptides->Count; j++)
	  {
			{
				printf(" %s %s Scan = %d, Score = %5.1f, RetTime(s)= %5.1f %8.5f %s %s\n", ListPeptides[j]->Peptide,
					ListPeptides[j]->sPepID, ListPeptides[j]->nScan, ListPeptides[j]->dIonscore,
					ListPeptides[j]->dRetTime, ListPeptides[j]->dExpect,ListPeptides[j]->sProtId,
					ListPeptides[j]->sEvidence);
	  
			}
	  }*/

	  for(int i=0; i < ProteinResultList->Count; i++)
	  {
		  //printf("Proteins ::: %s\n", ProteinResultList[i]->accession);

		  

		  ProteinResultList[i]->PeptideList = gcnew List<PeptideEntry^>;

		  for(int k=0; k < ProteinResultList[i]->PeptideEvidence->Count; k++)
		  {
			  //if(k < 3)
				  //printf("Evidence Protein %s\n", ProteinResultList[i]->PeptideEvidence[k]);

			for(int j= 0; j < ListPeptides->Count; j++)
			{
				//if(j < 2)
					//printf("Evidence PeptideList = %s\n", ListPeptides[j]->sEvidence);
				//if the peptide evidences are the same add the peptides to the Protein
				if(ListPeptides[j]->sEvidence->Equals(ProteinResultList[i]->PeptideEvidence[k]) )
				{
					/*printf("PROTEIND = %s %s %s\n", ProteinResultList[i]->accession, 
						ProteinResultList[i]->PeptideEvidence[k],ListPeptides[j]->Peptide );*/

					ProteinResultList[i]->PeptideList->Add(ListPeptides[j]);

					//break;    //this is new. Remove this comment, if works smooth
				}
			}
		  }
	  }
	  
	  delete(ListPeptides);



      if ( reader != nullptr )
            reader->Close();
   }

}

/*
* a routine to read the scan number
* sScanLine is one of two possible items holding the scan number
*  mzMLid=controllerType=0 controllerNumber=1 scan=15671
*  controllerType=0 controllerNumber=1 scan=15671
* 
*/
int ReadScanNumber(String ^sScanLine)
{
	int i = 0, j = 0, nScanPosition, nScan = -1;
					
	char szScan[1024];

	szScan[0] = '\0';

	nScanPosition = sScanLine->LastIndexOf("scan=");

	j = 0;

	for(i = nScanPosition + 5; i < sScanLine->Length; i++)
	{
		if(sScanLine[i] == '&' || sScanLine[i] == '"')
		{
			break;
		}

		szScan[j] = sScanLine[i];

		j++;
	}

	szScan[j] = '\0';

	nScan = atoi(szScan);

	return nScan;
}
/*
* This is the original method to read mzIdentML files
* Neo Med mzid filse as of 01.09.2017 work with this
* routine only. They have a different format from the 
* mzid files that we generate locally. For locally generated files
*  use the ReadProteins_mzIdent_2()
* read proteins
* num_tot_proteins - total number of proteins
* TempProtCollection  temporaroly holds  the protein description,
*     length and other information, which is stored at the start of mzid
*
*/

void mzIdentML::ReadProteins_mzIdent()
{
	String ^ nvalue, ^sPeptideEvidence;

	//String ^ protein, ^peptide, ^sTemp, ^sTempProtID;
	String  ^protein, ^sTemp, ^sTempProtID;
	
	String ^rank1_peptide, ^rank1_sTemPortID;

	List <String ^> ^ListEvidence;

	int num_tot_proteins = 0, iRank = 0, iCharge, iTemp, nScan;

	int nSpectralCount= 0, nDistinctSeqs = 0, nModLocation = 0;
	
	double dProteinScore, dSeqCoverage = 0;

	double dSeqMass, dExpMass, dRetTime, dModMass= 0.0;

	Boolean bProteinAmb = false, bRankOne = false, bSpectumId = false, bSpecResult = false;

	bool bDBSequence = false, bPeptideEvidence = false;


	Boolean sportAddCheck = false; //Check if the db sequence can be added


	PeptideAndID ^sTempPeptideID = gcnew PeptideAndID();

	List <PeptideAndID ^> ^lTempPeptides = gcnew List <PeptideAndID ^>;

	ProteinDescrip ^tempProtDes = gcnew ProteinDescrip();

	PeptideEntry ^tempPeptide = gcnew PeptideEntry();

	List <PeptideEntry ^> ^ListPeptides = gcnew List <PeptideEntry ^>;

	System::Collections::Generic::List<ProteinDescrip^> ^ TempProtCollection = 
		gcnew System::Collections::Generic::List<ProteinDescrip^>;

	ProteinEntry ^currProtein;

	ProteinResultList = gcnew System::Collections::Generic::List<ProteinEntry^>;
	
	XmlTextReader^ reader = nullptr;
	
	try
	{
      reader = gcnew XmlTextReader(smzIdentML);

      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
		 {
         case XmlNodeType::Element:  

			//section to read protein description, accession, start of mzid

			if("DBSequence" == reader->Name)
			{
				bDBSequence = true;

				tempProtDes = gcnew ProteinDescrip();

				//tempProtDes->PeptideList = gcnew List<PeptideEntry ^>;

				nvalue = ReadNode(reader, "id");

				if(String::IsNullOrEmpty(nvalue))
				{
					printf("FOUND IT\n");
				}

				tempProtDes->sId = nvalue;

				nvalue = ReadNode(reader, "length");

				tempProtDes->nSeqLength = int::Parse(nvalue);

				nvalue = ReadNode(reader,"accession");

				tempProtDes->accession = nvalue;
			}

			if(bDBSequence)
			{
				if("cvParam" == reader->Name)
				{
					nvalue = ReadNode(reader,"accession");

					if("MS:1001088" == nvalue)   //for protein description
					{
						tempProtDes->sDescription = ReadNode(reader, "value");
						//store the protein infor
						//TempProtCollection->Add(tempProtDes);
					}
				}
			}
			//section to read protein description, accession

			//read and save peptides
			if("Peptide" == reader->Name)
			{
				nvalue = ReadNode(reader, "id");

				sTempPeptideID = gcnew PeptideAndID();

				sTempPeptideID->ModLocations = gcnew List <int>;

				sTempPeptideID->dModMasses  = gcnew List <double>;
				
				sTempPeptideID->PeptideID = nvalue;

				//printf("Temp %s\n", sTempPeptideID->PeptideID);
			}
			else if("peptideSequence" == reader->Name || "PeptideSequence" == reader->Name)
			{
				
				nvalue = reader->ReadString();

				sTempPeptideID->sPeptide = nvalue;

				//printf("PeptideSequence and nvalue %s\n",nvalue);
				//lTempPeptides->Add(sTempPeptideID);    //MAY NEED TO REINSTATE THIS
			}
			
			else if("Modification" == reader->Name)
			{
				
				//printf("Modification\n");

				nvalue = ReadNode(reader, "location");

				nModLocation = int::Parse(nvalue);

				//printf("nModLocation %d\n",nModLocation);

				sTempPeptideID->ModLocations->Add(nModLocation);

				nvalue = ReadNode(reader, "monoisotopicMassDelta");

				dModMass = float::Parse(nvalue);

				//printf("monoisotopicMassDelta %f\n",dModMass);

				sTempPeptideID->dModMasses->Add(dModMass);

				//printf("ModLoc = %d %s\n", int::Parse(nvalue), sTempPeptideID->sPeptide);

			}

			//end of peptide reading
			
			//if("PeptideEvidence" == reader->Name) //Modified by Mahbubur Rahman to capture peptide evidence and sportid; date: 10/06/2015 
			//{

			//	nvalue = ReadNode(reader, "id");

			//	sPeptideEvidence = nvalue;

			//	//nvalue = ReadNode(reader, "DBSequence_Ref");

			//	////tempPeptide->sProtId = nvalue;
			//	//sTempProtID = nvalue;

			//	//tempPeptide->sEvidence = sPeptideEvidence;
			//	//ListPeptides->Add(tempPeptide); //Adding peptide

			//	//ListPeptides[ListPeptides->Count-1]->sProtId = sTempProtID; //adding db sequence


			//}



			//section to read peptide entries;
            if("SpectrumIdentificationResult" == reader->Name)
			{
				bSpecResult = true;
			}

			if("SpectrumIdentificationItem" == reader->Name)
			{
				bSpectumId = true;
				
				tempPeptide = gcnew PeptideEntry();

				nvalue = ReadNode(reader, "calculatedMassToCharge");

				dSeqMass = float::Parse(nvalue);

				nvalue = ReadNode(reader, "chargeState");

				iCharge = int::Parse(nvalue);

				nvalue = ReadNode(reader, "experimentalMassToCharge");
				
				dExpMass = float::Parse(nvalue);

				nvalue = ReadNode(reader, "Peptide_ref");

				tempPeptide->sPepID = nvalue;

				nvalue = ReadNode(reader, "rank");

				iRank = int::Parse(nvalue);

				if(1 == iRank)
				{
					tempPeptide->nScan   = nScan;

					tempPeptide->SeqMass = dSeqMass;

					tempPeptide->nCharge = iCharge;

					tempPeptide->SpecMass = dExpMass;

					tempPeptide->nRank = iRank;

					bRankOne = true;


					//////////////modification 10/09/2015////////////////
					
					
				    sPeptideEvidence = gcnew String("NONE");

					sTempProtID = gcnew String("NONE");      //some times the is no evidence reported for a peptide in mzid    change : 10/06/2015
				}
				else
				{
					bRankOne = false;
				}
			}

			if(bSpectumId && "cvParam" == reader->Name && bRankOne)
			{
				nvalue = ReadNode(reader,"accession");
				
				if("MS:1001171" == nvalue)   //scores, Mascot
				{
					tempPeptide->dIonscore = float::Parse(ReadNode(reader, "value"));
				}
				else if("MS:1001172" == nvalue)
				{
					tempPeptide->dExpect = double::Parse(ReadNode(reader, "value"));
				}

				//printf("MASCOT SCORES %10.5f\n", dIonScore);

				bPeptideEvidence = false;   //read about bPeptideEvidence
			}

			
			if(bSpecResult && "cvParam" == reader->Name)
			{
				nvalue = ReadNode(reader,"accession");

				if("MS:1001114" == nvalue || "MS:1000894" == nvalue)   //retention time <cvParam accession="MS:1001114" name="retention time(s)"  cvRef="PSI-MS" value="655.774" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" 
				{
					dRetTime  = float::Parse(ReadNode(reader, "value"));
				}
				//else if("MS:1000797" == nvalue)   //cvParam accession="MS:1000797" name="peak list scans"  cvRef="PSI-MS" value="1194" 
				else if("MS:1000796" == nvalue)     //name="spectrum title"  cvRef="PSI-MS" value="15apr2604.21098.21098.5 File:&quot;15apr2604.raw&quot;, 
                {                                    //NativeID:&quot;controllerType=0 controllerNumber=1 scan=21098&quot;" 
					//nvalue = ReadNode(reader, "value");
					String^ stemp;

					stemp = ReadNode(reader, "value");

					nScan = ReadScanNumber(stemp);
				}
			}
			//read the peptide evidence information;

			// <PeptideEvidence id="PE_2_4_LV106_HUMAN_0_47_51" start="47" end="51" pre="K" post="D" missedCleavages="0" 
			// added bPeptideEvidence, so that to read the evidence only once. It happens that in some proteins, the same
			// sequence is encountered multiple times, and then its peptide evidence is multiplexed as well, e.g.,
			//  <PeptideEvidence id="PE_481_1_SPA_STAA8_0_62_69" start="62" end="69" pre="R" post="D" missedCleavages="0" isDecoy="false" DBSequence_Ref="DBSeq_1_SPA_STAA8" />
			//  <PeptideEvidence id="PE_481_1_SPA_STAA8_0_123_130" start="123" end="130" pre="R" post="D" missedCleavages="0" isDecoy="false" DBSequence_Ref="DBSeq_1_SPA_STAA8" />
						
			
			if(bRankOne && "PeptideEvidence" == reader->Name && false == bPeptideEvidence)   //in Mascot 2.5.1 outputs do not reach here date:02.18.2016
			{
				nvalue = ReadNode(reader, "id");

				sPeptideEvidence = nvalue;

				nvalue = ReadNode(reader, "DBSequence_Ref");

				//tempPeptide->sProtId = nvalue;
				sTempProtID = nvalue;

				bPeptideEvidence = true;
			}

			if(bRankOne && "PeptideEvidenceRef" == reader->Name && false == bPeptideEvidence) //Mahbubur Rahman
			{

				nvalue = ReadNode(reader, "peptideEvidence_ref");

				//printf("peptideEvidence_ref %s\n", nvalue);

				sPeptideEvidence = nvalue;

				rank1_peptide = nvalue;

				//nvalue = ReadNode(reader, "DBSequence_Ref");

				//tempPeptide->sProtId = nvalue;
				//sTempProtID = nvalue;

				bPeptideEvidence = true;
				sportAddCheck = true;
			}
			//end of section for peptide entries

			//read protein information and scores
			if("ProteinAmbiguityGroup" == reader->Name)
			{

				nvalue = ReadNode(reader, "id");

				nSpectralCount = 0;

				bProteinAmb = true;

				ListEvidence = gcnew List <String ^>;
			}

			//read the protein ambiguity node
			if(bProteinAmb)
			{	
				if("cvParam" == reader->Name)
				{
					nvalue=ReadNode(reader,"accession");
				
					if (nvalue=="MS:1001171") // protein score
					{
						dProteinScore = double::Parse(ReadNode(reader, "value"));
					}
					else if (nvalue == "MS:1001093")    //sequence coverage
					{
						dSeqCoverage = double::Parse(ReadNode(reader, "value"));
					}
					else if (nvalue == "MS:1001097")    //number of distinct sequences
					{
						nDistinctSeqs = int::Parse(ReadNode(reader, "value"));
					}
					
				} // end of if("cvparams" == reader->name ...)
				else if("PeptideHypothesis" == reader->Name ) 
				{
					nvalue = ReadNode(reader, "PeptideEvidence_Ref");

					if(rank1_peptide == nvalue)
					{
						sTempProtID = rank1_sTemPortID; //Adding port id 10-12-2015

						/*if(sportAddCheck==true) //10/12/2015
						{ */ 
				            ListPeptides[ListPeptides->Count-1]->sProtId = sTempProtID;  
							//sportAddCheck = false;
						//}

					}

					ListEvidence->Add(nvalue);

					nSpectralCount++;
				}
				else if("ProteinDetectionHypothesis" == reader->Name)
				{
					protein = ReadNode(reader,"id");
					
					rank1_sTemPortID = ReadNode(reader, "dBSequence_ref"); //keeping Db_sequence
						//sTempProtID = nvalue;
						//sportAddCheck = false;
					
					//printf("rank1_sTemp = %s\n", rank1_sTemPortID);	      
				}
			}// if(bProteinAmb)

			break;


		 case XmlNodeType::EndElement:


			 if("ProteinDetectionHypothesis" == reader->Name)
			 {
				//date:02.18.2016 Commented out RGS. Because some protein ambiguity groups
				 //contain several proteins. In that case, the option only reads
				 //the last protein of several contained in ProteinAmbiguity
				if(nSpectralCount > 0)                           
				{
					currProtein = gcnew ProteinEntry();

					currProtein->nSpectralCount = nSpectralCount;

					currProtein->accession = protein;

					currProtein->ProteinScore = dProteinScore;

					currProtein->SeqCoverage = dSeqCoverage;

					currProtein->nDistinctSequences = nDistinctSeqs;

					currProtein->PeptideEvidence = ListEvidence;

					ProteinResultList->Add(currProtein);

                    /*printf("Protein %s SC %d ProtScore %10.2f Cov = %10.1f %d\n", 
							protein, nSpectralCount, dProteinScore, dSeqCoverage, num_tot_proteins);*/

					num_tot_proteins++;
				}
			 }
			 else if("ProteinAmbiguityGroup" == reader->Name)
			 {
				 //ReadNode(reader, "id");

				 //printf("Protein AmbigID %s\n", reader->Name);

				bProteinAmb = false;
				

				/*if(nSpectralCount > 0)      //date:02.18.2016 Commented out RGS. Because some protein ambiguity groups
				{									//contain several proteins. In that case, the option only reads
					currProtein = gcnew ProteinEntry();                //the last protein of several contained in ProteinAmbiguity

					currProtein->nSpectralCount = nSpectralCount;

					currProtein->accession = protein;

					currProtein->ProteinScore = dProteinScore;

					currProtein->SeqCoverage = dSeqCoverage;

					currProtein->nDistinctSequences = nDistinctSeqs;

					currProtein->PeptideEvidence = ListEvidence;

					ProteinResultList->Add(currProtein);

                    printf("Protein %s SC %d ProtScore %10.2f Cov = %10.1f %d\n", 
							protein, nSpectralCount, dProteinScore, dSeqCoverage, num_tot_proteins);

					num_tot_proteins++;
				}*/
			 }
			 else if("SpectrumIdentificationResult" == reader->Name)
			 {
				 bSpecResult = false;
                //update the ret time, as it is per spectrum Identification Result
                if(ListPeptides->Count > 0) //added on 05/16/2016
				{
					ListPeptides[ListPeptides->Count-1]->dRetTime = dRetTime;

					ListPeptides[ListPeptides->Count-1]->sProtId = sTempProtID;  
				
					ListPeptides[ListPeptides->Count-1]->nScan   = nScan;
				}
		
				//delete (tempPeptide);
			 }
			 else if("SpectrumIdentificationItem" == reader->Name && bRankOne)
			 {
				tempPeptide->sEvidence = sPeptideEvidence;

				ListPeptides->Add(tempPeptide);

				delete (tempPeptide);

				bRankOne = false;

				bSpectumId = false;
			 }
			 else if("DBSequence" == reader->Name)
			 {
				 bDBSequence = false;

				 //store the protein infor
				TempProtCollection->Add(tempProtDes);
			 }
			 else if("Peptide" == reader->Name)    //finished reading the peptide information 
			 {                                    // add it to the temporary list
				 //printf("End of peptide\n");
				 lTempPeptides->Add(sTempPeptideID);
			 }

            break; 
		 } 

} 
	   while (reader->Read());    

	   if(ProteinResultList->Count != TempProtCollection->Count)
	   {
			printf("Error in the mzid file, %s, different number of proteins %d %d\n", smzIdentML,
				ProteinResultList->Count, TempProtCollection->Count);

			puts("Exitting ...\n");

			exit (1);
	   }

	   printf("Number of Proteins in %s is %d\n", smzIdentML, num_tot_proteins);

	   //combine results from temporal protein holder and protein depository
	   // some results are in protein holder and some are in protein depository
	   // that is why it is needed to combine them

	   printf("Combining Protein Results\n");

	   for(int i=0; i < TempProtCollection->Count; i++)
	   {
			if(ProteinResultList[i]->accession->IndexOf(TempProtCollection[i]->accession) == -1)
			{
				printf("Problem for %s %s\n", ProteinResultList[i]->accession,
						TempProtCollection[i]->accession);

				exit(1);
			}
			else
			{
				ProteinResultList[i]->accession = TempProtCollection[i]->accession;

				ProteinResultList[i]->description = TempProtCollection[i]->sDescription;

				ProteinResultList[i]->nSeqLength = TempProtCollection[i]->nSeqLength;
			}
	   }

	   printf("Finished Combining Protein Results\n");

	   delete (TempProtCollection);

return ;

   }
   finally
   { 

	   printf("Combining Peptide Results\n");

      for(int i=0; i < lTempPeptides->Count; i++)
	  {
		  for(int j=0; j < ListPeptides->Count; j++)
		  {
			if(lTempPeptides[i]->PeptideID->Equals(ListPeptides[j]->sPepID) )
			{
				ListPeptides[j]->Peptide = lTempPeptides[i]->sPeptide;

				ListPeptides[j]->dModMasses = gcnew List <double>;

				ListPeptides[j]->dModMasses = lTempPeptides[i]->dModMasses;

				ListPeptides[j]->ModLocations = gcnew List <int>;

				ListPeptides[j]->ModLocations = lTempPeptides[i]->ModLocations;

				break;
			}
		  }
	  }

	  printf("Finished Combining Peptide Results\n");
	   

	  delete(lTempPeptides);


	/*  for(int j=0; j < ListPeptides->Count; j++)
	  {
			{
				printf(" %s %s Scan = %d, Score = %5.1f, RetTime(s)= %5.1f %8.5f %s %s\n", ListPeptides[j]->Peptide,
					ListPeptides[j]->sPepID, ListPeptides[j]->nScan, ListPeptides[j]->dIonscore,
					ListPeptides[j]->dRetTime, ListPeptides[j]->dExpect,ListPeptides[j]->sProtId,
					ListPeptides[j]->sEvidence);
	  
			}
	  }*/

	  for(int i=0; i < ProteinResultList->Count; i++)
	  {
		  //printf("Proteins ::: %s\n", ProteinResultList[i]->accession);

		  

		  ProteinResultList[i]->PeptideList = gcnew List<PeptideEntry^>;

		  for(int k=0; k < ProteinResultList[i]->PeptideEvidence->Count; k++)
		  {
			  //if(k < 3)
				  //printf("Evidence Protein %s\n", ProteinResultList[i]->PeptideEvidence[k]);

			for(int j= 0; j < ListPeptides->Count; j++)
			{
				//if(j < 2)
					//printf("Evidence PeptideList = %s\n", ListPeptides[j]->sEvidence);
				//if the peptide evidences are the same add the peptides to the Protein
				if(ListPeptides[j]->sEvidence->Equals(ProteinResultList[i]->PeptideEvidence[k]) )
				{
					/*printf("PROTEIND = %s %s %s\n", ProteinResultList[i]->accession, 
						ProteinResultList[i]->PeptideEvidence[k],ListPeptides[j]->Peptide );*/

					ProteinResultList[i]->PeptideList->Add(ListPeptides[j]);

					//break;    //this is new. Remove this comment, if works smooth
				}
			}
		  }
	  }
	  
	  delete(ListPeptides);



      if ( reader != nullptr )
            reader->Close();
   }

}