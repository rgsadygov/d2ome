// XmlMzML.h

/*
*	A program to read spectra from mzML files
*	scan_offset_array is an array that holds byte off-sets
*	of spectra in the mzML (by default mzML creates the indices)
*	The first line of mzML file indicates if it has been created with
*	default settings (create indices). The off-sets themselves are written
*	at the end of the mzML file.
*	A FileStream object,  inStream, is used because it is not possible to 
*	re-wind TextXMLReader. But using the stream allows one to manupilate the
*	position of the file off-set.
*
*
*   Rovshan Sadygov, Ph.D.
*   Department of Biochemistry and Molecular Biology
*   The University of Texas Medical Branch
*   Galveston, TX
*/

#pragma once

using namespace System;
using namespace System::IO;
using namespace System::Xml;
using namespace System::Collections::Generic;

namespace XmlMzML {


	public ref class FullScan
	{
	public:
		array<float, 1> ^ Intensity;

		array <double, 1> ^ moverz;

		float RetTime;

		long ScanNumber;
	};

	//public ref class FullScanChromatogram
	//{
	//public:
	//	List <int> ^ ScanList;

	//	List <FullScan ^> ^ FullScans;
	//};

	public ref class MzML
	{
// TODO: Add your methods for this class here.

		public:

			array<double, 2> ^ Spectrum, ^ TempSpectrum;

			System::Xml::XmlTextReader^ reader;

			array<long long> ^ scan_offset_array;

			long nTotalSpectra;

			String ^ sFilemzml;

			bool bOrbiTrap;

			float fCurrentRetTime;

			Boolean indexedmzml;

			//array<double, 2> ^ TempSpectrum;

			FileStream ^ inStream;

			List <FullScan ^> ^FullScanChromatogram;

			MzML(String ^ sInputFile)
			{

				bOrbiTrap = false;

				indexedmzml = false;

				sFilemzml = sInputFile;

				inStream = File::OpenRead(sFilemzml);

				reader = 
					gcnew System::Xml::XmlTextReader(sFilemzml);

				reader->WhitespaceHandling = 
					System::Xml::WhitespaceHandling::None;

				//read the general information pertinent to the 
					//whole raw file

				while (reader->Read())
				{
					if (reader->Name == "indexedmzML") 
					{
						indexedmzml=true;
					}

					//there are two if blocks to determine if an instrument is Ortibtrap
					// this is because PSI has a definition for Orbtrap but not for
					// Orbitrap Velos at this time. 07.08.2010
					if("cvParam" == reader->Name)
					{
						while (reader->MoveToNextAttribute()) 
						{
							if (reader->Name=="accession") 
							{
								//Orbitrap or FT
								if(reader->Value == "MS:1000449" ||
									reader->Value == "MS:1000448")
								{
									bOrbiTrap = true;

									break;
								}
							}
						}
					}

					if("userParam" == reader->Name)
					{
						while (reader->MoveToNextAttribute()) 
						{
							if ("value" == reader->Name) 
							{
								if(reader->Value->IndexOf("LTQ Orbitrap") != -1)
								{
									bOrbiTrap = true;

									break;
								}
							}
						}
					}

					if("spectrumList" == reader->Name )
					{
						while (reader->MoveToNextAttribute()) 
						{
						  if (reader->Name=="count") 
							{
								nTotalSpectra=long::Parse(reader->Value);

								scan_offset_array = gcnew array<long long> (nTotalSpectra+1);

								break;
							}
						}

						break;
					}
				}
			}

			int ReadZoomScanSequentially(long nScan);

			int ReadAScanSequentially(long nScan);

			int ReadFullScan(long nScan);

			int ReadAllFullScans();

			int ReadFullScan(long nScan, double dPrecursorMass, double dInterval);

			int RetrieveMass(double dPrecursorMass, double dInterval);

			int mzXMLScan(long int nScan);
			
			int ReadIndexedFullScan(long nScan, unsigned short *nMSLevel);

			int ReadIndexedScan(long nScan, unsigned short *nMSLevel);

			int ReadSpectrumByteOffset();

			int ReadAScanNoIndex(long nScan, XmlTextReader ^ xmzML);

			int PrecursorMassAndCharge();

			void IndexedChromatogramBuild();

			void SequentialChromatogramBuild();

			bool CloseFiles();

			Boolean ReadOffset();
	};

		//a class to read mzXML files.
	// tested this only on the MALDI files
	public ref class MzXML
	{
		public:

			array<double, 2> ^ Spectrum;

			System::Xml::XmlTextReader^ reader;

			long nTotalSpectra;


			String ^ sFilemzml;

			MzXML(String ^ sInputFile)
			{
				sFilemzml = sInputFile;

					reader = 
						gcnew System::Xml::XmlTextReader(sFilemzml);
			}

			int ReadmzXML();

			int ReadScan(long nScan);
	};

}


