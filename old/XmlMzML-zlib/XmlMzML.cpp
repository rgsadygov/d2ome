// This is the main DLL file.
/*
*   Rovshan Sadygov, Ph.D.
*   Department of Biochemistry and Molecular Biology
*   The University of Texas Medical Branch
*   Galveston, TX
*/

#include "stdafx.h"
#include <stdio.h>
#include "XmlMzML.h"
#include "stdlib.h"
#include "string.h"

#using <System.Xml.dll>

using namespace System;
using namespace System::IO;
using namespace System::Xml;
using namespace XmlMzML;
using namespace System::Collections::Generic;




String^ ReadNode(XmlTextReader^ reader,String^ nodeName)
{
	String^ nodeValue;

	while (reader->MoveToNextAttribute()) 
	{
	  if(reader->Name == nodeName) 
		{ 
			nodeValue = reader->Value;

			break;
		}
	}

	return nodeValue; 
}

int Base64 (XmlTextReader^ reader, int nSpect, int nb, array<double,2> ^numbers)
{ 
		int base64len = 0;
		
		array<Byte>^base64 = gcnew array<Byte>(nSpect*8+8);

		base64len = reader->ReadBase64( base64, 0, nSpect*nb );

		if ( nb==8 ) 
		{
			for (int k=0; k<nSpect*8; k=k+8)
			{
				numbers[0,k/8] = BitConverter::ToDouble(base64, k );

				//printf("FIRST VALUES HERE %10.5f\n", numbers[0,k/8]);
			}
		}
		else
		{
			for (int k=0; k<nSpect*4; k=k+4)
			{
				numbers[1,k/4] = BitConverter::ToSingle(base64, k );

				//printf("VALUES HERE %10.5f\n", numbers[1,k/4]);
			}
		}

	return 0;
}

/*
* converts the m/z and intensity values
*/

int Base64_Unified (XmlTextReader^ reader, int nSpect, int nb, array<double,2> ^numbers, bool isMZ, bool izIntens)
{ 
		int base64len = 0;
		
		array<Byte>^base64 = gcnew array<Byte>(nSpect*8+8);

		base64len = reader->ReadBase64( base64, 0, nSpect*nb );

		if(isMZ && !izIntens)
		{
			for (int k = 0; k < nSpect * nb; k = k + nb)
			{
				numbers[0, k/nb] = BitConverter::ToDouble(base64, k );

				//printf("m/z %10.5f\n", numbers[0, k/nb]);
			}
		}
		else if(!isMZ && izIntens)
		{
			for (int k=0; k < nSpect*nb; k = k + nb)
			{
				if(nb == 8)
				{
					numbers[1, k/nb] = BitConverter::ToDouble(base64, k );
				}
				else if(nb == 4)
				{
					numbers[1, k/nb] = BitConverter::ToSingle(base64, k );
				}
				else
				{
					printf("Error: Did not recognize the encoding: %d\n", nb);

					exit (1);
				}

				//printf("Abund %10.5f\n", numbers[1,k/nb]);
			}
		}
		else
		{
			printf("Error: Wrong input into Base64\n");

			exit(1);
		}


	return 0;
}


//*************

// a method returns the lenght of the spectrum
// (scan number) nScan

int MzML::ReadZoomScanSequentially(long nScan)

{
int nbit;
int i=0; 
int j=0, nEncodedLength = 0;
int nSpectLength = 0; 
long scanNumber, maxscanNumber;
Boolean isZoom=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
String^ nvalue;

	XmlTextReader^ reader = nullptr;
//   Console::WriteLine(L"ReadScan XmlTextReader...");
	try
	{
      reader = gcnew XmlTextReader(sFilemzml);

      reader->WhitespaceHandling =  WhitespaceHandling::None;

		do {
         switch (reader->NodeType)
         {
         case XmlNodeType::Element:

			if("spectrumList" == reader->Name )
			{
				nvalue=ReadNode(reader,"count");
				maxscanNumber=long::Parse(nvalue);
 			}
            
			if("spectrum" == reader->Name )
			{
				isZoom=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);
				if (scanNumber==nScan)
				{
					spectrumNumberOK=true;
					Spectrum = gcnew array<double, 2> (2,nSpectLength+1);

				}
	
			}
            


			if("cvParam" == reader->Name ) 
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1000497")
				{
					isZoom=true;
				}


				if (nvalue=="MS:1000523" && isZoom && spectrumNumberOK)
				{
					set64bit=true;

					nbit=8;
				}

				if (nvalue=="MS:1000521" && isZoom && spectrumNumberOK)
				{
					set32bit=true;

					nbit=4;
				}
				if (nvalue=="MS:1000574"  && isZoom && spectrumNumberOK)
				{
					zlib_compression=true;		
					printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
					exit(EXIT_FAILURE);
				}			

				if (nvalue=="MS:1000514" && isZoom && spectrumNumberOK)
				{
					ismoz=true;
				}
				if (nvalue=="MS:1000515" && isZoom && spectrumNumberOK)
				{
					isIntensity=true;
				}

 			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}

			if("binary" == reader->Name )
			{
				if (isZoom && spectrumNumberOK) 
				{	
						i = 0;
					if (ismoz)
					{
							Base64 (reader,nSpectLength,nbit,Spectrum);

						i++;
						ismoz=false;
						ismozOK=true;  		
					}

					if (isIntensity)
					{
						j = 0;

							Base64   (reader, nSpectLength, nbit, Spectrum);
						j++;
						isIntensity=false;

				if (ismozOK) return nSpectLength;
					}	
				}
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... ");
//           Console::Write(reader->Value);
            break;
         case XmlNodeType::EndElement:
//			 Console::Write(" endelement... ");
//           Console::Write("{0}", reader->Name);
            break;
		 } 

} 
	   while (reader->Read());    

return 0;

   }
   finally
   { 
	   if (!spectrumNumberOK)
		{	
			Console::WriteLine("");
			Console::Write("There is no zoom scan with number ");
			Console::Write(nScan); Console::WriteLine(" !!!");
		}
//	  Console::WriteLine("nSpectLength={0}",nSpectLength);
      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilemzml);
      if ( reader != nullptr )
            reader->Close();
   }
	return nSpectLength;
}

//*************

// a method for sequential reading of scan.
// the scan number is as argument.
// returns the lenght of the spectrum
// (scan number) nScan

int MzML::ReadAScanSequentially(long nScan)

{
int nbit;
int i=0; 
int j=0, nEncodedLength = 0;
int nSpectLength = 0; 
long scanNumber, maxscanNumber;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
String^ nvalue;

	XmlTextReader^ reader = nullptr;
//   Console::WriteLine(L"ReadScan XmlTextReader...");
	try
	{
      reader = gcnew XmlTextReader(sFilemzml);

      reader->WhitespaceHandling =  WhitespaceHandling::None;

		do {
         switch (reader->NodeType)
         {
         case XmlNodeType::Element:

			if("spectrumList" == reader->Name )
			{
				nvalue=ReadNode(reader,"count");
				maxscanNumber=long::Parse(nvalue);
 			}
            
			if("spectrum" == reader->Name )
			{
				spectrumNumberOK=false;
				
				nvalue=ReadNode(reader,"index");
				
				scanNumber=long::Parse(nvalue) + 1;
				
				reader->MoveToAttribute(0);
				

				nvalue=ReadNode(reader,"defaultArrayLength");
				
				nSpectLength=int::Parse(nvalue);

				if (scanNumber==nScan)
				{
					spectrumNumberOK=true;

					Spectrum = gcnew array<double, 2> (2,nSpectLength+1);

					
					if(0 == nSpectLength)
					{
						printf("Empty Scan\n");

						return 0;
					}

				}
	
			}
            


			if("cvParam" == reader->Name ) 
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1000523" && spectrumNumberOK)
				{
					set64bit=true;

					nbit=8;
				}

				if (nvalue=="MS:1000521" && spectrumNumberOK)
				{
					set32bit=true;

					nbit=4;
				}
				if (nvalue=="MS:1000574"   && spectrumNumberOK)
				{
					zlib_compression=true;	
					printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
					exit(EXIT_FAILURE);
				}			

				if (nvalue=="MS:1000514"  && spectrumNumberOK)
				{
					ismoz=true;
				}
				if (nvalue=="MS:1000515" && spectrumNumberOK)
				{
					isIntensity=true;
				}

 			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}

			if("binary" == reader->Name )
			{
				if (spectrumNumberOK) 
				{	
						i = 0;
					if (ismoz)
					{

							Base64 (reader,nSpectLength,nbit,Spectrum);

						i++;
						ismoz=false;
						ismozOK=true;  		
					}

					if (isIntensity)
					{
						j = 0;

							Base64   (reader, nSpectLength, nbit, Spectrum);
						j++;
						isIntensity=false;

				if (ismozOK)
				{
	
					TempSpectrum = gcnew array<double, 2> (2, nSpectLength+1);

					for(int ii=0; ii < nSpectLength +1; ii++)
					{
						TempSpectrum[0,ii] = Spectrum[0,ii];
						
						TempSpectrum[1,ii] = Spectrum[1,ii];
					}

					return nSpectLength;
				}
					}	
				}
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... ");
//           Console::Write(reader->Value);
            break;
         case XmlNodeType::EndElement:
//			 Console::Write(" endelement... ");
//           Console::Write("{0}", reader->Name);
            break;
		 } 

} 
	   while (reader->Read());    

return 0;

   }
   finally
   { 
	   if (!spectrumNumberOK)
		{	
			Console::WriteLine("");
			Console::Write("There is no Scan with number ");
			Console::Write(nScan); Console::WriteLine(" !!!");
		}
//	  Console::WriteLine("nSpectLength={0}",nSpectLength);
      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilemzml);
      if ( reader != nullptr )
            reader->Close();
   }
	return nSpectLength;
}


int MzML::ReadAScanNoIndex(long nScan, XmlTextReader ^ xmzML)
{
	int nbit;
	int i=0; 
	int j=0, nEncodedLength = 0;
	int nSpectLength = 0; 
	long scanNumber, maxscanNumber;
	Boolean set32bit=false;
	Boolean set64bit=false;
	Boolean spectrumList=false;
	Boolean spectrumNumberOK=false;
	Boolean ismoz=false;
	Boolean ismozOK=false;
	Boolean isIntensity=false;
	Boolean zlib_compression=false;
	String^ nvalue; 
	

		try
	{
      //reader = gcnew XmlTextReader(sFilemzml);

      //reader->WhitespaceHandling =  WhitespaceHandling::None;

		do {

			/*printf("BEFORE NODE %d %d %s %s %d\n", reader->NodeType, XmlNodeType::Element,
				reader->Name->ToString(), xmzML->Name->ToString(),
				nScan);*/

         switch (xmzML->NodeType)
         {
         case XmlNodeType::Element:

			if("spectrumList" == xmzML->Name )
			{
				nvalue=ReadNode(xmzML,"count");
				maxscanNumber=long::Parse(nvalue);
 			}
            
			if("spectrum" == xmzML->Name )
			{
				spectrumNumberOK=false;
				
				nvalue=ReadNode(xmzML,"index");
				
				scanNumber=long::Parse(nvalue) + 1;
				
				xmzML->MoveToAttribute(0);
				
				nvalue=ReadNode(xmzML,"defaultArrayLength");
				
				nSpectLength=int::Parse(nvalue);

				if (scanNumber==nScan)
				{
					spectrumNumberOK=true;

					Spectrum = gcnew array<double, 2> (2,nSpectLength+1);

					
					if(0 == nSpectLength)
					{
						printf("Empty Scan\n");

						return 0;
					}

				}
	
			}
            

			if("cvParam" == xmzML->Name ) 
			{
				nvalue=ReadNode(xmzML,"accession");

				if (nvalue=="MS:1000523" && spectrumNumberOK)
				{
					set64bit=true;

					nbit=8;
				}

				if (nvalue=="MS:1000521" && spectrumNumberOK)
				{
					set32bit=true;

					nbit=4;
				}
				if (nvalue=="MS:1000574"   && spectrumNumberOK)
				{
					zlib_compression=true;		
					printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
					exit(EXIT_FAILURE);
				}			

				if (nvalue=="MS:1000514"  && spectrumNumberOK)
				{
					ismoz=true;
				}
				if (nvalue=="MS:1000515" && spectrumNumberOK)
				{
					isIntensity=true;
				}

 			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == xmzML->Name)
			{
				nvalue = ReadNode(xmzML,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}

			if("binary" == xmzML->Name )
			{
				if (spectrumNumberOK) 
				{	
						i = 0;
					if (ismoz)
					{

							Base64 (xmzML,nSpectLength,nbit,Spectrum);

						i++;
						ismoz=false;
						ismozOK=true;  		
					}

					if (isIntensity)
					{
						j = 0;

							Base64   (xmzML, nSpectLength, nbit, Spectrum);
						j++;
						isIntensity=false;

				if (ismozOK)
				{
	
					TempSpectrum = gcnew array<double, 2> (2, nSpectLength+1);

					for(int ii=0; ii < nSpectLength +1; ii++)
					{
						TempSpectrum[0,ii] = Spectrum[0,ii];
						
						TempSpectrum[1,ii] = Spectrum[1,ii];
					}

					return nSpectLength;
				}
					}	
				}
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... ");
//           Console::Write(reader->Value);
            break;
         case XmlNodeType::EndElement:
//			 Console::Write(" endelement... ");
//           Console::Write("{0}", reader->Name);
            break;
		 } 

} 
	while (xmzML->Read());   

return 0;

   }
   finally
   { 
	   if (!spectrumNumberOK)
		{	
			Console::WriteLine("");
			Console::Write("There is no Scan with number ");
			Console::Write(nScan); Console::WriteLine(" !!!");
		}
//	  Console::WriteLine("nSpectLength={0}",nSpectLength);
      //Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilemzml);
      /*if ( xmzML != nullptr )
            xmzML->Close();*/
   }

	return nSpectLength;
}


//*************

// a method that reads the last Full Scan
// that comes before the scan number nScan
//i.e., if the nScan is 4056 and it is MS/MS
// then the read scan will be the last full
// scan before 4056. If 4056 is a full scan
// then it will be returned
// cv params type 1000579 is the full scan

int MzML::ReadFullScan(long nScan)
{
int nbitmz, nbition, nFullScan, nFullLength;
int i=0; 
int j=0;
int nSpectLength = 0, nEncodedLength; 
long scanNumber;
Boolean isFull=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false; 
String^ nvalue;

nFullLength = 0;

	do {
		if(reader->NodeType == XmlNodeType::Element)
		{
			if("spectrum" == reader->Name )
			{
				isFull=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);

				if (scanNumber == nScan - 1)
				{
					spectrumNumberOK = true;
				}
	
			}

			if((nScan - scanNumber) <= 10 && nScan >= scanNumber) 
			{


				if ("cvParam" == reader->Name)
				{
					nvalue=ReadNode(reader,"accession");

					if (nvalue=="MS:1000579")
					{
						isFull=true;

						Spectrum = gcnew array<double, 2> (2,nSpectLength);

						nFullScan = scanNumber;

						nFullLength = nSpectLength;
					}

					if(isFull)
					{
						if (nvalue=="MS:1000523" )
						{
							set64bit=true;

							nbitmz=8;
						}

						if (nvalue=="MS:1000521")
						{
							set32bit=true;

							nbition=4;
						}

						if (nvalue=="MS:1000574")
						{
							zlib_compression=true;
							printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
							exit(EXIT_FAILURE);
						}

						if (nvalue=="MS:1000514")
						{
							ismoz=true;

							isIntensity = false;
						}
						
						if (nvalue=="MS:1000515")
						{
							isIntensity=true;

							ismoz = false;
						}
					}

				}

 			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}


			if("binary" == reader->Name && isFull)
			{
				if (ismoz)
				{

						Base64 (reader,nSpectLength,nbitmz,Spectrum);
					

					ismoz=false;

					ismozOK=true;
				}

				if (isIntensity)
				{

						Base64   (reader, nSpectLength, nbition, Spectrum);
					

					isIntensity=false;

					isFull = false;


					if(spectrumNumberOK)
					{
						printf("Full Scan Number %d\n", (nFullScan+1));

						return nFullLength;
					}

				}	

			}

			if(scanNumber > nScan-1)
			{
				printf("Full Scan Number %d\n", (nFullScan+1));

				return nFullLength;
			}
		 } 

} 
	while (reader->Read());


	if (!spectrumNumberOK)
	{	
		Console::WriteLine("");
		Console::Write("There is no Full with number ");
		Console::Write(nScan); Console::WriteLine(" !!!");

		return 0;
	}

	return nSpectLength;

}






//*************

// a method that reads All Full Scans
// that comes before the scan number nScan
//i.e., if the nScan is 4056 and it is MS/MS
// then the read scan will be the last full
// scan before 4056. If 4056 is a full scan
// then it will be returned
// cv params type 1000579 is the full scan

int MzML::ReadAllFullScans()
{
int nbitmz, nbition;
int i=0; 
int j=0;
int nSpectLength = 0, nEncodedLength = 0; 
long scanNumber;
Boolean isFull=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
String^ nvalue;


	do {
		if(reader->NodeType == XmlNodeType::Element)
		{
			if("spectrum" == reader->Name )
			{
				isFull=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);
	
			}


			if ("cvParam" == reader->Name)
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1000579")
				{
					isFull=true;

					Spectrum = gcnew array<double, 2> (2,nSpectLength);
				}

				if(isFull)
				{
				
					if (nvalue=="MS:1000523" )
					{
						set64bit=true;

						nbitmz=8;
					}

					if (nvalue=="MS:1000521")
					{
						set32bit=true;

						nbition=4;
					}

					if (nvalue=="MS:1000574")
					{
						zlib_compression=true;
						printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
						exit(EXIT_FAILURE);
					}

					if (nvalue=="MS:1000514")
					{
						ismoz=true;

						isIntensity = false;
					}

					if (nvalue=="MS:1000515")
					{
						isIntensity=true;

						ismoz = false;
					}
				}

			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}


			if("binary" == reader->Name && isFull)
			{
				if (ismoz)
				{

						Base64 (reader,nSpectLength,nbitmz,Spectrum);
					
				}

				if (isIntensity)
				{

						Base64  (reader, nSpectLength, nbition, Spectrum);
					

					if(scanNumber > 2000)
					{
						printf("Scan = %d\n", scanNumber);

						for(i=0; i < nSpectLength && i < 1000; i++)
						{
							printf("Mass[%d] = %10.5f %10.5f\n", i, Spectrum[0,i],
								Spectrum[1,i]);
						}
					}
				}
			}
		}
	}
	while (reader->Read());


	if (!spectrumNumberOK)
	{	
		Console::WriteLine("");
		Console::Write("There is no Full with number ");
		//Console::Write(nScan); Console::WriteLine(" !!!");

		return 0;
	}

	return nSpectLength;

}

// a method that reads the last Full Scan
// that comes before the scan number nScan
//i.e., if the nScan is 4056 and it is MS/MS
// then the read scan will be the last full
// scan before 4056. If 4056 is a full scan
// then it will be returned
// cv params type 1000579 is the full scan
// will only read the full scans around the mass
// value of dPrecursorMass for following analysis
// by O16/O18 software

int MzML::ReadFullScan(long nScan, double dPrecursorMass, double dInterval)
{
int nbitmz, nbition, nFullScan, nFullLength;
int i=0; 
int j=0;
int nSpectLength = 0, nEncodedLength = 0; 
long scanNumber;
Boolean isFull=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
String^ nvalue;


nFullScan = -1;

nFullLength = 0;

	do {
		if(reader->NodeType == XmlNodeType::Element)
		{
			if("spectrum" == reader->Name )
			{
				isFull=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);

				if (scanNumber == nScan - 1)
				{
					spectrumNumberOK = true;
				}
				else if(scanNumber > nScan - 1)
				{
					return 0;
				}
			}

			if((nScan - scanNumber) <= 10 && nScan >= scanNumber) 
			{
				if ("cvParam" == reader->Name)
				{
					nvalue=ReadNode(reader,"accession");

					if (nvalue=="MS:1000579")
					{
						isFull=true;

						Spectrum = gcnew array<double, 2> (2,nSpectLength);

						nFullScan = scanNumber;

						nFullLength = nSpectLength;
					}

					if(isFull)
					{
						if (nvalue=="MS:1000523" )
						{
							set64bit=true;

							nbitmz=8;
						}

						if (nvalue=="MS:1000521")
						{
							set32bit=true;

							nbition=4;
						}

						if (nvalue=="MS:1000574")
						{
							zlib_compression=true;
							printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
							exit(EXIT_FAILURE);
						}

						if (nvalue=="MS:1000514")
						{
							ismoz=true;

							isIntensity = false;
						}
						
						if (nvalue=="MS:1000515")
						{
							isIntensity=true;

							ismoz = false;
						}
					}

				}

 			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}


			if("binary" == reader->Name && isFull)
			{
				if (ismoz)
				{

						Base64 (reader,nSpectLength,nbitmz,Spectrum);
					

					ismoz=false;

					ismozOK=true;
				}

				if (isIntensity)
				{

						Base64   (reader, nSpectLength, nbition, Spectrum);
					

					isIntensity=false;

					isFull = false;


					if(spectrumNumberOK)
					{
						j = 0;

						TempSpectrum = gcnew array<double, 2> (2,nFullLength);

						for(i=0; i < nFullLength; i++)
						{
							if( (dPrecursorMass - Spectrum[0,i]) <= 0.1 &&
								(Spectrum[0,i] - dPrecursorMass - dInterval) <= 0.1 )
							{
								j++;
							}

							TempSpectrum[0,i] = Spectrum[0,i];

							TempSpectrum[1,i] = Spectrum[1,i];
						}

						Spectrum = gcnew array<double, 2> (2,j);
						

						j = 0;

						//save only portion of the spectrum
						for(i=0; i < nFullLength; i++)
						{
							if( (dPrecursorMass - TempSpectrum[0,i]) <= 0.1 &&
								(TempSpectrum[0,i] - dPrecursorMass - dInterval) <= 0.1 )
							{
								Spectrum[0,j] = TempSpectrum[0,i];

								Spectrum[1,j] = TempSpectrum[1,i];
								
								j++;
							}
						}
						
						nFullLength = j;

						return nFullLength;
					}

				}	

			}

			if(scanNumber == nScan-1)
			{
				//printf("The Corresponding Full Scan Number %d\n", (nFullScan+1));

				//printf("scanNumber = %d\n", (scanNumber + 1));

				if((nFullScan + 1) == 0)
				{
					return 0;
				}

				j = 0;

				TempSpectrum = gcnew array<double, 2> (2,nFullLength);

				for(i=0; i < nFullLength; i++)
				{
					/*if( (dPrecursorMass - Spectrum[0,i]) <= 0.1 &&
						(Spectrum[0,i] - dPrecursorMass - dInterval) <= 0.1 )*/
					if( Spectrum[0,i] >= (dPrecursorMass - 0.1) &&
						Spectrum[0,i] <= (dPrecursorMass + dInterval + 0.1) )
					{
						j++;
					}

					TempSpectrum[0,i] = Spectrum[0,i];

					TempSpectrum[1,i] = Spectrum[1,i];
				}

				Spectrum = gcnew array<double, 2> (2,j);
				

				j = 0;

				//save only portion of the spectrum
				for(i=0; i < nFullLength; i++)
				{
					/*if( (dPrecursorMass - TempSpectrum[0,i]) <= 0.1 &&
						(TempSpectrum[0,i] - dPrecursorMass - dInterval) <= 0.1 )*/
					if( TempSpectrum[0,i] > (dPrecursorMass - 0.1) &&
						TempSpectrum[0,i] <= (dPrecursorMass + dInterval + 0.1) )
					{
						Spectrum[0,j] = TempSpectrum[0,i];

						Spectrum[1,j] = TempSpectrum[1,i];
						
						j++;
					}
				}
				
				nFullLength = j;
			

				return nFullLength;
			}
		 } 

} 
	while (reader->Read());


	if (!spectrumNumberOK)
	{	
		Console::WriteLine("");
		Console::Write("There is no Full Scan with number ");
		Console::Write(nScan); Console::WriteLine(" !!!");

		return 0;
	}

	return nSpectLength;

}

int MzML::RetrieveMass(double dPrecursorMass, double dInterval)
{
	int i, j, nFullLength;

	j = 0;

	for(i=0; i < TempSpectrum->Length/2; i++)
	{
		if( (dPrecursorMass - TempSpectrum[0,i]) <= 0.1 &&
			(TempSpectrum[0,i] - dPrecursorMass - dInterval) <= 0.1 )
		{
			j++;
		}
	}

	Spectrum = gcnew array<double, 2> (2,j);

	j = 0;

	//save only portion of the spectrum
	for(i=0; i < TempSpectrum->Length/2; i++)
	{
		if( (dPrecursorMass - TempSpectrum[0,i]) <= 0.1 &&
			(TempSpectrum[0,i] - dPrecursorMass - dInterval) <= 0.1 )
		{
			Spectrum[0,j] = TempSpectrum[0,i];

			Spectrum[1,j] = TempSpectrum[1,i];
			
			j++;
		}
	}
	
	nFullLength = j;

	return nFullLength;
}

/*
*	read mzXML only for the MALDI file test.
*/
int MzXML::ReadmzXML()
{
int i=0; 
int j=0;
int nSpectLength = 0; 
Boolean isZoom=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
Boolean mzXML=false;
int peaksCount; 
String^ nvalue;

do {
         switch (reader->NodeType)
         {
  		 case XmlNodeType::Element:

			 if("mzXML" == reader->Name )
			{
				mzXML=true;
 			}


			 if("scan" == reader->Name && mzXML )
			{
				nvalue=ReadNode(reader,"peaksCount");
				peaksCount = int::Parse(nvalue);
				Console::WriteLine("peaksCount={0} nvalue={1}",peaksCount, nvalue);
 			}

            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... \n"); 

			 if (mzXML)
			 {
					array<Byte>^binaryData = gcnew array<Byte>(peaksCount*8);
					array<Byte>^base64x = gcnew array<Byte>(4);
					Spectrum = gcnew array<double, 2> (2,peaksCount*2);
					nvalue=reader->Value;
					binaryData = Convert::FromBase64String( nvalue ); 
					
					for (int k=0; k<peaksCount*8; k+=8)
					{
				
					Array::Copy(binaryData, k,base64x,0,4);
					Array::Reverse(base64x);
					Spectrum[0,k/8] = BitConverter::ToSingle( base64x, 0 );
					Array::Copy(binaryData, k+4,base64x,0,4);
					Array::Reverse(base64x);
					Spectrum[1,k/8] = BitConverter::ToSingle( base64x, 0 );

					}	
					return 0;
			 }
            break;
         case XmlNodeType::EndElement:
//			 Console::Write(" endelement... ");
//           Console::Write("{0}", reader->Name);
          break;
		 } 

} 
	   while (reader->Read());    

return 0;

}

int MzML::mzXMLScan(long nScan)

{
int nbit;
int i=0; 
int j=0;
int nSpectLength = 0, nEncodedLength = 0; 
long scanNumber, maxscanNumber;
Boolean isZoom=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
Boolean mzXML=false; 
int peaksCount; 
String^ nvalue;

	XmlTextReader^ reader = nullptr;
//   Console::WriteLine(L"ReadScan XmlTextReader...");
//		int base64len = 0;
//		array<Byte>^base64 = gcnew array<Byte>(87121*8);
		array<Byte>^base64x = gcnew array<Byte>(4);
//		array<Byte>^binaryData = gcnew array<Byte>(87121*8);
//		Spectrum = gcnew array<double, 2> (2,87121*2);

	try
	{
      reader = gcnew XmlTextReader(sFilemzml);
      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
         {
  		 case XmlNodeType::Element:

			 if("mzXML" == reader->Name )
			{
				mzXML=true;
 			}


			 if("scan" == reader->Name && mzXML )
			{
				nvalue=ReadNode(reader,"peaksCount");
				peaksCount = int::Parse(nvalue);
				Console::WriteLine("peaksCount={0} nvalue={1}",peaksCount, nvalue);

 			}


			if("spectrumList" == reader->Name )
			{
				nvalue=ReadNode(reader,"count");
				maxscanNumber=long::Parse(nvalue);
 			}
            
			if("spectrum" == reader->Name )
			{
				isZoom=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);
				if (scanNumber==nScan)
				{
					spectrumNumberOK=true;
					Spectrum = gcnew array<double, 2> (2,nSpectLength+1);

				}
	
			}
             
			if("cvParam" == reader->Name ) 
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1000497")
				{
					isZoom=true;
				}

				if (nvalue=="MS:1000497")
				{
					isZoom=true;
				}


				if (nvalue=="MS:1000523" && isZoom && spectrumNumberOK)
				{
					set64bit=true;
					nbit=8;
				}

				if (nvalue=="MS:1000521" && isZoom && spectrumNumberOK)
				{
					set32bit=true;
					nbit=4;
				}
				if (nvalue=="MS:1000574"  && isZoom && spectrumNumberOK)
				{
					zlib_compression=true;			
					printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
					exit(EXIT_FAILURE);
				}			

				if (nvalue=="MS:1000514" && isZoom && spectrumNumberOK)
				{
					ismoz=true;
				}
				if (nvalue=="MS:1000515" && isZoom && spectrumNumberOK)
				{
					isIntensity=true;
				}

 			}
			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}


				if("binary" == reader->Name )
			{
if (isZoom && spectrumNumberOK) 
{	
		i = 0;
	if (ismoz)
	{

			Base64 (reader,nSpectLength,nbit,Spectrum);

		i++;
		ismoz=false;
		ismozOK=true;  		
	}

	if (isIntensity) 
	{
		j = 0;

			Base64   (reader, nSpectLength, nbit, Spectrum);
		j++;
		isIntensity=false;

if (ismozOK) return nSpectLength;
	}	
}
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... \n"); 

			 if (mzXML)
			 {
					array<Byte>^binaryData = gcnew array<Byte>(peaksCount*8);
					Spectrum = gcnew array<double, 2> (2,peaksCount*2);
					nvalue=reader->Value;
					binaryData = Convert::FromBase64String( nvalue ); 
					
					for (int k=0; k<peaksCount*8; k+=8)
					{
				
					Array::Copy(binaryData, k,base64x,0,4);
					Array::Reverse(base64x);
					Spectrum[0,k/8] = BitConverter::ToSingle( base64x, 0 );
					Array::Copy(binaryData, k+4,base64x,0,4);
					Array::Reverse(base64x);
					Spectrum[1,k/8] = BitConverter::ToSingle( base64x, 0 );

					}	

					TempSpectrum = gcnew array<double, 2> (2,peaksCount*2);

					for(int ii=0; ii < peaksCount*2; ii++)
					{
						TempSpectrum [0,ii] = Spectrum[0, ii];

						TempSpectrum [1,ii] = Spectrum[1, ii];
					}

					return 0;
			 }
            break;
         case XmlNodeType::EndElement:
//			 Console::Write(" endelement... ");
//           Console::Write("{0}", reader->Name);
          break;
		 } 

} 
	   while (reader->Read());    

return 0;

   }
   finally
   { 
	   if (!spectrumNumberOK && !mzXML)
		{	
			Console::WriteLine("");
			Console::Write("There is no zoomScan with number ");
			Console::Write(nScan); Console::WriteLine(" !!!");
		}
//	  Console::WriteLine("nSpectLength={0}",nSpectLength);
      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilemzml);
      if ( reader != nullptr )
            reader->Close();
   }
	return nSpectLength;
}




//a routine to mainly read 
// mzXML 

int MzXML::ReadScan(long nScan)

{
int nbit;
int i=0; 
int j=0;
int nSpectLength = 0, nEncodedLength = 0; 
long scanNumber, maxscanNumber;
Boolean isZoom=false;
Boolean set32bit=false;
Boolean set64bit=false;
Boolean spectrumList=false;
Boolean spectrumNumberOK=false;
Boolean ismoz=false;
Boolean ismozOK=false;
Boolean isIntensity=false;
Boolean zlib_compression=false;
Boolean mzXML=false;
int peaksCount; 
String^ nvalue;

	XmlTextReader^ reader = nullptr;
//   Console::WriteLine(L"ReadScan XmlTextReader...");
//		int base64len = 0;
//		array<Byte>^base64 = gcnew array<Byte>(87121*8);
		array<Byte>^base64x = gcnew array<Byte>(4);
//		array<Byte>^binaryData = gcnew array<Byte>(87121*8);
//		Spectrum = gcnew array<double, 2> (2,87121*2);

	try
	{
      reader = gcnew XmlTextReader(sFilemzml);
      reader->WhitespaceHandling =  WhitespaceHandling::None;

do {
         switch (reader->NodeType)
         {
  		 case XmlNodeType::Element:

			 if("mzXML" == reader->Name )
			{
				mzXML=true;
 			}


			 if("scan" == reader->Name && mzXML )
			{
				nvalue=ReadNode(reader,"peaksCount");
				peaksCount = int::Parse(nvalue);
				Console::WriteLine("peaksCount={0} nvalue={1}",peaksCount, nvalue);


 			}


			if("spectrumList" == reader->Name )
			{
				nvalue=ReadNode(reader,"count");
				maxscanNumber=long::Parse(nvalue);
 			}
            
			if("spectrum" == reader->Name )
			{
				isZoom=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);
				if (scanNumber==nScan)
				{
					spectrumNumberOK=true;
					Spectrum = gcnew array<double, 2> (2,nSpectLength+1);

				}
	
			}
             
			if("cvParam" == reader->Name ) 
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1000497")
				{
					isZoom=true;
				}

				if (nvalue=="MS:1000497")
				{
					isZoom=true;
				}


				if (nvalue=="MS:1000523" && isZoom && spectrumNumberOK)
				{
					set64bit=true;
					nbit=8;
				}

				if (nvalue=="MS:1000521" && isZoom && spectrumNumberOK)
				{
					set32bit=true;
					nbit=4;
				}
				if (nvalue=="MS:1000574"  && isZoom && spectrumNumberOK)
				{
					zlib_compression=true;		
					printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
					exit(EXIT_FAILURE);
				}			

				if (nvalue=="MS:1000514" && isZoom && spectrumNumberOK)
				{
					ismoz=true;
				}
				if (nvalue=="MS:1000515" && isZoom && spectrumNumberOK)
				{
					isIntensity=true;
				}

 			}
			
			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}

				if("binary" == reader->Name )
			{
if (isZoom && spectrumNumberOK) 
{	
		i = 0;
	if (ismoz)
	{

			Base64 (reader,nSpectLength,nbit,Spectrum);

		i++;
		ismoz=false;
		ismozOK=true;  		
	}

	if (isIntensity) 
	{
		j = 0;

			Base64   (reader, nSpectLength, nbit, Spectrum);
		j++;
		isIntensity=false;

if (ismozOK) return nSpectLength;
	}	
}
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... \n"); 

			 if (mzXML)
			 {
					array<Byte>^binaryData = gcnew array<Byte>(peaksCount*8);
					Spectrum = gcnew array<double, 2> (2,peaksCount*2);
					nvalue=reader->Value;
					binaryData = Convert::FromBase64String( nvalue ); 
					
					for (int k=0; k<peaksCount*8; k+=8)
					{
				
					Array::Copy(binaryData, k,base64x,0,4);
					Array::Reverse(base64x);
					Spectrum[0,k/8] = BitConverter::ToSingle( base64x, 0 );
					Array::Copy(binaryData, k+4,base64x,0,4);
					Array::Reverse(base64x);
					Spectrum[1,k/8] = BitConverter::ToSingle( base64x, 0 );

					}	
					return 0;
			 }
            break;
         case XmlNodeType::EndElement:
//			 Console::Write(" endelement... ");
//           Console::Write("{0}", reader->Name);
          break;
		 } 

} 
	   while (reader->Read());    

return 0;

   }
   finally
   { 
	   if (!spectrumNumberOK && !mzXML)
		{	
			Console::WriteLine("");
			Console::Write("There is no zoomScan with number ");
			Console::Write(nScan); Console::WriteLine(" !!!");
		}
//	  Console::WriteLine("nSpectLength={0}",nSpectLength);
      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilemzml);
      if ( reader != nullptr )
            reader->Close();
   }
	return nSpectLength;
}

/*
*	The function used in reading off-sets at the end of the mzML file
*	The first line of mzML file indicates if it has been created with
*	default settings (create indices). The off-sets themselves are written
*	at the end of the mzML file.
*/
//int ReadSpectrumByteOffset(FileStream^ inStream, array<long> ^ scan_offset_array)
int MzML::ReadSpectrumByteOffset()
{
	Boolean indexList=false;
	Boolean stopread=false;
	Boolean isoffset=false;
	String^ nvalue;
	String^ tmpstr;

	//unsigned long scan,scan_offset, index_offset;

	long scan;

	long long scan_offset;

	long long index_offset;

	Console::WriteLine("ReadSpectrumByteOffset...");

//**************************************************************
	array<Byte>^readbuffer = gcnew array<Byte>(512);

	try
	{
		inStream->Seek(-512,SeekOrigin::End);	
	}
	catch (IOException^ e ) 
	{
	   Console::WriteLine( "{0}: The write operation could not "
	   "be performed because the specified "
	   "part of the file is locked.", e->GetType()->Name );
	}
	catch (ObjectDisposedException^ e ) 
	{
      Console::WriteLine( "Caught: {0}", e->Message );
    }
	catch (NotSupportedException ^ ne)
	{
		Console::WriteLine("Caught: {0}", ne->Message);
	}
	catch (ArgumentException ^ae )
	{
		Console::WriteLine("Argument Exception: {0}:", ae->Message);
	}

	inStream->Read(readbuffer,0,512); 
			  
	MemoryStream^ ms = gcnew MemoryStream; 
				
	ms->Write(readbuffer,0,512); 
	
	ms->Position =0;
				
	StreamReader^ rrr = gcnew StreamReader(ms); 
		
	do 
	{
		nvalue= rrr->ReadLine(); 
				
		if (nvalue->Contains("<indexListOffset>"))
		{
			nvalue = nvalue->Substring(nvalue->IndexOf(">")+1);

			nvalue = nvalue->Substring(0, nvalue->IndexOf("<"));

			index_offset = System::Int64::Parse(nvalue);

			//index_offset  = System::Int64::Parse(nvalue);

			//printf("Values %s %d %f\n", nvalue, nvalue->Length, (double)index_offset);

			scan_offset_array[0] = index_offset;
		}

	} 
	while (!rrr->EndOfStream);
					
	inStream->Seek(index_offset,SeekOrigin::Begin);	
					
	array<Byte>^buffer = gcnew array<Byte>((inStream->Length)-index_offset);

	//printf("Buffer Size %d %d\n", buffer->Length, (inStream->Length)-index_offset);
					
	inStream->Read(buffer,0, (inStream->Length) - index_offset); 

	
	//for(int i=0; i < buffer->Length; i++)
	//{
	//	//printf("%c", (char)buffer[i]);
	//}

	printf("inStream size: %.0f %.0f\n", (double)inStream->Length, (double)index_offset);
					
	ms->Position =0;
	
	long long ltemp;
	
	ltemp = inStream->Length-index_offset;

	//printf("\nLLLLLLTEMP: %f %f %f\n", (double)inStream->Length, (double)index_offset, (double)ltemp);

	//ms->Write(buffer, 0, inStream->Length-index_offset); ms->Position =0;
	ms->Write(buffer, 0, ltemp); ms->Position =0;
					
	StreamReader^ sss = gcnew StreamReader(ms); 
					
	ms->Position =0; 
					
	XmlTextReader^ offsetreader = nullptr;
					
	offsetreader = gcnew XmlTextReader(sss); 


	/*printf("PASSEDE strem = %f %f %s %d\n", (double)(inStream->Length-index_offset), (double)ltemp, offsetreader->ReadString()->ToString(), offsetreader->Depth);

	printf("Streeeeem %s %s\n", offsetreader->LocalName, offsetreader->ReadString());*/
					
	offsetreader->Read();

	//printf("Starting node type %s\n", offsetreader->NodeType);
//**************************************************************
	do 
	{
         switch (offsetreader->NodeType)
         {
         case XmlNodeType::Element:			 
		 if (offsetreader->Name == "indexList")
		 {
			 indexList = true; 
		 }
		 if (indexList) 
		 {
			  if (offsetreader->Name == "offset") 
			  {
				  isoffset=true;
				  nvalue=ReadNode(offsetreader,"idRef");
				  tmpstr =  nvalue->Substring(nvalue->IndexOf("scan=")+5);
				  scan= long::Parse(tmpstr);
			  }
		  } 
		  break;
         case XmlNodeType::Text:
			 if(isoffset)
			 {
				 scan_offset = System::Int64::Parse(offsetreader->Value);

				 scan_offset_array[scan] = scan_offset;
				 
				 isoffset=false;
			 }
            break;
         case XmlNodeType::EndElement:
			 
			 nvalue = offsetreader->Name;
			 
			 if (nvalue=="index") 
				 return 0;
            break;
		 } 
	  } 
	   while (offsetreader->Read() );    

	  return 0;
}


/*
*	The method will read the off-sets (byte positions in the mzML file)
*	for every spectrum in the mzML file.
*
*
*/

Boolean MzML::ReadOffset()
{
	
	//FileStream^ inStream = File::OpenRead(sFilemzml);

	//ReadSpectrumByteOffset(inStream, scan_offset_array); 
	ReadSpectrumByteOffset(); 

	return true;
}


/*
*	 A method to read an ms spectrum (scan) from mzML file
*	 that has been created using indexing (default for msconvert)
*	 Indices (byte off-sets in the file) are recorded at the end of
*	 the mzML file. Each index points to the start of the spectrum node
*	 in the mzML file 
*	 the reader handle to the mzML file is instantiated in the ctor of the
*	 class in .h file. Also the total number of spectra are read there.
*	nMSLevel - indicates the MS level, 1 means Full MS, 2 means MS/MS, 3 means MS3
* //by Dr. Mirek Gilski
*	the indices are read in by the ReadSpectrumByteOffset method
*
*/
int MzML::ReadIndexedScan(long nScan, unsigned short *nMSLevel)
{
	int nbit, i=0, j=0, nSpectLength = 0; 

	unsigned short iMSLevel = 0;

	long scanNumber, specbytesize;

	int nencodedLength = 0;

	float fRetTime;

	Boolean set32bit=false, set64bit=false;

	Boolean  spectrumNumberOK=false;

	Boolean ismoz=false, ismozOK=false, isIntensity=false;

	Boolean zlib_compression=false;
	 
	String^ nvalue;
	
	if (nScan < nTotalSpectra)
	{
		specbytesize = scan_offset_array[nScan+1]-scan_offset_array[nScan];
	}				
	else 
	{
		specbytesize = scan_offset_array[0]-scan_offset_array[nScan];
	}


	if(specbytesize > 0) //Mahbubur Rahman; 10/26/2015
	{    

		array<Byte>^readbuffer = gcnew array<Byte>(specbytesize);

		inStream->Seek(scan_offset_array[nScan],SeekOrigin::Begin);
	
		inStream->Read(readbuffer,0,specbytesize); 
  	
		MemoryStream^ ms = gcnew MemoryStream; 
	
		ms->Write(readbuffer,0,specbytesize); 
	
		ms->Position = 0;

		XmlTextReader^ index_reader = nullptr;
	
		index_reader = gcnew XmlTextReader(ms);
	
		index_reader->WhitespaceHandling =  WhitespaceHandling::None;

		if (!indexedmzml) 
		{
			index_reader = reader;
		}


	

	do {
         switch (index_reader->NodeType)
         {
         case XmlNodeType::Element:
            
			 //printf("Element %s\n",index_reader->Name);

			if("spectrum" == index_reader->Name )
			{
				spectrumNumberOK=false;
			
				nvalue=ReadNode(index_reader,"index");

				scanNumber=long::Parse(nvalue)+1;

				index_reader->MoveToAttribute(0);

				nvalue=ReadNode(index_reader,"defaultArrayLength");

				nSpectLength=int::Parse(nvalue);

				if(0 == nSpectLength)
				{
					return 0;
				}

				if (scanNumber == nScan)
				{
					spectrumNumberOK = true;
					
					Spectrum = gcnew array<double, 2> (2, nSpectLength+1);
				}
				else
				{
					printf("Wrong Scan\n");
				}
	
			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == index_reader->Name)
			{
				nvalue = ReadNode(index_reader,"encodedLength");

				nencodedLength = int::Parse(nvalue);
			}


			if("cvParam" == index_reader->Name) 
			{			
				if (spectrumNumberOK)
				{
			
					nvalue=ReadNode(index_reader,"accession");

					if(nvalue == "MS:1000511")
					{
						String ^value;

						value = ReadNode(index_reader,"value");

						iMSLevel = short(int::Parse(value));

						*nMSLevel = iMSLevel;

						//printf("nMSLevel = %d\n", (*nMSLevel));
					}
					else if("MS:1000016" == nvalue)
					{
						String ^value;

						value = ReadNode(index_reader, "value");

						fRetTime = float::Parse(value);	

						fCurrentRetTime = fRetTime;  //fRetTimeCurrent is a public class variable. It is accessible via the class handle.
					}
					else if (nvalue=="MS:1000514")
					{
						ismoz=true;
					}
					else if (nvalue=="MS:1000515")
					{
						isIntensity=true;
					}

					if (nvalue=="MS:1000523")
					{
						set64bit=true;

						nbit=8;
					}

					if (nvalue=="MS:1000521")
					{
						set32bit=true;

						nbit=4;
					}
					if (nvalue=="MS:1000574")
					{
						zlib_compression=true;		
						printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
						exit(EXIT_FAILURE);
					}			

				}
				break; 
			}
		
			if("binary" == index_reader->Name)
			{
					if (spectrumNumberOK) 
					{	
							if (ismoz)
							{

									//printf("m/z Reading, nbiy = %d\n", nbit);

									//Base64 (index_reader,nSpectLength,nbit,Spectrum);

									Base64_Unified (index_reader,nSpectLength,nbit,Spectrum, true, false);
								

								ismoz=false;
								
								ismozOK=true;  		
							}

							if (isIntensity)
							{

									//printf("Intensity Reading, nbiy = %d\n", nbit);

									
									//Base64(index_reader, nSpectLength, nbit, Spectrum);
									
									Base64_Unified(index_reader, nSpectLength, nbit, Spectrum, false, true);
								

								isIntensity=false;

								//this where the program
								//ends successfull
								//it does not reach the terminul read
								if (ismozOK) 
								{
									/*for(i=0; i < nSpectLength; i++)
									{
										printf("%10.5f %10.5f \n", Spectrum[0,i], Spectrum[1,i]);
									}*/

									TempSpectrum = gcnew array<double, 2> (2, nSpectLength+1);

									for(int ii=0; ii < nSpectLength +1; ii++)
									{
										TempSpectrum[0,ii] = Spectrum[0,ii];
										
										TempSpectrum[1,ii] = Spectrum[1,ii];
									}

									index_reader->Close();

									return nSpectLength;
								}
							}	
					}  
					
					break; 
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... ");
//           Console::Write(reader->Value);
            break;
         case XmlNodeType::EndElement:
			 nvalue = index_reader->Name;
			 if (nvalue=="spectrum" && indexedmzml) 
				 return 0;
            break;
		 } 
} 
	   while (index_reader->Read()); 

}

    Console::WriteLine( "\nProcessing of the file {0} complete.xxx\n", sFilemzml);

	return nSpectLength;

} //The End !!!




/*
*	 A method to read a FULL ms spectrum (scan) from mzML file
*	 that has been created using indexing (default for msconvert)
*	 Indices (byte off-sets in the file) are recorded at the end of
*	 the mzML file. Each index points to the start of the spectrum node
*	 in the mzML file 
*	 the reader handle to the mzML file is instantiated in the ctor of the
*	 class in .h file. Also the total number of spectra are read there.
*	 nMSLevel is equal 1 then it returns the full spectrum,
*    otherwise it does not read the spectrum and return length 0 spectrum;
*
*/
int MzML::ReadIndexedFullScan(long nScan, unsigned short *nMSLevel)
{
	int nbit, i=0, j=0, nSpectLength = 0; 

	unsigned short iMSLevel = 0;

	long scanNumber, specbytesize;

	int nencodedLength = 0;

	float fRetTime;

	Boolean set32bit=false, set64bit=false;

	Boolean  spectrumNumberOK=false;

	Boolean ismoz=false, ismozOK=false, isIntensity=false;

	Boolean zlib_compression=false;
	 
	String^ nvalue;
	
	if (nScan < nTotalSpectra)
	{
		specbytesize = scan_offset_array[nScan+1]-scan_offset_array[nScan];
	}				
	else 
	{
		specbytesize = scan_offset_array[0]-scan_offset_array[nScan];
	}


	if(specbytesize > 0) //Mahbubur Rahman; 10/26/2015
	{    

		array<Byte>^readbuffer = gcnew array<Byte>(specbytesize);

		inStream->Seek(scan_offset_array[nScan],SeekOrigin::Begin);
	
		inStream->Read(readbuffer,0,specbytesize); 
  	
		MemoryStream^ ms = gcnew MemoryStream; 
	
		ms->Write(readbuffer,0,specbytesize); 
	
		ms->Position = 0;

		XmlTextReader^ index_reader = nullptr;
	
		index_reader = gcnew XmlTextReader(ms);
	
		index_reader->WhitespaceHandling =  WhitespaceHandling::None;

		if (!indexedmzml) 
		{
			index_reader = reader;
		}


	

	do {
         switch (index_reader->NodeType)
         {
         case XmlNodeType::Element:
            
			 //printf("Element %s\n",index_reader->Name);

			if("spectrum" == index_reader->Name )
			{
				spectrumNumberOK=false;
			
				nvalue=ReadNode(index_reader,"index");

				scanNumber=long::Parse(nvalue)+1;

				index_reader->MoveToAttribute(0);

				nvalue=ReadNode(index_reader,"defaultArrayLength");

				nSpectLength=int::Parse(nvalue);

				if(0 == nSpectLength)
				{
					return 0;
				}

				if (scanNumber == nScan)
				{
					spectrumNumberOK = true;
					
					Spectrum = gcnew array<double, 2> (2, nSpectLength+1);
				}
				else
				{
					printf("Wrong Scan\n");
				}
	
			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == index_reader->Name)
			{
				nvalue = ReadNode(index_reader,"encodedLength");

				nencodedLength = int::Parse(nvalue);
			}


			if("cvParam" == index_reader->Name) 
			{			
				if (spectrumNumberOK)
				{
			
					nvalue=ReadNode(index_reader,"accession");

					if(nvalue == "MS:1000511")
					{
						String ^value;

						value = ReadNode(index_reader,"value");

						iMSLevel = short(int::Parse(value));

						*nMSLevel = iMSLevel;

						if(1 != iMSLevel)
						{
							return 0;
						}
						//printf("nMSLevel = %d\n", (*nMSLevel));
					}
					else if("MS:1000016" == nvalue)
					{
						String ^value;

						value = ReadNode(index_reader, "value");

						fRetTime = float::Parse(value);	

						fCurrentRetTime = fRetTime;  //fRetTimeCurrent is a public class variable. It is accessible via the class handle.
					}
					else if (nvalue=="MS:1000514")
					{
						ismoz=true;
					}
					else if (nvalue=="MS:1000515")
					{
						isIntensity=true;
					}

					if (nvalue=="MS:1000523")
					{
						set64bit=true;

						nbit=8;
					}

					if (nvalue=="MS:1000521")
					{
						set32bit=true;

						nbit=4;
					}
					if (nvalue=="MS:1000574")
					{
						zlib_compression=true;			
						printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
						exit(EXIT_FAILURE);
					}			

				}
				break; 
			}
		
			if("binary" == index_reader->Name)
			{
					if (spectrumNumberOK) 
					{	
							if (ismoz)
							{

									//printf("m/z Reading, nbiy = %d\n", nbit);

									//Base64 (index_reader,nSpectLength,nbit,Spectrum);

									Base64_Unified (index_reader,nSpectLength,nbit,Spectrum, true, false);
								

								ismoz=false;
								
								ismozOK=true;  		
							}

							if (isIntensity)
							{

									//printf("Intensity Reading, nbiy = %d\n", nbit);

									
									//Base64(index_reader, nSpectLength, nbit, Spectrum);
									
									Base64_Unified(index_reader, nSpectLength, nbit, Spectrum, false, true);
								

								isIntensity=false;

								//this where the program
								//ends successfull
								//it does not reach the terminul read
								if (ismozOK) 
								{
									/*for(i=0; i < nSpectLength; i++)
									{
										printf("%10.5f %10.5f \n", Spectrum[0,i], Spectrum[1,i]);
									}*/

									TempSpectrum = gcnew array<double, 2> (2, nSpectLength+1);

									for(int ii=0; ii < nSpectLength +1; ii++)
									{
										TempSpectrum[0,ii] = Spectrum[0,ii];
										
										TempSpectrum[1,ii] = Spectrum[1,ii];
									}

									index_reader->Close();

									return nSpectLength;
								}
							}	
					}  
					
					break; 
			}
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... ");
//           Console::Write(reader->Value);
            break;
         case XmlNodeType::EndElement:
			 nvalue = index_reader->Name;
			 if (nvalue=="spectrum" && indexedmzml) 
				 return 0;
            break;
		 } 
} 
	   while (index_reader->Read()); 

}

    Console::WriteLine( "\nProcessing of the file {0} complete.xxx\n", sFilemzml);

	return nSpectLength;

} //The End !!!


/*
*	 A method to read precursor mass and charge
*    state from an mzML file
*
*/
int MzML::PrecursorMassAndCharge()
{
	String ^nvalue, ^sTemp;

	bool bPrecursor = false;

	bool bIon = false, bSpectrum = false;

	double m_z, dPrecursor, mProton, dRetTime;

	int nCharge, nScan, nTotal = 0, i;

	char szPrecursorFile[1024];

	FILE *fpPrs;

	i = sFilemzml->LastIndexOf(".mzML");

	sTemp = sFilemzml->Substring(0, i) + ".Prs";

	for(i=0; i < sTemp->Length; i++)
	{
		szPrecursorFile[i] = sTemp[i];
	}

	szPrecursorFile[i] = '\0';

	fpPrs = fopen(szPrecursorFile, "w");

	if(NULL == fpPrs)
	{
		printf("Cannot write to %s\n", szPrecursorFile);

		exit (1);
	}
	else
	{
		fprintf(fpPrs, "Scan\tCharge  \tMass\n");
	}

	mProton = 1.00727642;

	XmlTextReader^ reader = nullptr;
//   Console::WriteLine(L"ReadScan XmlTextReader...");
	
	try
	{
      reader = gcnew XmlTextReader(sFilemzml);

      reader->WhitespaceHandling =  WhitespaceHandling::None;

		do {
         switch (reader->NodeType)
         {
         case XmlNodeType::Element:

			if("spectrum" == reader->Name)
			{
				nvalue = ReadNode(reader, "id");

				bSpectrum = true;

				i = nvalue->LastIndexOf("=") + 1;

				nScan = int::Parse(nvalue->Substring(i));
			}
			else if("scan" == reader->Name)
			{

			}
			else if("precursor" == reader->Name )
			{
				nvalue=ReadNode(reader,"count");

				bPrecursor = true;
 			}
            
			if(bPrecursor)
			{
				if("SelectedIon" == reader->Name)
				{
					bIon = true;
				}
			}


		   if(false == bPrecursor && "cvParam" == reader->Name)
			{
				nvalue = ReadNode(reader, "accession");

				if("MS:1000016" == nvalue)
				{
					nvalue = ReadNode(reader, "value");

					dRetTime = float::Parse(nvalue);				  
				}
			}

			if(bPrecursor)
			{
				if("cvParam" == reader->Name)
				{
				  nvalue = ReadNode(reader, "accession");

				  if("MS:1000744" == nvalue)  // m/z
				  {
					 nvalue = ReadNode(reader, "value");

					 m_z = float::Parse(nvalue);

					 nTotal++;
				  }
				  else if("MS:1000041" == nvalue)
				  {
					nvalue = ReadNode(reader, "value");

					nCharge = int::Parse(nvalue);

					dPrecursor = m_z*(double)nCharge - mProton*(double)(nCharge - 1);

					fprintf(fpPrs, "%6d %5d %15.5f  %6.2f\n", nScan, nCharge, dPrecursor, dRetTime);
				  }
				  
				}
			}
			
            break;
         case XmlNodeType::Text:
//			 Console::Write(" mgtext... ");
//           Console::Write(reader->Value);
            break;
         case XmlNodeType::EndElement:
			if("precursor" == reader->Name )
			{
				nvalue=ReadNode(reader,"count");

				bPrecursor = false;
 			}
			else if("SelectedIon" == reader->Name)
			{
				bIon = true;
			}
            
			break;
		 } 

} 
	   while (reader->Read());  

	   fclose(fpPrs);

return 0;

   }
   finally
   { 

	  printf ("total of  %d spectra\n", nTotal);

      Console::WriteLine( "\nProcessing of the file {0} complete.\n", sFilemzml);
      if ( reader != nullptr )
            reader->Close();
   }

	return 0;
}

/*
*	the method closes all file streams, XML streams
*    and if necessary releases all memory (assigned in
*	global).
*/

bool MzML::CloseFiles()
{
	reader->Close();

	inStream->Close();

	delete Spectrum, TempSpectrum, scan_offset_array;

	int i, l;

	l = FullScanChromatogram->Count;

	if(l > 0)
	{
		for(i = 0; i < l; i++)
		{
			delete FullScanChromatogram[i]->Intensity;

			delete FullScanChromatogram[i]->moverz;
		}
	}

	delete FullScanChromatogram;     //good to delete this!!!

	return true;
}

/*
*  a method to build a chromatogram from full scans of all 
*  elution times. It uses the indexed address of each scan. The addresses are 
*  provided at the end of the indexed mzML file. The data are stored in the Lists
*  The chromatogram is build into the two-dimensional List called FullScanChromatogram
*  FullScanChromatogram is a publically available variable of the class MzML
*  Each element of the List is an MS1, scan number and elution time, a variable (also a List) defined
*  as FullScan.
*/

void MzML::IndexedChromatogramBuild()
{
	FullScan ^fullScan = gcnew FullScan;

	//FullScanChromatogram = gcnew List <FullScan>;

	long nScan;

	int nPeaks, i, k;

	unsigned short mslevel;

	FullScanChromatogram = gcnew List <FullScan ^>;

	
	i = sFilemzml->LastIndexOf("\\") + 1;

	if(i >= 1)
	{
		printf("Reading MS1 Chromatogram for %s\n", sFilemzml->Substring(i)); 
	}
	else
	{
		printf("Reading MS1 Chromatogram for %s\n", sFilemzml);
	}

	k = 0;

	for(nScan = 1; nScan <= nTotalSpectra; nScan++)
	{
		nPeaks = ReadIndexedFullScan(nScan, &mslevel);

		if(nPeaks > 0)
		{
			fullScan = gcnew FullScan;

			fullScan->Intensity = gcnew array <float> (Spectrum->Length/2);

			fullScan->moverz = gcnew array <double> (Spectrum->Length/2);


			fullScan->ScanNumber = nScan;

			fullScan->RetTime    = fCurrentRetTime;

			for(i = 0; i < Spectrum->Length/2; i++)
			{
				fullScan->moverz[i] = Spectrum[0, i];

			    fullScan->Intensity[i] = (unsigned int)Spectrum[1, i];
			}

			////add the newly read full scan into the chromatogram


			FullScanChromatogram->Add(fullScan);
			

			k++;

			//printf("nScan = %d RetTimd = %10.5f SpectLe = %d, k = %d\n", fullScan->ScanNumber, fullScan->RetTime, Spectrum->Length/2, k);

			delete fullScan;
		}
	}

	i = sFilemzml->LastIndexOf("\\") + 1;

	if(i >= 1)
	{
		printf("Finished reading and storing %d MS1 scans from %s\n",  k, sFilemzml->Substring(i)); 
	}
	else
	{
		printf("Finished reading and storing %d MS1 scans from %s\n",  k, sFilemzml);
	}
}



/*
*  a method to build a chromatogram from full scans of all 
*  elution times. It each scan sequentially. DOES NOT need file addresses that are 
*  provided at the end of the indexed mzML file. Can be used with non-indexed mzML files.
*  The data are stored in the Lists
*  The chromatogram is build into the two-dimensional List called FullScanChromatogram
*  FullScanChromatogram is a publically available variable of the class MzML
*  Each element of the List is an MS1, scan number and elution time, a variable (also a List) defined
*  as FullScan.
*/

void MzML::SequentialChromatogramBuild()
{
	FullScan ^fullScan = gcnew FullScan;

	//FullScanChromatogram = gcnew List <FullScan>;

	long nScan;

	int nPeaks, k, i;

	unsigned short mslevel;

	FullScanChromatogram = gcnew List <FullScan ^>;
	

	i = sFilemzml->LastIndexOf("\\") + 1;

	if(i >= 1)
	{
		printf("Started reading MS1 Chromatogram for %s\n", sFilemzml->Substring(i)); 
	}
	else
	{
		printf("Started reading MS1 Chromatogram for %s\n", sFilemzml);
	}


	int nbitmz, nbition;

	int j=0;
	int nSpectLength = 0, nEncodedLength = 0; 
	long scanNumber;
	Boolean isFull=false;
	Boolean set32bit=false;
	Boolean set64bit=false;
	Boolean spectrumList=false;
	Boolean spectrumNumberOK=false;
	Boolean ismoz=false;
	Boolean ismozOK=false;
	Boolean isIntensity=false;
	Boolean zlib_compression=false;
	String^ nvalue;

    k = 0;

	do {
		if(reader->NodeType == XmlNodeType::Element)
		{
			if("spectrum" == reader->Name )
			{
				isFull=false;
				spectrumNumberOK=false;
				nvalue=ReadNode(reader,"index");
				scanNumber=long::Parse(nvalue);
				reader->MoveToAttribute(0);
				nvalue=ReadNode(reader,"defaultArrayLength");
				nSpectLength=int::Parse(nvalue);
	
			}


			if ("cvParam" == reader->Name)
			{
				nvalue=ReadNode(reader,"accession");

				if (nvalue=="MS:1000579")
				{
					isFull=true;

					Spectrum = gcnew array<double, 2> (2,nSpectLength);
				}

				if(isFull)
				{
				
					if (nvalue=="MS:1000523" )
					{
						set64bit=true;

						nbitmz=8;
					}

					if (nvalue=="MS:1000521")
					{
						set32bit=true;

						nbition=4;
					}

					if (nvalue=="MS:1000574")
					{
						zlib_compression=true;
						printf("\n\n==============\nWe do not support compressed mz files. Exiting ...\n\n");
						exit(EXIT_FAILURE);
					}

					if (nvalue=="MS:1000514")
					{
						ismoz=true;

						isIntensity = false;
					}

					if (nvalue=="MS:1000515")
					{
						isIntensity=true;

						ismoz = false;
					}
				}

			}

			//read the encodedLength of the
			// binary data
			if("binaryDataArray" == reader->Name)
			{
				nvalue = ReadNode(reader,"encodedLength");

				nEncodedLength = int::Parse(nvalue);
			}


			if("binary" == reader->Name && isFull)
			{
				if (ismoz)
				{


						Base64 (reader,nSpectLength,nbitmz,Spectrum);
					
				}

				if (isIntensity)
				{

						Base64  (reader, nSpectLength, nbition, Spectrum);
					

					fullScan = gcnew FullScan;

					fullScan->Intensity = gcnew array <float> (Spectrum->Length/2);

					fullScan->moverz = gcnew array <double> (Spectrum->Length/2);


					//fullScan->ScanNumber = nScan;

					//fullScan->RetTime    = fCurrentRetTime;

					for(i = 0; i < Spectrum->Length/2; i++)
					{
						fullScan->moverz[i] = Spectrum[0, i];

						fullScan->Intensity[i] = (unsigned int)Spectrum[1, i];
					}

					////add the newly read full scan into the chromatogram


					FullScanChromatogram->Add(fullScan);
			

					k++;

					//printf("nScan = %d RetTimd = %10.5f SpectLe = %d, k = %d\n", fullScan->ScanNumber, fullScan->RetTime, Spectrum->Length/2, k);

					delete fullScan;

				}
			}
		}
	}
	while (reader->Read());


	if (!spectrumNumberOK)
	{	
		Console::WriteLine("");
		Console::Write("There is no Full with number ");
		//Console::Write(nScan); Console::WriteLine(" !!!");

		return;
	}



     k = 0;

	for(nScan = 1; nScan <= nTotalSpectra; nScan++)
	{
		nPeaks = ReadIndexedFullScan(nScan, &mslevel);

		if(nPeaks > 0)
		{
			fullScan = gcnew FullScan;

			fullScan->Intensity = gcnew array <float> (Spectrum->Length/2);

			fullScan->moverz = gcnew array <double> (Spectrum->Length/2);


			fullScan->ScanNumber = nScan;

			fullScan->RetTime    = fCurrentRetTime;

			for(i = 0; i < Spectrum->Length/2; i++)
			{
				fullScan->moverz[i] = Spectrum[0, i];

			    fullScan->Intensity[i] = (unsigned int)Spectrum[1, i];
			}

			////add the newly read full scan into the chromatogram


			FullScanChromatogram->Add(fullScan);
			

			k++;

			//printf("nScan = %d RetTimd = %10.5f SpectLe = %d, k = %d\n", fullScan->ScanNumber, fullScan->RetTime, Spectrum->Length/2, k);

			delete fullScan;
		}
	}

	i = sFilemzml->LastIndexOf("\\") + 1;

	if(i >= 1)
	{
		printf("Finished reading and storing %d MS1 scans from %s\n",  k, sFilemzml->Substring(i)); 
	}
	else
	{
		printf("Finished reading and storing %d MS1 scans from %s\n",  k, sFilemzml);
	}
}
