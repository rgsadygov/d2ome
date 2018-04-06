// This is the main DLL file.

#include "stdafx.h"
#include "stdio.h"
#include "XmlMzML.h"

#using <System.Xml.dll>

using namespace System;
using namespace System::IO;
using namespace System::Xml;
using namespace XmlMzML;
using namespace ManagedZLib;


String^ ReadNode(XmlTextReader^ reader,String^ nodeName)
{
	String^ nodeValue;
	while (reader->MoveToNextAttribute()) 
	{
	  if (reader->Name==nodeName) 
		{ 
			nodeValue=reader->Value;
		}
	}
	return nodeValue; 
}

int Base64 (XmlTextReader^ reader,int nSpect, int nb, array<double,2>^numbers)
{ 
		double value;
		int base64len = 0;
		array<Byte>^base64 = gcnew array<Byte>(nSpect*8+8);
		base64len = reader->ReadBase64( base64, 0, nSpect*nb );

		if ( nb==8 ) 
		{
					for (int k=0; k<nSpect*8; k=k+8)
					{
						numbers[0,k/8] = BitConverter::ToDouble( base64, k );
					}

		}
			else 
					for (int k=0; k<nSpect*4; k=k+4)
					{
						numbers[1,k/4] = BitConverter::ToSingle( base64, k );
					}

	return 0;
}

int ReadZlib (XmlTextReader^ reader,int nSpect, int nb, array<double,2>^numbers ) 
{ 
		//double value;
		int readlen,k;
		int base64len = 0, zliblen=0;
		array<Byte>^base64 = gcnew array<Byte>(nSpect*8+8);
		array<Byte>^zlib = gcnew array<Byte>(nSpect*8);
		array<Byte>^unzlib = gcnew array<Byte>(nSpect*8);
		CompressionStream^ zlibStream = nullptr;
		MemoryStream^ ms = gcnew MemoryStream; 
		zliblen = reader->ReadBase64( zlib, 0, nSpect*nb );
		ms->Write(zlib,0,zliblen);
		ms->Position = 0;
		zlibStream = gcnew CompressionStream (ms, CompressionOptions::Decompress );
		readlen = zlibStream->Read( unzlib, 0, nSpect*nb);
		//ms->Position = 0;

		if ( nb == 8 ) 
		{
						for ( k = 0; k < nSpect*8; k += 8 )
						{
							numbers[0,k/8] = BitConverter::ToDouble( unzlib, k );
						}
		}
			else 
						for ( k = 0; k < nSpect*4; k += 4 ) 
						{
							numbers[1,k/4] = BitConverter::ToSingle( unzlib, k );
						//	value = BitConverter::ToSingle( unzlib, 0 );
						}
	return 0;
}



//*************

// a method returns the lenght of the spectrum
// (scan number) nScan

int MzML::ReadScan(long nScan)

{
int nbit;
int i=0; 
int j=0;
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
double moz, Intensity; 
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

				if("binary" == reader->Name )
			{
if (isZoom && spectrumNumberOK) 
{	
		i = 0;
	if (ismoz)
	{
		if (zlib_compression) 
			
			ReadZlib (reader, nSpectLength, nbit, Spectrum);
		else
			Base64 (reader,nSpectLength,nbit,Spectrum);

		i++;
		ismoz=false;
		ismozOK=true;  		
	}

	if (isIntensity)
	{
		j = 0;
		if (zlib_compression) 
			ReadZlib (reader, nSpectLength, nbit, Spectrum);
		else
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


