// PeakDetectionIntegration.h
/*
*  PeakDetectionIntegration accepts
*  a file name - which is mzML file
*  and initialize a class to read it
*  Methods of this class will detect/integrate
*  peaks in m/z dimention and time dimensions
*
*/

#pragma once

using namespace System;
using namespace System::Collections::Generic;
using namespace XmlMzML; 
using namespace System::IO;

namespace PeakDetectionIntegration {

//Added by JA 2016.10.26
	public ref class RetTimeRange
	{
	public:
		float Intensity;

		int RTi;

		bool bRT;
	};

	public ref class PeakProfile
	{
	public:
		float Intensity;

		float RetTime;

		long ScanNumber;
	};

	public ref class PeakIndex    //three integers that show the position, start and end of a peak in m/z domain (MS1)
	{							  // iMaxPosition is the index (in the current MS1) of the peak
	public:						  // iLeftEnd is the index of the left end (start) of the peak
		int iMaxPosition;         // iRightEnd is the index of the right end (end) of the peak
								  // the chromatographic position is determined by the FullScanNumber
		int iLeftEnd;	          // iChromNumbe is the number in chromatogram, it is different from the FullScannumber
									// because the Full Scans are added one at a time into the chromatogram
		int iRightEnd;

		float PeakArea;

		float RetTime;

		long FullScanNumber;

		long iChromNumber; 

	};

	public ref class DetectAndIntegratePeak
	{
		public:

		MzML ^mzml;

		double dMassWindow, MassAccuracyPPM;

		double setMassWindow;

		double mProton = 1.00727642, mHydrogen = 1.007825035, deltaC13 = 13.003354838 - 12.0;

		double mDeuterium = 2.014102;

		float fElutionTimeWindow;    // difference from the time that MS2 was triggered

		float fPeakWidth, fPeakArea, fPeakHeight, fPeakStart, fPeakEnd;

		float fPeakThreshold;   //threshold value from the Apex intensity in 
								// elution profile to stop peak start and end
		StreamWriter ^fpPeakOutPut;

		DetectAndIntegratePeak(String ^file_mzML, double dWindow, float fTimeWindow)
		{
			mzml = gcnew MzML(file_mzML);

			dMassWindow = dWindow;

			setMassWindow = dWindow;

			fElutionTimeWindow = fTimeWindow;

			fPeakWidth = fPeakArea = fPeakHeight = fPeakStart = fPeakEnd = 0.;
			//read offset so that to 
			// use indexed reading
			mzml->ReadOffset();
		}

		

		//keep this till the previous definition works smoothly, then detelete.
		DetectAndIntegratePeak(String ^file_mzML, double dWindow, float fTimeWindow, double PPMMassAccuracy)
		{
			mzml = gcnew MzML(file_mzML);

			dMassWindow = dWindow;

			setMassWindow = dWindow;

			MassAccuracyPPM = PPMMassAccuracy;

			fElutionTimeWindow = fTimeWindow;

			fPeakWidth = fPeakArea = fPeakHeight = fPeakStart = fPeakEnd = 0.;
			//read offset so that to 
			// use indexed reading
			mzml->ReadOffset();
		}

		DetectAndIntegratePeak(String ^file_mzML, double dWindow, float fTimeWindow, double PPMMassAccuracy, StreamWriter ^fpOutPut)
		{
			mzml = gcnew MzML(file_mzML);

			dMassWindow = dWindow;

			setMassWindow = dWindow;

			MassAccuracyPPM = PPMMassAccuracy;

			fElutionTimeWindow = fTimeWindow;

			fPeakWidth = fPeakArea = fPeakHeight = fPeakStart = fPeakEnd = 0.;

			fpPeakOutPut = fpOutPut;

			fPeakThreshold = 0.2;

			//read offset so that to 
			// use indexed reading
			mzml->ReadOffset();
		}

		// TODO: Add your methods for this class here.

		void Close_mzML_file();

		void ReadIsotopes(int nScan, double dMonoPrecur, int nCharge, array <double, 2> ^Isotopes);

		void ReadIsotopes(double dMonoPrecur, int nCharge, array<double, 2> ^ Spectrum, array <double, 2> ^Isotopes);

		int BasePeakApex(int nScan, double dMonoPrecur, int nCharge, int iDirection);

		void IntegratePeakTimeMz(float RetTime, int nScan, double dMonoPrecur, int nCharge, float *fElutionStart,
			float *fElutionEnd, array <double, 2> ^Isotopes, array <double, 1> ^w);

		void FindPeakStartEnd(float RetTime, double dMonoPrecur, int nCharge, float *fStart, float *fEnd, int *jPBStart, int *jPBEnd, float *mzPBStart, float *mzPBEnd, float *LW, float *RW, float *RetentionTimePeak,array <double, 1> ^w);

		void FindPeakStartEnd_With_MassAccuracy_Charge(float RetTime, double dMonoPrecur, int nCharge, float *fStart, float *fEnd, int *jPBStart,
			int *jPBEnd, float *mzPBStart, float *mzPBEnd, float *LW, float *RW, float *RetentionTimePeak, array <double, 1> ^w);

		//Added by JA 2016.10.31
		void FindmzPeakWidth(float RetTime, double dMonoPrecur, int nCharge, float *fStart, float *fEnd, int *jPBStart, int *jPBEnd, float *mzPBStart, float *mzPBEnd, float *LW, float *RW, float *RetentionTimePeak,array <double, 1> ^w, int kPeak);

		void FindPeakPosition(float RetTime, double dMonoPrecur, int nCharge);

		void BuildChromatogram(short unsigned int iMethod);

		PeakIndex ^ MoverZ_Peak_Position_Width(FullScan ^MS1_Scan, double dTheoreticalM_z, double dMAccuracy);

		float QSelect_Float(List <float> ^lArray);

		//int CompareByIntensity(float ^x, float ^y);
	};
}
