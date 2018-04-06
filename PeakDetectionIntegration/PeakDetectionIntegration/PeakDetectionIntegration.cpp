// This is the main DLL file.

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"
#include <vector>

#include "PeakDetectionIntegration.h"

using namespace System;
using namespace PeakDetectionIntegration;
using namespace System::Collections::Generic;
//using namespace std;


double mProton = 1.00727642, mHydrogen = 1.007825035, deltaC13 = 13.003354838 - 12.0;

double mDeuterium = 2.014102;

double dMassAccuracy = 0.05, dNoiseFactor = 0.0001;

unsigned short int mass_accuracy_unit = 0; 

float PeakWidthMinutes = 4.0;


typedef struct _intensity
{
	float fIntensity, RetTime;
} Intensity;

typedef std::vector <Intensity> IonIntensity;

typedef IonIntensity::iterator nIonIntensity;

/*
int compare_float (const void * a, const void * b)
{
	return ( *(float*)a - *(float*)b );
} */

/*
*   The generic idiom for comparing two numerical values a and b
*   for qsort looks as (a > b) - (a < b).
*
*/
int compare_float(const void * a, const void * b)
{
	float fa = *(const float*)a;

	float fb = *(const float*)b;

	return (fa > fb) - (fa < fb);
}



void DetectAndIntegratePeak::Close_mzML_file()
{
	mzml->CloseFiles();
}

/*
*
* read isotopes from an mzML file
* and store them
* The peak is characterized by its scan number (elution time)
* and mass-to-charge ratio of the monoisotopic precursor
* Intensity holds all intensities in the precursor's elution profile
* Noise holds the intensity in +- 5.0 m/z interval. It is used to
* determine if the precursor's intensity falls below the noise level
* at which point the while loop is exited.
* mzml pointer is a public pointer in the DetectAndIntegratePeak class
* is defined and initiated in the class initiation.
* The integration results are returned in array <double, 2> ^Isotopes.
*/

void DetectAndIntegratePeak::ReadIsotopes(int nScan, double dMonoPrecur, int nCharge, array <double, 2> ^Isotopes)
{
	unsigned short mslevel = 2;  //2 is MS2, 1 is MS1 - full scan

	int i, i0, i1, i2, i3, i4, i5;

	int nCurrentScan, nPeakIndex = 0;

	int nStartElution, nEndElution;

	double dtemp, dtemp1, dtemp2, dtemp3, dtemp4, dtemp5;

	double d0, d1, d2, d3, d4, d5;

	double dSeqMZ, dPeak = 0.0, dtempPeak=0.0;

	double fmedian, dMonoArea, dFirstArea, dSecondArea, dThirdArea;

	double dFourthArea, dFifthArea;

	double dFirstIsotopeShift, dSecondIsotopeShift, dThirdIsotopeShift;

	double dFourthIsotopeShift, dFifthIsotopeShift, dBrakeMass;

	double dThreeElutionAverage, dPreviousElution, dTwoScansBefore;

	List <float> ^Intensity = gcnew List <float>;

	List <float> ^Noise;

	bool bStartElution = false, bEndElution = false;

	//set the mass accuracy

	//if(0 == mass_accuracy_unit)
	//{
	//	dMassWindow = 0.000001*dMassAccuracy*(dMonoPrecur*nCharge - nCharge*mProton);
	//}
	//else if(1 == mass_accuracy_unit)
	//{
	//	dMassWindow = dMassAccuracy;
	//}

	//printf("MassAccuracy = %10.5f\n", dMassWindow);
	nCurrentScan = nScan;

	dSeqMZ = dMonoPrecur;

	dMonoArea = dFirstArea = dSecondArea = dThirdArea = dFourthArea = dFifthArea = 0.0;

	//set this values equal to zero only one time
	//they are the sum of isotope intensities before
	//and after the peak detection

	d0 = d1 = d2 = d3 = d4 = d5 = fmedian = 0.0;

	dFirstIsotopeShift = dSeqMZ + deltaC13/(double)nCharge; 

	dSecondIsotopeShift = dSeqMZ + 2.*deltaC13/(double)nCharge;

	dThirdIsotopeShift = dSeqMZ + 3.*deltaC13/(double)nCharge;

	dFourthIsotopeShift = dSeqMZ + 4.*deltaC13/(double)nCharge; 

	dFifthIsotopeShift  = dSeqMZ + 5.*deltaC13/(double)nCharge;

	dBrakeMass          = dFifthIsotopeShift + 2*dMassWindow;  // for mass values larger than this stop integrattion routine

	dPreviousElution  = dTwoScansBefore = 0.0;

	//read hundred scans back
	//while(nCurrentScan >= nScan - 100 && nCurrentScan >= 1)
	while(nCurrentScan >= 1)
	{
		mzml->ReadIndexedScan(nCurrentScan, &mslevel);

		if(1 == mslevel)
		{
			dtemp = dtemp1 = dtemp2 = dtemp3 = dtemp4 = dtemp5 = fmedian = 0.0; 

			i0 = i1 = i2 = i3 = i4 = i5 = 0;

			Noise = gcnew List <float>;

			// for every MS1 file compute the isotope peak areas
			// in the m/z domain for everyone of the isotopes
			dtempPeak = 0.0;

			for(i=0; i < mzml->Spectrum->Length/2; i++)
			{
				if(Math::Abs(mzml->Spectrum[0, i] - dSeqMZ) <= 5.0 && mzml->Spectrum[1, i] > 1.)
				{
					Noise->Add(mzml->Spectrum[1,i]);
				}

				//do nothing for small m/z values
				if((dSeqMZ - mzml->Spectrum[0, i]) > dMassWindow )
				{
					continue;
				}

				//monoisotope
				if(Math::Abs(mzml->Spectrum[0, i] - dSeqMZ) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp)
					{
						dtemp = mzml->Spectrum[1, i];

						i0 = i;
					}

					d0 += mzml->Spectrum[1, i]; 

					dtempPeak = dtempPeak + mzml->Spectrum[1, i];
				} 

				//first isotope
				//else if(Math::Abs(mzml->Spectrum[0, i] - (dSeqMZ + mProton/(double)nCharge) ) <= dMassAccuracy)
				else if(Math::Abs(mzml->Spectrum[0, i] - dFirstIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp1)
					{
						dtemp1 = mzml->Spectrum[1, i];

						i1 = i;
					}

					d1 += mzml->Spectrum[1, i];
				}
				//second isotope
				//else if( Math::Abs(mzml->Spectrum[0, i] - (dSeqMZ + 2.*mProton/(double)nCharge) ) <= dMassAccuracy)
				else if( Math::Abs(mzml->Spectrum[0, i] - dSecondIsotopeShift ) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp2)
					{
						dtemp2 = mzml->Spectrum[1, i];

						i2 = i;
					}

					d2 += mzml->Spectrum[1, i];
				}
				//thirdd isotope
				//else if( Math::Abs(mzml->Spectrum[0, i] - (dSeqMZ + 3.*mProton/(double)nCharge) ) <= dMassAccuracy)
				else if( Math::Abs(mzml->Spectrum[0, i] - dThirdIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp3)
					{
						dtemp3 = mzml->Spectrum[1, i];

						i3 = i;
					}

					d3 += mzml->Spectrum[1, i];

				}
				//fourth isotope
				//else if( Math::Abs(mzml->Spectrum[0, i] - (dSeqMZ + 4.*mProton/(double)nCharge) ) <= dMassAccuracy)
				else if( Math::Abs(mzml->Spectrum[0, i] - dFourthIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp4)
					{
						dtemp4 = mzml->Spectrum[1, i];

						i4 = i;
					}

					d4 += mzml->Spectrum[1, i];
				}
				//fifth isotope, it will be present only if the 4th isotope is nonzero
				//else if( Math::Abs(mzml->Spectrum[0, i] - (dSeqMZ + 5.*mProton/(double)nCharge) ) <= dMassAccuracy)
				else if( Math::Abs(mzml->Spectrum[0, i] - dFifthIsotopeShift) <= dMassWindow &&
					i4 > 0)
				{
					if(mzml->Spectrum[1, i] > dtemp5)
					{
						dtemp5 = mzml->Spectrum[1, i];

						i5 = i;
					}

					d5 += mzml->Spectrum[1, i];
				}
				//break at exactly 5. Th, because the noise is determined in +- 5 Th interval
				else if(mzml->Spectrum[0, i] - dSeqMZ >= dBrakeMass)
				{
					break;
				}
			}  //for(i=0; i < mzml->Spectrum->Length/2; i++)

			if(Noise->Count > 0)
			{
				fmedian =  QSelect_Float(Noise);
			}

			if(mzml->Spectrum[1,i0] >= fmedian)
			{
				if(i0 > 0)
				{
					Intensity->Add(mzml->Spectrum[1,i0]);

					dMonoArea += mzml->Spectrum[1,i0];
				}

				if(i1 > 0) dFirstArea += mzml->Spectrum[1,i1];

				if(i2 > 0) dSecondArea += mzml->Spectrum[1,i2];

				if(i3 > 0) dThirdArea += mzml->Spectrum[1,i3];

				if(i4 > 0) dFourthArea += mzml->Spectrum[1,i4];

				//if the 4th isotope is small, the fift should be non-zero
				if(i5 > 0 && i4 > 0) dFifthArea += mzml->Spectrum[1,i5];

			}
			else if(mzml->Spectrum[1,i0] < fmedian * dNoiseFactor)
			{
				bStartElution = true;


				/*printf("Scan NUmber at brake = %d, fmedian = %10.5f mzML = %10.5f, i0 = %d\n",
				nCurrentScan, fmedian, mzml->Spectrum[1,i0], i0);  */ 

				//break;     //breaks out of the while loop
			}


			if(dtempPeak > dPeak)
			{

				dPeak = dtempPeak;

				nPeakIndex = nCurrentScan;
			}

			dThreeElutionAverage =  (dtempPeak + dPreviousElution + dTwoScansBefore)/3.0;

			dTwoScansBefore = dPreviousElution;

			dPreviousElution = dtempPeak;

			//if(dtempPeak < 0.001*dPeak)

			if(dThreeElutionAverage < 0.001*dPeak)
			{
				//printf("Braking here\n");

				if(nScan - nCurrentScan> 5)
				{
					break;
				}
			}

			//printf("nCurrent: = %d\n", nCurrentScan);
		}// (1 == mslevel)

		nCurrentScan--;

		if(0 == nCurrentScan)
		{
			printf("Reached the First Scan, error\n");

			break;

			//exit(1);
		}
	}

	nStartElution = nCurrentScan;

	//read hundred scans forward

	nCurrentScan = nScan;

	dPreviousElution  = dTwoScansBefore = 0.0;

	//while(nCurrentScan <= nScan + 100 && nCurrentScan < mzml->nTotalSpectra - 1)
	while(nCurrentScan < mzml->nTotalSpectra - 1)
	{
		mzml->ReadIndexedScan(nCurrentScan, &mslevel);

		if(1 == mslevel)
		{
			dtemp = dtemp1 = dtemp2 = dtemp3 = dtemp4 = dtemp5 = fmedian = 0.0; 

			i0 = i1 = i2 = i3 = i4 = i5 = 0;

			Noise = gcnew List <float>;

			//will keep peak information
			dtempPeak = 0.0;

			for(i=0; i < mzml->Spectrum->Length/2; i++)
			{
				//monoisotope
				if(Math::Abs(mzml->Spectrum[0, i] - dSeqMZ) <= 5.0 && mzml->Spectrum[1, i] > 1.)
				{
					Noise->Add(mzml->Spectrum[1,i]);
				}


				//do nothing for small m/z values
				if((dSeqMZ - mzml->Spectrum[0, i]) > dMassWindow )
				{
					continue;
				}

				//monoisotope
				if(Math::Abs(mzml->Spectrum[0, i] - dSeqMZ) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp)
					{
						dtemp = mzml->Spectrum[1, i];

						i0 = i;
					}

					d0 += mzml->Spectrum[1, i];

					dtempPeak = dtempPeak + mzml->Spectrum[1, i];
				}
				//first isotope
				else if( Math::Abs(mzml->Spectrum[0, i] - dFirstIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp1)
					{
						dtemp1 = mzml->Spectrum[1, i];

						i1 = i;
					}

					d1 += mzml->Spectrum[1, i];
				}
				//second isotope
				else if( Math::Abs(mzml->Spectrum[0, i] - dSecondIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp2)
					{
						dtemp2 = mzml->Spectrum[1, i];

						i2 = i;
					}

					d2 += mzml->Spectrum[1, i];
				}
				//third isotope
				else if( Math::Abs(mzml->Spectrum[0, i] - dThirdIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp3)
					{
						dtemp3 = mzml->Spectrum[1, i];

						i3 = i;
					}

					d3 += mzml->Spectrum[1, i];
				}
				//fourth isotope
				else if( Math::Abs(mzml->Spectrum[0, i] - dFourthIsotopeShift) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtemp4)
					{
						dtemp4 = mzml->Spectrum[1, i];

						i4 = i;
					}

					d4 += mzml->Spectrum[1, i];
				}
				//fifth isotope, it will be present only if the 4th isotope is non-zero
				else if( Math::Abs(mzml->Spectrum[0, i] - dFifthIsotopeShift) <= dMassWindow &&
					i4 > 0)
				{
					if(mzml->Spectrum[1, i] > dtemp5)
					{
						dtemp5 = mzml->Spectrum[1, i];

						i5 = i;
					}

					d5 += mzml->Spectrum[1, i];
				}
				//break at exactly 5. Th, because the noise is determined in +- 5 Th interval
				else if(mzml->Spectrum[0, i] - dSeqMZ >=  dBrakeMass)
				{
					break;
				}
			}

			if(Noise->Count > 0)
			{
				fmedian =  QSelect_Float(Noise);
			}

			if(mzml->Spectrum[1,i0] >= fmedian)
			{
				if(i0 > 0)
				{
					Intensity->Add(mzml->Spectrum[1,i0]);

					dMonoArea += mzml->Spectrum[1,i0];
				}

				if(i1 > 0) dFirstArea += mzml->Spectrum[1,i1];

				if(i2 > 0) dSecondArea += mzml->Spectrum[1,i2];

				if(i3 > 0) dThirdArea += mzml->Spectrum[1,i3];

				if(i4 > 0) dFourthArea += mzml->Spectrum[1,i4];

				if(i5 > 0 && i4 > 0) 
				{
					dFifthArea += mzml->Spectrum[1,i5];
				}
			}


			if(mzml->Spectrum[1,i0] < fmedian *dNoiseFactor)
			{
				bEndElution = true;

				//break;       //breaks out of the while loop
			}

			//store the peak or
			// exit if the peak is too small.

			if(dtempPeak > dPeak)
			{

				dPeak = dtempPeak;

				nPeakIndex = nCurrentScan;
			}

			dThreeElutionAverage =  (dtempPeak + dPreviousElution + dTwoScansBefore)/3.0;

			dTwoScansBefore = dPreviousElution;

			dPreviousElution = dtempPeak;

			//if(dtempPeak < 0.001*dPeak)
			if(dThreeElutionAverage < 0.001*dPeak)
			{
				//printf("%10.5f %10.5f\n", dtempPeak, dPeak);
				if(nCurrentScan - nScan > 5)
				{
					break;
				}
			}
		}// (1 == mslevel)

		nCurrentScan++;

		if(mzml->nTotalSpectra - 1 <= nCurrentScan)
		{
			printf("Reached the Last Scan, error\n");

			break;

			//exit(1);
		}

		/*if(bEndElution)
		{
		break;
		}*/
	}

	nEndElution = nCurrentScan;

	delete(Intensity); delete(Noise);

	/*printf("nS: %d m/z: %8.5f Isotopes: %10.f %10.f %10.f %10.f\n", 
	nScan, dMonoPrecur, dMonoArea, dFirstArea, dSecondArea, dThirdArea);*/

	/*printf("\tm/z: %8.5f Isotopes: %10.f %10.f %10.f %10.f %10.f %10.f\n", 
	dMonoPrecur, d0, d1, d2, d3, d4, d5);*/

	printf("dMonoPrecurosr = %10.5f nScan = %d, dStartScan = %d dEndScan = %d, Diff = %d, dMassWind = %10.5f\n", 
		dMonoPrecur, nScan, nStartElution, nEndElution, (nEndElution - nStartElution), dMassWindow);

	fpPeakOutPut->Write("dMonoPrecurosr = {0:10.5f} nScan = {1:%d}, dStartScan = {2} dEndScan = {3}, Diff = {4}, dMassWind = {5}\n",
		dMonoPrecur, nScan, nStartElution, nEndElution, (nEndElution - nStartElution), dMassWindow);

	Isotopes[1, 0] = d0; Isotopes[1, 1] = d1; Isotopes[1, 2] = d2; Isotopes[1, 3] = d3; Isotopes[1, 4] = d4; Isotopes[1, 5] = d5;


	/*if(nEndElution - nStartElution < 50 || nScan == 41749) 
	{
	int iil;

	printf("Elution length less than 50 scans ...\n");

	printf("Start = %d End Elution = %d\n", nStartElution, nEndElution);

	printf("The Peak Scan Number is %d\n", nPeakIndex);
	scanf("%d", &iil);
	}*/


	return;
}


/*
*
* integrate isotopes from a spectrum in the m/z domain.
* The spectrum is assumed to provide full scan peaks.
* Also in the input, is the dMonoPrecursor - monoisotopic
* precursor mass. The algorithm will simply move within the
* precursor mass tolerance (a class property), find the peak apex,
* and integrate it.
* 
*/

void DetectAndIntegratePeak::ReadIsotopes(double dMonoPrecur, int nCharge, array<double, 2> ^ Spectrum, array <double, 2> ^Isotopes)
{
	int i, i0, i1, i2, i3, i4, i5;

	double dtemp, dtemp1, dtemp2, dtemp3, dtemp4, dtemp5;

	double d0, d1, d2, d3, d4, d5;

	double dSeqMZ, dPeak = 0.0, dtempPeak=0.0;

	double  dMonoArea, dFirstArea, dSecondArea, dThirdArea;

	double dFourthArea, dFifthArea;

	double dFirstIsotopeShift, dSecondIsotopeShift, dThirdIsotopeShift;

	double dFourthIsotopeShift, dFifthIsotopeShift, dBrakeMass;

	double  dPreviousElution, dTwoScansBefore;


	//set the mass accuracy

	//if(0 == mass_accuracy_unit)
	//{
	//	dMassWindow = 0.000001*dMassAccuracy*(dMonoPrecur*nCharge - nCharge*mProton);
	//}
	//else if(1 == mass_accuracy_unit)
	//{
	//	dMassWindow = dMassAccuracy;
	//}

	//printf("MassAccuracy = %10.5f\n", dMassWindow);

	dSeqMZ = dMonoPrecur;

	dMonoArea = dFirstArea = dSecondArea = dThirdArea = dFourthArea = dFifthArea = 0.0;

	//set this values equal to zero only one time
	//they are the sum of isotope intensities before
	//and after the peak detection

	d0 = d1 = d2 = d3 = d4 = d5 = 0.0;

	dtemp = dtemp1 = dtemp2 = dtemp3 = dtemp4 = dtemp5 = 0;

	dFirstIsotopeShift = dSeqMZ + deltaC13/(double)nCharge; 

	dSecondIsotopeShift = dSeqMZ + 2.*deltaC13/(double)nCharge;

	dThirdIsotopeShift = dSeqMZ + 3.*deltaC13/(double)nCharge;

	dFourthIsotopeShift = dSeqMZ + 4.*deltaC13/(double)nCharge; 

	dFifthIsotopeShift  = dSeqMZ + 5.*deltaC13/(double)nCharge;

	dPreviousElution  = dTwoScansBefore = 0.0;

	dBrakeMass          = dFifthIsotopeShift + 2*dMassWindow;  // for mass values larger than this stop integrattion routine

	for(i=0; i < Spectrum->Length/2; i++)
	{

		//do nothing for small m/z values
		if((dSeqMZ - Spectrum[0, i]) > dMassWindow )
		{
			continue;
		}

		//monoisotope
		if(Math::Abs(Spectrum[0, i] - dSeqMZ) <= dMassWindow)
		{
			if(Spectrum[1, i] > dtemp)
			{
				dtemp = Spectrum[1, i];

				i0 = i;
			}

			d0 += Spectrum[1, i];

			dtempPeak = dtempPeak + Spectrum[1, i];
		}
		//first isotope
		else if( Math::Abs(Spectrum[0, i] - dFirstIsotopeShift) <= dMassWindow)
		{
			if(Spectrum[1, i] > dtemp1)
			{
				dtemp1 = mzml->Spectrum[1, i];

				i1 = i;
			}

			d1 += Spectrum[1, i];
		}
		//second isotope
		else if( Math::Abs(Spectrum[0, i] - dSecondIsotopeShift) <= dMassWindow)
		{
			if(Spectrum[1, i] > dtemp2)
			{
				dtemp2 = Spectrum[1, i];

				i2 = i;
			}

			d2 += Spectrum[1, i];
		}
		//third isotope
		else if( Math::Abs(Spectrum[0, i] - dThirdIsotopeShift) <= dMassWindow)
		{
			if(Spectrum[1, i] > dtemp3)
			{
				dtemp3 = mzml->Spectrum[1, i];

				i3 = i;
			}

			d3 += Spectrum[1, i];
		}
		//fourth isotope
		else if( Math::Abs(Spectrum[0, i] - dFourthIsotopeShift) <= dMassWindow)
		{
			if(Spectrum[1, i] > dtemp4)
			{
				dtemp4 = mzml->Spectrum[1, i];

				i4 = i;
			}

			d4 += Spectrum[1, i];
		}
		//fifth isotope, it will be present only if the 4th isotope is non-zero
		else if( Math::Abs(Spectrum[0, i] - dFifthIsotopeShift) <= dMassWindow &&
			i4 > 0)
		{
			if(Spectrum[1, i] > dtemp5)
			{
				dtemp5 = Spectrum[1, i];

				i5 = i;
			}

			d5 += mzml->Spectrum[1, i];
		}
		//break at exactly 5. Th, because the noise is determined in +- 5 Th interval
		else if(Spectrum[0, i] - dSeqMZ >=  dBrakeMass)
		{
			break;
		}
	}


	return ;
}

/*
*
*  The method will determine the peak apex based on the maximum peak 
*  height information only. It uses monisotopic m/z of the peak,
*  mass accuracy, and start scan information. If iDirection is -1,
*  the direction is to the left of the current Scan, if +1 it is to
*  the right of the current Scan;
*/

int DetectAndIntegratePeak::BasePeakApex(int nScan, double dMonoPrecur, int nCharge, int iDirection)
{
	unsigned short mslevel = 2;  //2 is MS2, 1 is MS1 - full scan

	int i, iMaxPeakHeight, nFullScanReturn;

	int nCurrentScan;

	double dMaxPeakHeight;

	float fRetTime, fPeakRetTime, fRetTimeMS2;

	double dSeqMZ, dPeak, dtempPeak;

	double  dBrakeMass;

	bool bStartElution = false, bEndElution = false;


	IonIntensity PeakProfile;

	Intensity currentPeak;

	//set the mass accuracy

	//if(0 == mass_accuracy_unit)
	//{
	//	dMassWindow = 0.000001*dMassAccuracy*(dMonoPrecur*nCharge - nCharge*mProton);
	//}
	//else if(1 == mass_accuracy_unit)
	//{
	//	dMassWindow = dMassAccuracy;
	//}

	//printf("MassAccuracy = %10.5f\n", dMassWindow);

	fPeakRetTime = fPeakRetTime = fRetTimeMS2 = 0.;

	dtempPeak = dPeak = 0.0;

	nCurrentScan = nScan;

	dSeqMZ = dMonoPrecur;

	dBrakeMass = dSeqMZ + 2*dMassWindow;  // for mass values larger than this stop integrattion routine

	dMaxPeakHeight = 0.0;

	iMaxPeakHeight = 0;

	mzml->ReadIndexedScan(nScan, &mslevel);

	fRetTimeMS2 = mzml->fCurrentRetTime;

	//read hundred scans back
	//while(nCurrentScan >= nScan - 100 && nCurrentScan >= 1)
	while(nCurrentScan >= 1)
	{
		nFullScanReturn = mzml->ReadIndexedFullScan(nCurrentScan, &mslevel);

		//do nothing if it is an MS2
		if(0 == nFullScanReturn)
		{
			nCurrentScan = nCurrentScan + iDirection;

			continue;
		}

		if(1 == mslevel)
		{ 
			fRetTime = mzml->fCurrentRetTime;

			// for every MS1 file compute the isotope peak areas
			// in the m/z domain for everyone of the isotopes
			dtempPeak = 0.0;

			for(i=0; i < mzml->Spectrum->Length/2; i++)
			{

				//low intensity peaks - ignore them.
				if(mzml->Spectrum[1, i] < 2.0)
				{
					continue;
				}
				//do nothing for small m/z values
				if((dSeqMZ - mzml->Spectrum[0, i]) > dMassWindow )
				{
					continue;
				}

				//monoisotope
				if(Math::Abs(mzml->Spectrum[0, i] - dSeqMZ) <= dMassWindow)
				{
					if(mzml->Spectrum[1, i] > dtempPeak)
					{
						dtempPeak = mzml->Spectrum[1, i]; 
					}
				} 
				else if(mzml->Spectrum[0, i] - dSeqMZ >= dBrakeMass)
				{
					break;
				}
			}  //for(i=0; i < mzml->Spectrum->Length/2; i++)

			currentPeak.RetTime = fRetTime;

			currentPeak.fIntensity = dtempPeak;

			PeakProfile.push_back(currentPeak);

			if(dtempPeak > dMaxPeakHeight)
			{

				dMaxPeakHeight = dtempPeak;

				iMaxPeakHeight = nCurrentScan;

				fPeakRetTime = fRetTime;
			}

			//printf("nCurrent: = %d\n", nCurrentScan);
		}// (1 == mslevel)

		if(1 == nCurrentScan)
		{
			printf("Reached the First Scan, error\n");

			break;

			//exit(1);
		}

		if(abs(fRetTime - fRetTimeMS2) > PeakWidthMinutes)
		{
			break;
		}

		nCurrentScan = nCurrentScan + iDirection;
	}


	printf("mz/ = %10.5f Max RetTime = %10.5f MS2time = %10.5f Peak = %d\n", dSeqMZ, fPeakRetTime, fRetTimeMS2,iMaxPeakHeight);

	for(i = 0; i < PeakProfile.size(); i++)
	{
		printf("Peak = %10.5f %10.5f\n", PeakProfile[i].RetTime, PeakProfile[i].fIntensity);
	}


	PeakProfile.empty();

	return nCurrentScan;
}



/*
* A recursive algorithm that returns a median
* of a list of float numbers
*/
float DetectAndIntegratePeak::QSelect_Float(List <float> ^lArray)
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
*
*  calls chromatogram builder from XmlMzML
*  the data from XmlMzML (FullChromatogramList)
*  then becomes available to DetedAndIntegrate class
*   if iMethod = 0 then build the Chromatogram using
*   Sequential Read 
*   if iMethod = 1 build the chromatogram using 
*    indexed read.
*  THE SEQUENTIAL READ DOES NOT WORK
*
*/
void  DetectAndIntegratePeak::BuildChromatogram(short unsigned int iMethod)
{

	if(0 == iMethod)
	{
		mzml->SequentialChromatogramBuild();
	}
	else if(1 == iMethod)
	{
		mzml->IndexedChromatogramBuild();
	}

	/*for(int i = 0; i < mzml->FullScanChromatogram->Count; i++)
	{
	printf("Time point %10.5f and Scan = %d\n", mzml->FullScanChromatogram[i]->RetTime, mzml->FullScanChromatogram[i]->ScanNumber);
	}

	exit (1);*/
}



/*
*
*  determines peaks and integrates them. Uses already 
*  built chromatograms.
* The peak is characterized by its scan number, elution time,
* and mass-to-charge ratio of the monoisotopic precursor
* Intensity holds all intensities in the precursor's elution profile
* Noise holds the intensity in +- 5.0 m/z interval. It is used to
* determine if the precursor's intensity falls below the noise level
* at which point the while loop is exited.
* mzml pointer is a public pointer in the DetectAndIntegratePeak class
* is defined and initiated in the class initiation.
* The integration results are returned in array <double, 2> ^Isotopes.
*/

void DetectAndIntegratePeak::IntegratePeakTimeMz(float RetTime, int nScan, double dMonoPrecur, int nCharge, float *fStartElution,
	float *fEndElution, array <double, 2> ^Isotopes,array <double, 1> ^w)
{
	unsigned short mslevel = 2;  //2 is MS2, 1 is MS1 - full scan

	int i, j;

	//Added by JA 2016.10.21
	bool b1 = true;
	bool b2 = true;

	double d0, d1, d2, d3, d4, d5;

	double dSeqMZ;

	double dFirstIsotopeShift, dSecondIsotopeShift, dThirdIsotopeShift;

	double dFourthIsotopeShift, dFifthIsotopeShift, dBrakeMass;

	List <float> ^Intensity = gcnew List <float>;

	float fStartTime, fEndTime;

	//Added by JA 2016.09.08
	int jPBStart, jPBEnd;

	//Added by JA 2016.10.21
	int jMaxPeak;
	float mzMaxPeak = 0.0f;

	float mzPBStart;
	float mzPBEnd;

	//Storing left and right peak widths
	float LW;
	float RW;

	//Storing Retention Time Peak
	float RetentionTimePeak;

	dSeqMZ = dMonoPrecur;

	//set this values equal to zero only one time
	//they are the sum of isotope intensities before
	//and after the peak detection

	d0 = d1 = d2 = d3 = d4 = d5 = 0.0;


	dFirstIsotopeShift = dSeqMZ + deltaC13/(double)nCharge; 

	dSecondIsotopeShift = dSeqMZ + 2.*deltaC13/(double)nCharge;

	dThirdIsotopeShift = dSeqMZ + 3.*deltaC13/(double)nCharge;

	dFourthIsotopeShift = dSeqMZ + 4.*deltaC13/(double)nCharge; 

	dFifthIsotopeShift  = dSeqMZ + 5.*deltaC13/(double)nCharge;

	// determine the start and end of the Peak Elution
	// passes the candidate start and end position from identifications of MS2

	fStartTime = (*fStartElution)/60.0;

	fEndTime   = (*fEndElution)/60.0;

	jPBStart = 0;

	jPBEnd = mzml->FullScanChromatogram[0]->moverz->Length;

	mzPBStart = 0.0f;

	mzPBEnd = 0.0f;

	LW = 0.0f;

	RW = 0.0f;


	//Passing nScan in order to check whether it lies in the Retention Time Window
	//w[10] = nScan;

	//FindPeakStartEnd(RetTime/60., dMonoPrecur,  nCharge, &fStartTime, &fEndTime, &jPBStart, &jPBEnd, &mzPBStart, &mzPBEnd, &LW, &RW,&RetentionTimePeak,w);

	FindPeakStartEnd_With_MassAccuracy_Charge(RetTime / 60., dMonoPrecur, nCharge, &fStartTime, &fEndTime, &jPBStart, &jPBEnd, &mzPBStart, &mzPBEnd, &LW, &RW, &RetentionTimePeak, w);
	
	if(dMassWindow >2 || dMassWindow <= 0)//||dMassWindow==setMassWindow)
	{
		dMassWindow = setMassWindow;
	}
	else if(LW <= RW && LW >= 0)
	{
		dMassWindow = LW;
	}
	else if(RW >= 0) 
	{
		dMassWindow = RW;
	}
	else
	{
		dMassWindow = setMassWindow;
	}
	//dMassWindow = Math::Abs(mzPBEnd - mzPBStart)/2;

	if (dMassWindow < 0.0000001)
	{
		printf("Defaulfing to Mass Width\n");

		dMassWindow = setMassWindow;;
	}
	
	dBrakeMass = dFifthIsotopeShift + 1.5*dMassWindow;


	mzMaxPeak = w[9];
	jMaxPeak = w[10];

	if(mzMaxPeak==0 || jMaxPeak==-1)
	{
		mzMaxPeak = dMonoPrecur;
	}

	//integrate the peaks between the start and end elution times: fStartTime, fEndTime

	if(fStartTime < 0. || fEndTime < 0.0)
	{
		fpPeakOutPut->Write("Untedectable: dMonoPrecursor = {0} nScan = {1}, StartELution = {2} dEndElution = {3}, ElutWidth = {4}, dMassWind = {5}\n",			dMonoPrecur, nScan, fStartTime, fEndTime, (fEndTime - fStartTime), dMassWindow);

		Isotopes[1, 0] = 0.0; Isotopes[1, 1] = 0.0; Isotopes[1, 2] = 0.0; Isotopes[1, 3] = 0.0; Isotopes[1, 4] = 0.0; Isotopes[1, 5] = 0.0;

		return;
	}

	//***********************
	//Peak integration is done in the following for loop
	for(i=0; i < mzml->FullScanChromatogram->Count; i++)
	{
		//Forcing start and end elution times to specific values for a given peptide
		if(mzml->FullScanChromatogram[i]->RetTime < fStartTime)
		{
			continue;
		}

		if(mzml->FullScanChromatogram[i]->RetTime >= fStartTime &&
			mzml->FullScanChromatogram[i]->RetTime <= fEndTime)
		{

			//Add to the left
			for(j=0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
			{
				//do nothing for small m/z values
				if((dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) > dMassWindow )
				{
					continue;
				}
				//monoisotope
				else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) <= dMassWindow)
				{
					d0 = d0 + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(Math::Abs(mzml->FullScanChromatogram[i]->moverz[j] - dFirstIsotopeShift) <= dMassWindow)
				{
					d1 = d1 + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(Math::Abs(mzml->FullScanChromatogram[i]->moverz[j] - dSecondIsotopeShift) <= dMassWindow)
				{
					d2 = d2 + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(Math::Abs(mzml->FullScanChromatogram[i]->moverz[j] - dThirdIsotopeShift) <= dMassWindow)
				{
					d3 = d3 + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(Math::Abs(mzml->FullScanChromatogram[i]->moverz[j] - dFourthIsotopeShift) <= dMassWindow)
				{
					d4 = d4 + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(Math::Abs(mzml->FullScanChromatogram[i]->moverz[j] - dFifthIsotopeShift) <= dMassWindow)
				{
					d5 = d5 + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(mzml->FullScanChromatogram[i]->moverz[j] > dBrakeMass)
				{
					break;
				}

			}//end for(j=0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
		}//end if(mzml->FullScanChromatogram[i]->RetTime >= fStartTime && mzml->FullScanChromatogram[i]->RetTime <= fEndTime) 
	}//end for(i=0; i < mzml->FullScanChromatogram->Count; i++)

	delete(Intensity); 

	*fStartElution = fStartTime;

	*fEndElution   = fEndTime;

	printf("dMonoPrecursor = %10.5f nScan = %d, StartELution = %10.5f dEndElution = %10.5f, ElutWidth = %10.5f, dMassWind = %10.5f\n", 
		dMonoPrecur, nScan, fStartTime, fEndTime, (fEndTime - fStartTime), dMassWindow);


	Isotopes[1, 0] = d0;Isotopes[1, 1] = d1; Isotopes[1, 2] = d2; Isotopes[1, 3] = d3; Isotopes[1, 4] = d4; Isotopes[1, 5] = d5;

	w[0] = dMassWindow;
	w[1] = LW;
	w[2] = RW;
	w[3] = RetentionTimePeak;

	return;
}


/*
*
* The method finds the peak position, start and end, given
* the monoisotopic mass, and retention time
*  the elutions times are assumed to be in minutes;
*  The Window that is used for elution time is 10 minutes!!!
*   It means a peak is searched for 10 minutes before and after
*   the MS2 has been triggered.
*/
void DetectAndIntegratePeak::FindPeakStartEnd_With_MassAccuracy_Charge(float RetTime, double dMonoPrecur, int nCharge, float *fStart, float *fEnd, int *jPBStart,
			int *jPBEnd, float *mzPBStart, float *mzPBEnd, float *LW, float *RW, float *RetentionTimePeak, array <double, 1> ^w)
{
	int i, j, k, nPeak, kPeak, increments, i1, kPeak1;

	float fHighEnd, fLowEnd, fPeakRetTime;

	float fEndTime, fStartTime, fThreshold;

	float fTimeWindow, fPeak, fSummedPeak, fSummedPeak0, dAreaUnderPeak;

	float fMaxLocalPeak;
	
	double dMoverZLocalPeak = 0., dtemp, dI0, dI1, dI2, dI3, dI4;

	double dtemp2;

	int iMaxLocalPeak, iMaxPeakMZ;

	fTimeWindow = fElutionTimeWindow;       //fElutionTimeWindow is a class variable, set at the time of declaration/initiation.

	if (false) //this is original start and end times
	{
		fHighEnd = RetTime + fTimeWindow;

		fLowEnd = RetTime - fTimeWindow;
	}
	else
	{
		fHighEnd = (*fEnd) + fTimeWindow;

		fLowEnd = (*fStart) - fTimeWindow;

		//printf("Initial Start %10.5f and End %10.5f Elution Times\n", fLowEnd, fHighEnd);
	}

	fPeak = fSummedPeak = 0.0;

	fEndTime = fStartTime = -1.;

	fThreshold = 0.01f;  kPeak = -1; kPeak1 = -1;

	//Added by JA 2016.09.08
	//Following variables store intensity values
	float fMaxPeak = 0.0f;
	float fPB = 0.0f;

	//Following variable stores incremental factor to dMassWindow
	double wincrement = 1.00;
	//Added by JA 2016.10.27
	double tempMassWindow = setMassWindow;//0.001;

										  //Following variables store m/z values
	float mzMaxPeak = 0.0f;
	float mzPBLow = 0.0f;
	float mzPBHigh = 0.0f;
	float mzPBLowT = 0.0f;
	float mzPBLowC = 0.0f;
	float mzPBHighT = 0.0f;
	float mzPBHighC = 0.0f;

	//Following variables store indices
	int jMaxPeak = 0;
	int jPBLow = 0;
	int jPBHigh = 0;
	int jPBLowT = 0;
	int jPBLowC = 0;
	int jPBHighT = 0;
	int jPBHighC = 0;
	//

	bool bWidthCalculated = true;

	RetTimeRange ^tempRTR;
	List <RetTimeRange ^> ^RTR;

	List <PeakIndex ^> ^PeakElution = gcnew List <PeakIndex ^>;

	PeakIndex ^currentPeak = gcnew PeakIndex;

	int iStartTime, iEndTime;

	double dUpperMass, dLowerMass;

	//Calculate Retention Time based Peak first with a masswindow that is incremented until we find a Retention Time Peak
	//wincrement = 0.001;


	dUpperMass = dMonoPrecur + MassAccuracyPPM * dMonoPrecur / 1000000.0;

	dLowerMass = dMonoPrecur - MassAccuracyPPM * dMonoPrecur / 1000000.0;

	iMaxPeakMZ = 0;

	for (i = 0; i < mzml->FullScanChromatogram->Count; i++)
	{
		if (mzml->FullScanChromatogram[i]->RetTime < fLowEnd)
		{
			continue;
		}

		if (mzml->FullScanChromatogram[i]->RetTime >= fLowEnd &&
			mzml->FullScanChromatogram[i]->RetTime <= fHighEnd)
		{
			fSummedPeak = 0.0;

			fMaxLocalPeak = 0.;  iMaxLocalPeak = 0; dMoverZLocalPeak = 0.0;


			for (j = 0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
			{
				//do nothing for small m/z values
				if ((dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) > setMassWindow)
				{
					continue;
				}
				//monoisotope
				else if (Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) <= setMassWindow)
				{
					fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[i]->Intensity[j];

					if (mzml->FullScanChromatogram[i]->Intensity[j] > fMaxLocalPeak)
					{
						
						fMaxLocalPeak = mzml->FullScanChromatogram[i]->Intensity[j];

						iMaxLocalPeak = j;

						dMoverZLocalPeak = mzml->FullScanChromatogram[i]->moverz[j];
					}
				}
				else if (mzml->FullScanChromatogram[i]->moverz[j] > dMonoPrecur + setMassWindow)
				{
					break;
				}
			}

			

			//if (fSummedPeak > fPeak)

			if(dMoverZLocalPeak >= dLowerMass && dMoverZLocalPeak <= dUpperMass)
			{


				/*fpPeakOutPut->Write("Trying this peak,  the peak = {0}, time = {1}\n",
					dMoverZLocalPeak, mzml->FullScanChromatogram[i]->RetTime);  */
				//added on 05/12/2017 RGS, Check to make sure that the peak is within the mass tolerance of the
				// precursor

				/*if (dMoverZLocalPeak > dMonoPrecur)
				{
					dtemp = dMoverZLocalPeak - dMonoPrecur;
				}
				else
				{
					dtemp = dMonoPrecur - dMoverZLocalPeak;
				} 

				dtemp = (dtemp / dMonoPrecur) * 1000000.;   // ppm mass accuracy

				*/
				//if the maximum of the peak is at required ppm mass accuracy from the theoretical
				// peptide mass, then the peak a candidate for this peak
				
				//if (dtemp <= MassAccuracyPPM)
				if(true)
				{
					
					//check the charge state for the peak

					//first determine the start of the peak, and make sure that the
					//iMaxLocalPeak is actually the Peak!!!

					int nStart, nEnd, iTempMax;

					nStart = nEnd = iTempMax = -1;

					dtemp2 = 0.0;

					dAreaUnderPeak = mzml->FullScanChromatogram[i]->Intensity[iMaxLocalPeak];

					//determine the start of the peak
					for (k = iMaxLocalPeak - 1; k >= 1;  k--)
					{
						dAreaUnderPeak = dAreaUnderPeak + mzml->FullScanChromatogram[i]->Intensity[k];

						if (mzml->FullScanChromatogram[i]->Intensity[k] > dtemp2)
						{
							dtemp2 = mzml->FullScanChromatogram[i]->Intensity[k];

							iTempMax = k;
						}

						if (mzml->FullScanChromatogram[i]->Intensity[k] < 10 ||
							mzml->FullScanChromatogram[i]->Intensity[k] < mzml->FullScanChromatogram[i]->Intensity[k-1] ||
							mzml->FullScanChromatogram[i]->moverz[j] < dMonoPrecur - setMassWindow)
						{

							nStart = k;

							break;
						}
					}

					//determine the end of the peak
					for (k = iMaxLocalPeak + 1; k < mzml->FullScanChromatogram[i]->moverz->Length; k++)
					{
						dAreaUnderPeak = dAreaUnderPeak + mzml->FullScanChromatogram[i]->Intensity[k];

						if (mzml->FullScanChromatogram[i]->Intensity[k] > dtemp2)
						{
							dtemp2 = mzml->FullScanChromatogram[i]->Intensity[k];

							iTempMax = k;
						}

						if (mzml->FullScanChromatogram[i]->Intensity[k] < 10 ||
							mzml->FullScanChromatogram[i]->Intensity[k] > mzml->FullScanChromatogram[i]->Intensity[k - 1] ||
							mzml->FullScanChromatogram[i]->moverz[j] > dMonoPrecur + setMassWindow)
						{
							nEnd = k;

							break;
						}
					}

					if (nEnd == -1 || nStart == -1)    //this does not seem to register often - may be a computational time waste
					{
						//fprintf(fpPeakOutPut, "This is not a Peak. It is a shoulder, m/z %15.7f  the peak = %15.7f, time = %f\n",
						//mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak],
						//mzml->FullScanChromatogram[i]->moverz[iTempMax], mzml->FullScanChromatogram[i]->RetTime);
						
						fpPeakOutPut->Write("This is not a Peak. It is a shoulder, m/z {0}  the peak = {1}, time = {2}\n", 
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak],
							mzml->FullScanChromatogram[i]->moverz[iTempMax], mzml->FullScanChromatogram[i]->RetTime);
					}

					if (dtemp2 > mzml->FullScanChromatogram[i]->Intensity[iMaxLocalPeak] && iTempMax > 0) // this does not seem to register often - may be a computational time waste
					{
						//fprintf(fpPeakOutPut, "SHOULDER CASE: the shoulder_m/z %15.7f  the peak = %15.7f, time = %f\n", mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak],
						//mzml->FullScanChromatogram[i]->moverz[iTempMax], mzml->FullScanChromatogram[i]->RetTime);	
						fpPeakOutPut->Write("SHOULDER CASE: the shoulder_m/z {0}  the peak = {1}, time = {2}\n", mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak],
							mzml->FullScanChromatogram[i]->moverz[iTempMax], mzml->FullScanChromatogram[i]->RetTime);
/*first continue*/		continue; /// very important!!! The peak position is not on the iMaxLocalPeak, this means this peaks shoulder is taller than valume of the true peak 
					}

					//second check:: go backward up to 1 moverz and see if this peak is a member of another
					// isotope envelop


					dI0 = dI1 = dI2 = dI3 = 0.0;

					for (k = iMaxLocalPeak; k >= 0; k--)
					{
						if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - 
							mzml->FullScanChromatogram[i]->moverz[k] - deltaC13/3.) <  setMassWindow)
						{
							dI0 = dI0 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - 
							mzml->FullScanChromatogram[i]->moverz[k] - deltaC13 / 2.) <  setMassWindow)
						{
							dI1 = dI1 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - 
							mzml->FullScanChromatogram[i]->moverz[k] - 2.*deltaC13/3.) <  setMassWindow)
						{
							dI2 = dI2 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - 
							mzml->FullScanChromatogram[i]->moverz[k] - deltaC13) <  setMassWindow)
						{
							dI3 = dI3 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if ((mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] -
							mzml->FullScanChromatogram[i]->moverz[k] - deltaC13) > setMassWindow )
						{
							break;
						}
					}

					if (fSummedPeak <= dI0 || fSummedPeak <= dI1)   //dI0 and dI1 are back peaks
					{
						/*fpPeakOutPut->Write("It is a part of a different isotope:  the precursor = {0}, time = {1}\n",
							dMonoPrecur, mzml->FullScanChromatogram[i]->RetTime);  */
				    
/* */	    			continue;     // this is a part of envolop of another isotope peak
					}

					//third check:: Determine if the peak has a correct charge state representation in its isotope distribution
					// go forward up to 1 moverz and record the +2 and +3 peaks, +4

					dI0 = dI1 = dI2 = dI3 = dI4 = 0.0;

					//printf("IMaxLog = %d, length = %d\n", iMaxLocalPeak, mzml->FullScanChromatogram[i]->moverz->Length);

					for (k = iMaxLocalPeak; k < mzml->FullScanChromatogram[i]->moverz->Length; k++)
					{
						if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[k] - 
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13 / 3.) <  setMassWindow)
						{
							dI3 = dI3 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[k] - 
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13 / 4.) <  setMassWindow)
						{
							dI4 = dI4 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[k] -
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13 / 2.) <  setMassWindow)
						{
							dI2 = dI2 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
					/*	else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[k] - 
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - 2.*deltaC13 / 3.) <  setMassWindow)
						{
							dI2 = dI2 + mzml->FullScanChromatogram[i]->Intensity[k];
						}
						else if (Math::Abs(mzml->FullScanChromatogram[i]->moverz[k] -
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] -  deltaC13) <  setMassWindow)
						{
							dI3 = dI3 + mzml->FullScanChromatogram[i]->Intensity[k];
						}   */
						else if ((mzml->FullScanChromatogram[i]->moverz[k] - 
							mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13) > setMassWindow)
						{
							break;
						}
					}

					if (nCharge == 2 && dI2 <= 0.1* fSummedPeak) 
					{
						//fprintf(fpPeakOutPut, "Not Correct Charge State CASE:  Peptide Charge = %d,   the precursor = %10.5f, time = %f\n",
							fpPeakOutPut->Write("Missing Charge State:  Peptide Charge = {0},   the precursor = {1}, time = {2}\n",
							nCharge, //(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13 / 3.),
								dMonoPrecur, mzml->FullScanChromatogram[i]->RetTime);
						
/* */	    			continue;     // this is a part of envolop of another isotope peak
					}
					else if (nCharge == 3 && dI3 <= 0.1* fSummedPeak)
					{
					//	fprintf(fpPeakOutPut, "Not Correct Charge State CASE:  Peptide Charge = %d the precursor = %10.5f, time = %f\n",
						fpPeakOutPut->Write("Missing first Isotope for Charge State:  Peptide Charge = {0} the precursor = {1}, time = {2}\n",
							nCharge, //(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13 / 3.),
							dMonoPrecur, mzml->FullScanChromatogram[i]->RetTime);

/* */	    			continue;     // this is a part of envolop of another isotope peak
					}
					else if (nCharge == 4 && dI4 <= 0.1* fSummedPeak)
					{
						//	fprintf(fpPeakOutPut, "Not Correct Charge State CASE:  Peptide Charge = %d the precursor = %10.5f, time = %f\n",
						fpPeakOutPut->Write("Missing First Isotope Charge State:  Peptide Charge = {0} the precursor = {1}, time = {2}\n",
							nCharge, //(mzml->FullScanChromatogram[i]->moverz[iMaxLocalPeak] - deltaC13 / 3.),
							dMonoPrecur, mzml->FullScanChromatogram[i]->RetTime);

						/* */	    			continue;     // this is a part of envolop of another isotope peak
					}

					// all filters are satisfied. Update the peak information with the new point

					currentPeak = gcnew PeakIndex;

					currentPeak->PeakArea = dAreaUnderPeak;

					currentPeak->iLeftEnd = nStart;

					currentPeak->iRightEnd = nEnd;

					currentPeak->iMaxPosition = iMaxLocalPeak;

					currentPeak->RetTime = mzml->FullScanChromatogram[i]->RetTime;

					currentPeak->FullScanNumber = mzml->FullScanChromatogram[i]->ScanNumber;

					currentPeak->iChromNumber = i;

					PeakElution->Add(currentPeak);

					delete(currentPeak);

				/* fpPeakOutPut->Write("Entering this peak,  the peak = {0}, time = {1} {2}\n",
						 PeakElution[PeakElution->Count - 1]->RetTime,
						PeakElution[PeakElution->Count -1]->FullScanNumber,
						mzml->FullScanChromatogram[i]->moverz[PeakElution[PeakElution->Count - 1]->iMaxPosition]);  */

					if (fSummedPeak > fPeak)
					{
						fPeak = fSummedPeak;
						nPeak = mzml->FullScanChromatogram[i]->ScanNumber;
						fPeakRetTime = mzml->FullScanChromatogram[i]->RetTime;
						kPeak = i;  //this values is used for referencing peak in the List FullScanChromatogram

						*RetentionTimePeak = mzml->FullScanChromatogram[kPeak]->RetTime; //This value is passed back to IntegratePeakTimeMz

						iMaxPeakMZ = iMaxLocalPeak;
					}
					
				}
				/*else 
				{
					printf("Mass Accuracy: This should not be a peak, dMono = %15.7f dPeakFromSpectra = %15.7f time = %f\n", 
						dMonoPrecur, dMoverZLocalPeak, mzml->FullScanChromatogram[i]->RetTime);
				} */

			}

			if (mzml->FullScanChromatogram[i]->RetTime > fHighEnd)
			{
				break;
			}

		}
	}

	fpPeakOutPut->Write("NumberOFPeaks = {0}\n", PeakElution->Count);

	if(kPeak > 0 )
	 //fprintf(fpPeakOutPut, "The peak m/z = %15.7f, time = %f\n", mzml->FullScanChromatogram[kPeak]->moverz[iMaxPeakMZ], mzml->FullScanChromatogram[kPeak]->RetTime);
	fpPeakOutPut->Write("The peak m/z = {0}, time = {1}\n", mzml->FullScanChromatogram[kPeak]->moverz[iMaxPeakMZ], mzml->FullScanChromatogram[kPeak]->RetTime);
	
	if (kPeak == -1)
	{
	//	fprintf(fpPeakOutPut, "Could not detect any signal\n");
		fpPeakOutPut->Write("Could not detect any signal\n");
		//fprintf(fpPeakOutPut, "Need to figure out why there is no detectable signal\n");

		fpPeakOutPut->Write("Need to figure out why there is no detectable signal\n");
		dMassWindow = setMassWindow;
	}
	else
	{
		//Initially use default Elution Time Window to set fStartTime and fEndTime
		for (i = kPeak - 3; i > 0; i--)
		{
			if (mzml->FullScanChromatogram[i]->RetTime < mzml->FullScanChromatogram[kPeak]->RetTime - fTimeWindow)
			{
				break;
			}
		}

		if (i>0)
		{
			fStartTime = mzml->FullScanChromatogram[i]->RetTime;

			iStartTime = i;
		}
		else //if(mzml->FullScanChromatogram[kPeak]->RetTime > 5.0)           //if the chromagraphic elution time is less than 5 minutes, a peak could reach the start.
		{																	// Otherwise it is a problem that the begining of elution has been reached.			
			fStartTime = mzml->FullScanChromatogram[0]->RetTime;
			
			iStartTime = 0;

			if (mzml->FullScanChromatogram[kPeak]->RetTime > 5.0)
			{
				//fprintf(fpPeakOutPut, "The Elution Profile reached the Start of Chromatogram, Start Scan number = %d\n", i);

				fpPeakOutPut->Write("The Elution Profile reached the Start of Chromatogram, Start Scan number = {0}\n", i);
			}
		}

		for (i = kPeak + 3; i < mzml->FullScanChromatogram->Count; i++)
		{
			if (mzml->FullScanChromatogram[i]->RetTime > mzml->FullScanChromatogram[kPeak]->RetTime + fTimeWindow)
			{
				break;
			}
		}

		if (i<mzml->FullScanChromatogram->Count)
		{
			fEndTime = mzml->FullScanChromatogram[i]->RetTime;
			
			iEndTime = i;
		}                                                                   //if the elution is at the end of chromatographic time, then the end can be reached. 
		else //if (mzml->FullScanChromatogram[mzml->FullScanChromatogram->Count - 1]->RetTime - mzml->FullScanChromatogram[kPeak]->RetTime < 5.0)     //Otherwise it is a problem that the end
		{																																	// elution has been reached.			
			fEndTime = mzml->FullScanChromatogram[mzml->FullScanChromatogram->Count - 1]->RetTime;
			
			iEndTime = mzml->FullScanChromatogram->Count - 1;

			if (mzml->FullScanChromatogram[mzml->FullScanChromatogram->Count - 1]->RetTime - mzml->FullScanChromatogram[kPeak]->RetTime < 5.0)
			{
				//fprintf(fpPeakOutPut, "The Elution Profile reached the End of Chromatogram, End Scan Number = %d\n", i);

				fpPeakOutPut->Write("The Elution Profile reached the End of Chromatogram, End Scan Number = {0}\n", i);
			}
		}

		//Calculate mz Peak Width
		FindmzPeakWidth(RetTime, dMonoPrecur, nCharge, fStart, fEnd, jPBStart, jPBEnd, mzPBStart, mzPBEnd, LW, RW, RetentionTimePeak, w, kPeak);


		bWidthCalculated = true;
		if (dMassWindow>2 || dMassWindow <= 0)
		{
			dMassWindow = setMassWindow;
			bWidthCalculated = false;
		}

		//Gather all Retention Time Range values of Retention Time index and corresponding Intensity
		RTR = gcnew List <RetTimeRange ^>;

		//Look for all peak values within the range of Retention Time values calculated above (i.e. between fStartTime and fEndTime)
		for (i = iStartTime; i <= iEndTime; i++)
		{
			tempRTR = gcnew RetTimeRange;
			fSummedPeak = 0.0;
			//tempRTR = gcnew RetTimeRange;
			for (j = 0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
			{
				//do nothing for small m/z values
				if ((dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) > setMassWindow)
				{
					continue;
				}
				//monoisotope
				else if (Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) <= setMassWindow)
				{
					fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if (mzml->FullScanChromatogram[i]->moverz[j] > dMonoPrecur + setMassWindow)
				{
					break;
				}
			}
			tempRTR->RTi = i;
			tempRTR->Intensity = fSummedPeak;
			tempRTR->bRT = true;
			RTR->Add(tempRTR);
			delete tempRTR;
		}

		//RTR->Sort(gcnew Comparison<float ^>(CompareByIntensity));
		wincrement = (double)nCharge;
		kPeak1 = kPeak;
		//tempMassWindow = 0.002;

		/*
		while(wincrement<6)
		{
		fSummedPeak = 0.0;
		for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)
		{
		//do nothing for small m/z values
		if((dMonoPrecur-(deltaC13/wincrement) - mzml->FullScanChromatogram[kPeak1]->moverz[j]) > dMassWindow )
		{
		continue;
		}
		//monoisotope
		else if(Math::Abs(dMonoPrecur-(deltaC13/wincrement) - mzml->FullScanChromatogram[kPeak1]->moverz[j]) <= dMassWindow)
		{
		fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[kPeak1]->Intensity[j];
		}
		else if(mzml->FullScanChromatogram[kPeak1]->moverz[j] > dMonoPrecur-(deltaC13/wincrement) + dMassWindow)
		{
		break;
		}

		}//end for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)

		fSummedPeak0 = 0.0;
		for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)
		{
		//do nothing for small m/z values
		if((dMonoPrecur - mzml->FullScanChromatogram[kPeak1]->moverz[j]) > dMassWindow )
		{
		continue;
		}
		//monoisotope
		else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[kPeak1]->moverz[j]) <= dMassWindow)
		{
		fSummedPeak0 = fSummedPeak0 + mzml->FullScanChromatogram[kPeak1]->Intensity[j];
		}
		else if(mzml->FullScanChromatogram[kPeak1]->moverz[j] > dMonoPrecur + dMassWindow)
		{
		break;
		}

		}//end for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)

		//if(fSummedPeak>0.1*RTR[kPeak1-iStartTime]->Intensity && RTR[kPeak1-iStartTime]->bRT==true)//if I0 condition violated
		if(fSummedPeak>0.1*fSummedPeak0 && RTR[kPeak1-iStartTime]->bRT==true)//if I0 condition violated
		{
		printf("Another peptide's higher isotope conflicts with current peptide. Recalculating Peak Retention Time.\n");
		wincrement = (double)nCharge;
		RTR[kPeak1-iStartTime]->bRT = false;

		//To the left
		i=kPeak1-1;
		if(i>=iStartTime) //When the peak happens at the first value (very rare)
		{
		while(RTR[i-iStartTime]->Intensity>0.1*RTR[kPeak1-iStartTime]->Intensity)
		{
		RTR[i-iStartTime]->bRT = false;
		if(i==iStartTime)
		{
		break;
		}
		i--;
		}
		}

		//To the right
		i=kPeak1+1;
		if(i<=iEndTime) //When the peak happens at the last value (very rare)
		{
		while(RTR[i-iStartTime]->Intensity>0.1*RTR[kPeak1-iStartTime]->Intensity)
		{
		RTR[i-iStartTime]->bRT = false;
		if(i==iEndTime)
		{
		break;
		}
		i++;
		}
		}

		kPeak1 = iStartTime;
		//Calculate new kPeak1
		for(i=iStartTime;i<=iEndTime;i++)
		{
		if(RTR[i-iStartTime]->Intensity > RTR[kPeak1-iStartTime]->Intensity && RTR[i-iStartTime]->bRT==true)
		{
		kPeak1 = i;
		}
		}


		//Recalculate m/z Peak width with the correct Retention Time Peak value
		FindmzPeakWidth(RetTime, dMonoPrecur, nCharge, fStart, fEnd, jPBStart, jPBEnd, mzPBStart, mzPBEnd, LW, RW, RetentionTimePeak, w, kPeak1);
		bWidthCalculated = true;
		if(dMassWindow>2 || dMassWindow<=0)
		{dMassWindow = setMassWindow; bWidthCalculated = false;}

		}
		else //if I0 condition met
		{
		wincrement++;
		}

		} //end while(wincrement<6)
		*/



		if (kPeak1 != kPeak || bWidthCalculated == false)
		{
			//Recalculate m/z Peak width with the correct Retention Time Peak value
			FindmzPeakWidth(RetTime, dMonoPrecur, nCharge, fStart, fEnd, jPBStart, jPBEnd, mzPBStart, mzPBEnd, LW, RW, RetentionTimePeak, w, kPeak1);

			bWidthCalculated = true;
			if (dMassWindow>2 || dMassWindow <= 0)
			{
				dMassWindow = setMassWindow; bWidthCalculated = false;
			}
		}

		//Calculate the Elution Time Start and End based on the accumulated list RTR which holds the new Retention Time Peak and adjacent indices
		if (kPeak1>iStartTime)
		{
			i = kPeak1 - 1;
		}
		else
		{
			i = kPeak1;
		}

		while (RTR[i - iStartTime]->Intensity > 0.1 * RTR[kPeak1 - iStartTime]->Intensity)
		{
			if (i == iStartTime)
			{
				break;
			}
			i--;
		}
		fStartTime = mzml->FullScanChromatogram[i]->RetTime;


		if (kPeak1<iEndTime)
		{
			i = kPeak1 + 1;
		}
		else
		{
			i = kPeak1;
		}
		while (RTR[i - iStartTime]->Intensity > 0.1 * RTR[kPeak1 - iStartTime]->Intensity)
		{
			if (i == iEndTime)
			{
				break;
			}
			i++;
		}
		fEndTime = mzml->FullScanChromatogram[i]->RetTime;


	} ////end if(kPeak==-1) else

	// this is completely need to test!!!! if works delete the previous block

	if (PeakElution->Count > 3)
	{
		fPeakArea = 0.;
		
		array <float> ^fPeakAreas = gcnew array <float>(PeakElution->Count);

		//3 point smoothing
		for (i = 0; i < PeakElution->Count ; i++)
		{
			/*fpPeakOutPut->Write("PeaksHHHHHH Its time Following {0} {1} {2}\n", 
				PeakElution[i]->PeakArea,
				PeakElution[i]->RetTime, PeakElution[i]->FullScanNumber);  */

			if (i == 0 || i == PeakElution->Count - 1)
			{
				fPeakAreas[i] = PeakElution[i]->PeakArea;
			}
			else
			{
				fPeakAreas[i] = PeakElution[i]->PeakArea + PeakElution[i + 1]->PeakArea / 2.0 + PeakElution[i - 1]->PeakArea / 2.0;
			}
			
		}

		for (i = 0; i < PeakElution->Count; i++)
		{

			PeakElution[i]->PeakArea = fPeakAreas[i];

			/*fpPeakOutPut->Write("Peak and Its time Following {0} {1} {2}\n", PeakElution[i]->PeakArea,
				PeakElution[i]->RetTime, PeakElution[i]->FullScanNumber);  */
				
			
			if (PeakElution[i]->PeakArea > fPeakArea)
			{
				fPeakArea = PeakElution[i]->PeakArea;

				iMaxLocalPeak = i;
			}	
		}


		fPeakArea = PeakElution[iMaxLocalPeak]->PeakArea;

		//start of the peak
		for (i = iMaxLocalPeak - 1; i >= 0; i--)
		{
			fPeakArea = fPeakArea + PeakElution[i]->PeakArea;

			if (PeakElution[i]->PeakArea < fPeakThreshold * PeakElution[iMaxLocalPeak]->PeakArea)
			{
				fStartTime = PeakElution[i + 1]->RetTime;

				break;
			}
		}

		//end of the peak
		for (i = iMaxLocalPeak + 1; i < PeakElution->Count; i++)
		{
			if (PeakElution[i]->PeakArea < fPeakThreshold * PeakElution[iMaxLocalPeak]->PeakArea)
			{
				fPeakArea = fPeakArea + PeakElution[i]->PeakArea;

				fEndTime = PeakElution[i-1]->RetTime;

				//fpPeakOutPut->Write("Peak Start End from NewCode Peak Width Following {0} {1}\n", fStartTime, PeakElution[i]->RetTime);

				break;
			}
		}

		fpPeakOutPut->Write("Peak Start End Inside NewCode Peak Width Following {0} {1}\n", fStartTime, fEndTime);

		delete fPeakAreas;

	}  //if (PeakElution->Count > 3)

	// the above block of the code is totally in need to test!!!! if it works delete the previous block

	*fStart = fStartTime;
	*fEnd = fEndTime;

	if (bWidthCalculated == false)
	{
		//fprintf(fpPeakOutPut, "Problem Peptide");

		fpPeakOutPut->Write("Problem Peptide\n");
	}

	if (((*fEnd) - (*fStart)) > 1.)
	{
		//fprintf(fpPeakOutPut, "Long Elution: PeakScan  = %d, Elut Time Peak = %10.5f, dMonoPrecursor = %10.5f., Startelut = %10.5f, Endelut = %10.5f\n",
			//nPeak, fPeakRetTime, dMonoPrecur, fStartTime, fEndTime);

		fpPeakOutPut->Write("Long Elution: PeakScan  = {0}, Elut Time Peak = {1}, dMonoPrecursor = {2}, Startelut = {3}, Endelut = {4}\n",
			nPeak, fPeakRetTime, dMonoPrecur, fStartTime, fEndTime);
	}
	//

	delete PeakElution;

	delete currentPeak;

	//		fHighEnd = mzml->FullScanChromatogram[kPeak]->RetTime + fTimeWindow;
	//		fLowEnd  = mzml->FullScanChromatogram[kPeak]->RetTime - fTimeWindow;


	fLowEnd = fStartTime;
	fHighEnd = fEndTime;

	w[9] = mzMaxPeak;
	w[10] = jMaxPeak;

	return;
}


/*
* 
* The method finds the peak position, start and end, given
* the monoisotopic mass, and retention time
*  the elutions times are assumed to be in minutes;
*  The Window that is used for elution time is 10 minutes!!!
*   It means a peak is searched for 10 minutes before and after
*   the MS2 has been triggered.
*/
void DetectAndIntegratePeak::FindPeakStartEnd(float RetTime, double dMonoPrecur, int nCharge, float *fStart, float *fEnd, int *jPBStart, int *jPBEnd, float *mzPBStart, float *mzPBEnd, float *LW, float *RW, float *RetentionTimePeak,array <double, 1> ^w)
{
	int i, j, nPeak, kPeak, increments, i1, kPeak1;

	float fHighEnd, fLowEnd, fPeakRetTime;

	float fEndTime, fStartTime, fThreshold;

	float fTimeWindow, fPeak, fSummedPeak, fSummedPeak0;

	fTimeWindow = fElutionTimeWindow;       //fElutionTimeWindow is a class variable, set at the time of declaration/initiation.

	if(false) //this is original start and end times
	{
		fHighEnd = RetTime + fTimeWindow;

		fLowEnd  = RetTime - fTimeWindow;
	}
	else
	{
		fHighEnd = (*fEnd) + fTimeWindow;

		fLowEnd  = (*fStart) - fTimeWindow;

		//printf("Initial Start %10.5f and End %10.5f Elution Times\n", fLowEnd, fHighEnd);
	}

	fPeak = fSummedPeak = 0.0;

	fEndTime =  fStartTime = -1.;

	fThreshold = 0.01f;  kPeak = -1; kPeak1 = -1;

	//Added by JA 2016.09.08
	//Following variables store intensity values
	float fMaxPeak = 0.0f;
	float fPB = 0.0f;

	//Following variable stores incremental factor to dMassWindow
	double wincrement = 1.00;
	//Added by JA 2016.10.27
	double tempMassWindow = setMassWindow;//0.001;

	//Following variables store m/z values
	float mzMaxPeak = 0.0f; 
	float mzPBLow = 0.0f;
	float mzPBHigh = 0.0f;
	float mzPBLowT = 0.0f;
	float mzPBLowC = 0.0f;
	float mzPBHighT =0.0f;
	float mzPBHighC = 0.0f;

	//Following variables store indices
	int jMaxPeak = 0;
	int jPBLow = 0;
	int jPBHigh = 0;
	int jPBLowT = 0;
	int jPBLowC = 0;
	int jPBHighT = 0;
	int jPBHighC = 0;
	//

	//Added by JA 2016.10.31
	bool bWidthCalculated = true;

	//Added by JA 2016.10.26
	RetTimeRange ^tempRTR;
	List <RetTimeRange ^> ^RTR;
	//List <float> ^RetTimeRangeIntensities;

	//float *LeftRTIntensities = (float*)calloc(200,sizeof(float));
	//float *RightRTIntensities = (float*)calloc(200,sizeof(float));
	int iStartTime, iEndTime;

	//Calculate Retention Time based Peak first with a masswindow that is incremented until we find a Retention Time Peak
	//wincrement = 0.001;


	for(i=0; i < mzml->FullScanChromatogram->Count; i++)
	{
		if(mzml->FullScanChromatogram[i]->RetTime < fLowEnd)
		{
			continue;
		}

		if(mzml->FullScanChromatogram[i]->RetTime >= fLowEnd &&
			mzml->FullScanChromatogram[i]->RetTime <= fHighEnd)
		{
			fSummedPeak = 0.0;


			//for(j=*jPBStart; j < *jPBEnd; j++)
			for(j=0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
			{
				//do nothing for small m/z values
				if((dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) > setMassWindow )
				{
					continue;
				}
				//monoisotope
				else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) <= setMassWindow)
				{
					fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(mzml->FullScanChromatogram[i]->moverz[j] > dMonoPrecur + setMassWindow)
				{
					break;
				}
			}

			if( fSummedPeak > fPeak)
			{
				fPeak = fSummedPeak;
				nPeak = mzml->FullScanChromatogram[i]->ScanNumber;
				fPeakRetTime = mzml->FullScanChromatogram[i]->RetTime;
				kPeak = i;  //this values is used for referencing peak in the List FullScanChromatogram
				//Added by JA 2016.09.16
				*RetentionTimePeak = mzml->FullScanChromatogram[kPeak]->RetTime; //This value is passed back to IntegratePeakTimeMz
				//
			}

			if(mzml->FullScanChromatogram[i]->RetTime > fHighEnd)
			{
				break; 
			}

		}
	}



	if(kPeak == -1)
	{
		fpPeakOutPut->Write("Could not detect any signal\n");
		
		fpPeakOutPut->Write("Need to figure out why there is no detectable signal\n");
		
		dMassWindow = setMassWindow;
	}
	else
	{
		//Initially use default Elution Time Window to set fStartTime and fEndTime
	for(i = kPeak-3; i > 0; i--)
		{
			if(mzml->FullScanChromatogram[i]->RetTime < mzml->FullScanChromatogram[kPeak]->RetTime - fTimeWindow)
			{
				break;
			}
		}
if(i>0)
{
	fStartTime = mzml->FullScanChromatogram[i]->RetTime;
	iStartTime = i;
}
else 
{
	printf("Suspicious activity");
	fStartTime = mzml->FullScanChromatogram[0]->RetTime;
	iStartTime = 0;
}

	for(i = kPeak + 3; i < mzml->FullScanChromatogram->Count; i++)
		{
			if(mzml->FullScanChromatogram[i]->RetTime > mzml->FullScanChromatogram[kPeak]->RetTime + fTimeWindow)
			{
				break;
			}
		}

if(i<mzml->FullScanChromatogram->Count)
{
	fEndTime = mzml->FullScanChromatogram[i]->RetTime;
	iEndTime = i;
}
else
{
	printf("Suspicious activity");
	fEndTime = mzml->FullScanChromatogram[mzml->FullScanChromatogram->Count-1]->RetTime;
	iEndTime = mzml->FullScanChromatogram->Count-1;
}

//Calculate mz Peak Width
		FindmzPeakWidth(RetTime, dMonoPrecur, nCharge, fStart, fEnd, jPBStart, jPBEnd, mzPBStart, mzPBEnd, LW, RW, RetentionTimePeak, w, kPeak);


		bWidthCalculated = true;
		if(dMassWindow>2 || dMassWindow<=0)
		{
			dMassWindow = setMassWindow;
			bWidthCalculated = false;
		}

		//Gather all Retention Time Range values of Retention Time index and corresponding Intensity
		RTR = gcnew List <RetTimeRange ^>;

	//Look for all peak values within the range of Retention Time values calculated above (i.e. between fStartTime and fEndTime)
	for(i=iStartTime;i<=iEndTime;i++)
	{
				tempRTR = gcnew RetTimeRange;
				fSummedPeak = 0.0;
				//tempRTR = gcnew RetTimeRange;
				for(j=0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
				{
					//do nothing for small m/z values
					if((dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) > setMassWindow )
					{
						continue;
					}
					//monoisotope
					else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) <= setMassWindow)
					{
						fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[i]->Intensity[j];
					}
					else if(mzml->FullScanChromatogram[i]->moverz[j] > dMonoPrecur + setMassWindow)
					{
						break;
					}
				}
				tempRTR->RTi = i;
				tempRTR->Intensity = fSummedPeak;
				tempRTR->bRT = true;
				RTR->Add(tempRTR);
				delete tempRTR;
	}

	//RTR->Sort(gcnew Comparison<float ^>(CompareByIntensity));
	wincrement = (double)nCharge;
	kPeak1 = kPeak;
	//tempMassWindow = 0.002;

	/*
	while(wincrement<6)
	{
	fSummedPeak = 0.0;
	for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)
	{
					//do nothing for small m/z values
					if((dMonoPrecur-(deltaC13/wincrement) - mzml->FullScanChromatogram[kPeak1]->moverz[j]) > dMassWindow )
					{
						continue;
					}
					//monoisotope
					else if(Math::Abs(dMonoPrecur-(deltaC13/wincrement) - mzml->FullScanChromatogram[kPeak1]->moverz[j]) <= dMassWindow)
					{
						fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[kPeak1]->Intensity[j];
					}
					else if(mzml->FullScanChromatogram[kPeak1]->moverz[j] > dMonoPrecur-(deltaC13/wincrement) + dMassWindow)
					{
						break;
					}

	}//end for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)

	fSummedPeak0 = 0.0;
	for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)
	{
					//do nothing for small m/z values
					if((dMonoPrecur - mzml->FullScanChromatogram[kPeak1]->moverz[j]) > dMassWindow )
					{
						continue;
					}
					//monoisotope
					else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[kPeak1]->moverz[j]) <= dMassWindow)
					{
						fSummedPeak0 = fSummedPeak0 + mzml->FullScanChromatogram[kPeak1]->Intensity[j];
					}
					else if(mzml->FullScanChromatogram[kPeak1]->moverz[j] > dMonoPrecur + dMassWindow)
					{
						break;
					}

	}//end for(j=0;j<mzml->FullScanChromatogram[kPeak1]->moverz->Length; j++)

			//if(fSummedPeak>0.1*RTR[kPeak1-iStartTime]->Intensity && RTR[kPeak1-iStartTime]->bRT==true)//if I0 condition violated
			if(fSummedPeak>0.1*fSummedPeak0 && RTR[kPeak1-iStartTime]->bRT==true)//if I0 condition violated
			{
					printf("Another peptide's higher isotope conflicts with current peptide. Recalculating Peak Retention Time.\n");
					wincrement = (double)nCharge;
					RTR[kPeak1-iStartTime]->bRT = false;
					
					//To the left
					i=kPeak1-1;
					if(i>=iStartTime) //When the peak happens at the first value (very rare)
					{
					while(RTR[i-iStartTime]->Intensity>0.1*RTR[kPeak1-iStartTime]->Intensity)
					{
						RTR[i-iStartTime]->bRT = false;
						if(i==iStartTime)
						{
							break;
						}
						i--;
					}
					}

					//To the right
					i=kPeak1+1;
					if(i<=iEndTime) //When the peak happens at the last value (very rare)
					{
					while(RTR[i-iStartTime]->Intensity>0.1*RTR[kPeak1-iStartTime]->Intensity)
					{
						RTR[i-iStartTime]->bRT = false;
						if(i==iEndTime)
						{
							break;
						}
						i++;
					}
					}

					kPeak1 = iStartTime;
					//Calculate new kPeak1
					for(i=iStartTime;i<=iEndTime;i++)
					{
						if(RTR[i-iStartTime]->Intensity > RTR[kPeak1-iStartTime]->Intensity && RTR[i-iStartTime]->bRT==true)
						{
							kPeak1 = i;
						}
					}

					
					//Recalculate m/z Peak width with the correct Retention Time Peak value
					FindmzPeakWidth(RetTime, dMonoPrecur, nCharge, fStart, fEnd, jPBStart, jPBEnd, mzPBStart, mzPBEnd, LW, RW, RetentionTimePeak, w, kPeak1);
					bWidthCalculated = true;
					if(dMassWindow>2 || dMassWindow<=0)
					{dMassWindow = setMassWindow; bWidthCalculated = false;}

				}
			else //if I0 condition met
				{
					wincrement++;
				}

	} //end while(wincrement<6)
	*/



	if(kPeak1!=kPeak || bWidthCalculated == false)
	{
	//Recalculate m/z Peak width with the correct Retention Time Peak value
	FindmzPeakWidth(RetTime, dMonoPrecur, nCharge, fStart, fEnd, jPBStart, jPBEnd, mzPBStart, mzPBEnd, LW, RW, RetentionTimePeak, w, kPeak1);

		bWidthCalculated = true;
		if(dMassWindow>2 || dMassWindow<=0)
		{dMassWindow = setMassWindow; bWidthCalculated = false;}
	}

		//Calculate the Elution Time Start and End based on the accumulated list RTR which holds the new Retention Time Peak and adjacent indices
					if(kPeak1>iStartTime)
					{
						i=kPeak1-1;
					}
					else
					{
						i = kPeak1;
					}

					while(RTR[i-iStartTime]->Intensity>0.1*RTR[kPeak1-iStartTime]->Intensity)
					{
						if(i==iStartTime)
						{
							break;
						}
						i--;
					}
					fStartTime = mzml->FullScanChromatogram[i]->RetTime;
					

					if(kPeak1<iEndTime)
					{
						i=kPeak1+1;
					}
					else
					{
						i = kPeak1;
					}
					while(RTR[i-iStartTime]->Intensity>0.1*RTR[kPeak1-iStartTime]->Intensity)
					{
						if(i==iEndTime)
						{
							break;
						}
						i++;
					}
					fEndTime = mzml->FullScanChromatogram[i]->RetTime;
					

		} ////end if(kPeak==-1) else

					*fStart = fStartTime;
					*fEnd   = fEndTime;

					if(bWidthCalculated == false)
					{
						printf("Problem Peptide");
					}

		if(((*fEnd) - (*fStart) ) > 5.)
		{
			printf("Long Elution: PeakScan  = %d, Elut Time Peak = %10.5f, dMonoPrecursor = %10.5f., Startelut = %10.5f, Endelut = %10.5f\n", 
				nPeak, fPeakRetTime, dMonoPrecur, fStartTime, fEndTime);
		}
		//



		//		fHighEnd = mzml->FullScanChromatogram[kPeak]->RetTime + fTimeWindow;
		//		fLowEnd  = mzml->FullScanChromatogram[kPeak]->RetTime - fTimeWindow;


		fLowEnd = fStartTime;
		fHighEnd = fEndTime;

	w[9] = mzMaxPeak;
	w[10] = jMaxPeak;

	return;
}

/*
*   a method to find the peak position
*   it uses the monoisotopic mass of the peak,
*   RetTime is the retention time at which MS2 was
*   triggered.
*   fTimeWindow - is a global variable and set in the quant.state file
*/
/*
*   Read the following the following conditions to understand  what a peak means:
*   First, do the weight based averaging of the points to reduce noise contribution
*   With the weight averaged elution profile
A point is a peak if:
1. It is highest point amond the 25 points centered on this point
2. Its intensity is higher that 20% of the averaged intensity as determined
from the start and end of the peak

*/
void  DetectAndIntegratePeak::FindPeakPosition(float RetTime, double dMonoPrecur, int nCharge)
{
	int i, j, k = 0, nPointWidth;

	float fHighEnd, fLowEnd;

	float fThreshold;

	float fPeak, fSummedPeak, fNinePeaks;

	bool bPeak = false;

	List <PeakProfile ^> ^PeakElution, ^PeakElution2, ^PeakPairs;

	PeakProfile ^SinglePeak;

	fHighEnd = RetTime + fElutionTimeWindow;

	fLowEnd  = RetTime - fElutionTimeWindow;

	fPeak = fSummedPeak = 0.0;

	fThreshold = 0.01f;

	PeakElution = gcnew List <PeakProfile ^>;

	PeakElution2 = gcnew List <PeakProfile ^>;

	PeakPairs = gcnew List <PeakProfile ^>;

	nPointWidth = 7;

	//generate the profile of the peak

	for(i=0; i < mzml->FullScanChromatogram->Count; i++)
	{
		if(mzml->FullScanChromatogram[i]->RetTime < fLowEnd)
		{
			continue;
		}

		if(mzml->FullScanChromatogram[i]->RetTime >= fLowEnd &&
			mzml->FullScanChromatogram[i]->RetTime <= fHighEnd)
		{

			SinglePeak = gcnew PeakProfile;

			fSummedPeak = 0.0;

			for(j=0; j < mzml->FullScanChromatogram[i]->moverz->Length; j++)
			{
				//do nothing for small m/z values
				if((dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) > dMassWindow )
				{
					continue;
				}
				//monoisotope
				else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i]->moverz[j]) <= dMassWindow)
				{
					fSummedPeak = fSummedPeak + mzml->FullScanChromatogram[i]->Intensity[j];
				}
				else if(mzml->FullScanChromatogram[i]->moverz[j] > dMonoPrecur + dMassWindow)
				{
					break;
				}
			}

			SinglePeak->Intensity = fSummedPeak;

			SinglePeak->RetTime = mzml->FullScanChromatogram[i]->RetTime;

			SinglePeak->ScanNumber = mzml->FullScanChromatogram[i]->ScanNumber;

			PeakElution->Add(SinglePeak);

			delete SinglePeak;
		}

		if(mzml->FullScanChromatogram[i]->RetTime > fHighEnd)
		{
			break; 
		}
	}

	// use the profile of the monoisotope to determine the
	// peak position
	// to be a peak a profile must have 9 non-zero points
	// averaging and the peak should be 
	// 
	fSummedPeak  = fNinePeaks = fPeak = 0.0f;

	for(i=4; i < PeakElution->Count - 4; i++)
	{
		//do not change or average the zero points
		// zero should stay zero
		if(PeakElution[i]->Intensity < 2.0)
		{
			PeakElution2->Add(PeakElution[i]);

			continue;
		}
		fNinePeaks = PeakElution[i-4]->Intensity/16. +  PeakElution[i-3]->Intensity/8. + PeakElution[i-2]->Intensity/4.;

		fNinePeaks += PeakElution[i-1]->Intensity/2. + PeakElution[i]->Intensity + PeakElution[i+1]->Intensity/2.;

		fNinePeaks += PeakElution[i+2]->Intensity/4. + PeakElution[i+3]->Intensity/8. + PeakElution[i+4]->Intensity/16.;

		SinglePeak = gcnew PeakProfile;

		SinglePeak->Intensity = fNinePeaks;

		SinglePeak->RetTime = PeakElution[i]->RetTime;

		SinglePeak->ScanNumber = PeakElution[i]->ScanNumber;

		PeakElution2->Add(SinglePeak);

		delete SinglePeak;
	}

	// a Peak in elution profile is a true peak
	// if it the highest point among the 9 points centered at it
	// and the points are non-zero intensity.

	for(i = nPointWidth; i < PeakElution2->Count - nPointWidth - 1; i++)
	{
		if(PeakElution2->Count < 2.0)
		{
			continue;
		}
		bPeak = true;

		fPeak = PeakElution2[i]->Intensity;

		for(j = 1; j <= nPointWidth; j++)
		{
			//the idea is to have the peak to be highest point - among 9 points centered at the peak
			// also there should be a clear gradient - that should be a monotonic increase to the peak
			// none of the 9 points can be zero intensity - meaning that a peak should be at least as
			// wider as 9 points
			//

			//Left side of the peak
			if(PeakElution2[i - j]->Intensity < 2.0 || PeakElution2[i - j]->Intensity > fPeak)
			{
				bPeak = false;

				break;
			}
			else if (PeakElution2[i - j]->Intensity < PeakElution2[i - j - 1]->Intensity)
			{
				bPeak = false;

				break;
			}

			//Right side of the peak
			if(PeakElution2[i + j]->Intensity < 2.0 || PeakElution2[i + j]->Intensity > fPeak)
			{
				bPeak = false;

				break;
			}
			else if (PeakElution2[i+j]->Intensity > PeakElution2[i + j - 1]->Intensity)
			{
				bPeak = false;

				break;
			}

		}

		if(bPeak)
		{
			SinglePeak = gcnew PeakProfile;

			SinglePeak->Intensity = PeakElution2[i]->Intensity;

			SinglePeak->RetTime = PeakElution2[i]->RetTime;

			SinglePeak->ScanNumber = i;                               // important for accessing later!!!!

			delete SinglePeak;

			PeakPairs->Add(SinglePeak);
		}

		printf("%10.5f %10.5f\n", PeakElution2[i]->RetTime, PeakElution2[i]->Intensity);
	}


	//Remove the "peaks" that are too small
	//

	for(i = 0; i < PeakPairs->Count; i++)
	{
		printf("OFFICIAL PEAKS = %10.5f %10.5f %d\n", PeakPairs[i]->RetTime, PeakPairs[i]->Intensity, PeakPairs[i]->ScanNumber);

		fPeak = 0.0;

		for(j = 0; j < PeakElution2->Count; j++)
		{
			if(PeakPairs[i]->ScanNumber == j)
			{

				//printf("k = %d %d %d, j = %d\n", k, j-k, j+k, j);

				//exit(1);

				fPeak = 0.0;

				for(k=1; k < 25 && (j - k >=0) && (j+k < PeakElution2->Count); k++)
				{
					//printf("k = %d %d %d\n", k, j-k, j+k);

					fPeak = fPeak + PeakElution2[j - k]->Intensity + PeakElution2[j + k]->Intensity;     //take average of 40 values 

					printf("Gradient %10.5f %10.5f  %10.5f %10.5f\n", (PeakElution2[j - k]->Intensity - PeakElution2[j - k - 1]->Intensity)/PeakElution2[j]->Intensity,
						(PeakElution2[j + k]->Intensity - PeakElution2[j + k + 1]->Intensity)/PeakElution2[j]->Intensity, PeakElution2[j-k]->Intensity, PeakElution2[j- k - 1]->Intensity);
				}


				PeakPairs[i]->ScanNumber = PeakElution2[j]->ScanNumber;      //restore the scan number

				break;
			}
		}

		fPeak = fPeak/50.0;

		if(fPeak/PeakPairs[i]->Intensity > 0.2)
		{
			PeakPairs[i]->ScanNumber = -1;        //if not a peak set scan number to -1
		}

		printf("%10.5f %10.5f\n", PeakPairs[i]->Intensity, fPeak);

	}

	for(i = 0; i < PeakPairs->Count; i++)
	{
		if(PeakPairs[i]->ScanNumber != -1)
		{
			printf("Serviving PEAKS = %10.5f %10.5f\n", PeakPairs[i]->RetTime, PeakPairs[i]->Intensity);
		}
	}

	delete SinglePeak, PeakElution, PeakElution2, PeakPairs;
	//check the width of the peak
}


/*
*  A the function will return three indices: 1) peak position, 2) peak left start, 3) peak right start
*  It uses MS1 survey scan, theoretical m/z value, mass accuracy in ppm
*/

PeakIndex ^ DetectAndIntegratePeak::MoverZ_Peak_Position_Width(FullScan ^MS1_Scan, double dActualM_z, double dMAccuracy)
{
	PeakIndex ^ PeakPosit;

	int j, iMaxLocalPeak;

	double dMoverZLocalPeak;
	
	float fmaxIntensity, fMaxLocalPeak;;


	fMaxLocalPeak = 0.0;

	for (j = 0; j < MS1_Scan->moverz->Length; j++)
	{
		//do nothing for small m/z values
		if ((dActualM_z - MS1_Scan->moverz[j]) > dMAccuracy)
		{
			continue;
		}
		//monoisotope
		else if (Math::Abs(dActualM_z - MS1_Scan->moverz[j]) <= dMAccuracy)
		{
			//fSummedPeak = fSummedPeak + MS1_Scan->Intensity[j];

			if (MS1_Scan->Intensity[j] > fMaxLocalPeak)
			{

				fMaxLocalPeak = MS1_Scan->Intensity[j];

				iMaxLocalPeak = j;

				dMoverZLocalPeak = MS1_Scan->moverz[j];
			}
		}
		else if (MS1_Scan->moverz[j] > dActualM_z + dMAccuracy)
		{
			break;
		}
	}

	return PeakPosit;
}

//Added by JA 2016.09 - 2016.10
void DetectAndIntegratePeak::FindmzPeakWidth(float RetTime, double dMonoPrecur, int nCharge, float *fStart, float *fEnd, int *jPBStart, int *jPBEnd, float *mzPBStart, float *mzPBEnd, float *LW, float *RW, float *RetentionTimePeak,array <double, 1> ^w, int kPeak)
{
		
	//Added by JA 2016.09.08
	//Following variables store intensity values
	float fMaxPeak = 0.0f;
	float fPB = 0.0f;
	float fSummedPeak = 0.0f;

	//Following variable stores incremental factor to dMassWindow
	double wincrement = 1.00;
	//Added by JA 2016.10.27
	double tempMassWindow = setMassWindow;//0.001;

	//Following variables store m/z values
	float mzMaxPeak = 0.0f; 
	float mzPBLow = 0.0f;
	float mzPBHigh = 0.0f;
	float mzPBLowT = 0.0f;
	float mzPBLowC = 0.0f;
	float mzPBHighT =0.0f;
	float mzPBHighC = 0.0f;

	//Following variables store indices
	int i1;
	int j;
	int jMaxPeak = 0;
	int jPBLow = 0;
	int jPBHigh = 0;
	int jPBLowT = 0;
	int jPBLowC = 0;
	int jPBHighT = 0;
	int jPBHighC = 0;
	
	
	i1 = kPeak;

		//Find m/z based peak at the calculated Retention Time Peak

		fSummedPeak = 0.0f;
		//Added by JA 2016.09.08
		//Following variables store intensity values
		fMaxPeak = 0.0f;
		fPB = 0.0f;

		//Following variables store m/z values
		mzMaxPeak = 0.0f;
		*mzPBStart = 0.0f;
		*mzPBEnd = 0.0f;
		mzPBLow = 0.0f;
		mzPBHigh = 0.0f;
		mzPBLowT = 0.0f;
		mzPBLowC = 0.0f;
		mzPBHighT = 0.0f;
		mzPBHighC = 0.0f;
		//

		//Following variables store indices
		jMaxPeak = 0;
		jPBLow = 0;
		jPBHigh = 0;
		jPBLowT = 0;
		jPBLowC = 0;
		jPBHighT = 0;
		jPBHighC = 0;
		*jPBStart = 0;
		*jPBEnd  = 0;
		//C and T suffixes above stand for Contaminant and Threshold respectively; 
		//i.e.if the PBLow or High is calculated dependent on the Threshold 0.01*fMaxPeak, then relevant info is stored in jPBLowT, etc.
		//if the PBLow or High is calculated dependent on the Contaminating peak entering the picture, the relevant info is stored in jPBLowC, etc.
		tempMassWindow = wincrement*setMassWindow;
		wincrement = 1.00;
		//Determine True Peak Intensity and m/z
		while(wincrement<=2)
		{
			for(j=0; j < mzml->FullScanChromatogram[i1]->moverz->Length; j++)
			{
				//do nothing for small m/z values
				if((dMonoPrecur - mzml->FullScanChromatogram[i1]->moverz[j]) > wincrement*setMassWindow )
				{
					continue;
				}//monoisotope
				else if(Math::Abs(dMonoPrecur - mzml->FullScanChromatogram[i1]->moverz[j]) <= wincrement*setMassWindow)
				{
					if(fMaxPeak <= mzml->FullScanChromatogram[i1]->Intensity[j])
					{
						fMaxPeak = mzml->FullScanChromatogram[i1]->Intensity[j];
						mzMaxPeak = mzml->FullScanChromatogram[i1]->moverz[j];
						jMaxPeak = j;
					}
				}
				else if(mzml->FullScanChromatogram[i1]->moverz[j] > dMonoPrecur + wincrement*setMassWindow)
				{
					break;
				}
			}


			if(jMaxPeak==0 || fMaxPeak==0)
			{
				wincrement = wincrement + 0.1;
			}
			else
			{
				wincrement = 2.5;
			}
		}
		dMassWindow = 0;
		if(fMaxPeak==0)
		{
			printf("Peak not found");
		}
		else
		{
		//Start Low
		//Determine the low end of the True Peak Base along the mz-axis
		if(wincrement==2.5)
			{
				wincrement = 1.00;
			}
		tempMassWindow = wincrement*setMassWindow;

		while(mzPBHigh==0 || mzPBLow==0 || wincrement<5)
		{
		for(j=jMaxPeak; j>=0; j--)
		{
			if(mzml->FullScanChromatogram[i1]->moverz[j] > mzMaxPeak + wincrement*setMassWindow)
			{
				continue;
			}
			//monoisotope
			else if(Math::Abs(mzMaxPeak - mzml->FullScanChromatogram[i1]->moverz[j]) <= wincrement*setMassWindow)
			{
				if(mzPBLowC==0 && jPBLowC==0)
				{
					if(mzml->FullScanChromatogram[i1]->Intensity[j] >= mzml->FullScanChromatogram[i1]->Intensity[j+1] && mzml->FullScanChromatogram[i1]->Intensity[j+1] < mzml->FullScanChromatogram[i1]->Intensity[j+2])
					{
						mzPBLowC = mzml->FullScanChromatogram[i1]->moverz[j+1];
						jPBLowC = j+1;
					}
				}

				//j<=jMaxPeak-2 && 
				if(mzPBLowT==0 && jPBLowT==0)
				{
					if(mzml->FullScanChromatogram[i1]->Intensity[j]<0.001*fMaxPeak && mzml->FullScanChromatogram[i1]->Intensity[j+1]>=0.001*fMaxPeak)
					{
						mzPBLowT = mzml->FullScanChromatogram[i1]->moverz[j+1];
						jPBLowT = j+1;
					}
				}
			}
			else if((mzMaxPeak - mzml->FullScanChromatogram[i1]->moverz[j]) > wincrement*setMassWindow )
			{
				break;
			}
		}
		//End Low


		//Start High
		//Determine the high end of the True Peak Base along the mz-axis
		for(j=jMaxPeak; j<=mzml->FullScanChromatogram[i1]->moverz->Length; j++)
		{
			if((mzMaxPeak - mzml->FullScanChromatogram[i1]->moverz[j]) > wincrement*setMassWindow )
			{
				continue;
			}
			//monoisotope
			else if(Math::Abs(mzMaxPeak - mzml->FullScanChromatogram[i1]->moverz[j]) <= wincrement*setMassWindow)
			{
				if(mzPBHighC==0 && jPBHighC==0)
				{
					if(mzml->FullScanChromatogram[i1]->Intensity[j] >= mzml->FullScanChromatogram[i1]->Intensity[j-1] && mzml->FullScanChromatogram[i1]->Intensity[j-1] < mzml->FullScanChromatogram[i1]->Intensity[j-2])
					{
						mzPBHighC = mzml->FullScanChromatogram[i1]->moverz[j-1];
						jPBHighC = j-1;
					}
				}
				if(mzPBHighT==0 && jPBHighT==0)
				{
					if(mzml->FullScanChromatogram[i1]->Intensity[j]<0.001*fMaxPeak && mzml->FullScanChromatogram[i1]->Intensity[j-1]>=0.001*fMaxPeak)
					{
						mzPBHighT = mzml->FullScanChromatogram[i1]->moverz[j-1];
						jPBHighT = j-1;
					}
				}
			}
			else if(mzml->FullScanChromatogram[i1]->moverz[j] > mzMaxPeak + wincrement*setMassWindow)
			{
				break;
			}
		}
		//End High

		/*
		printf("Right Intensity T = %f\n",mzml->FullScanChromatogram[i1]->Intensity[jPBHighT]);
		printf("Right m/z T = %f\n",mzPBHighT);
		printf("Right index T  = %d\n",jPBHighT);
		*/

		if(jPBLowT<=jPBLowC) {jPBLow = jPBLowC; mzPBLow = mzPBLowC;}
		if(jPBLowT>jPBLowC) {jPBLow = jPBLowT; mzPBLow = mzPBLowT;}

		if(jPBHighT<=jPBHighC && jPBHighT>0) 
		{jPBHigh = jPBHighT; mzPBHigh = mzPBHighT;}
		else if(jPBHighT==0 && jPBHighC>0)
		{jPBHigh = jPBHighC; mzPBHigh = mzPBHighC;}


		if(jPBHighT>=jPBHighC && jPBHighC>0) 
		{jPBHigh = jPBHighC; mzPBHigh = mzPBHighC;}
		else if(jPBHighC==0 && jPBHighT>0)
		{jPBHigh = jPBHighT; mzPBHigh = mzPBHighT;}

		wincrement = wincrement + 0.1;
		} //End while(mzPBHigh==0 || mzPBLow==0)

		*RW = mzPBHigh - mzMaxPeak;
		*LW = mzMaxPeak - mzPBLow;

		if(*LW <= *RW)
		{dMassWindow = *LW;}
		else
		{dMassWindow = *RW;}


		} //end if(fMaxPeak!=0) else
}



/*
int CompareByIntensity(float ^x, float ^y)
{
    if(x->Intensity < y->Intensity)
	{return y;}
	else
	{return x;}
};
*/