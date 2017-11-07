#include <iostream>
#include <string>
#include <fstream>
#include <wchar.h>
#include <time.h>
#include <numeric>

#include "getQuant.h"

void getQuant(	std::pair<std::vector<double>,std::vector<double> > &allQuants,
				std::pair<std::vector<std::vector<double> >, std::vector<std::vector<double> > > &allFWHMs,
				std::pair<std::vector<float>,std::vector<float> > &allChiSquareds,
				std::pair<std::vector<int>,std::vector<int> > &allDegreesOfFreedoms,
				std::pair<std::vector<float>,std::vector<float> > &allIntensities,
				std::pair<std::vector<float>,std::vector<float> > &allBestPeakRts,
				std::string &fullFileName,
				peptideInfo &newPeptideInfo){

	int numPeptides = newPeptideInfo.charges.size();

	fullFileName = "F:\\data\\CTAM\\selectedRuns\\RAW\\"+fullFileName+".raw";

	for(int sample=0; sample<2; ++sample){

		IXRawfile5* pISL;

		std::wstring wRawFile(fullFileName.begin(),fullFileName.end());
		wRawFile.append(1, 0);
		openRawFile(wRawFile, pISL);

		pISL->SetCurrentController ( 0, 1 ); // first MS controller

		BSTR filter = NULL;
		long scanNumber = 1;
		long nRet = pISL->GetFilterForScanNum(scanNumber,&filter);

		std::vector<double> peakAreas;
		std::vector<std::vector<double> > FWHMs;
		std::vector<float> chiSquareds;
		std::vector<int> degressOfFreedoms;
		std::vector<float> intensities;
		std::vector<float> bestPeakRts;

		for(int peptide=0; peptide<numPeptides; ++peptide){

			if(peptide==0) std::cout <<"Peptide " << peptide << " / " << numPeptides;		//
			if(peptide%200==0){																//
				std::cout << '\r';															//
				std::cout << "Peptide " << peptide << " / " << numPeptides;					//	Print to console the progress of quantitation,
			}																				//	updtating every 200 peptides.
			if(peptide+1==numPeptides){														//
				std::cout << '\r';															//
				std::cout << "Peptide " << peptide+1 << " / " << numPeptides << std::endl;	//
			}																				//

			double mz;
			float rt;
			int charge = newPeptideInfo.charges[peptide];

			if(sample==0){
				mz = newPeptideInfo.massToChargesMascot[peptide];
				rt = newPeptideInfo.retentionTimesMascot[peptide];
			} else{
				mz = newPeptideInfo.massToChargesXTandem[peptide];
				rt = newPeptideInfo.retentionTimesXTandem[peptide];
			}

			double mzWindow = 14; // This should be user define as a parameter
			float rtWindow = 1.5; // This should be user define as a parameter

			std::vector< std::pair<double, double> > timeIntPairs;
			std::vector< std::vector< std::pair<double, double> > > timeIntIsotopes;

			// First isotope
			double mzLower = mz - (mz * mzWindow) / 1000000;
			double mzUpper = mz + (mz * mzWindow) / 1000000;

			std::wstring mzRange = std::to_wstring(mzLower);
			mzRange=mzRange+std::wstring(L"-");
			mzRange+=std::to_wstring(mzUpper);
			mzRange.append(1, 0);
			float rtLower = rt - rtWindow;
			float rtUpper = rt + rtWindow;

			getTimeIntPairs(timeIntPairs, filter, mzRange, rtLower, rtUpper, pISL);
			timeIntIsotopes.push_back(timeIntPairs);
			timeIntPairs.clear();

			// Second isotope
			mz = mz + (1.0022 / charge);
			mzLower = mz - (mz * mzWindow) / 1000000;
			mzUpper = mz + (mz * mzWindow) / 1000000;

			mzRange = std::to_wstring(mzLower);
			mzRange=mzRange+std::wstring(L"-");
			mzRange+=std::to_wstring(mzUpper);
			mzRange.append(1, 0);
			getTimeIntPairs(timeIntPairs, filter, mzRange, rtLower, rtUpper, pISL);

			timeIntIsotopes.push_back(timeIntPairs);
			timeIntPairs.clear();

			// Third isotope
			mz = mz + (1.0022 / charge);
			mzLower = mz - (mz * mzWindow) / 1000000;
			mzUpper = mz + (mz * mzWindow) / 1000000;

			mzRange = std::to_wstring(mzLower);
			mzRange=mzRange+std::wstring(L"-");
			mzRange+=std::to_wstring(mzUpper);
			mzRange.append(1, 0);

			getTimeIntPairs(timeIntPairs, filter, mzRange, rtLower, rtUpper, pISL);

			timeIntIsotopes.push_back(timeIntPairs);

			peakInfo bestPeak;
			if(timeIntIsotopes[0].size()==0){
			    bestPeak.intensity=0;
			    bestPeak.area=0;
			    bestPeak.index=0;
			    std::vector<double> vec(3);
			    bestPeak.FWHMs=vec;
			    bestPeak.chiSquared=std::numeric_limits<float>::max();
			} else{
				peakHeight(timeIntIsotopes, bestPeak);
			}
			peakAreas.push_back(bestPeak.area);
			FWHMs.push_back(bestPeak.FWHMs);
			chiSquareds.push_back(bestPeak.chiSquared);
			degressOfFreedoms.push_back(bestPeak.degressOfFreedom);
			intensities.push_back(bestPeak.intensity);
			bestPeakRts.push_back(bestPeak.rt);
		}

		if(sample==0){
			allQuants.first = peakAreas;
			allFWHMs.first = FWHMs;
			allChiSquareds.first = chiSquareds;
			allDegreesOfFreedoms.first = degressOfFreedoms;
			allIntensities.first = intensities;
			allBestPeakRts.first = bestPeakRts;

		} else{
			allQuants.second = peakAreas;
			allFWHMs.second = FWHMs;
			allChiSquareds.second = chiSquareds;
			allDegreesOfFreedoms.second = degressOfFreedoms;
			allIntensities.second = intensities;
			allBestPeakRts.second = bestPeakRts;
		}

		pISL->Close();
		pISL->Release();

		CoUninitialize();
	}
}

void openRawFile(std::wstring wRawFile, IXRawfile5* &pISL){

	HRESULT hr;

	SetLastError(0);

	wchar_t* const clsid_MSFileReader_str = (wchar_t *)L"{1D23188D-53FE-4C25-B032-DC70ACDBDC02}";
	wchar_t* const iid_IXRawfile_str = (wchar_t *)L"{11B488A0-69B1-41FC-A660-FE8DF2A31F5B}";

	CLSID clsid_MSFileReader;
	CLSID clsid_IXRawfile;

	CLSIDFromString(clsid_MSFileReader_str, &clsid_MSFileReader);
	CLSIDFromString(iid_IXRawfile_str, &clsid_IXRawfile);

	DWORD dw = GetLastError();
	if(dw!=0){
		std::cout << "CLSIDFromString error: " << dw << ". ";
		getError();
		SetLastError(0);
		dw = 0;
	}

	//hr = CoInitializeEx(NULL, COINIT_MULTITHREADED);
	hr=CoInitialize(NULL);

	if(SUCCEEDED(hr)){
		std::cout << "CoInitialize success" << std::endl;

		hr = CoCreateInstance ( clsid_MSFileReader,		// CLSID of coclass
								NULL,					// not used - aggregation
								CLSCTX_INPROC_SERVER,	// type of server
								clsid_IXRawfile,		// IID of interface
								(void**) &pISL );		// Pointer to our interface pointer

		if ( SUCCEEDED ( hr ) ){

			std::cout << "CoCreateInstance success" << std::endl;

			wchar_t* wcharRawFile = (wchar_t*)wRawFile.c_str();
			std::wcout << wcharRawFile << std::endl;
			long nRet = pISL->Open( wcharRawFile );

			if(nRet != 0){
				std::cout << "Open error : " << nRet << ". " << std::endl;
				dw = GetLastError();
				if(dw!=0){
					getError();
					SetLastError(0);
					dw = 0;
				}
			}

			std::cout << "Opened RAW file" << std::endl;
		}
		else{
			std::cout << "CoCreateInstance failure" << std::endl;
			dw = GetLastError();
			if(dw!=0){
				std::cout << "CoCreateInstance error: " << dw << ". ";
				getError();
				SetLastError(0);
				dw = 0;
			}
			// Couldn't create the COM object.  hr holds the error code.
		}

	} else {
		std::cout << "CoInitialize failure" << std::endl;

		dw = GetLastError();
		if(dw!=0){
			std::cout << "CoInitialize error: " << dw << ". ";
			getError();
			SetLastError(0);
			dw = 0;
		}
	}
}

void getTimeIntPairs(std::vector<std::pair<double, double> > &timeIntPairs, BSTR &filter, std::wstring mzRange, float rtLower, float rtUpper, IXRawfile5* &pISL){

	wchar_t* wcharMzRange = &mzRange[0];
	VARIANT varChroData;
	VariantInit(&varChroData);
	VARIANT varPeakFlags;
	VariantInit(&varPeakFlags);
	long nArraySize = 0;
	double dStartTime = double(rtLower);
	double dEndTime = double(rtUpper);

	long nRet = pISL->GetChroData ( 0,
									0,
									0,
									filter,
									wcharMzRange,
									NULL,
									0.0,
									&dStartTime,
									&dEndTime,
									1, //Boxcar chromatogram smoothing
									9, //Smoothing value between 3-15 (must be odd)
									&varChroData,
									&varPeakFlags,
									&nArraySize );

	if( nArraySize ){
		// Get a pointer to the SafeArray
		SAFEARRAY FAR* psa = varChroData.parray;

		double* spectrumData;

		SafeArrayAccessData( psa, (void**)(&spectrumData) );
		std::pair<double, double> retrievedPair;
		for( long j=0; j<2*nArraySize; j=j+2 ){
			retrievedPair.first = spectrumData[j];
			retrievedPair.second = spectrumData[j+1];

			timeIntPairs.push_back(retrievedPair);
		}
		// Release the data handle
		SafeArrayUnaccessData( psa );

		SafeArrayDestroy(psa);
		psa = NULL;
	}

}

void peakHeight(std::vector< std::vector<std::pair<double, double> > > &timeIntIsotopes, peakInfo &bestPeak){
	std::vector<std::pair<double, double> > &isotope1 = timeIntIsotopes[0];

	int numIndices = isotope1.size();

	std::vector<int> listOfPeaks;

	for(int slidingWindowIndex=6; slidingWindowIndex<numIndices-6; slidingWindowIndex++){
		//First let's do a basic check that this point is at least at a local maximum about 3 points

		if(!(	(isotope1[slidingWindowIndex].second>=isotope1[slidingWindowIndex-1].second) &&
				(isotope1[slidingWindowIndex].second>=isotope1[slidingWindowIndex+1].second))){
			continue;
		}

		int higherThan = 0;
		for(int index=slidingWindowIndex-6; index<slidingWindowIndex; index++){
			if(isotope1[slidingWindowIndex].second > isotope1[index].second){
				higherThan++;
			}
		}
		if(higherThan>2){
			higherThan=0;
			for(int index=slidingWindowIndex+1; index<=slidingWindowIndex+6; index++){
				if(isotope1[slidingWindowIndex].second > isotope1[index].second){
					higherThan++;
				}
			}
			if(higherThan>2) listOfPeaks.push_back(slidingWindowIndex);
		}
	}

	int numPeaks = listOfPeaks.size();
	double bestArea = 0;
	double bestIntensity = 0;
	int bestIndex = 0;
	std::vector<double> bestFWHM = {0,0,0};
	double bestFWHMMins=0;

	for(int peakIndex=0; peakIndex<numPeaks; peakIndex++){

		std::vector<double> peakAreas;
		std::vector<double> FWHMs;
		double FWHMMins;

		for(int isotopeIndex=0; isotopeIndex<3; isotopeIndex++){

			int intensityIndex = listOfPeaks[peakIndex];
			double peakIntensity = timeIntIsotopes[isotopeIndex][intensityIndex].second;
			while(timeIntIsotopes[isotopeIndex][intensityIndex].second>peakIntensity/2 && intensityIndex>0){
				intensityIndex--;
			}
			int lowerIndex = intensityIndex;

			intensityIndex = listOfPeaks[peakIndex];
			while(timeIntIsotopes[isotopeIndex][intensityIndex].second>peakIntensity/2 && intensityIndex<numIndices){
				intensityIndex++;
			}
			int upperIndex = intensityIndex;

			double fwhm=0;

			double peakArea = getPeakArea(upperIndex, lowerIndex, numIndices, peakIntensity, listOfPeaks[peakIndex], fwhm, timeIntIsotopes, isotopeIndex);

			peakAreas.push_back(peakArea);
			FWHMs.push_back(fwhm);

			// Let's get the FWHM in minutes for the first isotope
			if(isotopeIndex==0){
				FWHMMins = fwhm;
			}
		}

		double totalArea = peakAreas[0] + peakAreas[1] + peakAreas[2];
		if(totalArea>bestArea && (peakAreas[1] > (peakAreas[0]*0.2)) ){
			bestArea=totalArea;
			bestIntensity = isotope1[listOfPeaks[peakIndex]].second;
			bestIndex = listOfPeaks[peakIndex];
			bestFWHM = FWHMs;
			bestFWHMMins = FWHMMins;
		}
	}

	double rtOfBestPeak = timeIntIsotopes[0][bestIndex].first;
	int degressOfFreedom=0;
	float chiSquared = calculateGaussPearson(rtOfBestPeak, bestFWHMMins, bestIntensity, timeIntIsotopes, degressOfFreedom);

	bestPeak.chiSquared = chiSquared;
	bestPeak.degressOfFreedom = degressOfFreedom;
	bestPeak.intensity = bestIntensity;
	bestPeak.area = bestArea;
	bestPeak.index = bestIndex;
	bestPeak.rt = rtOfBestPeak;
	bestPeak.FWHMs = bestFWHM;
}

void outputTimeIntPairs(std::vector<std::pair<double, double> > &timeIntPairs, int isotopeNum){
		std::ofstream datFile;
		std::string csvFileName = "timeInt"+std::to_string(isotopeNum)+".csv";
		const char *csvFileNameChar = csvFileName.c_str();
		datFile.open (csvFileNameChar);

		datFile << "time,int\n";
		for(int i=0; i<timeIntPairs.size(); i++){
			datFile << timeIntPairs[i].first << "," << timeIntPairs[i].second << "\n";
		}
		datFile.close();
}

void getError(){

	//Returns the last Win32 error, in string format. Returns an empty string if there is no error.

	//Get the error message, if any.
	DWORD errorMessageID = ::GetLastError();
	if(errorMessageID == 0)
		std::cout << "No error" << std::endl;

	LPSTR messageBuffer = NULL;
	size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
								 NULL, errorMessageID, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);

	std::string message(messageBuffer, size);

	//Free the buffer.
	LocalFree(messageBuffer);

	std::cout << message << std::endl;

}

float calculateGaussPearson(double rtOfBestPeak, double bestFWHMMins, double bestIntensity, std::vector< std::vector<std::pair<double, double> > > &timeIntIsotopes, int &degressOfFreedom){
	// This function calculates the Pearson's Chi-Squared value for the peptide's first isotopic peak's closeness to a Gaussian function

	// We will calculate the closeness-of-fit for each intensity value +/- 2-sigma of the mean
	double sigma = bestFWHMMins / (2*sqrt(2*log(2)));
	double minimumRt = rtOfBestPeak - 2*sigma;
	double maximumRt = rtOfBestPeak + 2*sigma;

	degressOfFreedom=0;
	float chiSquaredSum=0;
	int numPoints = timeIntIsotopes[0].size();
	for(int index=0; (index<numPoints); ++index){
		if(timeIntIsotopes[0][index].first > minimumRt && timeIntIsotopes[0][index].first < maximumRt){
			double observed = timeIntIsotopes[0][index].second;
			double expected = bestIntensity*exp(-(pow(timeIntIsotopes[0][index].first-rtOfBestPeak,2))/2*pow(sigma,2));
			chiSquaredSum += pow(observed - expected,2) / expected;
			++degressOfFreedom;
		}
	}

	return chiSquaredSum;
}

float calculateShapiroWilk(double rtOfBestPeak, double bestFWHMMins, double bestIntensity, std::vector< std::vector<std::pair<double, double> > > &timeIntIsotopes, int &degressOfFreedom){
	// This function calculates the Shaprio-Wilk test statistic for goodness of fit. It returns a float which is a p-value


}

double getPeakArea(int upperIndex, int lowerIndex, int numIndices, double peakIntensity, int peakIndex, double &fwhm, std::vector< std::vector<std::pair<double, double> > > &timeIntIsotopes, int isotopeIndex){

	double peakWidth=0, peakArea=0;
	//If the half maximum can't be found to the left, then we will use the half maximum to the right and multiple by 2, and visa-versa
	//If the half maximum can't be found to the left AND the right, then this suggests a peak width wider than our XIC window - we can treat this as a bad peak and report its areas as zero
	if(lowerIndex==0){
		if(upperIndex==numIndices-1){
			return peakArea;
		} else{
			peakWidth = (upperIndex-peakIndex)*2*60;
			fwhm = (timeIntIsotopes[isotopeIndex][upperIndex].first - timeIntIsotopes[isotopeIndex][peakIndex].first)*2;
		}
	} else if(upperIndex==numIndices-1){
		peakWidth = (peakIndex-lowerIndex)*2*60;
		fwhm = (timeIntIsotopes[isotopeIndex][peakIndex].first - timeIntIsotopes[isotopeIndex][lowerIndex].first)*2;
	} else{
		peakWidth = (upperIndex-lowerIndex)*60;
		fwhm = timeIntIsotopes[isotopeIndex][upperIndex].first - timeIntIsotopes[isotopeIndex][lowerIndex].first;
	}

	peakArea = peakWidth*peakIntensity;

	return peakArea;
}