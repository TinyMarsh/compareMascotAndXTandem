#ifndef PESCAL_QUANTITATION_H
#define PESCAL_QUANTITATION_H

#include <vector>
#include <string>
#include <tchar.h>
#include <Windows.h>
#include <cmath>
#include <limits>
#include "getPeptideInfo.h"

#include "XRawfile2_x64.h"

struct peakInfo{
    double intensity;
    double area;
    int index;
    float rt;
    std::vector<double> FWHMs;
    float chiSquared;
    int degressOfFreedom;
};

void getQuant(	std::pair<std::vector<double>,std::vector<double> > &allQuants,
				std::pair<std::vector<std::vector<double> >,std::vector<std::vector<double> > > &allFWHMs,
				std::pair<std::vector<float>,std::vector<float> > &allChiSquareds,
				std::pair<std::vector<int>,std::vector<int> > &allDegreesOfFreedoms,
                std::pair<std::vector<float>,std::vector<float> > &allIntensities,
                std::pair<std::vector<float>,std::vector<float> > &allBestPeakRts,
				std::string &fullFileName,
				peptideInfo &newPeptideInfo);

void openRawFile(std::wstring wRawFile, IXRawfile5* &pISL);

void getTimeIntPairs(std::vector<std::pair<double, double> > &timeIntPairs, BSTR &filter, std::wstring mzRange, float rtLower, float rtUpper, IXRawfile5* &pISL);

void peakHeight(std::vector <std::vector<std::pair<double, double> > > &timeIntIsotopes, peakInfo &bestPeak);

void outputTimeIntPairs(std::vector<std::pair<double, double> > &timeIntPairs, int isotopeNum);

void getError();

float calculateGaussPearson(double rtOfBestPeak, double bestFWHM, double bestIntensity, std::vector< std::vector<std::pair<double, double> > > &timeIntIsotopes, int &degressOfFreedom);

double getPeakArea(int upperIndex, int lowerIndex, int numIndices, double peakIntensity, int peakIndex, double &fwhm, std::vector< std::vector<std::pair<double, double> > > &timeIntIsotopes, int isotopeIndex);

#endif //PESCAL_QUANTITATION_H
