#include "getQuant.h"
#include "getPeptideInfo.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

/*
This program gets the XICs for the three isotopes of each peptide and outputs the quant information to a file, peakAreas.csv. It does this for
*/

int main(int argc, char* argv[]) {

	std::string rawFileName = argv[1];

	peptideInfo newPeptideInfo = getPeptideInfo(rawFileName);

	//output peptideInfo
	std::ofstream peptideInfofile("peptideInfo.csv");

	int numPeptides = newPeptideInfo.dbIDs.size();
	for(int peptide=0; peptide<numPeptides; ++peptide){
		peptideInfofile << newPeptideInfo.dbIDs[peptide].first << "," << newPeptideInfo.dbIDs[peptide].second << "," << newPeptideInfo.massToChargesMascot[peptide] << "," << newPeptideInfo.massToChargesXTandem[peptide] << "," << newPeptideInfo.retentionTimesMascot[peptide] << "," << newPeptideInfo.retentionTimesXTandem[peptide] << "," << newPeptideInfo.charges[peptide] << "\n";
	}
	peptideInfofile.close();

	std::pair<std::vector<double>,std::vector<double> > allQuants;
	std::pair<std::vector<std::vector<double> >,std::vector<std::vector<double> > > allFWHMs;
	std::pair<std::vector<float>,std::vector<float> > allChiSquareds;
	std::pair<std::vector<int>,std::vector<int> > allDegreesOfFreedoms;
	std::pair<std::vector<float>,std::vector<float> > allIntensities;
	std::pair<std::vector<float>,std::vector<float> > allBestPeakRts;
	getQuant(allQuants, allFWHMs, allChiSquareds, allDegreesOfFreedoms, allIntensities, allBestPeakRts, rawFileName, newPeptideInfo);

	//output quantitation info to "peakAreas.csv"
	std::ofstream file("peakAreas.csv");
	file << "Mascot" << "," << "X!Tandem" << "\n";
	for(int peptide=0; peptide<numPeptides; ++peptide){
		file << allQuants.first[peptide] << "," << allQuants.second[peptide] << "\n";
	}
	file.close();

	//output FWHMs to "FWHMs.csv"
	std::ofstream fileFWHM("FWHMs.csv");
	fileFWHM << "Mascot" << "," << "X!Tandem" << "\n";
	for(int peptide=0; peptide<numPeptides; ++peptide){
		fileFWHM << allFWHMs.first[peptide][0] << "," << allFWHMs.first[peptide][1] << "," << allFWHMs.first[peptide][2] << "," << allFWHMs.second[peptide][0] << "," << allFWHMs.second[peptide][1] << "," << allFWHMs.second[peptide][2] << "," << "\n";
	}
	fileFWHM.close();

	//output chi-squareds to "ChiSquared.csv"
	std::ofstream fileChiSquared("ChiSquared.csv");
	fileChiSquared << "Mascot" << "," << "X!Tandem" << "\n";
	for(int peptide=0; peptide<numPeptides; ++peptide){
		fileChiSquared <<	allChiSquareds.first[peptide] << "," <<
							allChiSquareds.second[peptide] << ", " <<
							allDegreesOfFreedoms.first[peptide] << ", " <<
							allDegreesOfFreedoms.second[peptide] << "\n";
	}
	fileChiSquared.close();

	//output intensities to "intensities.csv"
	std::ofstream fileIntensities("intensities.csv");
	fileIntensities << "Mascot" << "," << "X!Tandem" << "\n";
	for(int peptide=0; peptide<numPeptides; ++peptide){
		fileIntensities <<	allIntensities.first[peptide] << "," <<
							allIntensities.second[peptide] << "\n";
	}
	fileIntensities.close();

	//output best peak retention times to "bestPeakRts.csv"
	std::ofstream fileBestRts("bestPeakRts.csv");
	fileBestRts << "Mascot" << "," << "X!Tandem" << "\n";
	for(int peptide=0; peptide<numPeptides; ++peptide){
		fileBestRts <<	allBestPeakRts.first[peptide] << "," <<
						allBestPeakRts.second[peptide] << "\n";
	}
	fileBestRts.close();

	return 0;
}