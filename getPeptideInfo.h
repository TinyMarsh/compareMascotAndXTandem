#ifndef PESCAL_GETPEPTIDEINFO_H
#define PESCAL_GETPEPTIDEINFO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

struct peptideInfo{
	std::vector<double> massToChargesMascot;
	std::vector<double> massToChargesXTandem;
	std::vector<float> retentionTimesMascot;
	std::vector<float> retentionTimesXTandem;
	std::vector<int> charges;
	std::vector<std::pair<int, int> > dbIDs;
};

peptideInfo getPeptideInfo(std::string &rawFileName);

int getMascot(std::vector<std::string> &mascotSequences, std::vector<std::string> &mascotModPos, std::vector<std::string> &mascotCharges,std::vector<std::string> &mascotMz, std::string &rawFileName);
int getXTandem(std::vector<std::string> &xTandemSequences, std::vector<std::string> &xTandemModPos, std::vector<std::string> &xTandemCharges, std::vector<std::string> &xTandemMz, std::string &rawFileName);

void getMascotRTs(std::vector<double> &mascotRTs, int mascotRawFileIndex);
void getXTandemRTs(std::vector<double> &xTandemRTs, int xTandemRawFileIndex);

#endif //PESCAL_GETPEPTIDEINFO_H