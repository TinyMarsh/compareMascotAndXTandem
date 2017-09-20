#include "getPeptideInfo.h"

peptideInfo getPeptideInfo(std::string &rawFileName){

	//Create our struct which is going to store all our info
	peptideInfo newPeptideInfo;
	
	//We will fetch all the sequence, modification position, and charge for each peptide and store them in vectors

	//...said vectors...
	std::vector<std::string> mascotSequences;
	std::vector<std::string> mascotModPos;
	std::vector<std::string> mascotCharges;
	std::vector<std::string> mascotMz;

	//this function fills said vectors
	int mascotRawFileIndex = getMascot(mascotSequences, mascotModPos, mascotCharges, mascotMz, rawFileName);

	//more of said vectors, but these will hold the X!Tandem peptides
	std::vector<std::string> xTandemSequences;
	std::vector<std::string> xTandemModPos;
	std::vector<std::string> xTandemCharges;
	std::vector<std::string> xTandemMz;

	//this is another function filling said vectors, but this time for X!Tandem results
	int xTandemRawFileIndex = getXTandem(xTandemSequences, xTandemModPos, xTandemCharges, xTandemMz, rawFileName);

	//Let's get the retention times of each matched peptide
	std::vector<double> mascotRTs;
	std::vector<double> xTandemRTs;

	getMascotRTs(mascotRTs, mascotRawFileIndex);
	getXTandemRTs(xTandemRTs, xTandemRawFileIndex);

	//We'll store the peptide db IDs of peptides which are identical between each identification software in this handy vector of type pair
	std::vector<std::pair<int, int> > dbIDs; 

	//Now we will fill the above handy vector. For every mascot sequence, modpos, and charge, check if there exists an xTandem equivalent. If there is, store the db ID from Mascot and X!Tandem in an int pair
	//While we're here we will fill our newPeptideInfo struct
	int xTandemIndex=0;
	for(int i=0, size=mascotSequences.size(); i<size; ++i){
		int compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
		while(compare>0 && xTandemIndex<xTandemSequences.size()){
			++xTandemIndex;
			compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
		}
		bool found=0;
		int tempIndex=xTandemIndex;
		while(compare==0&&!found){
			if(mascotModPos[i]==xTandemModPos[tempIndex] && mascotCharges[i]==xTandemCharges[tempIndex]){
				found=1;
				std::pair<int, int> indexes;
				indexes.first = i+1;
				indexes.second = tempIndex+1;
				dbIDs.push_back(indexes);

				//this is the bit where we store the info in our struct
				newPeptideInfo.charges.push_back(std::stoi(mascotCharges[i]));
				newPeptideInfo.massToChargesMascot.push_back(std::stod(mascotMz[i]));
				newPeptideInfo.massToChargesXTandem.push_back(std::stod(xTandemMz[tempIndex]));
				newPeptideInfo.retentionTimesMascot.push_back(mascotRTs[i]);
				newPeptideInfo.retentionTimesXTandem.push_back(xTandemRTs[tempIndex]);
			}
			if(tempIndex<xTandemSequences.size()){
				++tempIndex;
				compare = mascotSequences[i].compare(xTandemSequences[tempIndex]);
			} else{
				break;
			}
		}
	}	

	newPeptideInfo.dbIDs = dbIDs;

	return newPeptideInfo;
}

int getMascot(std::vector<std::string> &mascotSequences, std::vector<std::string> &mascotModPos, std::vector<std::string> &mascotCharges, std::vector<std::string> &mascotMz, std::string &rawFileName){

	std::ifstream file("F:\\data\\CTAM\\analysis\\mascotPosOfMod\\combiPeptData.csv");

	std::string value;

	//We need to find the index of the raw file name in the column headers so we know where to grab info from later
	bool found=0;
	int index=0;
	while(!found){
		std::getline(file, value, ',');
		if(value.compare(rawFileName)==0){
			found=1;
		}
		++index;
	}
	index=index-35;

	//first getline just moves the cursor past the first row, since we don't want the headers
	std::getline(file, value, '\n');

	for(int peptide=0; peptide<45187; ++peptide){

		//skip the first five quatation marks so that we get the sequence
		for(int i=0; i<5; ++i){
			std::getline(file, value, '\"');
		}

		//get the sequence
		std::getline(file, value, '\"');
		mascotSequences.push_back(value);

		//skip the next 3 quatation marks to get the m/z
		for(int i=0; i<3; ++i){
			std::getline(file, value, '\"');
		}

		std::getline(file, value, '\"');
		mascotMz.push_back(value);

		//skip the next 3 quatation marks to get the charge
		for(int i=0; i<3; ++i){
			std::getline(file, value, '\"');
		}

		//get the charge
		std::getline(file, value, '\"');
		mascotCharges.push_back(value);

		//skip the next quatation mark to get the modPos string
		std::getline(file, value, '\"');

		//get the modPos string
		std::getline(file, value, '\"');
		mascotModPos.push_back(value);

		std::getline(file, value, '\n');
	}

	return index;
}

int getXTandem(std::vector<std::string> &xTandemSequences, std::vector<std::string> &xTandemModPos, std::vector<std::string> &xTandemCharges, std::vector<std::string> &xTandemMz, std::string &rawFileName){


	std::ifstream file("F:\\data\\CTAM\\analysis\\peptideShakerPosOfMod\\combiPeptData.csv");

	std::string value;

	//We need to find the index of the raw file name in the column headers so we know where to grab info from later
	bool found=0;
	int index=0;
	while(!found){
		std::getline(file, value, ',');
		if(value.compare(rawFileName)==0){
			found=1;
		}
		++index;
	}
	index=index-35;

	//first getline just moves the cursor past the first row, since we don't want the headers
	std::getline(file, value, '\n');

	for(int peptide=0; peptide<59410; ++peptide){

		//skip the first five quatation marks so that we get the sequence
		for(int i=0; i<5; ++i){
			std::getline(file, value, '\"');
		}

		//get the sequence
		std::getline(file, value, '\"');
		xTandemSequences.push_back(value);

		//skip the next 3 quatation marks to get the m/z
		for(int i=0; i<3; ++i){
			std::getline(file, value, '\"');
		}

		std::getline(file, value, '\"');
		xTandemMz.push_back(value);

		//skip the next 3 quatation marks to get the charge
		for(int i=0; i<3; ++i){
			std::getline(file, value, '\"');
		}

		//get the charge
		std::getline(file, value, '\"');
		xTandemCharges.push_back(value);

		//skip the next quatation mark to get the modPos string
		std::getline(file, value, '\"');

		//get the modPos string
		std::getline(file, value, '\"');

		std::size_t found=0;
		while(found!=std::string::npos){
			found = value.find("4", found+1);
			if(found!=std::string::npos){
				value[found] = '3';
			}
		}

		xTandemModPos.push_back(value);

		std::getline(file, value, '\n');
	}

	return index;
}

void getMascotRTs(std::vector<double> &mascotRTs, int mascotRawFileIndex){
	std::ifstream file("F:\\data\\CTAM\\analysis\\mascotPosOfMod\\calculatedRTs.csv");

	std::string value;

	//not interested in first line
	std::getline(file, value, '\n');

	for(int i=0; i<45187; ++i){
		for(int j=0; j<mascotRawFileIndex; ++j){
			std::getline(file, value, '\"');
			std::getline(file, value, '\"');
		}

		//get quant value
		std::getline(file, value, '\"');
		std::getline(file, value, '\"');

		mascotRTs.push_back(std::stod(value));

		std::getline(file, value, '\n');
	}
}

void getXTandemRTs(std::vector<double> &xTandemRTs, int xTandemRawFileIndex){
	std::ifstream file("F:\\data\\CTAM\\analysis\\peptideShakerPosOfMod\\calculatedRTs.csv");

	std::string value;

	//not interested in first line
	std::getline(file, value, '\n');

	for(int i=0; i<59410; ++i){
		for(int j=0; j<xTandemRawFileIndex; ++j){
			std::getline(file, value, '\"');
			std::getline(file, value, '\"');
		}

		//get quant value
		std::getline(file, value, '\"');
		std::getline(file, value, '\"');

		xTandemRTs.push_back(std::stod(value));

		std::getline(file, value, '\n');
	}
}