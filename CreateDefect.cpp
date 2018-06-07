/*
 *  CreateDefect.cpp
 *  
 *
 *  Created by Yifei Shi on 1/30/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */
#include<fstream>
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<vector>
#include <sstream>
#include "MersenneTwister.h"

int main(int argc, char* argv[]) {
	int N_, tDf_;
	
	if(argc != 3) {
		std::cout << "Wrong usage!\n"; }
	
	sscanf(argv[1], "%d", &N_);
	sscanf(argv[2], "%d", &tDf_);
	
	std::vector<int> AvlbSites;
	
	std::vector<int> Defects_;
	
	MTRand ran;		//random number generator
	
	int NAsites = 4*N_*N_; //# of available sites
	
	AvlbSites.resize(4*N_*N_); // list of av sites
	
	Defects_.resize(tDf_);
	
	for(int i = 0; i < 4*N_*N_; i ++)
		AvlbSites[i] = i;
	
	for(int j = 0; j < tDf_; j ++) {
		int randidx = ran.randDblExc()*NAsites;
		
		Defects_[j] = AvlbSites[randidx];
		
		AvlbSites.erase(AvlbSites.begin() + randidx);
		
		NAsites --;
	}
	
	FILE* output;
	output = fopen("Defects.dat", "w");
	
	for(int i = 0; i < tDf_; i ++)
		fprintf(output, "%d\n", Defects_[i]);
	
	fclose(output);
}
