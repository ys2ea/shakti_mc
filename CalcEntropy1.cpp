/*
 *  CalcEntropy.cpp
 *  
 *  Do the integral to calculate entropy, using simpson's rule
 * 
 *  Created by Yifei Shi on 11/5/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */


#include<fstream>
#include<iostream>
#include<cmath>
#include<sstream>
#include<stdlib.h>
#include<vector>


int main(int argc, char* argv[]) {
	const std::string agmt = argv[1];
	const char* inname = agmt.c_str();
	
	FILE* input;
	
	std::vector<double> Scurve;
	
	input = fopen(inname,"r");
	
	double T, AcRate, C, V4a, V4b, V4c, V4d, V3A, V3B, V3C, V3D, V3E, prevC, prevT, ppC, ppT, integral;
	
	double n5, n6;
	
	prevC = 0.;
	prevT = 0.;
	
	ppC = 0.;
	ppT = 0.;
	
	integral = 0.;
	
	int Count = 0;
	
	while(fscanf(input, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				  &T, &AcRate, &C, &V4a, &V4b, &V4c, &V4d, &V3A, &V3B, &V3C, &V3D, &V3E, &n5, &n6) != EOF) {
		
	
		 if(ppT != 0.) {
			if(ppT - prevT == prevT - T) {
				if(Count == 2) {
					integral += (ppT - T)*(ppC/ppT/ppT/ppT+C/T/T/T+4.*prevC/prevT/prevT/prevT)/6.; 
					Count = 0;
					Scurve.push_back(T);
					Scurve.push_back(integral);
				}

			}
			
			else {
				if(Count == 2) {
					integral += (ppT - prevT)*(ppC/ppT/ppT/ppT+prevC/prevT/prevT/prevT)/2.;
					Count = 1;
					Scurve.push_back(prevT);
					Scurve.push_back(integral);

				}
				
			}
		}  
		
if(integral > 99999999) std::cout << T << "\n";
		//if(prevT != 0.)
		//integral += (prevT - T)*(prevC/prevT/prevT/prevT+C/T/T/T)/2.;
		
		ppT = prevT;
		ppC = prevC;
		
		prevT = T;
		prevC = C;	
		
		Count ++;
		
	}
	std::cout << integral << "\n";
	
	FILE* output;
	output = fopen("S_Curve.dat", "w");
	int Ndata = Scurve.size()/2;
	for(int i = 0; i < Ndata; i ++)
		fprintf(output, "%lf\t%lf\n", Scurve[2*i], Scurve[2*i+1]);
	fclose(output);
	
}
