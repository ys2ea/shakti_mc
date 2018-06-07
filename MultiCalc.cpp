/*
 *  MultiCalc.cpp
 *  
 *
 *  Created by Yifei Shi on 3/3/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include<fstream>
#include<iostream>
#include<cmath>
#include<sstream>
#include<stdlib.h>
#include<vector>


int main(int argc, char* argv[]) {

	FILE* input;
	FILE* output;
	
	int Ndef, Nnode, Neffnode = 0;
	int sec=0;
	if(argc == 3) {
		sscanf(argv[1], "%d", &Ndef);
		sscanf(argv[2], "%d", &Nnode);
	}
	else if(argc == 4) {
		sscanf(argv[1], "%d", &Ndef);
		sscanf(argv[2], "%d", &Nnode);
		sscanf(argv[3], "%d", &sec);
	}
	
	else {
		std::cout << "Wrong usage!\n";
		exit(1); }
	
	double* Sarray = (double*) malloc(sizeof(double)*Nnode);
	
	for(int i = 0; i < Nnode; i ++)
		Sarray[i] = 0.0;
	
	std::stringstream output_name;
	output_name << "S_D=" << Ndef << ".dat";
	const std::string &temp1=output_name.str();
	const char *name1=temp1.c_str();
	output = fopen(name1, "w");		

	double Saverage = 0.0;
	
	for(int rank = 1; rank < Nnode +1; rank ++) {
		std::stringstream input_name;
		if(argc == 3)
			input_name << "SKT_D=" << Ndef << "node" << rank << ".dat";
		else 
			input_name << "SKT_D=" << Ndef << "sec" << sec << "node" << rank << ".dat";

		const std::string &temp=input_name.str();
		const char *name=temp.c_str();
		input = fopen(name, "r");		
		
		//std::vector<double> Scurve;
				
		double T, AcRate, C, V4a, V4b, V4c, V4d, V3A, V3B, V3C, V3D, V3E, prevC, prevT, ppC, ppT;
		
		double n5, n6;
		
		prevC = 0.;
		prevT = 0.;
		
		ppC = 0.;
		ppT = 0.;
				
		int Count = 0;
		
		while(fscanf(input, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					 &T, &AcRate, &C, &V4a, &V4b, &V4c, &V4d, &V3A, &V3B, &V3C, &V3D, &V3E, &n5, &n6) != EOF) {
			
			
			if(ppT != 0.) {
				if(ppT - prevT == prevT - T) {
					if(Count == 2) {
						Sarray[rank-1] += (ppT - T)*(ppC/ppT/ppT/ppT+C/T/T/T+4.*prevC/prevT/prevT/prevT)/6.; 
						Count = 0;
						//Scurve.push_back(T);
						//Scurve.push_back(integral);
					}
					
				}
				
				else {
					if(Count == 2) {
						Sarray[rank-1] += (ppT - prevT)*(ppC/ppT/ppT/ppT+prevC/prevT/prevT/prevT)/2.;
						Count = 1;
						//Scurve.push_back(prevT);
						//Scurve.push_back(integral);
						
					}
					
				}
			}  
			
			//if(prevT != 0.)
			//integral += (prevT - T)*(prevC/prevT/prevT/prevT+C/T/T/T)/2.;
			
			ppT = prevT;
			ppC = prevC;
			
			prevT = T;
			prevC = C;	
			
			Count ++;
			
		}
		
		fprintf(output, "%d\t%lf\n", rank, Sarray[rank-1]);
		n6 = 0.;
		if(T<=0.13) 
		{ Saverage += Sarray[rank-1]+n6*log(2.); Neffnode += 1;}

	/*	FILE* output;
		output = fopen("S_Curve.dat", "w");
		int Ndata = Scurve.size()/2;
		for(int i = 0; i < Ndata; i ++)
			fprintf(output, "%lf\t%lf\n", Scurve[2*i], Scurve[2*i+1]);
		fclose(output); */
	} //end for rank	
	
	Saverage /= double(Neffnode);
	
	std::cout << "Average Entropy for " << Ndef << " defects:" << log(2.)-Saverage/double(1280-Ndef) << "\n";
}

