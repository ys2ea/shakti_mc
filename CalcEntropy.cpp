/*
 *  CalcEntropy.cpp
 *  
 *
 *  Created by Yifei Shi on 11/5/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */


#include<fstream>
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<vector>


int main() {
	FILE* input;
	
	input = fopen("V1.dat","r");
	
	double T, AcRate, C, V4a, V4b, V4c, V4d, V3A, V3B, V3C, prevC, prevT, integral;
	
	prevC = 0.;
	prevT = 0.;
	
	integral = 0.;
	
	while(fscanf(input, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				  &T, &AcRate, &C, &V4a, &V4b, &V4c, &V4d, &V3A, &V3B, &V3C) != EOF) {
		
		
		if(T >= 0.3) {
			integral += prevC/T/T/T*(T-prevT);
			printf("%lf\t", T); 
		}

		prevT = T;
		prevC = C;
				
		
	}
	
	std::cout << integral << "\n";
	
}
