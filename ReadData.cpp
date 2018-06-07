/*
 *  ReadDate.cpp
 *  
 *
 *  Created by Yifei Shi on 3/3/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

int main(int argc, char* argv[]) {
	
	int Nnode, Ndef;
	
	sscanf(argv[1], "%d", &Ndef);
	sscanf(argv[2], "%d", &Nnode);
	
	FILE* input1;
	
	FILE* input2;
	
	FILE* input3;
	
	FILE* output;
		

	//int ndata1 = 0, ndata2 = 0;
	
	for(int node = 1; node <= Nnode; node ++) {
		std::stringstream output_name4;
		output_name4 << "SKT_D=" << Ndef << "sec" << 7 << "node" << node << ".dat";
		const std::string &temp4=output_name4.str();
		const char *name4=temp4.c_str();
		input1 = fopen(name4, "r");
		
		
		if(input1 == NULL) {				
			printf("Cannot open File: D=%dnode%d.dat!\n", Ndef, node);
			exit(1);			
		}
		
		std::stringstream output_name1;
		output_name1 << "SKT_D=" << Ndef << "sec" << 8 << "node" << node << ".dat";
		const std::string &temp1=output_name1.str();
		const char *name1=temp1.c_str();
		input2 = fopen(name1, "r");
		
		
		if(input2 == NULL) {				
			printf("Cannot open File: D=%dnode%d.dat!\n", Ndef, node);
			exit(1);			
		}
		
	/*	std::stringstream output_name2;
		output_name2 << "SKT_D=" << Ndef << "sec" << 2 << "node" << node << ".dat";
		const std::string &temp2=output_name2.str();
		const char *name2=temp2.c_str();
		input3 = fopen(name2, "r"); 
		
		if(input3 == NULL) {				
			printf("Cannot open File: D=%dnode%d.dat!\n", Ndef, node);
			exit(1);			
		}  */
		
		std::stringstream output_name3;
		output_name3 << "SKT_D=" << Ndef << "MSnode" << node << ".dat";
		const std::string &temp3=output_name3.str();
		const char *name3=temp3.c_str();
		output = fopen(name3, "w");
		
		double dt[14];
		
		while(fscanf(input1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					 &dt[0], &dt[1], &dt[2], &dt[3], &dt[4], &dt[5], &dt[6], &dt[7], &dt[8], &dt[9], &dt[10], &dt[11], &dt[12], &dt[13]) != EOF ) {
			
			fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					dt[0], dt[1], dt[2], dt[3], dt[4], dt[5], dt[6], dt[7], dt[8], dt[9], dt[10], dt[11], dt[12], dt[13]);	
		}
		
		while((fscanf(input2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					 &dt[0], &dt[1], &dt[2], &dt[3], &dt[4], &dt[5], &dt[6], &dt[7], &dt[8], &dt[9], &dt[10], &dt[11], &dt[12], &dt[13]) != EOF )) {
			if(dt[0]<0.52)
			fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				   dt[0], dt[1], dt[2], dt[3], dt[4], dt[5], dt[6], dt[7], dt[8], dt[9], dt[10], dt[11], dt[12], dt[13]);	
		} 
				  
			  
	/*	while(fscanf(input1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					  &dt[0], &dt[1], &dt[2], &dt[3], &dt[4], &dt[5], &dt[6], &dt[7], &dt[8], &dt[9], &dt[10], &dt[11], &dt[12], &dt[13]) != EOF ) {
			if(dt[0] < 0.2)
			fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					dt[0], dt[1], dt[2], dt[3], dt[4], dt[5], dt[6], dt[7], dt[8], dt[9], dt[10], dt[11], dt[12], dt[13]);	
		}  */ 
		
		fclose(input1);  fclose(output);
	}
}
