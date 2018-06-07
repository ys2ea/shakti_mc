#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

int main(int argc, char * argv[]){
	int npo, N, D, sec, Ncol=15;
	
	int Nrow = 200;
	

	if(argc != 5) {
		printf("wrong useage!\n");
		//exit(1);
	}
	sscanf(argv[1], "%d", &npo);
	sscanf(argv[2], "%d", &N);
	sscanf(argv[3], "%d", &D);
	sscanf(argv[4], "%d", &sec);
	
	FILE * input;		
	        
	double DataArray[Ncol][Nrow];
	int Ntrue[Nrow];

	for(int i = 0; i< Nrow; i++) {
		Ntrue[i] = 0; 
		for(int k = 0; k < Ncol; k ++)
		DataArray[k][i] = 0.; }
	

	int np = npo;
	for(int i = 0; i < npo; i++)  {
		int ndata = 0;
		std::stringstream input_name;				
		input_name << "SKT_N=" << N << "D=" << D << "sec" << sec << "node" << i+1 << ".dat";					
		const string &temp1=input_name.str();
		const char *inname=temp1.c_str();				
		printf("filename: %s\n", inname);		        
		input = fopen(inname, "r");			
			
		if(input == NULL) {				
			printf("Cannot open File: D=%dnode%d.dat!\n", D, i+1);
			np --;
			
		}
		else {
		double dt[Ncol];

		while(fscanf(input, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf%lf\t%lf\n",
					 &dt[0], &dt[1], &dt[2], &dt[3], &dt[4], &dt[5], &dt[6], &dt[7], &dt[8], &dt[9], &dt[10], &dt[11], &dt[12], &dt[13], &dt[14]) != EOF) {

			if(DataArray[0][ndata] == 0.) DataArray[0][ndata] = dt[0];
			
			for(int k = 1; k < Ncol; k++) {
				DataArray[k][ndata] += dt[k];
			}
			Ntrue[ndata++] ++;
		}
			
		fclose(input);}
		
	}

		
	FILE* output;

	std::stringstream output_name;
	output_name << "DataN=" << N << "D=" << D << "sec=" << sec <<".dat";
	const string &temp=output_name.str();
	const char *name=temp.c_str();
	output = fopen(name, "w");
	//printf("file%d\n", L);
		//fprintf(output, "#\n");

	for(int i = 0; i < Nrow; i++) {
		int dn = Ntrue[i];
		if(dn != 0)
		fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				DataArray[0][i], DataArray[1][i]/dn,  DataArray[2][i]/dn,  DataArray[3][i]/dn,  DataArray[4][i]/dn,  DataArray[5][i]/dn,  
				DataArray[6][i]/dn, DataArray[7][i]/dn,  DataArray[8][i]/dn,  DataArray[9][i]/dn,  DataArray[10][i]/dn,  DataArray[11][i]/dn, DataArray[12][i]/dn, DataArray[13][i]/dn, DataArray[14][i]/dn);
	}
	fclose(output);
}
