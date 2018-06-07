/*
 *  DfShakti.cpp
 *  
 *
 *  Created by Yifei Shi on 10/28/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "DfShakti.h"
#include "Functions.h"

int main(int argc, char* argv[]) {

	int N, Ndef, Nthermal, Nmeasure, Nthermal1, Nmeasure1, CorrLength, rank, sec;
	double Ez4[4];
	double Ez3[5];

	if(argc != 3) {
		std::cout << "Wrong usage!\n"; }
	
	sscanf(argv[1], "%d", &sec);
	sscanf(argv[2], "%d", &rank);
	
	double MaxT, MinT, DeltaT;
		
	// parameters form file
	std::ifstream input;
	input.open("param.dat");
	if(input==NULL) { 
		std::cout << "FEED ME THE FILE!!\n";
		exit(1); }
	
	input >> N;
	input >> Ndef;
	input >> CorrLength;
	input >> MaxT;
	input >> MinT;
	input >> Nthermal;
	input >> Nmeasure;
	input >> Ez4[0];
	input >> Ez4[1];
	input >> Ez4[2];
	input >> Ez4[3];
	input >> Ez3[0];
	input >> Ez3[1];
	input >> Ez3[2];
	input >> Ez3[3];
	input >> Ez3[4];
	
	FILE* output;
	
	std::stringstream output_name;
	output_name << "SKT_N=" << N << "D=" << Ndef << "sec" << sec << "node" << rank << ".dat";
	const std::string &temp=output_name.str();
	const char *name=temp.c_str();
	output = fopen(name, "w");
		
//	fprintf(output, "#T\t ExceptRatio\t SpecificHeat\t V4a\t V4b\t V4c\t V4d\t V3A\t V3B\t V3C\n");
	
	if(MaxT < MinT) {
		std::cout << "T input error!!!!!!\n";
		exit(1); 
	}
	
	DfShakti SpinIce(N, Ndef, rank, MinT, Ez4[0], Ez4[1], Ez4[2], Ez4[3], Ez3[0], Ez3[1], Ez3[2], Ez3[3], Ez3[4]);
	
	double Eground = Ez4[0]*4.*N*N + Ez3[0]*4.*N*N + Ez3[1]*4.*N*N;	//ground state energy
	
	DeltaT = 0.5;
	
	for(double T = MaxT; T >= MinT; T -= DeltaT) {
		
		double vertex3[3] = {0., 0., 0.};			//averge vertex count
		double vertex4[4] = {0., 0., 0., 0.};
		double Isv1 = 0., Isv2 = 0.;
	
		if(T>30.) {DeltaT = 1; Nthermal1 = Nthermal; Nmeasure1 = Nmeasure; }
		else if(T>5.) {DeltaT = 0.5; Nthermal1 = Nthermal; Nmeasure1 = Nmeasure; }
		else if(T>1.6) {DeltaT = 0.1; Nthermal1 = Nthermal; Nmeasure1 = Nmeasure; CorrLength = 20; }
		else if(T>0.6) {DeltaT = 0.05; Nthermal1 = 2*Nthermal; Nmeasure1 = 2*Nmeasure; CorrLength = 32;}
		else {DeltaT = 0.03; Nthermal1 = 2*Nthermal; Nmeasure1 = 2*Nmeasure; CorrLength = 100; }
		
		double EAverage=0., EsqAverage=0., ExceptRatio=0., Ej = 0.;	//energy, specific heat, and exceptance rate

		SpinIce.ResetTemp(T);
		
		//if(T != MinT) 
			//SpinIce.ReadVertex("VertexConfig.dat");
		
		for(int i = 0; i < Nthermal1; i++) {
			for(int k = 0; k < 4*N*N*CorrLength; k++) {
				SpinIce.Spinflip(0);
				SpinIce.Spinflip(0);
				SpinIce.Spinflip(0);
				SpinIce.Spinflip(0);
				SpinIce.DoubleLongflip();
				SpinIce.DoubleLongflip();
				SpinIce.Z3flip();
				//SpinIce.Z4flip();
				if(k < (4*N*N-Ndef)*CorrLength)
					SpinIce.Spinflip(1);
			}
			
		}//end thermalize loop
			
	
		for(int j = 0; j < Nmeasure1; j++) {
			Ej = 0.;
			for(int k =0 ; k < 4*N*N*CorrLength; k++) {
				SpinIce.Spinflip(0);
				SpinIce.Spinflip(0);
				SpinIce.Spinflip(0);
				SpinIce.Spinflip(0);
				SpinIce.Z3flip();
				//SpinIce.Z4flip();
				Ej += double(SpinIce.DoubleLongflip());
				Ej += double(SpinIce.DoubleLongflip());
				
				//ExceptRatio = ExceptRatio*double(k+j*4*N*N*CorrLength)/double(k+j*4*N*N*CorrLength+1) + SpinIce.Spinflip(0)/double(k+j*4*N*N*CorrLength+1);
				if(k < (4*N*N-Ndef)*CorrLength)
					SpinIce.Spinflip(1);
			}
			
			ExceptRatio = ExceptRatio*double(j)/double(j+1) + Ej/8./N/N/CorrLength/double(j+1);
		
			//measurements here
			for(int tp4 = 0; tp4 < 4; tp4++) {
				vertex4[tp4] = vertex4[tp4]*double(j)/double(j+1)+SpinIce.Showz4Vertex(tp4)/double(j+1);
			}
			for(int tp3 = 0; tp3 < 3; tp3++) {
				vertex3[tp3] = vertex3[tp3]*double(j)/double(j+1)+SpinIce.Showz3Vertex(tp3)/double(j+1);
			}
								
			int ism = SpinIce.IsingMv2();
			Isv1 = Isv1*double(j)/double(j+1) + ism*ism/double(j+1)/double(64*N*N*N*N);
			Isv2 = Isv2*double(j)/double(j+1) + abs(ism)/double(j+1)/double(8*N*N);

			EAverage = EAverage*double(j)/double(j+1) + (SpinIce.ShowEnergy()-Eground)/double(j+1);
			EsqAverage = EsqAverage*double(j)/double(j+1) + (SpinIce.ShowEnergy()-Eground)*(SpinIce.ShowEnergy()-Eground)/double(j+1);
		
		}//end measure loop 
		
		//SpinIce.SaveVertex("VertexConfig.dat");
		SpinIce.CountVertex();
		int * Dn = SpinIce.Dn56();
				
		fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n", \
				T, ExceptRatio, EsqAverage-EAverage*EAverage, vertex4[0], vertex4[1], vertex4[2], vertex4[3], vertex3[0], vertex3[1], vertex3[2], Isv1, Isv2, Dn[0], Dn[1]);
		
		
		fflush(output);
		
		std::cout << "D= " << Ndef << "sec" << sec << "T= " << T << " Finished " << ExceptRatio << "\n";
		fflush(stdout);
		
	}
	
	/* std::stringstream ss1, ss2;
	ss1 << "SS1node" << rank << ".dat";
	ss2 << "SS2node" << rank << ".dat";
	
	const std::string &s1name=ss1.str();
	const char *name1=s1name.c_str();
	
	const std::string &s2name=ss2.str();
	const char *name2=s2name.c_str();
	
	SpinIce.SnapShot(name1, name2); */
	
	fclose(output);
	
}
