/*
 *  Shakti.cpp
 *  This code used parallel tempering. This is a copy of the old version which donot use replica
 *
 *  Edited on 08/20/2016
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "DfShakti.h"
#include "Functions.h"

int main(int argc, char* argv[]) {

	int N, Ndef, Nthermal, Nmeasure, Nthermal1, Nmeasure1, CorrLength, rank, sec;
	int M = 10; //number of templates
	double Ez4[4];
	double Ez3[5];

	if(argc != 3) {
		std::cout << "Wrong usage!\n"; }
	
	sscanf(argv[1], "%d", &sec);
	sscanf(argv[2], "%d", &rank);
	
	double MaxT, MinT, DeltaT= 0.02;
		
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
	
	MTRand ran;
	
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
	
	std::vector<DfShakti> SIarray;
	
	std::vector<double> Betaarray;
	
	std::vector<int> Idx;
			
    //Temps from minT to MaxT, every M Ts are grouped together, the last group can have less than M Ts.
	double T = MaxT;
	bool Isfirst = 1;
	while(T > MinT) {
		Betaarray.clear();
		Idx.clear();
		
		for(int i = 0; i < M; i ++) {
			Betaarray.push_back(1./T);
			Idx.push_back(i);
			if(Isfirst == 1)
				SIarray.push_back(DfShakti(N, Ndef, rank, T, Ez4[0], Ez4[1], Ez4[2], Ez4[3], Ez3[0], Ez3[1], Ez3[2], Ez3[3], Ez3[4])); 
			else {
				if(i!=M-1) SIarray[i].Copy(SIarray[M-1]);
				SIarray[i].ResetTemp(T);
			}
			T -= DeltaT;
		}
		
		Isfirst = 0;
		double vertex3[5][10];			//averge vertex count
		double vertex4[4][10];
		
		double EAverage[10], EsqAverage[10], ExceptRatio[10], Ej[10];	//energy, specific heat, and exceptance rate
		for(int ti = 0; ti < 10; ti ++) {
			vertex3[0][ti] = 0;
			vertex3[1][ti] = 0;
			vertex3[2][ti] = 0;
			vertex3[3][ti] = 0;
			vertex3[4][ti] = 0;
			vertex4[0][ti] = 0;
			vertex4[1][ti] = 0;
			vertex4[2][ti] = 0;
			vertex4[3][ti] = 0;
			
			EAverage[ti] = 0.;
			EsqAverage[ti] = 0.;
			ExceptRatio[ti] = 0.;
			Ej[ti] = 0.;
		}
		
		for(int ti = 0; ti < M; ti ++) {
			for(int i = 0; i < Nthermal; i++) {
				for(int k = 0; k < 4*N*N*CorrLength; k++) {
					SIarray[ti].Spinflip(0);
					SIarray[ti].Spinflip(0);
					SIarray[ti].Spinflip(0);
					SIarray[ti].Spinflip(0);
					SIarray[ti].Doubleflip();
					SIarray[ti].Doubleflip();
					SIarray[ti].Z3flip();
					SIarray[ti].Z4flip();
					if(k < (4*N*N-Ndef)*CorrLength)
						SIarray[ti].Spinflip(1);
				}
				
			}//end thermalize loop
		}
		
		
		for(int j = 0; j < Nmeasure; j++) {
			for(int ti = 0; ti < M; ti ++) {
				for(int k =0 ; k < 4*N*N*CorrLength; k++) {
					SIarray[ti].Spinflip(0);
					SIarray[ti].Spinflip(0);
					SIarray[ti].Spinflip(0);
					SIarray[ti].Spinflip(0);
					
					SIarray[ti].Z4flip();
					SIarray[ti].Z3flip();
					
					//ExceptRatio = ExceptRatio*double(k+j*4*N*N*CorrLength)/double(k+j*4*N*N*CorrLength+1) + SpinIce.Spinflip(0)/double(k+j*4*N*N*CorrLength+1);
					if(k < (4*N*N-Ndef)*CorrLength)
						SIarray[ti].Spinflip(1);
				}
				
				//pt swap
				for(int p = 0; p < M-1; p ++) {
					double E1, E2, dt;
					E1 = SIarray[Idx[p]].ShowEnergy();
					E2 = SIarray[Idx[p+1]].ShowEnergy();
					
					dt = (Betaarray[p+1] - Betaarray[p])*(E1 - E2);
					
					double r = ran.randDblExc();
					if(dt < 0 || r < exp(-dt)) {
						
						SIarray[Idx[p]].ResetTemp(1./Betaarray[p+1]);
						SIarray[Idx[p+1]].ResetTemp(1./Betaarray[p]);
						
						int cp = Idx[p];
						Idx[p] = Idx[p+1];
						Idx[p+1] = cp;
						
						
					}
				}
				
			}
			for(int ti = 0; ti < M; ti ++) {
				
				//measurements here
				for(int tp4 = 0; tp4 < 4; tp4++) {
					vertex4[tp4][ti] = vertex4[tp4][ti]*double(j)/double(j+1)+SIarray[Idx[ti]].Showz4Vertex(tp4)/double(j+1);
				}
				for(int tp3 = 0; tp3 < 5; tp3++) {
					vertex3[tp3][ti] = vertex3[tp3][ti]*double(j)/double(j+1)+SIarray[Idx[ti]].Showz3Vertex(tp3)/double(j+1);
				}
				
				EAverage[ti] = EAverage[ti]*double(j)/double(j+1) + SIarray[Idx[ti]].ShowEnergy()/double(j+1);
				EsqAverage[ti] = EsqAverage[ti]*double(j)/double(j+1) + SIarray[Idx[ti]].ShowEnergy()*SIarray[Idx[ti]].ShowEnergy()/double(j+1);
				
			}
			
		}
		
		for(int ti = 0; ti < M; ti ++) {
			int * Dn = SIarray[Idx[ti]].Dn56();
			fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n", \
					1./Betaarray[ti], 0., EsqAverage[ti]-EAverage[ti]*EAverage[ti], vertex4[0][ti], vertex4[1][ti], vertex4[2][ti], vertex4[3][ti], \
					vertex3[0][ti], vertex3[1][ti], vertex3[2][ti], vertex3[3][ti], vertex3[4][ti], Dn[0], Dn[1]);
		}
		
	}
	fclose(output);
	
}
