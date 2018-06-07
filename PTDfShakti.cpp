/*
 *  Shakti.cpp
 *  This code used parallel tempering
 *
 *  Edited on 08/20/2016
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "DfShakti.h"
#include "Functions.h"

int main(int argc, char* argv[]) {

	int N, Ndef, Nthermal, Nmeasure, CorrLength, Nstep, rank, sec;
	int M = 30; //number of templates
	double Ez4[4];
	double Ez3[5];

	if(argc != 3) {
		std::cout << "Wrong usage!\n"; }
	
	sscanf(argv[1], "%d", &sec);
	sscanf(argv[2], "%d", &rank);
	
    double MaxT, MinT;
		
	// parameters form file
	std::ifstream input;
	input.open("param.dat");
	if(input.fail()) {
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
	
    Nstep = CorrLength;
    
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
	
	std::vector<DfShakti> SIarray1;
	
	std::vector<double> Betaarray1;
	
	std::vector<int> Idx1;
    
    std::vector<DfShakti> SIarray2;
    
    std::vector<double> Betaarray2;
    
    std::vector<int> Idx2;
    
	double T = MaxT, DeltaT= (MaxT-MinT)/M;
    
    
    for(int i = 0; i < M; i ++){
        Betaarray1.push_back(1./T);
        Idx1.push_back(i);
        
        Betaarray2.push_back(1./T);
        Idx2.push_back(i);
        
        SIarray1.push_back(DfShakti(N, Ndef, rank, T, Ez4[0], Ez4[1], Ez4[2], Ez4[3], Ez3[0], Ez3[1], Ez3[2], Ez3[3], Ez3[4]));
        SIarray2.push_back(DfShakti(N, Ndef, rank, T, Ez4[0], Ez4[1], Ez4[2], Ez4[3], Ez3[0], Ez3[1], Ez3[2], Ez3[3], Ez3[4]));
    
        T -= DeltaT;
    }

    double vertex3[5][M];			//averge vertex count
    double vertex4[4][M];
    
    double EAverage[M], EsqAverage[M], ExceptRatio[M], Ej[M];	//energy, specific heat, and exceptance rate
    double Ism[M], Ismsq[M], OvLap0[M], OvLap1[M], OvLap2[M];
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
        Ism[ti] = 0.;
        Ismsq[ti] = 0.;
        OvLap0[ti] = 0.;
        OvLap1[ti] = 0.;
        OvLap2[ti] = 0.;
        
    }
    
    int checkpoint[6] = {5000, 10000, 20000, 30000, 40000, 50000};  //check the convergence at these steps.
    int cpnt = 0;
    double sus = 0., sus1 = 0.;
    for(int i = 0; i < Nthermal; i++) {
        for(int ti = 0; ti < M; ti ++) {
            if(1./Betaarray1[ti] < 0.605)
                Nstep = CorrLength * 2;
            
            if(1./Betaarray1[ti] < 0.305)
                Nstep = CorrLength * 3;
            
            for(int k = 0; k < 4*N*N*Nstep; k++) {
                SIarray1[Idx1[ti]].Spinflip(0);
                SIarray1[Idx1[ti]].Spinflip(0);
                SIarray1[Idx1[ti]].Spinflip(0);
                SIarray1[Idx1[ti]].Spinflip(0);
                SIarray1[Idx1[ti]].Doubleflip();
                SIarray1[Idx1[ti]].Doubleflip();
                SIarray1[Idx1[ti]].Z3flip();
                //SIarray1[ti].Z4flip();
                if(k < (4*N*N-Ndef)*Nstep)
                    SIarray1[ti].Spinflip(1);
                
                SIarray2[Idx2[ti]].Spinflip(0);
                SIarray2[Idx2[ti]].Spinflip(0);
                SIarray2[Idx2[ti]].Spinflip(0);
                SIarray2[Idx2[ti]].Spinflip(0);
                SIarray2[Idx2[ti]].Doubleflip();
                SIarray2[Idx2[ti]].Doubleflip();
                SIarray2[Idx2[ti]].Z3flip();
                //SIarray2[ti].Z4flip();
                if(k < (4*N*N-Ndef)*Nstep)
                    SIarray2[ti].Spinflip(1);
                
            }
        }
        
        //pt swap
        //Betaarray always in order. The i th Beta correspond in Idx[i]th SIarray.
        for(int p = 0; p < M-1; p ++) {
            double E1, E2, dt;
            E1 = SIarray1[Idx1[p]].ShowEnergy();
            E2 = SIarray1[Idx1[p+1]].ShowEnergy();
            
            dt = (Betaarray1[p+1] - Betaarray1[p])*(E1 - E2);
            
            double r = ran.randDblExc();
            if(dt < 0 || r < exp(-dt)) {
                
                SIarray1[Idx1[p]].ResetTemp(1./Betaarray1[p+1]);
                SIarray1[Idx1[p+1]].ResetTemp(1./Betaarray1[p]);
                
                int cp = Idx1[p];
                Idx1[p] = Idx1[p+1];
                Idx1[p+1] = cp;
                
                
            }
        }
        
        for(int p = 0; p < M-1; p ++) {
            double E1, E2, dt;
            E1 = SIarray2[Idx2[p]].ShowEnergy();
            E2 = SIarray2[Idx2[p+1]].ShowEnergy();
            
            dt = (Betaarray2[p+1] - Betaarray2[p])*(E1 - E2);
            
            double r = ran.randDblExc();
            if(dt < 0 || r < exp(-dt)) {
                
                SIarray2[Idx2[p]].ResetTemp(1./Betaarray2[p+1]);
                SIarray2[Idx2[p+1]].ResetTemp(1./Betaarray2[p]);
                
                int cp = Idx2[p];
                Idx2[p] = Idx2[p+1];
                Idx2[p+1] = cp;
            }
        }
        
        
        //record convergence behavior at choosen steps
        if(i>=checkpoint[cpnt] && i<checkpoint[cpnt]+2000) {
            sus += SIarray1[Idx1[M-1]].OverLap(SIarray2[Idx2[M-1]], 0)/2000.;
            sus1 += SIarray1[Idx1[M-1]].OverLap(SIarray2[Idx2[M-1]], 1)/2000.;
        }
        
        if(i==checkpoint[cpnt]+2000) {
            printf("%d\t%lf\t%lf\n", checkpoint[cpnt], sus, sus1);
            cpnt ++;
            sus = 0.;
            sus1 = 0.;
        }
    }//end thermalize
    
    
    for(int j = 0; j < Nmeasure; j++) {
        for(int ti = 0; ti < M; ti ++) {
            if(1./Betaarray1[ti] < 0.6)
                Nstep = CorrLength * 2;
            
            if(1./Betaarray1[ti] < 0.3)
                Nstep = CorrLength * 3;
            for(int k =0 ; k < 4*N*N*Nstep; k++) {
                SIarray1[ti].Spinflip(0);
                SIarray1[ti].Spinflip(0);
                SIarray1[ti].Spinflip(0);
                SIarray1[ti].Spinflip(0);
                
                SIarray1[ti].Doubleflip();
                SIarray1[ti].Z3flip();
                
                //ExceptRatio = ExceptRatio*double(k+j*4*N*N*CorrLength)/double(k+j*4*N*N*CorrLength+1) + SpinIce.Spinflip(0)/double(k+j*4*N*N*CorrLength+1);
                if(k < (4*N*N-Ndef)*Nstep)
                    SIarray1[ti].Spinflip(1);
                
                SIarray2[ti].Spinflip(0);
                SIarray2[ti].Spinflip(0);
                SIarray2[ti].Spinflip(0);
                SIarray2[ti].Spinflip(0);
                
                SIarray2[ti].Doubleflip();
                SIarray2[ti].Z3flip();
                
                //ExceptRatio = ExceptRatio*double(k+j*4*N*N*CorrLength)/double(k+j*4*N*N*CorrLength+1) + SpinIce.Spinflip(0)/double(k+j*4*N*N*CorrLength+1);
                if(k < (4*N*N-Ndef)*Nstep)
                    SIarray2[ti].Spinflip(1);
                
            }
        }
        
        //pt swap
        for(int p = 0; p < M-1; p ++) {
            double E1, E2, dt;
            E1 = SIarray1[Idx1[p]].ShowEnergy();
            E2 = SIarray1[Idx1[p+1]].ShowEnergy();
            
            dt = (Betaarray1[p+1] - Betaarray1[p])*(E1 - E2);
            
            double r = ran.randDblExc();
            if(dt < 0 || r < exp(-dt)) {
                
                SIarray1[Idx1[p]].ResetTemp(1./Betaarray1[p+1]);
                SIarray1[Idx1[p+1]].ResetTemp(1./Betaarray1[p]);
                
                int cp = Idx1[p];
                Idx1[p] = Idx1[p+1];
                Idx1[p+1] = cp;
                
                
            }
        }
        
        for(int p = 0; p < M-1; p ++) {
            double E1, E2, dt;
            E1 = SIarray2[Idx2[p]].ShowEnergy();
            E2 = SIarray2[Idx2[p+1]].ShowEnergy();
            
            dt = (Betaarray2[p+1] - Betaarray2[p])*(E1 - E2);
            
            double r = ran.randDblExc();
            if(dt < 0 || r < exp(-dt)) {
                
                SIarray2[Idx2[p]].ResetTemp(1./Betaarray2[p+1]);
                SIarray2[Idx2[p+1]].ResetTemp(1./Betaarray2[p]);
                
                int cp = Idx2[p];
                Idx2[p] = Idx2[p+1];
                Idx2[p+1] = cp;
                
                ExceptRatio[p] = ExceptRatio[p]*double(j)/double(j+1) + 1./double(j+1);
                
            }
            
            else
                ExceptRatio[p] = ExceptRatio[p]*double(j)/double(j+1);
        }
        
        for(int ti = 0; ti < M; ti ++) {
            
            //measurements here
            for(int tp4 = 0; tp4 < 4; tp4++) {
                vertex4[tp4][ti] = vertex4[tp4][ti]*double(j)/double(j+1)+(SIarray1[Idx1[ti]].Showz4Vertex(tp4)+SIarray2[Idx2[ti]].Showz4Vertex(tp4))/double(j+1)/2.;
            }
            for(int tp3 = 0; tp3 < 5; tp3++) {
                vertex3[tp3][ti] = vertex3[tp3][ti]*double(j)/double(j+1)+(SIarray1[Idx1[ti]].Showz3Vertex(tp3)+SIarray2[Idx2[ti]].Showz3Vertex(tp3))/double(j+1)/2.;
            }
            
            EAverage[ti] = EAverage[ti]*double(j)/double(j+1) + (SIarray1[Idx1[ti]].ShowEnergy()+SIarray2[Idx2[ti]].ShowEnergy())/2./double(j+1);
            EsqAverage[ti] = EsqAverage[ti]*double(j)/double(j+1) + \
            (SIarray1[Idx1[ti]].ShowEnergy()*SIarray1[Idx1[ti]].ShowEnergy()+SIarray2[Idx2[ti]].ShowEnergy()*SIarray2[Idx2[ti]].ShowEnergy())/2./double(j+1);
            
            int ism1 = SIarray1[Idx1[ti]].IsingMv2(), ism2 = SIarray2[Idx2[ti]].IsingMv2();
            
            Ismsq[ti] = Ismsq[ti]*double(j)/double(j+1) + (ism1*ism1+ism2*ism2)/2./double(j+1)/(64*N*N*N*N);
            Ism[ti] = Ism[ti]*double(j)/double(j+1) + (abs(ism1)+abs(ism2))/2./double(j+1)/double(8*N*N);
            
            //total # of spins 20*N*N - NDef;
            double q0 =  SIarray1[Idx1[ti]].OverLap(SIarray2[Idx2[ti]], 0);
            double q1 =  SIarray1[Idx1[ti]].OverLap(SIarray2[Idx2[ti]], 1);
            double q2 =  SIarray1[Idx1[ti]].OverLap(SIarray2[Idx2[ti]], 2);
            OvLap0[ti] = OvLap0[ti]*double(j)/double(j+1) + q0/double(20*N*N-Ndef)/double(j+1);
            OvLap1[ti] = OvLap1[ti]*double(j)/double(j+1) + q1/double(20*N*N-Ndef)/double(j+1);
            OvLap2[ti] = OvLap2[ti]*double(j)/double(j+1) + q2/double(20*N*N-Ndef)/double(j+1);
            
        }
        
    }// end for j
    ExceptRatio[M-1] = ExceptRatio[M-2];
    
    //The data structure: T, Exceptance rate(not very accurate), Var of energy, # of 4 types of z4, # of 3 type of z3, Ising order |m^2|, Ising order |m|, replica overlap
    for(int ti = 0; ti < M; ti ++) {
        fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", \
                1./Betaarray1[ti], ExceptRatio[ti], EsqAverage[ti]-EAverage[ti]*EAverage[ti], vertex4[0][ti], vertex4[1][ti], vertex4[2][ti], vertex4[3][ti], \
                vertex3[0][ti], vertex3[1][ti], vertex3[2][ti], vertex3[3][ti], vertex3[4][ti], Ismsq[ti], Ism[ti], OvLap0[ti], OvLap1[ti], OvLap2[ti]);
    }
    
	fclose(output);
	
}
