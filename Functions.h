/*
 *  Functions.h Implement class functions.
 *  
 *
 *  Created by Yifei Shi on 11/1/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "DfShakti.h"

const double PI = 3.14159265;

inline void DfShakti::SetCurrentConfig_(const Link_* lk, int crt) {
	if(lk->lr == 0) 
		z4CurrentConfig_[4*(lk->lx+2*N_*lk->ly)+lk->ldir] = crt;
	
	else if(lk->lr == 1)
		z3r1CurrentConfig_[3*(lk->lx+2*N_*lk->ly)+lk->ldir] = crt;
	
	else if(lk->lr == 2)
		z3r2CurrentConfig_[3*(lk->lx+2*N_*lk->ly)+lk->ldir] = crt;
	
	else 
		std::cout << "SetCurrent Error!\n" ;
}

DfShakti::DfShakti(int N, int Ndef, int node, double T, double Ez4a, double Ez4b, double Ez4c, double Ez4d, double Ez3A, double Ez3B, double Ez3C, double Ez3D, double Ez3E)
:	Thermalized(0)
,	N_(N)
,	tDf_(Ndef)
,	Beta_(1./T)
,	Node_(node)
{
	Ez4_[0]=Ez4a;
	Ez4_[1]=Ez4b;
	Ez4_[2]=Ez4c;
	Ez4_[3]=Ez4d;
	Ez3_[0]=Ez3A;
	Ez3_[1]=Ez3B;
	Ez3_[2]=Ez3C;
	Ez3_[3]=Ez3D;
	Ez3_[4]=Ez3E;
	
	Link_ nlk;  //neighbor link

	z3r1CurrentConfig_.resize(4*N_*N_*3);
	z3r2CurrentConfig_.resize(4*N_*N_*3);
	z4CurrentConfig_.resize(4*N_*N_*4);
	
	z3r1VertexType_.resize(4*N_*N_);
	z3r2VertexType_.resize(4*N_*N_);
	z4VertexType_.resize(4*N_*N_);
	
	//init Neighbor structure, when there are no defects
	for(int y = 0; y < 2*N_; y++){
		for(int x=0; x < 2*N_; x++) {
			
			//z=4 sites
			nlk.lr = 1;
			nlk.lx = x;
			nlk.ly = y;
			nlk.ldir = 1;
			r0NeighborStruct_.push_back(nlk);

			nlk.lr = 2;
			nlk.lx = x;
			nlk.ly = y;
			nlk.ldir = 1;
			r0NeighborStruct_.push_back(nlk);
			
			nlk.lr = 1;
			nlk.lx = x==0? 2*N_-1 : x-1;
			nlk.ly = y;
			nlk.ldir = 0;
			r0NeighborStruct_.push_back(nlk);
			
			nlk.lr = 2;
			nlk.lx = x;
			nlk.ly = y==0? 2*N_-1 : y-1;
			nlk.ldir = 0;
			r0NeighborStruct_.push_back(nlk);
			
			
			//z=3 sites, more complicated than z=4 cases;
			
			nlk.lr = 0;
			nlk.lx = x==2*N_-1? 0 : x+1;
			nlk.ly = y;
			nlk.ldir = 2;
			r1NeighborStruct_.push_back(nlk);
			
			nlk.lr = 0;
			nlk.lx = x;
			nlk.ly = y;
			nlk.ldir = 0;
			r1NeighborStruct_.push_back(nlk);
			
			nlk.lr = 0;
			nlk.lx = x;
			nlk.ly = y==2*N_-1? 0 : y+1;
			nlk.ldir = 3;
			r2NeighborStruct_.push_back(nlk);
			
			nlk.lr = 0;
			nlk.lx = x;
			nlk.ly = y;
			nlk.ldir = 1;
			r2NeighborStruct_.push_back(nlk);			
			
			
			if((x+y)%2==0) {
				
				
				nlk.lr = 1;
				nlk.lx = x;
				nlk.ly = y==0? 2*N_-1 : y-1;
				nlk.ldir = 2;
				r1NeighborStruct_.push_back(nlk);
				
								
				nlk.lr = 2;
				nlk.lx = x==2*N_-1? 0 : x+1;
				nlk.ly = y;
				nlk.ldir = 2;
				r2NeighborStruct_.push_back(nlk);
				
			}
			
			else {
								
				nlk.lr = 1;
				nlk.lx = x;
				nlk.ly = y==2*N_-1? 0 : y+1;
				nlk.ldir = 2;
				r1NeighborStruct_.push_back(nlk);
								
				nlk.lr = 2;
				nlk.lx = x==0? 2*N_-1 : x-1;
				nlk.ly = y;
				nlk.ldir = 2;
				r2NeighborStruct_.push_back(nlk);
			}
			
			
		}
	}//end init neighbor no defects
	
	//add defects to the system
	
	Defects_.resize(tDf_);

	GenerateDefects();
			
	Listofr1_.resize(4*N_*N_);
	Listofr2_.resize(4*N_*N_);
	
	for(int i=0; i<4*N_*N_; i++) {
		Listofr1_[i] = i;
		Listofr2_[i] = i;
	}

	for(int d = 0; d < tDf_; d ++) {
		int dfx = Defects_[d]%(2*N_);
		int dfy = Defects_[d]/(2*N_);
		
		Link_ dnlk;
		
		//z=4 vertices
		int ndfx = dfx == 2*N_-1? 0 : dfx + 1;
		int ndfy = dfy == 2*N_-1? 0 : dfy + 1;
		
		if( (dfx + dfy)%2 == 0) {
										
			//z=3 vertices, they are taken out, by setting lr == -1;
			dnlk.lr = -1;
			dnlk.lx = -1;
			dnlk.ly = -1;
			dnlk.ldir = -1;
			
			Listofr2_[Defects_[d]] = -1;
			Listofr2_[ndfx + 2*N_*dfy] = -1;
			
			r2NeighborStruct_[3*Defects_[d]+2] = nlk;
			r2NeighborStruct_[3*(ndfx + 2*N_*dfy)+2] = nlk;
			z3r2CurrentConfig_[3*Defects_[d]+2] = 0;
			z3r2CurrentConfig_[3*(ndfx + 2*N_*dfy)+2] = 0;
			
		}  //end if ()%2 == 0
		
		else {
						
			//z=3 vertices, they are taken out, by setting lr == -1;
			dnlk.lr = -1;
			dnlk.lx = -1;
			dnlk.ly = -1;
			dnlk.ldir = -1;
			
			Listofr1_[Defects_[d]] = -1;
			Listofr1_[dfx + 2*N_*ndfy] = -1;

			r1NeighborStruct_[3*Defects_[d]+2] = nlk;
			r1NeighborStruct_[3*(dfx + 2*N_*ndfy)+2] = nlk;
			z3r1CurrentConfig_[3*Defects_[d]+2] = 0;
			z3r1CurrentConfig_[3*(dfx + 2*N_*2)+2] = 0;
			
		}	//end else ()%2 == 1
	}// end add defects

	//Pick out all the existing r1 and r2 sites
	int ndfr1 = 0, ndfr2 = 0;
	for(int i = 0; i < 4*N_*N_; i ++) {
		if(Listofr1_[i-ndfr1] == -1) {
			Listofr1_.erase(Listofr1_.begin() + i - ndfr1);
			ndfr1 ++;
		}
		
		if(Listofr2_[i-ndfr2] == -1) {
			Listofr2_.erase(Listofr2_.begin() + i - ndfr2);
			ndfr2 ++;
		}
	}

	if(ndfr1 + ndfr2 != 2*tDf_)
		std::cout << "Error! # of defects mismatch!\n";
	
	tr1Ild_ = 4*N_*N_ - ndfr1;
	tr2Ild_ = 4*N_*N_ - ndfr2;
	
	if(tr1Ild_ < 0 || tr2Ild_ < 0) 
		std::cout << "Error! # of island is negative!\n";
		
	//random configuration to start. To count the currents only once, only consider links connected to z3 sites, and only for (x+y)%2=0 consider the middle link
		
	int cc, drcn; //randomly generated current direction
		
	//init all the currents of z4 vertices
	for(int x = 0; x < 2*N_; x ++) {
		for(int y = 0; y < 2*N_; y ++) {
			for(drcn = 0; drcn < 4; drcn ++) {
				cc = 2*int(ran.randDblExc()*2)-1;
				z4CurrentConfig_[z4NeighborCount_(x,y,drcn)] = cc;
				
				Link_ nlk = r0NeighborStruct_[z4NeighborCount_(x,y,drcn)];
				
				if(nlk.lr == 1) 
					z3r1CurrentConfig_[z3NeighborCount_(nlk.lx,nlk.ly,nlk.ldir)] = -cc;
				
				else if(nlk.lr == 2)
					z3r2CurrentConfig_[z3NeighborCount_(nlk.lx,nlk.ly,nlk.ldir)] = -cc;
				
			}
		}
	}
		
	//init all island currents
	
	for(int dd = 0; dd < tr1Ild_; dd ++) {
		int coord1 = Listofr1_[dd];
		
		if((coord1%(2*N_)+coord1/(2*N_))%2 == 0) {
			cc = 2*int(ran.randDblExc()*2)-1;
			z3r1CurrentConfig_[coord1*3+2] = cc;
			Link_ nlk = r1NeighborStruct_[coord1*3+2];
			
			z3r1CurrentConfig_[z3NeighborCount_(nlk.lx,nlk.ly,nlk.ldir)] = -cc;
		}
	}
	
	for(int dd = 0; dd < tr2Ild_; dd ++) {
		int coord2 = Listofr2_[dd];
		
		if((coord2%(2*N_)+coord2/(2*N_))%2 == 0) {
			cc = 2*int(ran.randDblExc()*2)-1;
			z3r2CurrentConfig_[coord2*3+2] = cc;
			Link_ nlk = r2NeighborStruct_[coord2*3+2];
			
			z3r2CurrentConfig_[z3NeighborCount_(nlk.lx,nlk.ly,nlk.ldir)] = -cc;			
		}
	}	
		
	int c0, c1, c2, c3, ctotal;	//currents of z4 vertex, and total
	int cc0, cc1, cc2;	//same for z3 but no total
	int p3, p4;
	int ta = 0, tb = 0, tc = 0, td = 0, tA = 0, tB = 0, tC = 0, tD = 0, tE = 0;
	
	Energy_ = 0;
	//vertex type and energy init. Vertex types are 'a', 'b', 'c', 'd' for z=4 and 'A', 'B', 'C' for z=3;
	for(int x = 0; x < 2*N_; x++) {
		for(int y = 0; y < 2*N_; y++) {
			p4 = x + 2*N_*y;
			c0 = z4CurrentConfig_[4*p4];
			c1 = z4CurrentConfig_[4*p4+1];
			c2 = z4CurrentConfig_[4*p4+2];
			c3 = z4CurrentConfig_[4*p4+3];
			
			ctotal = c0+c1+c2+c3;
			
			if(ctotal==4 || ctotal==-4) {
				td++;
				z4VertexType_[p4] = 3;
				Energy_ += Ez4_[3]; }
			
			else if(ctotal==2 || ctotal==-2) {
				tc++;
				z4VertexType_[p4] = 2;
				Energy_ += Ez4_[2]; }
			
			else if(ctotal==0 && c0 == -c2) {
				tb++;
				z4VertexType_[p4] = 1;
				Energy_ += Ez4_[1]; }
			
			else if(ctotal==0 && c0 == c2) {
				ta++;
				z4VertexType_[p4] = 0;
				Energy_ += Ez4_[0]; }
			
			else 
				std::cout << "z4 current config error!!!!\n";
			
			p3 = x + 2*N_*y;
			cc0 = z3r1CurrentConfig_[3*p3];
			cc1 = z3r1CurrentConfig_[3*p3+1];
			cc2 = z3r1CurrentConfig_[3*p3+2];
			
			if(cc2 != 0) {
				if(cc0 != cc1) {
					tB ++;
					z3r1VertexType_[p3] = 1;
					Energy_ += Ez3_[1]; }
			
				else if(cc0==cc1 && cc0==cc2) {
					tC ++;
					z3r1VertexType_[p3] = 2;
					Energy_ += Ez3_[2]; }
			
				else if(cc0==cc1 && cc0!=cc2) {
					tA++;
					z3r1VertexType_[p3] = 0;
					Energy_ += Ez3_[0]; }
				else 
					std::cout << "z3 current config error!!!!\n";
				}
			
			else {
				if(cc0 != cc1) {
					tD ++;
					z3r1VertexType_[p3] = 3;
					Energy_ += Ez3_[3];
				}
				
				else {
					tE ++;
					z3r1VertexType_[p3] = 4;
					Energy_ += Ez3_[4];
				}
				
			}
			
			cc0 = z3r2CurrentConfig_[3*p3];
			cc1 = z3r2CurrentConfig_[3*p3+1];
			cc2 = z3r2CurrentConfig_[3*p3+2];
			
			if(cc2 != 0) {
				if(cc0 != cc1) {
					tB ++;
					z3r2VertexType_[p3] = 1;
					Energy_ += Ez3_[1]; }
			
				else if(cc0==cc1 && cc0==cc2) {
					tC ++;
					z3r2VertexType_[p3] = 2;
					Energy_ += Ez3_[2]; }
			
				else if(cc0==cc1 && cc0!=cc2) {
					tA++;
					z3r2VertexType_[p3] = 0;
					Energy_ += Ez3_[0]; }
				else 
					std::cout << "z3 current config error!!!!\n";
				}
			
			else {
				if(cc0 != cc1) {
					tD ++;
					z3r2VertexType_[p3] = 3;
					Energy_ += Ez3_[3];
				}
				
				else {
					tE ++;
					z3r2VertexType_[p3] = 4;
					Energy_ += Ez3_[4];
				}
				
			}
			
		}
	}  //end vertex type and energy init
		
	z4VertexCount_[0] = ta;
	z4VertexCount_[1] = tb;
	z4VertexCount_[2] = tc;
	z4VertexCount_[3] = td;
	
	z3VertexCount_[0] = tA;
	z3VertexCount_[1] = tB;
	z3VertexCount_[2] = tC;
	z3VertexCount_[3] = tD;
	z3VertexCount_[4] = tE;
	
	
	Thermalized = false;
	
}


inline int DfShakti::Decidez4VertexType_(int c0, int c1, int c2, int c3) {
	int ctotal = c0+c1+c2+c3;
	
	if(ctotal==4 || ctotal==-4) 
		return 3;
	
	else if(ctotal==2 || ctotal==-2) 
		return 2;
	
	else if(ctotal==0 && c0 == -c2) 
		return 1;
	
	else if(ctotal==0 && c0 == c2) 
		return 0;
	else {
		std::cout << "z4 current config error!!!!\n";
		return -1; }	
	
}

inline int DfShakti::Decidez3VertexType_(int cc0, int cc1, int cc2) {
	if(cc2 != 0) {
		if(cc0 != cc1) 
			return 1;
	
		else if(cc0==cc1 && cc0==cc2) 
			return 2;
	
		else if(cc0==cc1 && cc0!=cc2) 
			return 0;
		else {
			std::cout << "z3 current config error!!!!\n";
			return -1; }
	}
	
	else {
		if(cc0 != cc1)
			return 3;
		else return 4;
	}
	
}


//flip a spin. First specify the vertex, then specify which link of the vertex.
int DfShakti::Spinflip(int onwhichsite) {
	//whether the site is z4 or not
	if(onwhichsite==0) {
		int randomsite = ran.randDblExc()*4*N_*N_;
		int drcn = ran.randDblExc()*4;
		
		Link_ nlk = r0NeighborStruct_[randomsite*4+drcn];
		int nbsite = nlk.lx + nlk.ly * 2*N_;
				

		int c4[4];	//curents involved;
		int c3[3];
		
		c4[0] = z4CurrentConfig_[4*randomsite];
		c4[1] = z4CurrentConfig_[4*randomsite+1];
		c4[2] = z4CurrentConfig_[4*randomsite+2];
		c4[3] = z4CurrentConfig_[4*randomsite+3];
		
		if(nlk.lr == 1) {
			c3[0] = z3r1CurrentConfig_[3*nbsite];
			c3[1] = z3r1CurrentConfig_[3*nbsite+1];
			c3[2] = z3r1CurrentConfig_[3*nbsite+2];
		}
		
		else {
			c3[0] = z3r2CurrentConfig_[3*nbsite];
			c3[1] = z3r2CurrentConfig_[3*nbsite+1];
			c3[2] = z3r2CurrentConfig_[3*nbsite+2];
		}
		
		double Ebf, Eaf;		//energy before and after flip
		
		int vt3bf, vt4bf, vt3af, vt4af;		//vertex types before and after
		
		vt4bf = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
		vt3bf = Decidez3VertexType_(c3[0], c3[1], c3[2]);
		
		Ebf = Ez4_[vt4bf] + Ez3_[vt3bf];
		
		c4[drcn] *= -1;
		c3[nlk.ldir] *= -1;
		
		vt3af = Decidez3VertexType_(c3[0], c3[1], c3[2]);
		vt4af = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
		
		Eaf = Ez4_[vt4af] + Ez3_[vt3af];
		
		double prob = ran.randDblExc();
		
		if(Eaf<Ebf || prob < exp(-Beta_*(Eaf-Ebf))) {
			
			z4CurrentConfig_[4*randomsite+drcn] *= -1;
			z4VertexType_[randomsite] = vt4af;				
				
			z4VertexCount_[vt4bf] --;
			z4VertexCount_[vt4af] ++;
				
			if(nlk.lr == 1) {
				z3r1CurrentConfig_[3*nbsite+nlk.ldir] *= -1;
				z3r1VertexType_[nlk.lx+2*N_*nlk.ly] = vt3af; 
			}
			
			else {
				z3r2CurrentConfig_[3*nbsite+nlk.ldir] *= -1;
				z3r2VertexType_[nlk.lx+2*N_*nlk.ly] = vt3af; 
			}
			
			z3VertexCount_[vt3bf] --;
			z3VertexCount_[vt3af] ++;
			
			Energy_ += Eaf - Ebf;
						
			return 1;
			
		}
		else return 0;
	} //onwhichsite == 0
	
	else {
		int randidx = ran.randDblExc()*(tr1Ild_+tr2Ild_);
				
		if(randidx < tr1Ild_) {
		
		int randomsite = Listofr1_[randidx];
		int drcn = 2;
		
		Link_ nlk = r1NeighborStruct_[3*randomsite+drcn];
		
		int nbsite = nlk.lx + 2*N_*nlk.ly;
		
			//whether a middle(or island) link
			if(nlk.lr == 1) {
				int c3[3], nc3[3];
				
				c3[0] = z3r1CurrentConfig_[3*randomsite];
				c3[1] = z3r1CurrentConfig_[3*randomsite+1];
				c3[2] = z3r1CurrentConfig_[3*randomsite+2];
				
				nc3[0] = z3r1CurrentConfig_[3*nbsite];
				nc3[1] = z3r1CurrentConfig_[3*nbsite+1];
				nc3[2] = z3r1CurrentConfig_[3*nbsite+2];
				
				if(c3[2] == 0 || nc3[2] ==0 ) std::cout << "current 0!\n";
				
				double Ebf, Eaf;
				
				int vt3bf, vtn3bf, vt3af, vtn3af;
				
				vt3bf = Decidez3VertexType_(c3[0], c3[1], c3[2]);
				vtn3bf = Decidez3VertexType_(nc3[0], nc3[1], nc3[2]);
				
				Ebf = Ez3_[vt3bf] + Ez3_[vtn3bf];
				
				c3[drcn] *= -1;
				nc3[nlk.ldir] *= -1;
				
				vt3af = Decidez3VertexType_(c3[0], c3[1], c3[2]);
				vtn3af = Decidez3VertexType_(nc3[0], nc3[1], nc3[2]);
				
				Eaf = Ez3_[vt3af] + Ez3_[vtn3af];
				
				double prob = ran.randDblExc();
				
				if(Eaf<Ebf || prob < exp(-Beta_*(Eaf-Ebf))) {
					z3r1CurrentConfig_[3*randomsite+drcn] *= -1;
					z3r1CurrentConfig_[3*nbsite+nlk.ldir] *= -1;
					
					z3r1VertexType_[nbsite] = vtn3af;
					z3r1VertexType_[randomsite] = vt3af;
					
					z3VertexCount_[vt3bf] --;
					z3VertexCount_[vt3af] ++;
					
					z3VertexCount_[vtn3bf] --;
					z3VertexCount_[vtn3af] ++;
					
					Energy_ += Eaf - Ebf;
					
                    if(std::abs(Eaf-Ebf) > 0.001)
						return 1;
					else return 0;
				}
				else return 0;
			} // lr == 1
			
			else {
				std::cout << "r1 neighbor type wrong!\n";
				return -1; }		
		} //onwhichsitee == 1 (or randidx < tr1Ild_)
		
		else  {
			randidx -= tr1Ild_;
			int randomsite = Listofr2_[randidx];
			int drcn = 2;
			
			Link_ nlk = r2NeighborStruct_[3*randomsite+drcn];
			
			int nbsite = nlk.lx + 2*N_*nlk.ly;
			
			//whether a middle link
			if(nlk.lr == 2) {
				int c3[3], nc3[3];
				
				c3[0] = z3r2CurrentConfig_[3*randomsite];
				c3[1] = z3r2CurrentConfig_[3*randomsite+1];
				c3[2] = z3r2CurrentConfig_[3*randomsite+2];
				
				nc3[0] = z3r2CurrentConfig_[3*nbsite];
				nc3[1] = z3r2CurrentConfig_[3*nbsite+1];
				nc3[2] = z3r2CurrentConfig_[3*nbsite+2];
				
				if(c3[2] == 0 || nc3[2] ==0 ) std::cout << "current 0!\n";
				
				double Ebf, Eaf;
				
				int vt3bf, vtn3bf, vt3af, vtn3af;
				
				vt3bf = Decidez3VertexType_(c3[0], c3[1], c3[2]);
				vtn3bf = Decidez3VertexType_(nc3[0], nc3[1], nc3[2]);
				
				Ebf = Ez3_[vt3bf] + Ez3_[vtn3bf];
				
				c3[drcn] *= -1;
				nc3[nlk.ldir] *= -1;
				
				vt3af = Decidez3VertexType_(c3[0], c3[1], c3[2]);
				vtn3af = Decidez3VertexType_(nc3[0], nc3[1], nc3[2]);
				
				Eaf = Ez3_[vt3af] + Ez3_[vtn3af];
				
				double prob = ran.randDblExc();
				
				if(Eaf<Ebf || prob < exp(-Beta_*(Eaf-Ebf))) {
					z3r2CurrentConfig_[3*randomsite+drcn] *= -1;
					z3r2CurrentConfig_[3*nbsite+nlk.ldir] *= -1;
					
					z3r2VertexType_[nbsite] = vtn3af;
					z3r2VertexType_[randomsite] = vt3af;
					
					z3VertexCount_[vt3bf] --;
					z3VertexCount_[vt3af] ++;
					
					z3VertexCount_[vtn3bf] --;
					z3VertexCount_[vtn3af] ++;
					
					Energy_ += Eaf - Ebf;
					
                    if(std::abs(Eaf-Ebf) > 0.001)
						return 1;
					else return 0;
				}
				else return 0;
			}  //lr == 2 			
			
			else {
				std::cout << "r2 neighbor wrong!\n";
				return -1;				
			} 
			
		} //onwhichsitee == 2 (randidx >= tr1ild_)
	}
	
}

int DfShakti::Doubleflip() {
	int randomsite = ran.randDblExc()*4*N_*N_;
	int drcn1, drcn2;
	
	if(z4VertexType_[randomsite] == 0 || z4VertexType_[randomsite] == 1) {
		int c4[4];
		for(int i = 0; i < 4; i ++)
			c4[i] = z4CurrentConfig_[randomsite*4+i];	
		
		if(z4VertexType_[randomsite] == 0) {		
			int quadrant = ran.randDblExc()*4;
			
			if(quadrant == 0) { drcn1 = 3; drcn2 = 0;}
			
			else if(quadrant == 1) { drcn1 = 0; drcn2 = 1;}
			
			else if(quadrant == 2) { drcn1 = 1; drcn2 = 2;}
			
			else { drcn1 = 2; drcn2 = 3;}
		}
		
		else {
			int quadrant = ran.randDblExc()*2;
			
			if(c4[0] == c4[1]) {
				if(quadrant == 0) { drcn1 = 3; drcn2 = 0; }
				else { drcn1 = 1; drcn2 = 2; }
			}
			
			else {
				if(quadrant == 0) { drcn1 = 0; drcn2 = 1; }
				else { drcn1 = 2; drcn2 = 3; }
			}
		}
		
		Link_ nlk1 = r0NeighborStruct_[4*randomsite+drcn1];
		Link_ nlk2 = r0NeighborStruct_[4*randomsite+drcn2];
		int c3n1[3], c3n2[3];
		
		int nbsite1 = nlk1.lx + nlk1.ly * 2*N_;
		int nbsite2 = nlk2.lx + nlk2.ly * 2*N_;
		
		if(nlk1.lr == 1) {
			c3n1[0] = z3r1CurrentConfig_[3*nbsite1];
			c3n1[1] = z3r1CurrentConfig_[3*nbsite1+1];
			c3n1[2] = z3r1CurrentConfig_[3*nbsite1+2];
		}
		
		else {
			c3n1[0] = z3r2CurrentConfig_[3*nbsite1];
			c3n1[1] = z3r2CurrentConfig_[3*nbsite1+1];
			c3n1[2] = z3r2CurrentConfig_[3*nbsite1+2];
		}
		
		if(nlk2.lr == 1) {
			c3n2[0] = z3r1CurrentConfig_[3*nbsite2];
			c3n2[1] = z3r1CurrentConfig_[3*nbsite2+1];
			c3n2[2] = z3r1CurrentConfig_[3*nbsite2+2];
		}
		
		else {
			c3n2[0] = z3r2CurrentConfig_[3*nbsite2];
			c3n2[1] = z3r2CurrentConfig_[3*nbsite2+1];
			c3n2[2] = z3r2CurrentConfig_[3*nbsite2+2];
		}
		
		
		double Ebf, Eaf;		//energy before and after flip
		
		int vt3bf1, vt3bf2, vt4bf, vt3af1, vt3af2, vt4af;		//vertex types before and after
		
		vt4bf = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
		vt3bf1 = Decidez3VertexType_(c3n1[0], c3n1[1], c3n1[2]);
		vt3bf2 = Decidez3VertexType_(c3n2[0], c3n2[1], c3n2[2]);
		
		//int vcp4 = z4VertexType_[randomsite];
		//int vcp31 = nlk1.lr == 1 ? z3r1VertexType_[nbsite1] : z3r2VertexType_[nbsite1];
		//int vcp32 = nlk2.lr == 1 ? z3r1VertexType_[nbsite2] : z3r2VertexType_[nbsite2];
		
		Ebf = Ez4_[vt4bf] + Ez3_[vt3bf1] + Ez3_[vt3bf2];
		
		c4[drcn1] *= -1;
		c4[drcn2] *= -1;
		c3n1[nlk1.ldir] *= -1;
		c3n2[nlk2.ldir] *= -1;
		
		vt3af1 = Decidez3VertexType_(c3n1[0], c3n1[1], c3n1[2]);
		vt3af2 = Decidez3VertexType_(c3n2[0], c3n2[1], c3n2[2]);
		vt4af = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
				
		Eaf = Ez4_[vt4af] + Ez3_[vt3af1] + Ez3_[vt3af2];
		
		double DtE = vt4bf == 0? Eaf-Ebf-log(2.)/Beta_ : Ebf-Eaf-log(2.)/Beta_;
		
		double prob = ran.randDblExc();
		
		
		//trickiy when dealing the probabilities
		if((vt4bf == 0 && (DtE < 0. || prob < exp(-Beta_*DtE))) || (vt4bf == 1 && (DtE > 0. || prob < exp(Beta_*DtE)))) {
			
			z4CurrentConfig_[4*randomsite+drcn1] *= -1;
			z4CurrentConfig_[4*randomsite+drcn2] *= -1;
			
			z4VertexType_[randomsite] = vt4af;				
			
			z4VertexCount_[vt4bf] --;
			z4VertexCount_[vt4af] ++;
			
			if(nlk1.lr == 1) {
				z3r1CurrentConfig_[3*nbsite1+nlk1.ldir] *= -1;
				z3r1VertexType_[nbsite1] = vt3af1; 
			}
			
			else {
				z3r2CurrentConfig_[3*nbsite1+nlk1.ldir] *= -1;
				z3r2VertexType_[nbsite1] = vt3af1; 
			}
			
			z3VertexCount_[vt3bf1] --;
			z3VertexCount_[vt3af1] ++;
			
			if(nlk2.lr == 1) {
				z3r1CurrentConfig_[3*nbsite2+nlk2.ldir] *= -1;
				z3r1VertexType_[nbsite2] = vt3af2; 
			}
			
			else {
				z3r2CurrentConfig_[3*nbsite2+nlk2.ldir] *= -1;
				z3r2VertexType_[nbsite2] = vt3af2; 
			}
			
			z3VertexCount_[vt3bf2] --;
			z3VertexCount_[vt3af2] ++;
						
			Energy_ += Eaf - Ebf;
			
			return 1;
			
		}
		
		else return 0;
		
	}  //end if vtype == 0,1
	
	else {return 0;}
	
}

//A better verion of doubleflip, if the leg is connected to a Z3 tpye 3 vetex (looks like ->-> ), the flip the whole Z3 vertex also.
int DfShakti::DoubleLongflip() {
	int randomsite = ran.randDblExc()*4*N_*N_;
	int drcn1, drcn2;
	
	if(z4VertexType_[randomsite] == 0 || z4VertexType_[randomsite] == 1) {
		int c4[4];
		for(int i = 0; i < 4; i ++)
			c4[i] = z4CurrentConfig_[randomsite*4+i];	
		
		if(z4VertexType_[randomsite] == 0) {		
			int quadrant = ran.randDblExc()*4;
			
			if(quadrant == 0) { drcn1 = 3; drcn2 = 0;}
			
			else if(quadrant == 1) { drcn1 = 0; drcn2 = 1;}
			
			else if(quadrant == 2) { drcn1 = 1; drcn2 = 2;}
			
			else { drcn1 = 2; drcn2 = 3;}
		}
		
		else {
			int quadrant = ran.randDblExc()*2;
			
			if(c4[0] == c4[1]) {
				if(quadrant == 0) { drcn1 = 3; drcn2 = 0; }
				else { drcn1 = 1; drcn2 = 2; }
			}
			
			else {
				if(quadrant == 0) { drcn1 = 0; drcn2 = 1; }
				else { drcn1 = 2; drcn2 = 3; }
			}
		}
		
		Link_ nlk1 = r0NeighborStruct_[4*randomsite+drcn1];
		Link_ nlk2 = r0NeighborStruct_[4*randomsite+drcn2];
		
		int nbsite1 = nlk1.lx + nlk1.ly * 2*N_;
		int nbsite2 = nlk2.lx + nlk2.ly * 2*N_;
		
		if(z3MixVertexType_(nlk1.lr, nbsite1) == 3 && z3MixCurrentConfig_(nlk2.lr, 2, nbsite2) != 0) {
			Link_ nlk3 = z3MixNeighborStruct_( nlk1.lr, 1-nlk1.ldir, nbsite1);
			int nbsite3 = nlk3.lx + 2*N_*nlk3.ly;
			
			int cn4[4];
			cn4[0] = z4CurrentConfig_[4*nbsite3];
			cn4[1] = z4CurrentConfig_[4*nbsite3+1];
			cn4[2] = z4CurrentConfig_[4*nbsite3+2];
			cn4[3] = z4CurrentConfig_[4*nbsite3+3];
			
			int cn3[3];
			cn3[0] = z3MixCurrentConfig_(nlk2.lr, 0, nbsite2);
			cn3[1] = z3MixCurrentConfig_(nlk2.lr, 1, nbsite2);
			cn3[2] = z3MixCurrentConfig_(nlk2.lr, 2, nbsite2);
			
			int vt4bf1, vt4bf2, vt3bf, vt4af1, vt4af2, vt3af;
			
			double Ebf, Eaf;
			
			vt4bf1 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt4bf2 = Decidez4VertexType_(cn4[0], cn4[1], cn4[2], cn4[3]);
			vt3bf = Decidez3VertexType_(cn3[0], cn3[1], cn3[2]);
			
			Ebf = Ez4_[vt4bf1] + Ez4_[vt4bf2] + Ez3_[vt3bf];
			
			c4[drcn1] *= -1;
			c4[drcn2] *= -1;
			cn4[nlk3.ldir] *= -1;
			cn3[nlk2.ldir] *= -1;
			
			vt4af1 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt4af2 = Decidez4VertexType_(cn4[0], cn4[1], cn4[2], cn4[3]);
			vt3af = Decidez3VertexType_(cn3[0], cn3[1], cn3[2]);
			
			Eaf = Ez4_[vt4af1] + Ez4_[vt4af2] + Ez3_[vt3af];
			
			double DtE = vt4bf1 == 0? Eaf-Ebf-log(2.)/Beta_ : Ebf-Eaf-log(2.)/Beta_;
			
			double prob = ran.randDblExc();
						
			//trickiy when dealing the probabilities
			if((vt4bf1 == 0 && (DtE < 0. || prob < exp(-Beta_*DtE))) || (vt4bf1 == 1 && (DtE > 0. || prob < exp(Beta_*DtE)))) {
				z4CurrentConfig_[4*randomsite+drcn1] *= -1;
				z4CurrentConfig_[4*randomsite+drcn2] *= -1;
				z4VertexType_[randomsite] = vt4af1;
				
				if(nlk1.lr == 1) {
					z3r1CurrentConfig_[3*nbsite1] *= -1;
					z3r1CurrentConfig_[3*nbsite1+1] *= -1;
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite1] *= -1;
					z3r2CurrentConfig_[3*nbsite1+1] *= -1;
				}
				
				z4CurrentConfig_[nbsite3*4+nlk3.ldir] *= -1;
				z4VertexType_[nbsite3] = vt4af2;
				
				if(nlk2.lr ==1) {
					z3r1CurrentConfig_[3*nbsite2+nlk2.ldir] *= -1;
					z3r1VertexType_[nbsite2] = vt3af;
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite2+nlk2.ldir] *= -1;
					z3r2VertexType_[nbsite2] = vt3af;

				}	
				
				z4VertexCount_[vt4bf1] --;
				z4VertexCount_[vt4bf2] --;
				z4VertexCount_[vt4af1] ++;
				z4VertexCount_[vt4af2] ++;
				
				z3VertexCount_[vt3bf] --;
				z3VertexCount_[vt3af] ++;
				
				return 1;				
			}			
			else return 0;
			
		} // one leg long
		
		else if(z3MixVertexType_(nlk2.lr, nbsite2) == 3 && z3MixCurrentConfig_(nlk1.lr, 2, nbsite1) != 0) {
			Link_ nlk3 = z3MixNeighborStruct_( nlk2.lr, 1-nlk2.ldir, nbsite2);
			int nbsite3 = nlk3.lx + 2*N_*nlk3.ly;
			
			int cn4[4];
			cn4[0] = z4CurrentConfig_[4*nbsite3];
			cn4[1] = z4CurrentConfig_[4*nbsite3+1];
			cn4[2] = z4CurrentConfig_[4*nbsite3+2];
			cn4[3] = z4CurrentConfig_[4*nbsite3+3];
			
			int cn3[3];
			cn3[0] = z3MixCurrentConfig_(nlk1.lr, 0, nbsite1);
			cn3[1] = z3MixCurrentConfig_(nlk1.lr, 1, nbsite1);
			cn3[2] = z3MixCurrentConfig_(nlk1.lr, 2, nbsite1);
			
			int vt4bf1, vt4bf2, vt3bf, vt4af1, vt4af2, vt3af;
			
			double Ebf, Eaf;
			
			vt4bf1 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt4bf2 = Decidez4VertexType_(cn4[0], cn4[1], cn4[2], cn4[3]);
			vt3bf = Decidez3VertexType_(cn3[0], cn3[1], cn3[2]);
			
			Ebf = Ez4_[vt4bf1] + Ez4_[vt4bf2] + Ez3_[vt3bf];
			
			c4[drcn1] *= -1;
			c4[drcn2] *= -1;
			cn4[nlk3.ldir] *= -1;
			cn3[nlk1.ldir] *= -1;
			
			vt4af1 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt4af2 = Decidez4VertexType_(cn4[0], cn4[1], cn4[2], cn4[3]);
			vt3af = Decidez3VertexType_(cn3[0], cn3[1], cn3[2]);
			
			Eaf = Ez4_[vt4af1] + Ez4_[vt4af2] + Ez3_[vt3af];
			
			double DtE = vt4bf1 == 0? Eaf-Ebf-log(2.)/Beta_ : Ebf-Eaf-log(2.)/Beta_;
			
			double prob = ran.randDblExc();
			
			//trickiy when dealing the probabilities
			if((vt4bf1 == 0 && (DtE < 0. || prob < exp(-Beta_*DtE))) || (vt4bf1 == 1 && (DtE > 0. || prob < exp(Beta_*DtE)))) {
				z4CurrentConfig_[4*randomsite+drcn1] *= -1;
				z4CurrentConfig_[4*randomsite+drcn2] *= -1;
				z4VertexType_[randomsite] = vt4af1;
				
				if(nlk2.lr == 1) {
					z3r1CurrentConfig_[3*nbsite2] *= -1;
					z3r1CurrentConfig_[3*nbsite2+1] *= -1;
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite2] *= -1;
					z3r2CurrentConfig_[3*nbsite2+1] *= -1;
				}
				
				z4CurrentConfig_[nbsite3*4+nlk3.ldir] *= -1;
				z4VertexType_[nbsite3] = vt4af2;
				
				if(nlk1.lr ==1) {
					z3r1CurrentConfig_[3*nbsite1+nlk1.ldir] *= -1;
					z3r1VertexType_[nbsite1] = vt3af;
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite1+nlk1.ldir] *= -1;
					z3r2VertexType_[nbsite1] = vt3af;
					
				}		
				
				z4VertexCount_[vt4bf1] --;
				z4VertexCount_[vt4bf2] --;
				z4VertexCount_[vt4af1] ++;
				z4VertexCount_[vt4af2] ++;
				
				z3VertexCount_[vt3bf] --;
				z3VertexCount_[vt3af] ++;
				
				return 1;				
			}			
			else return 0;			
		} // other one leg long
		
		else if(z3MixVertexType_(nlk1.lr, nbsite1) == 3 && z3MixVertexType_(nlk2.lr, nbsite2) == 3) {
			int c4n1[4], c4n2[4];
			
			Link_ nlk3 = z3MixNeighborStruct_( nlk1.lr, 1-nlk1.ldir, nbsite1);
			Link_ nlk4 = z3MixNeighborStruct_( nlk2.lr, 1-nlk2.ldir, nbsite2);
			
			int nbsite3 = nlk3.lx + 2*N_*nlk3.ly;
			int nbsite4 = nlk4.lx + 2*N_*nlk4.ly;
			
			for(int i = 0; i < 4; i ++) {
				c4n1[i] = z4CurrentConfig_[4*nbsite3+i];
				c4n2[i] = z4CurrentConfig_[4*nbsite4+i];
			}
			
			int vt4bf1, vt4bf2, vt4bf3, vt4af1, vt4af2, vt4af3;
			
			double Ebf, Eaf;
			
			vt4bf1 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt4bf2 = Decidez4VertexType_(c4n1[0], c4n1[1], c4n1[2], c4n1[3]);
			vt4bf3 = Decidez4VertexType_(c4n2[0], c4n2[1], c4n2[2], c4n2[3]);
			
			c4[drcn1] *= -1;
			c4[drcn2] *= -1;
			c4n1[nlk3.ldir] *= -1;
			c4n2[nlk4.ldir] *= -1;
			
			vt4af1 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt4af2 = Decidez4VertexType_(c4n1[0], c4n1[1], c4n1[2], c4n1[3]);
			vt4af3 = Decidez4VertexType_(c4n2[0], c4n2[1], c4n2[2], c4n2[3]);
			
			Ebf = Ez4_[vt4bf1] + Ez4_[vt4bf2] + Ez4_[vt4bf3];
			
			Eaf = Ez4_[vt4af1] + Ez4_[vt4af2] + Ez4_[vt4af3];
				
			double DtE = vt4bf1 == 0? Eaf-Ebf-log(2.)/Beta_ : Ebf-Eaf-log(2.)/Beta_;
			
			double prob = ran.randDblExc();
			
			//trickiy when dealing the probabilities
			if((vt4bf1 == 0 && (DtE < 0. || prob < exp(-Beta_*DtE))) || (vt4bf1 == 1 && (DtE > 0. || prob < exp(Beta_*DtE)))) {
				z4CurrentConfig_[4*randomsite+drcn1] *= -1;
				z4CurrentConfig_[4*randomsite+drcn2] *= -1;
				z4VertexType_[randomsite] = vt4af1;
				
				if(nlk1.lr == 1) {
					z3r1CurrentConfig_[3*nbsite1] *= -1;
					z3r1CurrentConfig_[3*nbsite1+1] *= -1;
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite1] *= -1;
					z3r2CurrentConfig_[3*nbsite1+1] *= -1;
				}
					
				if(nlk2.lr == 1) {
					z3r1CurrentConfig_[3*nbsite2] *= -1;
					z3r1CurrentConfig_[3*nbsite2+1] *= -1;
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite2] *= -1;
					z3r2CurrentConfig_[3*nbsite2+1] *= -1;
				}
				
				z4CurrentConfig_[4*nbsite3+nlk3.ldir] *= -1;
				z4VertexType_[nbsite3] = vt4af2;
				
				z4CurrentConfig_[4*nbsite4+nlk4.ldir] *= -1;
				z4VertexType_[nbsite4] = vt4af3;
				
				z4VertexCount_[vt4bf1] --;
				z4VertexCount_[vt4bf2] --;
				z4VertexCount_[vt4bf3] --;
				
				z4VertexCount_[vt4af1] ++;
				z4VertexCount_[vt4af2] ++;
				z4VertexCount_[vt4af3] ++;
				
				return 1;
			}
			else return 0;
				
		} //both let long
		
		else if(z3MixCurrentConfig_(nlk1.lr, 2, nbsite1) != 0 && z3MixCurrentConfig_(nlk2.lr, 2, nbsite2) != 0){
			int c3n1[3], c3n2[3];

			if(nlk1.lr == 1) {
				c3n1[0] = z3r1CurrentConfig_[3*nbsite1];
				c3n1[1] = z3r1CurrentConfig_[3*nbsite1+1];
				c3n1[2] = z3r1CurrentConfig_[3*nbsite1+2];
			}
			
			else {
				c3n1[0] = z3r2CurrentConfig_[3*nbsite1];
				c3n1[1] = z3r2CurrentConfig_[3*nbsite1+1];
				c3n1[2] = z3r2CurrentConfig_[3*nbsite1+2];
			}
			
			if(nlk2.lr == 1) {
				c3n2[0] = z3r1CurrentConfig_[3*nbsite2];
				c3n2[1] = z3r1CurrentConfig_[3*nbsite2+1];
				c3n2[2] = z3r1CurrentConfig_[3*nbsite2+2];
			}
			
			else {
				c3n2[0] = z3r2CurrentConfig_[3*nbsite2];
				c3n2[1] = z3r2CurrentConfig_[3*nbsite2+1];
				c3n2[2] = z3r2CurrentConfig_[3*nbsite2+2];
			}
			
			
			double Ebf, Eaf;		//energy before and after flip
			
			int vt3bf1, vt3bf2, vt4bf, vt3af1, vt3af2, vt4af;		//vertex types before and after
			
			vt4bf = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			vt3bf1 = Decidez3VertexType_(c3n1[0], c3n1[1], c3n1[2]);
			vt3bf2 = Decidez3VertexType_(c3n2[0], c3n2[1], c3n2[2]);
			
			//int vcp4 = z4VertexType_[randomsite];
			//int vcp31 = nlk1.lr == 1 ? z3r1VertexType_[nbsite1] : z3r2VertexType_[nbsite1];
			//int vcp32 = nlk2.lr == 1 ? z3r1VertexType_[nbsite2] : z3r2VertexType_[nbsite2];
			
			Ebf = Ez4_[vt4bf] + Ez3_[vt3bf1] + Ez3_[vt3bf2];
			
			c4[drcn1] *= -1;
			c4[drcn2] *= -1;
			c3n1[nlk1.ldir] *= -1;
			c3n2[nlk2.ldir] *= -1;
			
			vt3af1 = Decidez3VertexType_(c3n1[0], c3n1[1], c3n1[2]);
			vt3af2 = Decidez3VertexType_(c3n2[0], c3n2[1], c3n2[2]);
			vt4af = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			
			Eaf = Ez4_[vt4af] + Ez3_[vt3af1] + Ez3_[vt3af2];
			
			double DtE = vt4bf == 0? Eaf-Ebf-log(2.)/Beta_ : Ebf-Eaf-log(2.)/Beta_;
			
			double prob = ran.randDblExc();
			
			
			//trickiy when dealing the probabilities
			if((vt4bf == 0 && (DtE < 0. || prob < exp(-Beta_*DtE))) || (vt4bf == 1 && (DtE > 0. || prob < exp(Beta_*DtE)))) {
				
				z4CurrentConfig_[4*randomsite+drcn1] *= -1;
				z4CurrentConfig_[4*randomsite+drcn2] *= -1;
				
				z4VertexType_[randomsite] = vt4af;				
				
				z4VertexCount_[vt4bf] --;
				z4VertexCount_[vt4af] ++;
				
				if(nlk1.lr == 1) {
					z3r1CurrentConfig_[3*nbsite1+nlk1.ldir] *= -1;
					z3r1VertexType_[nbsite1] = vt3af1; 
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite1+nlk1.ldir] *= -1;
					z3r2VertexType_[nbsite1] = vt3af1; 
				}
				
				z3VertexCount_[vt3bf1] --;
				z3VertexCount_[vt3af1] ++;
				
				if(nlk2.lr == 1) {
					z3r1CurrentConfig_[3*nbsite2+nlk2.ldir] *= -1;
					z3r1VertexType_[nbsite2] = vt3af2; 
				}
				
				else {
					z3r2CurrentConfig_[3*nbsite2+nlk2.ldir] *= -1;
					z3r2VertexType_[nbsite2] = vt3af2; 
				}
				
				z3VertexCount_[vt3bf2] --;
				z3VertexCount_[vt3af2] ++;
				
				Energy_ += Eaf - Ebf;
				
				return 1;
				
			}
			
			else return 0;		
		} //no long legs
		
		else return 0;
	} 
	
	else {return 0;}
	
}


//flip all current of a Z4 type 0 vertex
int DfShakti::Z4flip() {	
	int randomsite = ran.randDblExc()*4*N_*N_;
	
	if(z4VertexType_[randomsite] == 0) {
		Link_ nlk[4];
		for(int drcn = 0; drcn < 4; drcn ++)
			nlk[drcn] = r0NeighborStruct_[4*randomsite + drcn];
		
		int c4[4], nbsite[4];
		
		int c3[12];
		
		double Ebf, Eaf;
		
		for(int drcn = 0; drcn < 4; drcn ++) {
			c4[drcn] = z4CurrentConfig_[4*randomsite + drcn];
			nbsite[drcn] = nlk[drcn].lx + 2*N_*nlk[drcn].ly;
			
			for(int d3 = 0; d3 < 3; d3 ++) {
				c3[3*drcn + d3] = ( nlk[drcn].lr == 1 ? z3r1CurrentConfig_[nbsite[drcn]*3+d3] :  z3r2CurrentConfig_[nbsite[drcn]*3+d3]);
			}
		}
			
		int vt4bf, vt4af; 
		int vt3bf[4], vt3af[4];
		
		vt4bf = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
		
		for(int drcn = 0; drcn < 4; drcn ++) {
			vt3bf[drcn] = Decidez3VertexType_(c3[3*drcn], c3[3*drcn+1], c3[3*drcn+2]);
		}
		
		Ebf = Ez4_[vt4bf] + Ez3_[vt3bf[0]] + Ez3_[vt3bf[1]] + Ez3_[vt3bf[2]] + Ez3_[vt3bf[3]];
				
		c4[0] *= -1;
		c4[1] *= -1;
		c4[2] *= -1;
		c4[3] *= -1;
		
		for(int drcn = 0; drcn < 4; drcn ++)
			c3[3*drcn+nlk[drcn].ldir] *= -1;
		
		vt4af = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
		
		for(int drcn = 0; drcn < 4; drcn ++) {
			vt3af[drcn] = Decidez3VertexType_(c3[3*drcn], c3[3*drcn+1], c3[3*drcn+2]);
		}
		
		Eaf = Ez4_[vt4bf] + Ez3_[vt3af[0]] + Ez3_[vt3af[1]] + Ez3_[vt3af[2]] + Ez3_[vt3af[3]];
		
		double prob = ran.randDblExc();
		
		if(Eaf<Ebf || prob < exp(-Beta_*(Eaf-Ebf))) {
			for(int drcn = 0; drcn < 4; drcn ++) {
				z4CurrentConfig_[4*randomsite + drcn] *= -1;
				
				if( nlk[drcn].lr == 1 ) z3r1CurrentConfig_[nbsite[drcn]*3+nlk[drcn].ldir] *= -1; 
				else z3r2CurrentConfig_[nbsite[drcn]*3+nlk[drcn].ldir] *= -1;
			}
			
			z4VertexType_[randomsite] = vt4af;
			for(int drcn = 0; drcn < 4; drcn ++) {
				
				if( nlk[drcn].lr == 1 ) z3r1VertexType_[nbsite[drcn]] = vt3af[drcn]; 
				else z3r2VertexType_[nbsite[drcn]] = vt3af[drcn];
			
			}
			
			z4VertexCount_[vt4bf] --;
			z4VertexCount_[vt4af] ++;
			
			for(int drcn = 0; drcn < 4; drcn ++) {
				z3VertexCount_[vt3bf[drcn]] --;
				z3VertexCount_[vt3af[drcn]] ++;
			}
			
			Energy_ += Eaf - Ebf;
			
			return 1;
			
		} //if accepted
		
		else return 0;
	} // if tpye = 0
	
	else {
		return 0; }
	
}

int DfShakti::Z3flip() {
	int randomsite = ran.randDblExc()*4*N_*N_;
	
	int result = 0;
	
	if(z3r1VertexType_[randomsite] == 3) {
		Link_ nlk1 = r1NeighborStruct_[randomsite*3];
		Link_ nlk2 = r1NeighborStruct_[randomsite*3+1];
		int nbsite1 = nlk1.lx + 2*N_*nlk1.ly;
		int nbsite2 = nlk2.lx + 2*N_*nlk2.ly;
		
		int c3[2]; 
		int c4n1[4], c4n2[4];
		
		c3[0] = z3r1CurrentConfig_[randomsite*3];
		c3[1] = z3r1CurrentConfig_[randomsite*3+1];
		
		for(int dd = 0; dd < 4; dd ++) {
			c4n1[dd] = z4CurrentConfig_[nbsite1*4+dd];
			c4n2[dd] = z4CurrentConfig_[nbsite2*4+dd];
		}
		
		double Ebf, Eaf;
		
		int vt4bf1, vt4bf2, vt3bf, vt4af1, vt4af2, vt3af;
		
		vt3bf = Decidez3VertexType_(c3[0], c3[1], 0);
		vt4bf1 = Decidez4VertexType_(c4n1[0], c4n1[1], c4n1[2], c4n1[3]);
		vt4bf2 = Decidez4VertexType_(c4n2[0], c4n2[1], c4n2[2], c4n2[3]);
		
		Ebf = Ez4_[vt4bf1] + Ez4_[vt4bf2] + Ez3_[vt3bf];
		
		c3[0] *= -1;
		c3[1] *= -1;
		
		c4n1[nlk1.ldir] *= -1;
		c4n2[nlk2.ldir] *= -1;
		
		vt3af = Decidez3VertexType_(c3[0], c3[1], 0);
		vt4af1 = Decidez4VertexType_(c4n1[0], c4n1[1], c4n1[2], c4n1[3]);
		vt4af2 = Decidez4VertexType_(c4n2[0], c4n2[1], c4n2[2], c4n2[3]);
		
		Eaf = Ez4_[vt4af1] + Ez4_[vt4af2] + Ez3_[vt3af];
		
		double prob = ran.randDblExc();
		
		if(Eaf<Ebf || prob < exp(-Beta_*(Eaf-Ebf))) {
			z3r1CurrentConfig_[randomsite*3] *= -1;
			z3r1CurrentConfig_[randomsite*3+1] *= -1;
			
			z4CurrentConfig_[nbsite1*4+nlk1.ldir] *= -1;
			z4CurrentConfig_[nbsite2*4+nlk2.ldir] *= -1;
			
			z3r1VertexType_[randomsite] = vt3af;
			z4VertexType_[nbsite1] = vt4af1;
			z4VertexType_[nbsite2] = vt4af2;
			
			z4VertexCount_[vt4bf1] --;
			z4VertexCount_[vt4bf2] --;
			z4VertexCount_[vt4af1] ++;
			z4VertexCount_[vt4af2] ++;
			
			result = 1;
		} //if accepted(r1 case)
		
	} // if r1 case
		
	if(z3r2VertexType_[randomsite] == 3) {
		Link_ nlk1 = r2NeighborStruct_[randomsite*3];
		Link_ nlk2 = r2NeighborStruct_[randomsite*3+1];
		int nbsite1 = nlk1.lx + 2*N_*nlk1.ly;
		int nbsite2 = nlk2.lx + 2*N_*nlk2.ly;
		
		int c3[2]; 
		int c4n1[4], c4n2[4];
		
		c3[0] = z3r2CurrentConfig_[randomsite*3];
		c3[1] = z3r2CurrentConfig_[randomsite*3+1];
		
		for(int dd = 0; dd < 4; dd ++) {
			c4n1[dd] = z4CurrentConfig_[nbsite1*4+dd];
			c4n2[dd] = z4CurrentConfig_[nbsite2*4+dd];
		}
		
		double Ebf, Eaf;
		
		int vt4bf1, vt4bf2, vt3bf, vt4af1, vt4af2, vt3af;
		
		vt3bf = Decidez3VertexType_(c3[0], c3[1], 0);
		vt4bf1 = Decidez4VertexType_(c4n1[0], c4n1[1], c4n1[2], c4n1[3]);
		vt4bf2 = Decidez4VertexType_(c4n2[0], c4n2[1], c4n2[2], c4n2[3]);
		
		Ebf = Ez4_[vt4bf1] + Ez4_[vt4bf2] + Ez3_[vt3bf];
		
		c3[0] *= -1;
		c3[1] *= -1;
		
		c4n1[nlk1.ldir] *= -1;
		c4n2[nlk2.ldir] *= -1;
		
		vt3af = Decidez3VertexType_(c3[0], c3[1], 0);
		vt4af1 = Decidez4VertexType_(c4n1[0], c4n1[1], c4n1[2], c4n1[3]);
		vt4af2 = Decidez4VertexType_(c4n2[0], c4n2[1], c4n2[2], c4n2[3]);
		
		Eaf = Ez4_[vt4af1] + Ez4_[vt4af2] + Ez3_[vt3af];
		
		double prob = ran.randDblExc();
		
		if(Eaf<Ebf || prob < exp(-Beta_*(Eaf-Ebf))) {
			z3r2CurrentConfig_[randomsite*3] *= -1;
			z3r2CurrentConfig_[randomsite*3+1] *= -1;
			
			z4CurrentConfig_[nbsite1*4+nlk1.ldir] *= -1;
			z4CurrentConfig_[nbsite2*4+nlk2.ldir] *= -1;
			
			z3r2VertexType_[randomsite] = vt3af;
			z4VertexType_[nbsite1] = vt4af1;
			z4VertexType_[nbsite2] = vt4af2;
			
			z4VertexCount_[vt4bf1] --;
			z4VertexCount_[vt4bf2] --;
			z4VertexCount_[vt4af1] ++;
			z4VertexCount_[vt4af2] ++;
			
			result = 1;
		} //if accepted(r2 case)
		
	} // if r2 case
	return result;
}

double DfShakti::ShowEnergy() {
	double E = Ez4_[0]*z4VertexCount_[0] + Ez4_[1]*z4VertexCount_[1] + Ez4_[2]*z4VertexCount_[2] + Ez4_[3]*z4VertexCount_[3] + \
	Ez3_[0]*z3VertexCount_[0] + Ez3_[1]*z3VertexCount_[1] + Ez3_[2]*z3VertexCount_[2] + Ez3_[3]*z3VertexCount_[3] + Ez3_[4]*z3VertexCount_[4];
	Energy_ = E;
	
	return E;
}

int DfShakti::Showz3Vertex(int tp) {
	if(tp < 0 || tp > 5) {
		std::cout << "DON'T MESS WITH VERTEX TYPE!!!!\n";
		return -1;
	}
	
	else return z3VertexCount_[tp];
}

int DfShakti::Showz4Vertex(int tp) {
	if(tp < 0 || tp > 4) {
		std::cout << "DON'T MESS WITH VERTEX TYPE!!!!\n";
		return -1;
	}
	
	else return z4VertexCount_[tp];
}

 void DfShakti::SaveVertex(const char* filename) {
	FILE* oput;
	oput = fopen(filename,"w");
	
	for(int i=0; i<4*N_*N_; i++) {
		int x = i%(2*N_);
		int y = i/(2*N_);
		
		fprintf(oput, "(%d,%d)\t(%d,%d,%d,%d)%d\t(%d,%d,%d)%d\t(%d,%d,%d)%d\n",\
				x,y,z4CurrentConfig_[4*i],z4CurrentConfig_[4*i+1],z4CurrentConfig_[4*i+2],z4CurrentConfig_[4*i+3],z4VertexType_[i],\
				z3r1CurrentConfig_[3*i],z3r1CurrentConfig_[3*i+1],z3r1CurrentConfig_[3*i+2],z3r1VertexType_[i],\
				z3r2CurrentConfig_[3*i],z3r2CurrentConfig_[3*i+1],z3r2CurrentConfig_[3*i+2],z3r2VertexType_[i]);
	}
}

void DfShakti::ReadVertex(const char* filename) {
	FILE* input;
	input = fopen(filename, "r");
	
	int x ,y, vt4, vt3r1, vt3r2;
	int z4c[4];
	int z3r1c[3];
	int z3r2c[3];
	
	for(int i=0; i<4*N_*N_; i++) {
		fscanf(input, "(%d,%d)\t(%d,%d,%d,%d)%d\t(%d,%d,%d)%d\t(%d,%d,%d)%d\n",\
			   &x, &y, &z4c[0], &z4c[1], &z4c[2], &z4c[3], &vt4, \
			   &z3r1c[0], &z3r1c[1], &z3r1c[2], &vt3r1, \
			   &z3r2c[0], &z3r2c[1], &z3r2c[2], &vt3r2); 
		
		for(int j = 0; j < 4; j ++) {
			z4CurrentConfig_[4*i+j] = z4c[j];
			z4VertexType_[j] = vt4; 
		}
		
		for(int j = 0; j < 3; j++) {
			z3r1CurrentConfig_[3*i+j] = z3r1c[j];
			z3r2CurrentConfig_[3*i+j] = z3r2c[j];
			z3r1VertexType_[j] = vt3r1;
			z3r2VertexType_[j] = vt3r2;
		}
	}
	
	CountVertex();
}
	
void DfShakti::CountVertex() {
	int actz4v[4] = {0, 0, 0, 0};
	int actz3v[5] = {0, 0, 0, 0, 0};
	int vz4, vz3;
	
	int c4[4];
	int c3r1[3];
	int c3r2[3];
	for(int y = 0; y < 2*N_; y ++) {
		for(int x = 0; x < 2*N_; x ++) {
			int p4 = x + y*2*N_;
			for(int d = 0; d < 4; d++)
				c4[d] = z4CurrentConfig_[p4*4+d];
			
			vz4 = Decidez4VertexType_(c4[0], c4[1], c4[2], c4[3]);
			if(vz4 != z4VertexType_[p4]) std::cout << "z4 type mismatch!\n";
			actz4v[vz4] ++;
			
			for(int d = 0; d < 3; d++)
				c3r1[d] = z3r1CurrentConfig_[3*p4+d];
			
			vz3 = Decidez3VertexType_(c3r1[0], c3r1[1], c3r1[2]);
			if(vz3 != z3r1VertexType_[p4]) std::cout << "z3r1 type mismatch!\n";
			actz3v[vz3] ++;
			
			for(int d = 0; d < 3; d++)
				c3r2[d] = z3r2CurrentConfig_[3*p4+d];
			
			vz3 = Decidez3VertexType_(c3r2[0], c3r2[1], c3r2[2]);
			if(vz3 != z3r2VertexType_[p4]) std::cout << "z3r2 type mismatch!\n";
			actz3v[vz3] ++;
		}
	}
	
	for(int i = 0; i < 4; i ++) {
		
		if(z4VertexCount_[i] != actz4v[i]) {
			std::cout << "Initialize Vertex Count4\n";
			z4VertexCount_[i] = actz4v[i]; }
	}
	
	for(int i = 0; i < 5; i ++) {
		if(z3VertexCount_[i] != actz3v[i]) {
			std::cout << "Initialize Vertex Count3\n";
			z3VertexCount_[i] = actz3v[i]; }
	}
		
}	
	

	
void DfShakti::ResetTemp(double tp) {
	Beta_ = 1./tp; }

void DfShakti::SnapShot(const char* filename1, const char* filename2, const char* filename3, const char* filename4) {
	// file2 for spins on the islands
	FILE* output1;
	output1 = fopen(filename1, "w");
	FILE* output2;
	output2 = fopen(filename2, "w");
	FILE* output3;
	output3 = fopen(filename3, "w");
	FILE* output4;
	output4 = fopen(filename4, "w");
	int cc, drcn;
	
	for(int x = 0; x < 2*N_; x ++) {
		for(int y = 0; y < 2*N_; y ++) {
			
			if(z4VertexType_[x + 2*N_*y] == 0) {
				if((2*((x+y)%2)-1)*z4CurrentConfig_[4*(x+2*N_*y)] == 1) {
					fprintf(output3, "%d\t%d\n", 2*x, 2*y); }
				else {
					fprintf(output4, "%d\t%d\n", 2*x, 2*y); }
			}
			
			drcn = 0;
			cc = z4CurrentConfig_[4*(x+2*N_*y)+drcn];
				
			if(cc == 1) 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x , 2*y, 1, 0);
			else 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x+1, 2*y, -1, 0);
			
			drcn = 1;
			cc = z4CurrentConfig_[4*(x+2*N_*y)+drcn];
			
			if(cc == 1) 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x , 2*y, 0, 1);
			else 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x, 2*y+1, 0, -1);
			
			drcn = 2;
			cc = z4CurrentConfig_[4*(x+2*N_*y)+drcn];
			
			if(cc == 1) 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x , 2*y, -1, 0); 
			else 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x-1, 2*y, 1, 0);
			
			drcn = 3;
			cc = z4CurrentConfig_[4*(x+2*N_*y)+drcn];
			
			if(cc == 1) 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x , 2*y, 0, -1);
			else 
				fprintf(output1, "%d\t%d\t%d\t%d\n", 2*x, 2*y-1, 0, 1);
			
			drcn = 2;
			cc = z3r1CurrentConfig_[3*(x+2*N_*y)+drcn];
			
			if(cc == 1) {
				if((x+y)%2 == 0)
					fprintf(output2, "%d\t%d\t%d\t%d\n", 2*x+1, 2*y, 0, -2);
				else 
					fprintf(output2, "%d\t%d\t%d\t%d\n", 2*x+1, 2*y, 0, 2);
			}
			

			cc = z3r2CurrentConfig_[3*(x+2*N_*y)+drcn];
			
			if(cc == 1) {
				if((x+y)%2 == 0)
					fprintf(output2, "%d\t%d\t%d\t%d\n", 2*x, 2*y+1, 2, 0);
				else 
					fprintf(output2, "%d\t%d\t%d\t%d\n", 2*x, 2*y+1, -2, 0);
			}

		}
	}
	
	fclose(output1);
	fclose(output2);
	fclose(output3);
	fclose(output4);
}

int* DfShakti::Dn56() {
	int* result;
	result = new int[2];
	int n5 = 0; 
	int n6 = 0;
	int nx, ny;
	for(int x = 0; x < 2*N_; x ++) {
		for(int y = 0; y < 2*N_; y ++) {
			nx = x==2*N_-1 ? 0 : x+1;
			ny = y==2*N_-1 ? 0 : y+1;

			if((x+y)%2 == 0) {
				if(z3r2VertexType_[x+2*N_*y] == 1 && z3r2VertexType_[nx+2*N_*y] == 1)
					n6 ++;
				
				if(z3r1VertexType_[x+2*N_*y] == 1 && z3r1VertexType_[x+2*N_*ny] == 1)
					n5 ++;
			}
			
			if((x+y)%2 == 1) {
				if(z3r2VertexType_[x+2*N_*y] == 1 && z3r2VertexType_[nx+2*N_*y] == 1)
					n5 ++;
				
				if(z3r1VertexType_[x+2*N_*y] == 1 && z3r1VertexType_[x+2*N_*ny] == 1)
					n6 ++;
			}
		}
	}
	result[0] = n5;
	result[1] = n6;
	return result;
			
}

//find the z4 vertex Ising order parameter, only consider z3 type 0 vertex.
int DfShakti::IsingM(){
	int ism = 0, pos;
	for(int x = 0; x < 2*N_; x ++) {
		for(int y = 0; y < 2*N_; y ++) {
			pos = x + 2*N_*y;
			if(z4VertexType_[pos] == 0) {
				if(z4CurrentConfig_[4*pos]*(2*((x+y)%2)-1) == 1)
					ism ++;
				else ism --;
			}
		}
	}
	
	return ism;
}

//v2. Count all vertex types.
int DfShakti::IsingMv2(){
	int ism = 0, pos;
	for(int x = 0; x < 2*N_; x ++) {
		for(int y = 0; y < 2*N_; y ++) {
			pos = x + 2*N_*y;
			ism += z4CurrentConfig_[4*pos]*(2*((x+y)%2)-1);
			ism -= z4CurrentConfig_[4*pos+1]*(2*((x+y)%2)-1);

		}
	}
	
	return ism;
}

// |replica overlap|^2. k denotes the momentum. NOTE::because of the boundary condition, the periodicity is 2. So Q=2*k*\pi/L, ----k should be even in the code----
//The overlap now is complex, and this function returns the modulus square
double DfShakti::OverLap(const DfShakti &s, double k) {
    if(N_ != s.N_ || Beta_ != s.Beta_ || tDf_ != s.tDf_)
        std::cout << "OverLap mismatch!\n";
    
    double result1 = 0., result2 = 0.; int pos;
    for(int i = 0; i < 2*N_; ++i) {
        for(int j = 0; j < 2*N_; ++j) {
            pos = i + 2*N_*j;
            
/*            result1 += z4CurrentConfig_[4*pos] * s.z4CurrentConfig_[4*pos] * cos(PI*k/4./N_ + PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos] * s.z4CurrentConfig_[4*pos] * sin(PI*k/4./N_ + PI*k*i/N_);
            
            result1 += z4CurrentConfig_[4*pos+2] * s.z4CurrentConfig_[4*pos+2] * cos(PI*k*i/N_ - PI*k/4./N_);
            result2 += z4CurrentConfig_[4*pos+2] * s.z4CurrentConfig_[4*pos+2] * sin(PI*k*i/N_ - PI*k/4./N_);
            
            result1 += z4CurrentConfig_[4*pos+1] * s.z4CurrentConfig_[4*pos+1] * cos(PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos+1] * s.z4CurrentConfig_[4*pos+1] * sin(PI*k*i/N_);
            
            result1 += z4CurrentConfig_[4*pos+3] * s.z4CurrentConfig_[4*pos+3] * cos(PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos+3] * s.z4CurrentConfig_[4*pos+3] * sin(PI*k*i/N_);
            
            
            if((i+j)%2==0) {
                result1 += z3r1CurrentConfig_[3*pos+2] * s.z3r1CurrentConfig_[3*pos+2] * cos(PI*k*(i+0.5)/N_);
                result2 += z3r1CurrentConfig_[3*pos+2] * s.z3r1CurrentConfig_[3*pos+2] * sin(PI*k*(i+0.5)/N_);
                
                result1 += z3r2CurrentConfig_[3*pos+2] * s.z3r2CurrentConfig_[3*pos+2] * cos(PI*k*(i-0.5)/N_);
                result2 += z3r2CurrentConfig_[3*pos+2] * s.z3r2CurrentConfig_[3*pos+2] * sin(PI*k*(i-0.5)/N_);
            } */
            
            result1 += z4CurrentConfig_[4*pos] * s.z4CurrentConfig_[4*pos] * cos(PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos] * s.z4CurrentConfig_[4*pos] * sin(PI*k*i/N_);
            
            result1 += z4CurrentConfig_[4*pos+2] * s.z4CurrentConfig_[4*pos+2] * cos(PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos+2] * s.z4CurrentConfig_[4*pos+2] * sin(PI*k*i/N_);
            
            result1 += z4CurrentConfig_[4*pos+1] * s.z4CurrentConfig_[4*pos+1] * cos(PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos+1] * s.z4CurrentConfig_[4*pos+1] * sin(PI*k*i/N_);
            
            result1 += z4CurrentConfig_[4*pos+3] * s.z4CurrentConfig_[4*pos+3] * cos(PI*k*i/N_);
            result2 += z4CurrentConfig_[4*pos+3] * s.z4CurrentConfig_[4*pos+3] * sin(PI*k*i/N_);
            
            
            if((i+j)%2==0) {
                result1 += z3r1CurrentConfig_[3*pos+2] * s.z3r1CurrentConfig_[3*pos+2] * cos(PI*k*i/N_);
                result2 += z3r1CurrentConfig_[3*pos+2] * s.z3r1CurrentConfig_[3*pos+2] * sin(PI*k*i/N_);
                
                result1 += z3r2CurrentConfig_[3*pos+2] * s.z3r2CurrentConfig_[3*pos+2] * cos(PI*k*i/N_);
                result2 += z3r2CurrentConfig_[3*pos+2] * s.z3r2CurrentConfig_[3*pos+2] * sin(PI*k*i/N_);
            }            
        }
    }
    return (result1*result1 + result2*result2);
}


void DfShakti::GenerateDefects() {
	
	FILE* input;
	
	std::stringstream output_name;
	output_name << "DefectsN" << N_ << "T" << tDf_ << "node"  << Node_ << ".dat";
	const std::string &temp=output_name.str();
	const char *name=temp.c_str();
	input = fopen(name, "r");
	
	if(input == NULL) {
		
		std::vector<int> AvlbSites;
		
		int NAsites = 4*N_*N_; //# of available sites
		
		AvlbSites.resize(4*N_*N_); // list of av sites
		
		for(int i = 0; i < 4*N_*N_; i ++)
			AvlbSites[i] = i;
		
		for(int j = 0; j < tDf_; j ++) {
			int randidx = ran.randDblExc()*NAsites;
			
			Defects_[j] = AvlbSites[randidx];
			
			AvlbSites.erase(AvlbSites.begin() + randidx);
			
			NAsites --;
		}
		
		FILE* output;
				
		output = fopen(name, "w");
		
		for(int i = 0; i < tDf_; i ++)
			fprintf(output, "%d\n", Defects_[i]);
		
		fclose(output);
	}		
	
	else {
		for(int i = 0; i < tDf_; i++)
			if(fscanf(input, "%d\n", &Defects_[i])==EOF)
				std::cout << "Not Enough defects found!\n";
	}
	
}
	
void DfShakti::Copy(const DfShakti& s) {
	if(N_ != s.N_) {std::cout << "size mismatch!\n"; exit(1); }
	int pos;
	for(int y = 0; y < 2*N_; y ++) {
		for(int x = 0; x < 2*N_; x ++) {
			pos = x + 2*N_*y;
			z4CurrentConfig_[4*pos] = s.z4CurrentConfig_[4*pos];
			z4CurrentConfig_[4*pos+1] = s.z4CurrentConfig_[4*pos+1];
			z4CurrentConfig_[4*pos+2] = s.z4CurrentConfig_[4*pos+2];
			z4CurrentConfig_[4*pos+3] = s.z4CurrentConfig_[4*pos+3];
			
			z3r1CurrentConfig_[3*pos] = s.z3r1CurrentConfig_[3*pos];
			z3r1CurrentConfig_[3*pos+1] = s.z3r1CurrentConfig_[3*pos+1];
			z3r1CurrentConfig_[3*pos+2] = s.z3r1CurrentConfig_[3*pos+2];
			
			z3r2CurrentConfig_[3*pos] = s.z3r2CurrentConfig_[3*pos];
			z3r2CurrentConfig_[3*pos+1] = s.z3r2CurrentConfig_[3*pos+1];
			z3r2CurrentConfig_[3*pos+2] = s.z3r2CurrentConfig_[3*pos+2];
			
			z4VertexType_[pos] = s.z4VertexType_[pos];
			z3r1VertexType_[pos] = s.z3r1VertexType_[pos];
			z3r2VertexType_[pos] = s.z3r2VertexType_[pos];
		}
	}
	
	for(int i = 0; i < 4; i ++)
		z4VertexCount_[i] = s.z4VertexCount_[i];
	
	for(int j = 0; j < 5; j ++)
		z3VertexCount_[j] = s.z3VertexCount_[j];
	
	Energy_ = s.Energy_;
						
}
#endif
