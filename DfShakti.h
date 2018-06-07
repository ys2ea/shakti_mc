/*
 *  DfShakti.h
 *  Spin ice on a Shakti lattice. Df means defects are introduced.
 *  Simple Metropolis update. 
 *	Definition of class
 *
 *  Created by Yifei Shi on 10/28/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DfSHAKTI_H
#define DfSHAKTI_H

#include<fstream>
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<vector>
#include <sstream>
#include "MersenneTwister.h"

// Specify a link, using: r(relative position in a unit cell), x, y, dir(direction)
struct Link_ {
	int lr;
	int lx;
	int ly;
	int ldir;
} ;	


class DfShakti {
private:
	bool Thermalized;

	int N_;		//Number of Unit cells. Lattice size is 4*N_; site = vertex in the text.
	int tDf_;	//Number of defects. smaller than 4*N_*N_;
	int tr1Ild_, tr2Ild_;	//Number of Islands.
	
	double Beta_;
	
	int Node_;	//for parallel compute
	
	//Unit cell is 4x4, there are 12 sites and 20 links in one unit cell. But in practice we use a smaller 2x2 cell as our elementary cell, with 3 sites and 5 links.
	//A site(vertex) is labeled by: (x, y) which are the coordinates of the ele cell; and (a, b), the coordinates inside the unit cell. 0<=x,y<2*N, 0<=a,b<2;
	//If x+y==odd, then middle link vertical, otherwise middle link horizental
	
	//~~~~~~~IMPORTANT. Site notation in the code: for z=4 sites: coord: (2*x,2*y) denoted by x+y*2N.
	// z=3 sites, coord: (2*x+1, 2*y) and (2*x, 2*y+1) donoted by: 2(x+y*2N)+r
	
	//energy of different vertex types.
	double Ez4_[4];
	double Ez3_[5];
		
	//std::vector<bool> HasIsland_;		//1 if there is an island on the cell, 0 if not. 2Nx2N dimensional
	std::vector<int> Defects_;		//List of coordinates of defects.
	//std::vector<int> Islands_;		//List of Islands
	
	// Since some of the z3 vertices are taken out, this tells the new index in CurrentConfig_ of the z3 vertices 
	//std::vector<int> r1NewIndex_;
	//std::vector<int> r2NewIndex_;
	
	// List of remaining r1 and r2 sites (x + 2*N_*y);
	std::vector<int> Listofr1_;
	std::vector<int> Listofr2_;
		
	//Neighbor structure for three types of vertices.
	std::vector<Link_> r0NeighborStruct_;
	std::vector<Link_> r1NeighborStruct_;
	std::vector<Link_> r2NeighborStruct_;
		
	std::vector<int> z4CurrentConfig_;
	std::vector<int> z4VertexType_;

	//Same for z3 vertex.
	std::vector<int> z3r1CurrentConfig_;    //going out is positive
	std::vector<int> z3r2CurrentConfig_;
	std::vector<int> z3r1VertexType_;
	std::vector<int> z3r2VertexType_;

	//translate from r, x, y, dir, to neighborstruct index;
	inline int z4NeighborCount_(int, int, int);
	inline int z3NeighborCount_(int, int, int);
	
	inline int Decidez4VertexType_(int, int, int, int); //decide the current type given currents, must be ordered!	
	inline int Decidez3VertexType_(int, int, int);
	
	//so you don't deal with r1 and r2 differently
	inline int z3MixVertexType_(int, int);
	inline int z3MixCurrentConfig_(int, int, int);
	inline Link_ z3MixNeighborStruct_(int, int, int);
	
	inline void SetCurrentConfig_(const Link_* lk, int crt);
	
	//There are total of 16 type of sites for z=4 sites, and 8 for z=3 sites. So this number is 0~15 for z=4 site and 0~7 for z=3 site. 

	double Energy_;			//energy of system
	
	const char* Ofname;
	
	MTRand ran;		//random number generator
	
	
	int z4VertexCount_[4];
	int z3VertexCount_[5];	//number of each vertex, updated for each spin flip
	
public:
	DfShakti(int, int, int, double, double, double, double, double, double, double, double, double, double);
	
	void Copy(const DfShakti & s);
		
	int Spinflip(int); //bool tells whether on z4 site or z3 site. return whether or not update is accepted
	
	int Doubleflip();	//flip 2 spins on the corner of z4 vertex
	
	int DoubleLongflip();	//better verion of doubleflip
	
	int Z4flip(); //flip all currents of a z4 vertex
	
	int Z3flip();	//flip currents of z3(3) vertex
	
	double ShowEnergy();	//Calculate and reset the energy.
	
	int Showz3Vertex(int);
	
	int Showz4Vertex(int);
	
	void SaveVertex(const char*);		//Save all the informatin about the config
	
	void ReadVertex(const char*);
	
	void CheckConsist();	//Debuging, test whether VertexCount is consistant
	
	void CountVertex();
    
    double OverLap(const DfShakti &s, double k);         //Edward-Anderson oder parameter in the replica overlap version
	
	void GenerateDefects();
	
	void ResetTemp(double);
	
	void SnapShot(const char*, const char*,const char*, const char*);		//take a snapshot of current config
	
	void ReadDefects();
	
	int* Dn56();					//number of n5 and n6 verteces in the mapped F1 model
	
	int IsingM();					//Ising order parameter
	int IsingMv2();
		
};

inline int DfShakti::z4NeighborCount_(int x, int y, int dir) {
	return 4*(x+2*N_*y)+dir; }

inline int DfShakti::z3NeighborCount_(int x, int y, int dir) {
	return 3*(x+2*N_*y)+dir; }

inline int DfShakti::z3MixVertexType_(int r, int site) {
	if(r!=1 && r!=2) {std::cout << "r wrong!\n"; exit(1); }
	return r == 1 ? z3r1VertexType_[site] : z3r2VertexType_[site]; 
}

inline int DfShakti::z3MixCurrentConfig_(int r, int drcn, int site) {
	if(r!=1 && r!=2) {std::cout << "r wrong!\n"; exit(1); }
	return r == 1 ? z3r1CurrentConfig_[3*site+drcn] : z3r2CurrentConfig_[3*site+drcn]; 
}

inline Link_ DfShakti::z3MixNeighborStruct_(int r, int drcn, int site) {
	if(r!=1 && r!=2) {std::cout << "r wrong!\n"; exit(1); }
	return r == 1 ? r1NeighborStruct_[3*site+drcn] : r2NeighborStruct_[3*site+drcn];
}
	
#endif
	

		
		

	
