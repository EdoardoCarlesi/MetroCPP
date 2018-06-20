#include <math.h>
#include <malloc.h>
#include <iostream>
#include <string>

#include "general.h"
#include "Particle.h"
#include "Halo.h"

using namespace std;


Halo::Halo()
{
	mTot = 0.0; 	mGas = 0.0; 	mDM = 0.0; 	
	rVir = 0.0;	lambda = 0.0;
	ID = 0;		nPart = 0;	hostID = 0;
	isToken = false;
};


Halo::~Halo()
{
	//delete Part;
};


float Halo::Distance(float *Pos)
{
	int iX = 0;
	float dX = 0.0;

	for (iX = 0; iX < 3; iX++)
	{
		dX += pow(Pos[iX] - X[iX], 2);	
	}
	
	dX = sqrt(dX);
	
	return dX;
};


void Halo::ReadLineAHF(const char * lineRead)
{
	float dummy;
	
//ID(1)  hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)        lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)        SurfP(41)       Phi0(42)        cNFW(43)        n_gas(44)       M_gas(45)       lambda_gas(46)  lambdaE_gas(47) Lx_gas(48)      Ly_gas(49)      Lz_gas(50)      b_gas(51)       c_gas(52)       Eax_gas(53)     Eay_gas(54)     Eaz_gas(55)     Ebx_gas(56)     Eby_gas(57)     Ebz_gas(58)     Ecx_gas(59)     Ecy_gas(60)     Ecz_gas(61)     Ekin_gas(62)    Epot_gas(63)    n_star(64)      M_star(65)      lambda_star(66) lambdaE_star(67)        Lx_star(68)     Ly_star(69)     Lz_star(70)     b_star(71)      c_star(72)      Eax_star(73)    Eay_star(74)    Eaz_star(75)    Ebx_star(76)    Eby_star(77)    Ebz_star(78)    Ecx_star(79)    Ecy_star(80)    Ecz_star(81)    Ekin_star(82)   Epot_star(83)   
	// Col:		    1   2    3  4  5  6  7  8  9 10 11 
	sscanf(lineRead, "%llu %llu %d %f %d %f %f %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f\
			  %f   %f   %f", 
			&ID, &hostID, &nSub, &mTot, &nPart, &X[0], &X[1], &X[2], &V[0], &V[1], &V[2], 	// 11
			&rVir, &dummy, &dummy, &dummy, &dummy, &vMax, &dummy, &sigV, &lambda, &dummy, // 21
			&L[0], &L[1], &L[2]); 

};


void Halo::Info(void)
{
	cout << "Task " << locTask << " Halo ID " << ID << endl;
	printf("Mtot: %.3e, ID: %.llu, Npart: %d, X:(%.3f, %.3f, %.3f), V:(%.3f, %.3f, %.f) \n", 
		mTot, ID, nPart, X[0], X[1], X[2], V[0], V[1], V[2]);

};

