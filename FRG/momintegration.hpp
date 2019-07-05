/*
		File containing the integrands of the flow equations


		The most important thing in this file is to choose the values of
		LPAp and DH, so that the functions solves for the desired approximation.

			LPA:	LPAp = DH = 0
			LPAp: LPAp = 1, DH = '
			Full: LPAp = DH = 1'
*/

int flowequations(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {

	// Choosing approximation
	double	LPAp = 1;
	double	DH = 1;

	// Definition of rg-variables, momentum grid and hcubature external variables
	state_type x0(9);
	double P2[NL2];
	const double K = ((double*)fdata)[0];
	for(unsigned int i1 = 0; i1 < 9; i1++){   x0[i1] = ((double*)fdata)[1 + i1];   }

	// Definition of relevant quantities
	double expr10,expr11,expr12,expr13, expr14;
	double rk1, rk2,drk1, drk2;
	double rb1, rb2;
	double rf1, rf2;
  double rdf1 , rdb1;
	double sqrtP2,sqrtpPp,sqrtpPm,pPm,pPp,pPm2,pPp2, prodpP;
	double Dpsi1, Dp1, Ds1;
	double Dpsi2, Dp2,Ds2;
	double hq12, hq13, hq14, hq16, hq17, hq18;
	double hq21, hq22, hq23, hq24, hq25, hq26, hq27, hq28;
	double hq31, hq32, hq33, hq34, hq35, hq36, hq37, hq38;
	double hq41, hq42, hq43, hq44, hq45, hq46, hq47, hq48;
	double hq52, hq53, hq54, hq56, hq57, hq58;
	double hq61, hq62, hq63, hq64, hq65, hq66, hq67, hq68;
	double hq71, hq72, hq73, hq74, hq75, hq76, hq77, hq78;
	double hq81, hq82, hq83, hq84, hq85, hq86, hq87, hq88;
	double z12, z13, z14, z23, z24, z34;
	double aq1, aq2;
	double fP1,fP2;
	double fs1,fs2;


	// Assignation of hcubature integration variables and jacobian (multiplier)
	double p2 = x[0];
  double z = x[1];
	double sqrtp2 = sqrt(p2);
	double multiplier = (1.0/(pow(2.0*PI,3)))*sqrt(1.0-x[1]*x[1])*p2;

	// Mass terms in every RG-step
	double mp = x0[1];
	double ms = x0[1] + 2*Rho0*x0[2];

	// Interpolated dressings
	double fP0 = splineFP(p2);
	double fs0 = splineFS(p2);
	double aq0 = splineA(p2);
	double hq15 = YK(p2,p2,-1);

	// Computation of every regulator term
	double rk0 = Rk(p2,K);
	double drk0 = dRk(p2,K);
	double rb0 = p2*rk0;
	double rdb0 = p2*drk0;
	double rf0 = -1.0 + sqrt(1.0 + rk0);
	double rdf0 = drk0/(2.0*(1.0 + rf0));

	// Assignation of regulator term (in case they are not comparable)
	double rdbetabs = rdb0;
	double rdbetabp = rdb0;
	double rdfetaf = rdf0;

	// Computation of propagators
	double Dp0 = Meson_Propagator(fP0,p2,mp,rb0);
	double Ds0 = Meson_Propagator(fs0,p2,ms,rb0);
	double Dpsi0 = Fermion_Propagator(aq0,p2,hq15,rf0);


	// Evaluate the expressions of the expressions folder to compute potential terms

	double expr0=
	#include "Expressions/flowpot0.txt"			// Innecessary but left for completeness
	;

	double expr1 =
	#include "Expressions/flowpot1.txt"
	;

	double expr2 =
	#include "Expressions/flowpot2.txt"
	;

	double expr3 =
	#include "Expressions/flowpot3.txt"
	;

	double expr4 =
	#include "Expressions/flowpot4.txt"
	;

	double expr5 =
	#include "Expressions/flowpot5.txt"
	;

	double expr6 =
	#include "Expressions/flowpot6.txt"
	;

	double expr7 =
	#include "Expressions/flowpot7.txt"
	;

	double expr8 =
	#include "Expressions/flowpot8.txt"
	;


	// Assign final value of integrands
	fval[0] = multiplier*expr0;
	fval[1] = multiplier*expr1;
	fval[2] = multiplier*expr2;
	fval[3] = multiplier*expr3;
	fval[4] = multiplier*expr4;
	fval[5] = multiplier*expr5;
	fval[6] = multiplier*expr6;
	fval[7] = multiplier*expr7;
	fval[8] = multiplier*expr8;


	// Initialisation of parallelised procedure using openMP
	// In this loop one computes external momentum dependent dressings
	// Variables defined private in this way for optimisation's sake
	omp_set_num_threads(8);
	#pragma omp parallel for private(sqrtP2,prodpP,pPm,sqrtpPm,pPp,sqrtpPp,
		fP2,fs2,aq2,z12,z13,z14,z23,z24,z34,hq12,hq18,hq28,hq45,hq46,hq48,hq65,
		rk2,rb2,rf2,Dp2,Ds2,Dpsi2,expr10,expr11,expr12,expr13,expr14)
	for(unsigned int i1 = 0; i1 < NL2; i1++)
	{

		// External momentum dependent parameters
		P2[i1] = ((double*)fdata)[10 + i1];
		sqrtP2 = sqrt(P2[i1]);
		prodpP = sqrtp2*sqrtP2*x[1];
		pPm = p2 + P2[i1] - 2.0*prodpP;
		pPp = p2 + P2[i1] + 2.0*prodpP;
		sqrtpPp = sqrt(pPp);
		sqrtpPm = sqrt(pPm);

		// External momentum dependent dressing interpolation
		fP2 = splineFP(pPp);
		fs2 = splineFS(pPp);
		aq2 = splineA(pPp);

		// Calculation of relative angles (from Mathematica file)
		z12 = z;
		z13 = (sqrtp2 - sqrtP2*z)/sqrtpPm;
		z14 = (sqrtp2 + sqrtP2*z)/sqrtpPp;
		z23 = -(sqrtP2 - sqrtp2*z)/sqrtpPm;
		z24 = (sqrtP2 + sqrtp2*z)/sqrtpPp;
		z34 = (p2 - P2[i1])/(sqrtpPm*sqrtpPp);

		// Interpolation of the momentum dependent yukawas
		hq12 = YK(p2,P2[i1],z12);
		hq18 = YK(p2,pPp,-z14);
		hq28 = YK(P2[i1],pPp,-z24);
		hq45 = YK(pPp,p2,-z14);
		hq46 = YK(pPp,P2[i1],-z24);
		hq48 = YK(pPp,pPp,-1);
		hq65 = YK(P2[i1],p2,z12);

		// Calculation of the momentum dependent regulator terms
		rk2 = Rk(pPp,K);
		rb2 = pPp*rk2;
		rf2 = -1.0 + sqrt(1.0 + rk2);

		// Calculation of the momentum dependent propagator terms
		Dp2 = Meson_Propagator(fP2,pPp,mp,rb2);
		Ds2 = Meson_Propagator(fs2,pPp,ms,rb2);
		Dpsi2 = Fermion_Propagator(aq2,pPp,hq48,rf2);


		// Evaluation of the expressions for the yukawa, Zp, Zs, Zq and DH

		expr10 =
		#include "Expressions/quarkpropb.txt"
		;

		expr11 = LPAp *
		#include "Expressions/anomalousp0.txt"
		;

		expr12 = LPAp *
		#include "Expressions/anomalouss0.txt"
		;

		expr13 = LPAp *
		#include "Expressions/anomalousf0.txt"
		;

		expr14 = LPAp * DH *
		#include "Expressions/dynamical.txt"
		;


		// Assigning values and properly normalising (Mathematica)
		fval[9 + i1] = multiplier*expr10;
		fval[9 + i1 + NL2] = multiplier*expr11/P2[i1];
		fval[9 + i1 + 2*NL2] = multiplier*expr12/P2[i1];
		fval[9 + i1 + 3*NL2] = multiplier*expr13/P2[i1];
		fval[9 + i1 + 4*NL2] = multiplier*expr14;;
	}

return 0; // success
}
