#include "CohesiveSurface.h"
#include "OPS_Globals.h"
#include "OPS_Stream.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>


#include <string.h>
#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>




CohesiveSurface::CohesiveSurface(void)
	:Kalg(2, 2), Kpen(2, 2), sigma(2), s_di(2), s_di_nplus1(2) {
}


CohesiveSurface::CohesiveSurface(Matrix& Kpen, double c, double mu, double muR, double Gc, double dropDrift, bool elasticSolution)
	: Kalg(Kpen), Kpen(Kpen), KCommitted(Kpen), sigma(2), s_di(2), dD_ds(2), mu(mu), muR(muR), c(c), Gc(Gc), GfII(dropDrift), flagAlg(false), deltaLambda(0.0), s_di_nplus1(2), 
	  elasticSolution(elasticSolution), k(0.), redStiffness(1.),redStiffnessComm(1.), x(0.), xCommitted(0.) {
		
		D = 0.;
		Dnplus1 = D;
}

/*
CohesiveSurface::CohesiveSurface(Matrix& Kpen, double mu, double muR, double kStar, double ft, double sN0, double sNT0, double etaN, double etaNT, bool elasticSolution)
: Kalg(Kpen), Kpen(Kpen), KCommitted(Kpen), sigma(2), s_di(2), dD_ds(2), mu(mu), muR(muR), kStar(kStar), flagAlg(false), sN0(sN0), sNT0(sNT0), etaN(etaN), etaNT(etaNT), 
  deltaLambda(0.0), s_di_nplus1(2), elasticSolution(elasticSolution), k(0.), redStiffness(1.),redStiffnessComm(1.)  {

	  // change names for Sacco Alfano 2006  and my modification. 
	  	this->s01 = sN0;
		this->s02 = sNT0;

		this->sc1 = etaN;
		this->sc2 = etaNT;
		
		this->eta = 1.0 - s02/sc2;

		c = s02*Kpen(1,1);
		GfII = sc2*c/2.0;

		D = 0.;
		Dnplus1 = D;
}
*/

CohesiveSurface::~CohesiveSurface(void) {
}

double 
CohesiveSurface::signTau(double tau) {
	if (tau == 0.0)
		return 0.0;
	else if (tau > 0.0)
		return 1.0;
	else
		return -1.0;
}

double
CohesiveSurface::heaviside(double arg) {
	if (arg > 0.0)
		return 1.0;
	else
		return 0.0;
}

double 
CohesiveSurface::macauley(double arg) {
	if (arg > 0.0)
		return arg;
	else
		return 0.0;
}

double 
CohesiveSurface::getMu(double kk) {
      return muR; 
}

double 
CohesiveSurface::getdMu_dk(double kk) {
    return 0; 
}


void 
CohesiveSurface::updatePlasticStrain(Vector& sigmaTrial) {
	
	if (mu-muR>0) {
		// return mapping for constant muR, no iterations needed
		double yieldFunction = muR*sigmaTrial(0) + abs(sigmaTrial(1));
		deltaLambda = 0.0;
		s_di_nplus1 = s_di;

		if (yieldFunction > DBL_EPSILON) {
			deltaLambda = 1.0 / (Kpen(1,1)) * yieldFunction;
			s_di_nplus1(1) += deltaLambda * this->signTau(sigmaTrial(1));
			double sign = this->signTau(sigmaTrial(1));
			sigmaTrial(1) -= deltaLambda * Kpen(1,1) * sign;		
		} 


	} else {

    // return mapping for constant mu, no iterations needed
	double yieldFunction = mu*sigmaTrial(0) + abs(sigmaTrial(1));
	deltaLambda = 0.0;
	s_di_nplus1 = s_di;

	if (yieldFunction > DBL_EPSILON) {
		//opserr << "called plasticity\n";
		deltaLambda = 1.0 / Kpen(1, 1) * yieldFunction;
		s_di_nplus1(1) += deltaLambda * this->signTau(sigmaTrial(1));
		double sign = this->signTau(sigmaTrial(1));
		sigmaTrial(1) -= deltaLambda * Kpen(1,1) * sign;
		
	} 

	}

	return;
}

double 
CohesiveSurface::evolveDamage(const Vector& s) {
	// from Sacco 2006
	/*
	double beta = sqrt( pow(this->macauley(s(0))/s01, 2.0) + pow(s(1)/s02, 2.0) ) - 1.0;
	double Dtrial;
	if (beta == -1.0)
		Dtrial = -1.0;
	else
		Dtrial = 1.0 / eta *beta / (1.0 + beta);
	*/

	// my model
	/*
	//double muStar = (mu-muR)*Kpen(0,0)/Kpen(1,1);
	double muStar = (eta*(mu-2.0*muR)+muR+2*sqrt((eta-1)*eta*muR*(muR-mu))) *Kpen(0,0)/Kpen(1,1);
	double beta = sqrt( pow(this->macauley(s(0))/s01, 2.0) + pow(s(1)/(s02 + muStar*this->macauley(-s(0)) ), 2.0) ) - 1;
	double Dtrial;
	if (beta == -1.0)
		Dtrial = -1.0;
	else
		Dtrial = 1.0 / eta *beta / (1.0 + beta);
	
	
	if (Dtrial > 1.0) Dtrial = 1.0;

	
	opserr << "d = " << s(1) << ", ";
	opserr << "N = " << s(0)*Kpen(0,0) << ", ";
	opserr << "s2 = " << s02 << ", ";
	opserr << "D = " << Dtrial << ", ";
	opserr << "beta = " << beta << ", ";
	opserr << "eta = " << eta << ", ";
	opserr << "dstar= "  << s02 + muStar*this->macauley(-s(0)) << endln;
	*/

	double Dtrial;
	double n = min(0.0,Kpen(0,0)*s(0));
	
	double VMax = c - n*mu;  // mohr-coulomb yield criterion
	double dVMax_dn = -1.0*mu;

	double s0 = -n*muR/Kpen(1,1);
	double sMax = -muR*n/Kpen(1,1) + Gc*(VMax - (-n)*muR)/Kpen(1,1);
	double sDrop = GfII;  // constant drift model

	double zeta = 0.20; // capacity drop percentage
	double slope = zeta *(c-mu*n)/ (std::max(sDrop, 1.001*sMax) - sMax);

	if (abs(s(1))<s0) {
		x = 0.0;
		Dtrial = 0.0;
	} else {
		if (abs(s(1))<sMax) {
			x = (abs(s(1)) - (-muR*n/Kpen(1,1)))/(sMax - (-muR*n/Kpen(1,1)));
			Dtrial = (1.0 - 1.0/Gc)*pow(x, 1.0/(Gc - 1.0));
		} else {
			x = (abs(s(1)) - (-muR*n/Kpen(1,1)))/(sMax - (-muR*n/Kpen(1,1)));  // not used but calculated
			Dtrial = (Kpen(1,1)*Kpen(1,1)*abs(s(1)) + muR*n*slope + Kpen(1,1)*abs(s(1))*slope - muR*n*Gc*slope - Kpen(1,1)*VMax - Gc*slope*VMax) / (Kpen(1,1)*(muR*n + Kpen(1,1)*abs(s(1))));
		}
	}

	if (Dtrial < 0.0)   Dtrial = 0.0;
	if (Dtrial > 1.0)   Dtrial = 1.0;
	
	if (Dtrial == 1.0 || D > Dtrial) {
		flagAlg = false;
	    dD_ds.Zero();

	} else {
		flagAlg = true;
		dD_ds.Zero();

		// Sacco 2006
		/*
		dD_ds(0) = this->macauley(s(0))/pow(s01,2);
		dD_ds(1) = s(1)/pow(s02,2);
		dD_ds /= eta*pow(1+beta,3);
		*/

		// my model
		/*
		if (s(0)<0.0) {
			dD_ds(0) = muStar*abs(s(1)) / pow(s02 - muStar*s(0),2);
			dD_ds(1) = this->signTau(s(1)) / (s02 - muStar*s(0));
			dD_ds /= eta*pow(1+beta,2);

		} else {
			dD_ds(0) = this->macauley(s(0))/pow(s01,2);
			dD_ds(1) = s(1)/pow(s02,2);
			dD_ds /= eta*pow(1+beta,3);

		}
		*/

		/*
		double sign;
		if (s(1)>=0.0)
			sign = 1.0;
		else
			sign = -1.0;

	    if (abs(s(1))<=sMax && Dtrial>0.0) {
			dD_ds(0)  = Kpen(0,0) *( (c*muR*+ Kpen(1,1)*(mu - muR)*abs(s(1)))/(4.*pow(-(mu*n) + muR*n + c,2.0)) );
	    	dD_ds(1)  = sign*Kpen(1,1)/2.0 /(muR*n  +Kpen(1,1)*sMax);
	    } else {

			double dSlope_dN = (beta*Kpen(1,1)*(-(Kpen(1,1)*mu*sDrop) + c*muR))/pow(2*mu*n - muR*n + Kpen(1,1)*sDrop - 2*c, 2);

			if (sDrop>1.001*sMax) {

				dD_ds(1)  = -((mu*n - muR*n - c)*(Kpen(1,1) + 2.0*slope)* sign) / pow(muR*n + Kpen(1,1)*abs(s(1)), 2);

				dD_ds(0)  = ((c*muR + Kpen(1,1)*(mu-muR)*abs(s(1)))*(Kpen(1,1) + 2*slope) + (muR*n + Kpen(1,1)*abs(s(1)))*(2*mu*n - muR*n - 2*c + Kpen(1,1)*abs(s(1)))*dSlope_dN)/(Kpen(1,1)*pow(muR*n + Kpen(1,1)*abs(s(1)),2));
			    dD_ds(0)*= Kpen(0,0);

			}   else {
				dD_ds(0)  = ((muR*n + Kpen(1,1)*abs(s(1)))*(2*mu + 2000.*mu*beta - (2000.*Kpen(1,1)*(-2*mu + muR)* abs(s(1))*((-mu)*n + c)*beta)/pow(-2*mu*n + muR*n + 
                            2*c, 2.0) - (2000.*Kpen(1,1)*mu*abs(s(1))*beta)/(-2*mu*n + muR*n + 2*c)) - muR*(2*mu*n + 2*Kpen(1,1)*abs(s(1)) - 2*c + 2000.*(mu*n - c)*beta + 
                            (2000.*Kpen(1,1)*abs(s(1))*((-mu)*n + c)*beta)/(-2*mu*n + muR*n + 2*c)))/(2*pow(muR*n + Kpen(1,1)*abs(s(1)), 2.0));
			    dD_ds(0)*= Kpen(0,0);
				
				dD_ds(1)  = (Kpen(1,1)*(muR*muR*n*n + mu*n*(c*(-4.-4000.*beta) + muR*n*(-3.-2000.*beta)) + mu*mu*n*n*(2. + 2000.*beta) + c*c*(2.+2000.*beta) + 
                            c*muR*n*(3. + 2000.*beta)))/(pow(muR*n + Kpen(1,1)*abs(s(1)), 2.0)*(-2.*mu*n + muR*n + 2.*c));
	    }  

		*/

		double sign;
		if (s(1)>=0.0)
			sign = 1.0;
		else
			sign = -1.0;

	    if (abs(s(1))<=sMax && Dtrial>0.0) {

			dD_ds(0)  = -((pow((muR*n + Kpen(1,1)*abs(s(1)))/(Gc*muR*n + Gc*VMax),pow(-1 + Gc,-1))*(k*muR*abs(s(1)) - muR*VMax + (muR*n + Kpen(1,1)*abs(s(1)))* dVMax_dn))
				         /(Gc*(muR*n + Kpen(1,1)*abs(s(1)))*(muR*n + VMax))); 
			
			dD_ds(0)*= Kpen(0,0);

	    	dD_ds(1)  = sign*  (Kpen(1,1)*pow((muR*n + Kpen(1,1)*abs(s(1)))/(Gc*muR*n + Gc*VMax),pow(-1 + Gc,-1)))/(Gc*(muR*n + Kpen(1,1)*abs(s(1))));
	    } else {

			double dSlope_dN = (Kpen(1,1)*zeta*(muR*(-1 + Gc)*VMax + (Kpen(1,1)*sDrop + muR*(n - n*Gc))*dVMax_dn))/pow(-(Kpen(1,1)*sDrop) + muR*n*(-1 + Gc) + Gc*VMax, 2.0);

			if (sDrop<0.0) {

				dD_ds(1)  = -((mu*n - muR*n - c)*(Kpen(1,1) + 2.0*slope)* sign) / pow(muR*n + Kpen(1,1)*abs(s(1)), 2);

				dD_ds(0)  = ((c*muR + Kpen(1,1)*(mu-muR)*abs(s(1)))*(Kpen(1,1) + 2*slope) + (muR*n + Kpen(1,1)*abs(s(1)))*(2*mu*n - muR*n - 2*c + Kpen(1,1)*abs(s(1)))*dSlope_dN)/(Kpen(1,1)*pow(muR*n + Kpen(1,1)*abs(s(1)),2));
			    
				dD_ds(0)*= Kpen(0,0);

			}   else {

				dD_ds(0)  = (-((Kpen(1,1) + Gc*slope)*(Kpen(1,1)*muR*abs(s(1)) - muR*VMax + (muR*n + Kpen(1,1)*abs(s(1)))* dVMax_dn)) + 
                            (muR*n + Kpen(1,1)*abs(s(1)))*(Kpen(1,1)*abs(s(1)) + muR*(n - n*Gc) - Gc*VMax)* dSlope_dN)/(Kpen(1,1)*pow(muR*n + Kpen(1,1)*abs(s(1)),2.0));
			    
				dD_ds(0)*= Kpen(0,0);
				
								
				dD_ds(1)  = sign*((Kpen(1,1) + Gc*slope)*(muR*n + VMax))/pow(muR*n + Kpen(1,1)*abs(s(1)), 2.0);
	    }  

		}
		
	}

	if (x<xCommitted)
		x=xCommitted;

	if (Dtrial > D)
		return Dtrial;
	else
		return D;

}


int 
CohesiveSurface::setTrialDisplacement(const Vector& s) {

	if (elasticSolution) {

		sigma = Kpen * s;
		Kalg = Kpen;

	} else {

    x = xCommitted;
	Dnplus1 = this->evolveDamage(s);
	Matrix H = Kpen;
	H(0,0) *= 1.0 - this->heaviside(s(0));
	H(1,1) *= 1.000;
	Vector sigmaTrial = H*(s-s_di);
	
	this->updatePlasticStrain(sigmaTrial);

	sigma = ((1.0 - Dnplus1) * (Kpen*s)) + (Dnplus1 * sigmaTrial);	

	// update tangent matrix
	Matrix Ht(2, 2);
	if (s(0) <= -1.0*DBL_EPSILON) {
		if (deltaLambda > DBL_EPSILON) {
			  Ht(0, 0) = Kpen(0, 0);
			  Ht(1, 0) = -1.0*muR*Kpen(0,0) *this->signTau(sigma(1));
		}
		else {
			Ht = Kpen;
			Ht(1,1) *= 1.000;
		}
	}

	Kalg = (1.0 - Dnplus1)*Kpen + Dnplus1*Ht;
	
	if (flagAlg) {
	  if (Dnplus1>D) {
		Kalg -= Kpen* (s % dD_ds);
		Kalg += sigmaTrial % dD_ds;
	  }

	}
	}

	return 0;

}

double 
CohesiveSurface::getMode(Vector& s) {
	// not relevant
	return xCommitted;
}

Vector 
CohesiveSurface::getSigma() {
	return sigma;
}

Matrix 
CohesiveSurface::getElasticTangent() {
	return Kpen;
}

Matrix 
CohesiveSurface::getAlgorithmicTangent() {
	if (elasticSolution) {
		return Kpen;
	} else {
		return Kalg;
	}
}

Matrix
CohesiveSurface::getCommittedTangent() {
	return KCommitted;
}
	
	
void
CohesiveSurface::commit() {
	if (!elasticSolution) {
		D = Dnplus1;
		xCommitted = x;
		s_di = s_di_nplus1;
		KCommitted = Kalg;
		k += deltaLambda;
 
	}

	return;
}

CohesiveSurface* 
CohesiveSurface::getCopy() {
	CohesiveSurface* theCopy;
	theCopy = new CohesiveSurface(Kpen, c, mu, muR, Gc, GfII, elasticSolution);	
	return theCopy;
}