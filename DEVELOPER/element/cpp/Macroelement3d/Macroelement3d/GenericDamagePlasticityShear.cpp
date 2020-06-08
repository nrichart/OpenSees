/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

/*
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/CohesiveSurface.cpp
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/


#include "GenericDamagePlasticityShear.h"
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

GenericDamagePlasticityShear::GenericDamagePlasticityShear(void)
	:Kalg(2, 2), Kpen(2, 2), sigma(2), s_di(2), s_di_nplus1(2), eps(2) {
}

GenericDamagePlasticityShear::GenericDamagePlasticityShear(int tag, double G, int modelType, Vector matProp, double Gc, double dropDrift, double alpha, bool elasticSolution, double initialDamage)
	: NDMaterial(tag, 0), Kalg(2,2), Kpen(2,2), KCommitted(2,2), sigma(2), s_di(2), dD_ds(2), eps(2), modelType(modelType), props(matProp), Gc(Gc), dropDrift(dropDrift),
	  flagAlg(false), deltaLambda(0.0), s_di_nplus1(2), D(initialDamage), alpha(alpha),
	  elasticSolution(elasticSolution), k(0.), x(0.), xCommitted(0.) 
{			
	Dnplus1 = D;
	Kpen(0, 0) = Kpen(1, 1) = G;
	Kalg = Kpen;
	KCommitted = Kpen;
	
}

GenericDamagePlasticityShear::~GenericDamagePlasticityShear(void) {
}

double 
GenericDamagePlasticityShear::signTau(double tau) {
	if (tau == 0.0)
		return 0.0;
	else if (tau > 0.0)
		return 1.0;
	else
		return -1.0;
}

double
GenericDamagePlasticityShear::heaviside(double arg) {
	if (arg > 0.0)
		return 1.0;
	else
		return 0.0;
}

double 
GenericDamagePlasticityShear::macauley(double arg) {
	if (arg > 0.0)
		return arg;
	else
		return 0.0;
}


// general methods for force capacity
double 
GenericDamagePlasticityShear::Vmax(const Vector& s) {
	// params are H0, axial load ratio
	// axial force (negative in compression)
	double N = min(0.0, Kpen(0, 0)*s(0));
	if (modelType == 1) {
		// Mohr-Coulomb type criterion
		// props = [L, t, c, mu]
		return props(0)*props(1)*props(2) - props(3)*N;
	}
	else if (modelType==2) {
		// Turnsek-Cacovic type criterion
		// props = [L, t, b, ft]
		double Ltft = props(0)*props(1)*props(3);
		return Ltft / props(2) * sqrt(1 - (N / Ltft));
	}
	else if (modelType == 3) {
		// Mohr-Coulomb on the compressed section + condition on the brick cracking 0.065fm*lc*t
		// props = [L, t, c, mu, 0.065fm] + params(0) = H0
		double V1 = -N / 2.0*(3.*props(2)*props(0)*props(1) - 2.*props(3)*N );
		V1 /= 3.*props(2)*props(1)*params(0) - N;
		return V1;
	}
	else {
		return 0.0;
	}
}

double
GenericDamagePlasticityShear::dVmax_dN(const Vector& s) {
	// in s: sAxial, sShear, L'/L
	// axial force (negative in compression)
	double N = min(0.0, Kpen(0, 0)*s(0));
	if (modelType == 1) {
		// Mohr-Coulomb type criterion
		// props = [L, t, c, mu]
		if (s(0) > 0.0) {
			return 0.0;
		}
		else {
			return -1.0*props(3);
		}	
	}
	else if (modelType == 2) {
		// Turnsek-Cacovic type criterion
		// props = [L, t, b, ft]
		double Ltft = props(0)*props(1)*props(3);
		if (s(0) > 0.0) {
			return 0.0;
		}
		else {
			return -Kpen(0, 0)/(2.*props(2)*sqrt(1.-N/(Ltft)));
		}
	}
	else if (modelType == 3) {
		// Mohr-Coulomb on the compressed section + condition on the brick cracking 0.065fm*lc*t
		// props = [L, t, c, mu, 0.065fm] + params(0) = H0
		double Ltft = props(0)*props(1)*props(3);
		if (s(0) > 0.0) {
			return 0.0;
		}
		else {
			double V1 = N / 2.0*(3.*props(2)*props(0)*props(1) + 2.*props(3)*N);
			V1 /= 3.*props(2)*props(1)*params(0) + N;
			double dV1_dN = -2.*N*N*props(3) + 3.*props(2)*params(0)*props(1)*(-3.*props(2)*props(0)*props(1)+4.*N*props(2));
			dV1_dN /= 2.*pow(N - 3.*props(2)*params(0)*props(1), 2.);
			return dV1_dN;
		}
	}
	else {
		return 0.0;
	}
}


void 
GenericDamagePlasticityShear::updatePlasticStrain(Vector& sigmaTrial, double maxShearCapacity) {

	// return mapping for constant muR, no iterations needed
	double yieldFunction = abs(sigmaTrial(1)) - alpha*maxShearCapacity;
	deltaLambda = 0.0;
	s_di_nplus1 = s_di;

	if (yieldFunction > DBL_EPSILON) {
		deltaLambda = 1.0 / (Kpen(1, 1)) * yieldFunction;
		s_di_nplus1(1) += deltaLambda * this->signTau(sigmaTrial(1));
		double sign = this->signTau(sigmaTrial(1));
		sigmaTrial(1) -= deltaLambda * Kpen(1, 1) * sign;
	}
	return;
}


double 
GenericDamagePlasticityShear::evolveDamage(const Vector& s) {
	
	double Dtrial;
	
	// force-capacity model 
	double VMax = this->Vmax(s);          
	double dVMax_dn = this->dVmax_dN(s);       

	double s0 = alpha*VMax/Kpen(1,1);
	double sMax = s0 + (1.0-alpha)*VMax/Kpen(1, 1) * Gc;
	double sDrop = dropDrift;      

	double zeta = 0.20;               // capacity drop percentage
	double slope = zeta *VMax/ (std::max(sDrop, 1.001*sMax) - sMax);

	if (abs(s(1))<s0) {
		x = 0.0;
		Dtrial = 0.0;
	} else {
		if (abs(s(1))<sMax) {
			x = (abs(s(1)) - s0)/(sMax - s0);
			Dtrial = (1.0 - 1.0/Gc)*pow(x, 1.0/(Gc - 1.0));
		} else {
			x = (abs(s(1)) - s0)/(sMax - s0);  // not used but calculated for output
			Dtrial = (Kpen(1, 1)*abs(s(1)) + abs(s(1))*slope - sMax*slope - VMax) / (Kpen(1, 1)*abs(s(1)) - alpha*VMax);
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

		double sign;
		if (s(1)>=0.0)
			sign = 1.0;
		else
			sign = -1.0;

	    if (abs(s(1))<=sMax && Dtrial>0.0) {

			double dD_dx = pow(x, -1. + 1. / (Gc - 1.)) / Gc;

			dD_ds(0) = dD_dx * Kpen(1, 1)*abs(s(1)) * dVMax_dn / (Gc*(alpha - 1.)*VMax*VMax);		
			dD_ds(0)*= Kpen(0,0);

			dD_ds(1) = sign*  Kpen(1, 1) / (Gc*(alpha - 1.)*VMax) *dD_dx;

	    } else {

			double dSlope_dN = (Kpen(1, 1)*Kpen(1, 1) *sDrop*zeta * dVMax_dn) / pow(Kpen(1, 1)*sDrop + (Gc*(alpha - 1.0) - alpha)*VMax, 2.0);

			dD_ds(0) = (Kpen(1, 1)*abs(s(1)) + (Gc*(alpha - 1) - alpha)*VMax)*(Kpen(1, 1)*abs(s(1)) - alpha*VMax)* dSlope_dN
			        	+ Kpen(1, 1)*(alpha - 1.)*abs(s(1))*(Kpen(1, 1) + Gc*slope)*dVMax_dn;
			dD_ds(0) /= Kpen(1, 1)*pow(Kpen(1, 1)*abs(s(1)) - alpha*VMax, 2.0) / Kpen(0, 0);

			dD_ds(1) = ((alpha - 1)*(Kpen(1, 1) + Gc*slope)*VMax*sign) / pow(Kpen(1, 1)*abs(s(1)) - alpha*VMax, 2.0);
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
GenericDamagePlasticityShear::setTrialStrainExtraParams(const Vector& s, Vector params) {
	if (elasticSolution) {

		sigma = Kpen * s;
		Kalg = Kpen;

	} else {

    Vector increment = s - eps;
		if (increment.Norm() > DBL_EPSILON) {
			eps = s;
		}
		else {
			return 0;
		}

	this->params = params;
    x = xCommitted;
	Dnplus1 = this->evolveDamage(s);
	Matrix H = Kpen;
	H(0,0) *= 1.0 - this->heaviside(s(0));
	Vector sigmaTrial = H*(s-s_di);
	
	double VMax = this->Vmax(s);
	this->updatePlasticStrain(sigmaTrial, VMax);

	sigma = ((1.0 - Dnplus1) * (Kpen*s)) + (Dnplus1 * sigmaTrial);	

	// update tangent matrix
	Matrix Ht(2, 2);
	if (s(0) <= -1.0*DBL_EPSILON) {
		if (deltaLambda > DBL_EPSILON) {
			  Ht(0, 0) = Kpen(0, 0);
			  Ht(1, 0) = this->dVmax_dN(s) *Kpen(0,0) *this->signTau(sigma(1));
		}
		else {
			Ht = Kpen;
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


int
GenericDamagePlasticityShear::setTrialStrain(const Vector &s) {
	Vector empty(0);
	return this->setTrialStrainExtraParams(s, empty);
}


int 
GenericDamagePlasticityShear::setTrialStrain(const Vector &s, const Vector &rate) {
	return this->setTrialStrain(s);
}



double 
GenericDamagePlasticityShear::getMode() {
	return xCommitted;
}

const Vector &
GenericDamagePlasticityShear::getStress() {
	return sigma;
}

const Matrix & 
GenericDamagePlasticityShear::getInitialTangent() {
	return Kpen;
}


const Matrix&
GenericDamagePlasticityShear::getTangent()
{
	// provide a (fake) small stiffness if the zero tangent 
	if (abs(Kalg(1, 1))<DBL_EPSILON) {
		Matrix el(2, 2);
		Kalg(1, 1) = 0.00001*Kpen(1, 1);
	}
	return Kalg;
}


const Vector &
GenericDamagePlasticityShear::getStrain(void) {
	return eps;
}



const Matrix &
GenericDamagePlasticityShear::getCommittedTangent() {
	return KCommitted;
}
	

int
GenericDamagePlasticityShear::commitState() {
	if (!elasticSolution) {
		D = Dnplus1;
		xCommitted = x;
		s_di = s_di_nplus1;
		KCommitted = Kalg;
		k += deltaLambda;
	}
	return 0;
}

int
GenericDamagePlasticityShear::revertToLastCommit() {
	if (!elasticSolution) {
		Dnplus1 = D;
		x = xCommitted;
		s_di_nplus1 = s_di;
		Kalg = KCommitted;
		deltaLambda = 0.0;
	}
	return 0;
}

int
GenericDamagePlasticityShear::revertToStart() {
	if (!elasticSolution) {
		Dnplus1 = D = 0.0;
		x = xCommitted = 0.0;
		s_di_nplus1 = s_di = 0.*s_di;
		Kalg = KCommitted = Kpen;
		deltaLambda = 0.0;
	}
	return 0;
}




NDMaterial*
GenericDamagePlasticityShear::getCopy() {
	GenericDamagePlasticityShear* theCopy;
	theCopy = new GenericDamagePlasticityShear(this->getTag(), Kpen(1,1), modelType, props,  Gc,  dropDrift, alpha, elasticSolution, D);
	return theCopy;
}



NDMaterial *
GenericDamagePlasticityShear::getCopy(const char *type) {
	return this->getCopy();
}


const char*
GenericDamagePlasticityShear::getType(void) const
{
	return "BeamFiber";
}

int
GenericDamagePlasticityShear::getOrder(void) const
{
	return 2;
}


void 
GenericDamagePlasticityShear::Print(OPS_Stream &s, int flag) {
	s << endln;
	s << "DamageShearInterface : nothing to declare ";
	s << endln;
}


int
GenericDamagePlasticityShear::sendSelf(int commitTag, Channel &theChannel)
{
	opserr << "WARNING: DamageShearInterface has no method sendSelf implemented" << endln;
	return -1;
}

int
GenericDamagePlasticityShear::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	opserr << "WARNING: DamageShearInterface has no method recvSelf implemented" << endln;
	return -1;
}