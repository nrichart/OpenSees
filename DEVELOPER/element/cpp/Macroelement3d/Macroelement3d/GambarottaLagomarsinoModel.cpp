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
                                                                        
// $Revision: 1.7 $
// $Date: 2009/03/23 23:17:04 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/GambarottaLagomarsinoModel.cpp,v $
                                                                        
// Written: Francesco Vanin
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) GambarottaLagomarsinoModel.C, revA"

#include <elementAPI.h>
#include "GambarottaLagomarsinoModel.h"


#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Information.h>
#include <Channel.h>
#include <Message.h>
#include <NDMaterial.h>
#include <MaterialResponse.h>



#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numGambarottaLagomarsinoModel = 0;
static double tol = 1e-7;
static int maxIter = 100;

OPS_Export void *
OPS_GambarottaLagomarsinoModel()
{
  // print out some KUDO's
  if (numGambarottaLagomarsinoModel == 0) {
    opserr << "GambarottaLagomarsinoModel - Loaded from external library\n";
    numGambarottaLagomarsinoModel =1;
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[9];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag" << endln;
    return 0;
  }

  numData = 9;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid input parameters" << endln;
    return 0;	
  }

  bool _elastic= false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
	  const char* type = OPS_GetString();
	  if (strcmp(type, "-elastic") == 0 || strcmp(type, "-Elastic") == 0) {
		  _elastic = true;
	  }
  }
  
  // create a new material
  theMaterial = new GambarottaLagomarsinoModel(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],dData[7],dData[8], _elastic);

  if (theMaterial == 0) {
    opserr << "WARNING could not create NDMaterial of type GambarottaLagomarsinoModel\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

//0.000001
// full constructor

GambarottaLagomarsinoModel::GambarottaLagomarsinoModel(int tag, double _E, double _G, double _c, double _mu, double _ct, double _beta,  double _L, double _t, double _h, bool _elasticSolution)
						  :NDMaterial(tag, 0), stress(2), stressCommitted(2), u(2), uCommitted(2), Kpen(2,2), K(2,2), KCommitted(2,2), elasticSolution(_elasticSolution),
						  c(_c), mu(_mu), ct(_ct), beta(_beta), alpha(0.0000000), alphaCommitted(0.0000000), deltaLambda(0.0), deltaAlpha(0.0), s(0.0), sCommitted(0.0),
						  L(_L), t(_t), h(_h)
{ 
   Kpen(0,0) = _E*L*t/h;
   Kpen(1,1) = _G*L*t/h;

   K = Kpen;
   KCommitted = Kpen;
}

GambarottaLagomarsinoModel::GambarottaLagomarsinoModel()
:NDMaterial(0, 0), stress(2), stressCommitted(2), u(2), uCommitted(2), Kpen(2,2), K(2,2), KCommitted(2,2), elasticSolution(false),
			c(0.0), mu(0.0), ct(0.0), beta(0.0), alpha(0.0), alphaCommitted(0.0), deltaLambda(0.0), deltaAlpha(0.0), s(0.0), sCommitted(0.0),
			L(0.0), t(0.0), h(0.0)
{ 
}

GambarottaLagomarsinoModel::~GambarottaLagomarsinoModel() {
}


//internal methods for toughness function
double 
GambarottaLagomarsinoModel::toughnessFunction(double _alpha) {
	
	if (_alpha<=1.0) {
		return 0.5*c*c*ct*h*L*t*_alpha;
	} else {
		return 0.5*c*c*ct*h*L*t*pow(_alpha, -beta);
	}
}

double 
GambarottaLagomarsinoModel::derToughnessFunction(double _alpha) {

	if (_alpha<=1.0) {
		return 0.5*c*c*ct*h*L*t;
	} else {
		return 0.5*c*c*ct*h*L*t*(-beta)*pow(_alpha, -beta-1.0);
	}
}

double 
GambarottaLagomarsinoModel::sign(double _V) {

	if (_V>=0.0) {
		return 1.0;
	} else {
		return -1.0;
	}
}


// standard methods
int 
GambarottaLagomarsinoModel::setTrialStrain(const Vector &strain, const Vector &rate) {
	return this->setTrialStrain(strain);
}

int
GambarottaLagomarsinoModel::setTrialStrain(const Vector &strain) { 

	Vector increment = strain - u;
	if (increment.Norm() > DBL_EPSILON) {
		u = strain;
		stress(0) = Kpen(0,0)*u(0);
		stress(1) = Kpen(1,1)*(u(1) - sCommitted);
		Vector trialStress = stress;
		K = Kpen;

		if (elasticSolution)
			return 0;
		
		// initialise values of state variables to last committed
		alpha = alphaCommitted;
		s = sCommitted;
		deltaLambda = 0.0;
		deltaAlpha = 0.0;

		double csi = 0.;
		if (alphaCommitted>0.0) csi = sCommitted / alphaCommitted;

		double N = stress(0);
		double Vf;
		// kill it here if in tension? Yes
		if (N>-DBL_EPSILON) {
			stress(1) = 0.0;
			K.Zero();
			//opserr << "GambarottaLagomarsinoModel: killed in tension\n";
			return 0;
		}

		Vector residuals(2);
		// damage function	
		residuals(0) = 0.5*L*t / ( ct*h)*pow(csi,2) - this->toughnessFunction(alpha);
		Vf = Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ;
		residuals(1) = this->sign(Vf)* (Vf) + mu*N; 
		
		//opserr<< "RESIDUALS 0 =" << residuals;
		
		if (residuals(0)<=tol && residuals(1)<=tol ) {
			// no evolution of damage parameters, elastic step 
			
		} else {
			// evolution of damage parameters

			if (residuals(0)<=-tol && alpha>DBL_EPSILON) {   // -tol because otherwise we enter here even if there should be evolution of damage: for the first iteration tha damage function is equal to the lat converged value
				// try sliding without damage evolution (unloading/reloading)
				//opserr<<"sliding without damage evolution (unloading/reloading)\n";
				//opserr << "from " << uCommitted(1) << " to " << u(1) << ". csi0 = " << csi << ", sCommitted " << sCommitted << endln;
				
				deltaLambda = residuals(1) / (Kpen(1,1) + L*t/(ct*h*alpha));
				s += deltaLambda * this->sign(Vf);
								
				// stress update
				stress(1) = Kpen(1,1)*(u(1) - s);

				// algorithmic operator
				double df_dLambda = -L*t/(ct*h*alpha);

				Vector df_dsigma(2);
				df_dsigma(0) = mu;
				df_dsigma(1) = this->sign(Vf);

				Vector dg_dsigma(2);
				dg_dsigma(1) = this->sign(Vf);

				K.Zero();
				K.addMatrixTripleProduct(0.0, Kpen, (dg_dsigma % df_dsigma), Kpen, -1.0);
				K /= (df_dsigma^(Kpen*dg_dsigma)) - df_dLambda;
				K += Kpen;

				//opserr<< "alpha=" << alpha << ", csi=" << csi << ". s=" << s << ", deltaLambda = " << deltaLambda <<  "\n";

			} 
			else {
				// sliding with evolution of damage (new loading)
				// Newton Raphson scheme to solve the nonlinear system

				
				int iter = 0;
			
				Matrix inverseDerivative(2,2);
				Matrix M(2,2);
				Matrix preM(2,2);
				preM(0,0) = 1.;
				preM(1,1) = this->sign(Vf);

				Vector updates(2);
				double G = Kpen(1,1) / (L*t/h);

				//opserr << "from " << uCommitted(1) << " to " << u(1) << ". csi0 = " << csi << ", sCommitted " << sCommitted << endln;

				while ( (residuals.Norm() >= tol || (alpha<DBL_EPSILON && iter<=1)) && iter<=maxIter) {
					iter++;
					// jacobian (-jacobian, in fact)
					
					/*
					if (alpha>tol) {

					M(0,0) = L*t*(s*s)/(ct*h*pow(alpha,3)) + this->derToughnessFunction(alpha);
					M(0,1) = -L*t*s/(ct*h*pow(alpha,2))   * this->sign(Vf);
					M(1,0) = M(0,1);
					M(1,1) = Kpen(1,1) + L*t/(ct*h*alpha);

					} else {
						M(0,0) =  this->derToughnessFunction(alpha);
						M(0,1) = 0.0;
						M(1,0) = M(0,1);
						M(1,1) = Kpen(1,1);
					}
					*/

					M(0,0) = this->derToughnessFunction(alpha);
					M(0,1) = -L*t/(ct*h)*csi;
					M(1,0) = G*L*t/h*csi * this->sign(Vf);
					M(1,1) = this->sign(Vf)*(Kpen(1,1)*alpha + L*t/(ct*h));


					inverseDerivative(0,0) = M(1,1);
					inverseDerivative(0,1) = -M(0,1);
					inverseDerivative(1,0) = -M(1,0);
					inverseDerivative(1,1) = M(0,0);

					inverseDerivative /= M(0,0)*M(1,1) - M(0,1)*M(1,0);

					//updates = (preM*inverseDerivative)*residuals;
					updates = (inverseDerivative)*residuals;
				
					deltaAlpha  += updates(0);
					if (deltaAlpha<DBL_EPSILON)      deltaAlpha=0.0;
					alpha += updates(0);
					csi += updates(1);

					//deltaLambda += updates(1) * this->sign(Vf);
					deltaLambda = (csi*alpha - sCommitted) / this->sign(Vf);
					//s += updates(1);
					//s += updates(1)*alpha;
					s = csi*alpha;

					//opserr<< "updated: alpha=" << alpha << ", csi=" << csi << ". s=" << s << ", deltaLambda = " << deltaLambda <<  "\n";
									

					/*
					if (alpha>DBL_EPSILON)
						residuals(0) = 0.5*L*t / ( ct*h*pow(alpha,2) )*pow(s,2) - this->toughnessFunction(alpha);
					else
						residuals(0) = this->toughnessFunction(0.0);

 
					if (alpha>DBL_EPSILON)
						Vf = Kpen(1,1)*(u(1) - s)  - s*L*t/(ct*h*alpha);			
					else
						Vf = Kpen(1,1)*(u(1) - s);  

					residuals(1) = abs(Vf) + mu*N;  

					//preM(1,1) = this->sign(Vf);
					opserr<< "RESIDUALS old =" << residuals;
					*/

					

					residuals(0) = 0.5*L*t / ( ct*h)*pow(csi,2) - this->toughnessFunction(alpha);
					residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

					//opserr<< "RESIDUALS new =" << residuals;




					if (alpha<alphaCommitted)   {  // if we enter here: there should have been the update only for deltaLambda and we are here because the damage function is close to zero

						alpha=alphaCommitted;
						deltaAlpha = 0.0;
						csi = sCommitted / alphaCommitted;					
						s = sCommitted;

						residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

						csi +=  residuals(1) / ( this->sign(Vf)*(Kpen(1,1)*alpha + L*t/(ct*h)) );
						deltaLambda = (csi*alpha - sCommitted) / this->sign(Vf);
						s = csi*alpha;

						residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

						residuals(0) = 0.5*L*t / ( ct*h*pow(alpha,2) )*pow(s,2) - this->toughnessFunction(alpha);
						if (residuals(0)<0.0)
							residuals(0) = 0.0;
						
					}
					
				}  // end while (NR iterations)

				// stress update
				stress(1) = Kpen(1,1)*(u(1) - s);

				if (iter>maxIter) {
					//opserr << "GambarottaLagomarsinoModel: did not converge the internal Newton-Raphson scheme. N=" <<N << ". displ=" << s <<". Residuals " << residuals;
					// maybe insert here fatal error
					K=Kpen*0;
					stress(1) = -mu*N* this->sign(Vf);
					s = sCommitted;
					alpha = alphaCommitted;
					deltaLambda = 0.0;
					deltaAlpha = 0.0;
					/*
					
					opserr << "NOT converged in " << iter << " iterations. N=" <<N << ". displ=" << u(1) <<" from " << uCommitted(1) << ". Residuals " << residuals << "\n";

					// repeat iterations, verbose ------------------------------------------------------------------------------------
					

					csi = 0.;
					if (alphaCommitted>0.0) csi = sCommitted / alphaCommitted;

	
					residuals(0) = 0.5*L*t / ( ct*h)*pow(csi,2) - this->toughnessFunction(alpha);
					residuals(1) = abs( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N; 
					Vf = Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ;



				iter = 0;
			
				Matrix inverseDerivative(2,2);
				Matrix M(2,2);
				Matrix preM(2,2);
				preM(0,0) = 1.;
				preM(1,1) = this->sign(Vf);

				opserr << "from " << uCommitted(1) << " to " << u(1) << ". csi0 = " << csi << ", sCommitted " << sCommitted << ", alphaCommitted " << alphaCommitted << endln;



				while (residuals.Norm() >= tol && iter<=maxIter) {

					opserr<< "RESIDUALS new =" << residuals;


					iter++;
					// jacobian (-jacobian, in fact)
					
					M(0,0) = this->derToughnessFunction(alpha);
					M(0,1) = -L*t/(ct*h)*csi;
					M(1,0) = G*L*t/h*csi * this->sign(Vf);
					M(1,1) = this->sign(Vf)*(Kpen(1,1)*alpha + L*t/(ct*h));


					inverseDerivative(0,0) = M(1,1);
					inverseDerivative(0,1) = -M(0,1);
					inverseDerivative(1,0) = -M(1,0);
					inverseDerivative(1,1) = M(0,0);

					inverseDerivative /= M(0,0)*M(1,1) - M(0,1)*M(1,0);

					//updates = (preM*inverseDerivative)*residuals;
					updates = (inverseDerivative)*residuals;
				
					deltaAlpha  += updates(0);
					if (deltaAlpha<DBL_EPSILON)      deltaAlpha=0.0;
					alpha += updates(0);
					csi += updates(1);

					//deltaLambda += updates(1) * this->sign(Vf);
					deltaLambda = (csi*alpha - sCommitted) / this->sign(Vf);
					//s += updates(1);
					//s += updates(1)*alpha;
					s = csi*alpha;

					opserr<< "updated: alpha=" << alpha << ", csi=" << csi << ". s=" << s << ", deltaLambda = " << deltaLambda <<  "\n";
									
					residuals(0) = 0.5*L*t / ( ct*h)*pow(csi,2) - this->toughnessFunction(alpha);
					//residuals(1) = abs( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N; 
					residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

					
					

					if (alpha<alphaCommitted)   {  // if we enter here: there should have been the update only for deltaLambda and we are here because the damage function is close to zero
						alpha=alphaCommitted;
						deltaAlpha = 0.0;
						
						s = sCommitted;
						
						
					if (alpha>DBL_EPSILON)
						Vf = Kpen(1,1)*(u(1) - s)  - s*L*t/(ct*h*alpha);			
					else
						Vf = Kpen(1,1)*(u(1) - s);  

					residuals(1) = abs(Vf) + mu*N;  

						deltaLambda = residuals(1) / (Kpen(1,1) + L*t/(ct*h*alpha));
				        s += deltaLambda * this->sign(Vf);

						if (alpha>DBL_EPSILON)
							residuals(1) = abs( Kpen(1,1)*(u(1) - s)  - s*L*t/(ct*h*alpha) ) + mu*N;    
						else
							residuals(1) = abs( Kpen(1,1)*(u(1) - s)) + mu*N; 

						//residuals(0) = 0.0; // fake, just to exit the while loop.
						residuals(0) = 0.5*L*t / ( ct*h*pow(alpha,2) )*pow(s,2) - this->toughnessFunction(alpha);

						//opserr << "Corrected in the loop (iter " << iter << "). N/DBL_EPS = " << N/DBL_EPSILON << ". Residuals: " << residuals;

						if (residuals(0)<-tol)
							residuals(0) = 0.0;

						
						
					}
					
				}  // end while (NR iterations)



				*/
				



					return 0;

				}	else {

					//opserr << "converged in " << iter << " iterations: alphaCommitted " << alphaCommitted << ", alpha " << alpha << ", sCommitted " << sCommitted << ", s " << s << "\n";

				// tangent operator update
				double hh;				
				hh = -(L*t*alpha*s*this->sign(Vf)) / (L*t*s*s + ct*h*pow(alpha,3)* this->derToughnessFunction(alpha));

				if (deltaAlpha<= DBL_EPSILON)   hh = 0.0;

				double df_dalpha = (L*t*s*this->sign(Vf)) / (ct*h*pow(alpha,2));
				double df_dLambda = -L*t/(ct*h*alpha);

				Vector df_dsigma(2);
				df_dsigma(0) = mu;
				df_dsigma(1) = this->sign(Vf);

				Vector dg_dsigma(2);
				dg_dsigma(1) = this->sign(Vf);

				K.Zero();
				K.addMatrixTripleProduct(0.0, Kpen, (dg_dsigma % df_dsigma), Kpen, -1.0);
				K /= df_dalpha*hh + (df_dsigma^(Kpen*dg_dsigma)) - df_dLambda;
				K += Kpen;

				//opserr << "forces = " << stress << "tangent = " << K;
				//opserr << "elastic tangent = " << Kpen << endln;

				}


			}  // end second case

		}      // end evolution of damage parameters (2 cases)

	}          // end set new strain
    return 0;
}

int
GambarottaLagomarsinoModel::setTrialStrainIncr (const Vector &strain)
{
	opserr << "WARNING: GambarottaLagomarsinoModel has no method setTrialStrainIncr implemented" << endln; 
	return 0;
}

int
GambarottaLagomarsinoModel::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return this->setTrialStrainIncr(strain);
}

const Matrix& 
GambarottaLagomarsinoModel::getCommittedTangent( )  {
	return KCommitted;
}
	
const Matrix& 
GambarottaLagomarsinoModel::getTangent( ) 
{
	//opserr << "get tangent=" << Kpen;
	if (alpha>1.0 && alphaCommitted<1.0) {
		//opserr << "trying to pass from alpha " << alphaCommitted << " to " << alpha << endln;
		return KCommitted;
	} else {
		if (abs(K(0,1))>-DBL_EPSILON && abs(K(1,1))>-DBL_EPSILON)		
			return K;
		else
			return Kpen;
	}
}

const Matrix&
GambarottaLagomarsinoModel::getInitialTangent (void)
{ 
	return Kpen;
}

const Vector&
GambarottaLagomarsinoModel::getStress (void)
{
	return stress;
}

const Vector&
GambarottaLagomarsinoModel::getStrain (void)
{
  return u;
}

int
GambarottaLagomarsinoModel::commitState (void)
{
	// implement all commit variables
	alphaCommitted = alpha;
	sCommitted = s;
	KCommitted = K;

	uCommitted = u;
	stressCommitted = stress;

	return 0;
}

int
GambarottaLagomarsinoModel::revertToLastCommit (void)     
{ 
	alpha = alphaCommitted;
	s = sCommitted;
	K = KCommitted;

	u = uCommitted;
	stress = stressCommitted;

  return 0;
}

int
GambarottaLagomarsinoModel::revertToStart (void)
{ 
	alpha = 0;
	s = 0;
	K = Kpen;

	alphaCommitted = 0;
	sCommitted = 0;
	KCommitted = Kpen;

	u.Zero();
	stress.Zero();

  return 0;
}

double 
GambarottaLagomarsinoModel::getMode(void) {
	return alphaCommitted;
}

NDMaterial*
GambarottaLagomarsinoModel::getCopy (void)
{
	GambarottaLagomarsinoModel *theCopy = new GambarottaLagomarsinoModel(this->getTag(), Kpen(0,0)/(L*t/h), Kpen(1,1)/(L*t/h), c, mu, ct, beta, L, t, h, elasticSolution);

  return theCopy;
}

NDMaterial*
GambarottaLagomarsinoModel::getCopy (const char *type) {
  return this->getCopy();
}

const char*
GambarottaLagomarsinoModel::getType (void) const
{
  return "BeamFiber";
}

int
GambarottaLagomarsinoModel::getOrder (void) const
{
  return 2;
}

//print out material data
void 
GambarottaLagomarsinoModel::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "GambarottaLagomarsinoModel : " ; 
  s << this->getType( ) << endln ;
  s << endln ;
}

int 
GambarottaLagomarsinoModel::sendSelf(int commitTag, Channel &theChannel) 
{  return -1;  }

int 
GambarottaLagomarsinoModel::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{  return -1;  }


// methods to overwrite the standard implementation in NDMaterial for asking outputs
//----------------------------------------------------------------------------------

Response*
GambarottaLagomarsinoModel::setResponse (const char **argv, int argc, OPS_Stream &output) {

  Response *theResponse = 0;
  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    output.tag("ResponseType","sigma");
	output.tag("ResponseType","tau");   
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
	output.tag("ResponseType","epsilon");
	output.tag("ResponseType","gamma");
    theResponse =  new MaterialResponse(this, 2, this->getStrain());

  } else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"stiffness") == 0) {
	const Matrix &res = this->getCommittedTangent();
    int size = 4;
	output.tag("ResponseType","K11");
	output.tag("ResponseType","K12");
	output.tag("ResponseType","K21");
	output.tag("ResponseType","K22");
    theResponse =  new MaterialResponse(this, 3, this->getCommittedTangent());
  }

    else if (strcmp(argv[0], "stateVariables") == 0 || strcmp(argv[0], "state") == 0) {
	  double res = this->getMode();
	  int size = 1;
	  output.tag("ResponseType", "mode");
	  theResponse = new MaterialResponse(this, 4, this->getMode());
  }
  

  output.endTag(); // NdMaterialOutput

  return theResponse;
}


int 
GambarottaLagomarsinoModel::getResponse (int responseID, Information &matInfo)
{
  Vector ep_dot(2);
  double mode;

  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getStress());
    
  case 2:
    return matInfo.setVector(this->getStrain());

  case 3:
	  return matInfo.setMatrix(this->getCommittedTangent());

  case 4:
	  return matInfo.setDouble(this->getMode());

  default:
    return -1;
  }
}