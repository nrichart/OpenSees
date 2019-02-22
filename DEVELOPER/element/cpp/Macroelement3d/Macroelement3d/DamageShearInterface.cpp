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
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/DamageShearInterface.cpp,v $
                                                                        
// Written: Francesco Vanin
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) DamageShearInterface.C, revA"

#include <elementAPI.h>
#include "DamageShearInterface.h"


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

static int numDamageShearInterface = 0;


OPS_Export void *
OPS_DamageShearInterface()
{
  // print out some KUDO's
  if (numDamageShearInterface == 0) {
    opserr << "DamageShearInterface - Loaded from external library\n";
    numDamageShearInterface =1;
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[13];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag" << endln;
    return 0;
  }

  numData = 13;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid input parameters" << endln;
    return 0;	
  }

  // check for additional input parameters
  numData = OPS_GetNumRemainingInputArgs();
  double optData[2];
  if (numData>2) numData = 2;
  optData[0] = 0.50;
  optData[1] = 0.80;
  if (numData==2) {
	  OPS_GetDoubleInput(&numData, optData);
  }

  bool linearUnloading = false;
  bool explicitFormulation = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
	  const char* type = OPS_GetString();
	  if (strcmp(type, "-linearUnloading") == 0 || strcmp(type, "-linearUnloading") == 0) {
		  linearUnloading = true;
	  }

	  if (strcmp(type, "-expl") == 0 || strcmp(type, "-explicit")==0 || strcmp(type, "-Explicit") == 0) {
		  explicitFormulation = true;
	  }
  }
  
  // create a new material
  //theMaterial = new DamageShearInterface(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], optData[0], optData[1], linearUnloading, explicitFormulation);

  if (theMaterial == 0) {
    opserr << "WARNING could not create NDMaterial of type DamageShearInterface\n";
    return 0;
  }

  // return the material
  return theMaterial;
}


// full constructor
DamageShearInterface::DamageShearInterface(int tag, double _E, double _G, double c, double mu, double muR, double beta, double dropDrift, bool elasticSolution)
		             :NDMaterial(tag, 0), CommittedStrain(2), epsTrial(2), elasticSolution(elasticSolution)
{ 
	Matrix Kpen(2, 2);
	Kpen(0, 0) = _E;
	Kpen(1, 1) = _G;

	K = Kpen;

	if (mu==muR)
	     cohesiveSurface = new CohesiveSurface(Kpen, c, mu, muR, beta, dropDrift, elasticSolution); 
	else
		cohesiveSurface = new CohesiveSurface(Kpen, c, mu, muR, beta, dropDrift, elasticSolution); 
}

// Constructor for copies
DamageShearInterface::DamageShearInterface(int tag, CohesiveSurface* _copyCohesiveSurface)
	 :NDMaterial(tag, 0), CommittedStrain(2), epsTrial(2), cohesiveSurface(_copyCohesiveSurface) {
		 K = cohesiveSurface->getElasticTangent();
 }


DamageShearInterface::DamageShearInterface()
:NDMaterial(0, 0), CommittedStrain(2), epsTrial(2)
{
	// check if it is necessary to initialise the other members
}

DamageShearInterface::~DamageShearInterface() {
}

int 
DamageShearInterface::setTrialStrain(const Vector &strain, const Vector &rate) {
	return this->setTrialStrain(strain);
}

int
DamageShearInterface::setTrialStrain(const Vector &strain) { 
	//opserr << "setTrialDispl = " << strain;

	Vector increment = strain - epsTrial;
	if (increment.Norm() > DBL_EPSILON) {
		epsTrial = strain;
		cohesiveSurface->setTrialDisplacement(strain);
	}

    return 0;
}

int
DamageShearInterface::setTrialStrainIncr (const Vector &strain)
{
	opserr << "WARNING: DamageShearInterface has no method setTrialStrainIncr implemented" << endln; 
	return 0;
}

int
DamageShearInterface::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return this->setTrialStrainIncr(strain);
}

const Matrix& 
DamageShearInterface::getCommittedTangent( )  {
	K = cohesiveSurface->getCommittedTangent();
	return K;
}
	

const Matrix& 
DamageShearInterface::getTangent( ) 
{
	K = cohesiveSurface->getAlgorithmicTangent();

	if (abs(K(1,1))<DBL_EPSILON) {
		Matrix el(2,2);
		el = this->getInitialTangent();
		K(1,1) = 0.00001*el(1,1);
	}

	return K;
}

const Matrix&
DamageShearInterface::getInitialTangent (void)
{ 
	K = cohesiveSurface->getElasticTangent();
	return K;
}

const Vector&
DamageShearInterface::getStress (void)
{
	sigma = cohesiveSurface->getSigma();
	return sigma;
}

const Vector&
DamageShearInterface::getStrain (void)
{
  return epsTrial;
}

int
DamageShearInterface::commitState (void)
{
	cohesiveSurface->commit();
	return 0;
}


int
DamageShearInterface::revertToLastCommit (void)     
{ 
  return 0;
}

int
DamageShearInterface::revertToStart (void)
{ 
  return 0;
}

NDMaterial*
DamageShearInterface::getCopy (void)
{
	CohesiveSurface* copyCohesiveSurface = cohesiveSurface->getCopy();    
	DamageShearInterface *theCopy = new DamageShearInterface(this->getTag(), copyCohesiveSurface);

  return theCopy;
}

NDMaterial*
DamageShearInterface::getCopy (const char *type) {
  return this->getCopy();
}

const char*
DamageShearInterface::getType (void) const
{
  return "BeamFiber";
}

int
DamageShearInterface::getOrder (void) const
{
  return 2;
}

double
DamageShearInterface::getMode() {
		return cohesiveSurface->getMode(epsTrial);
}




//print out material data
void 
DamageShearInterface::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "DamageShearInterface : " ; 
  s << this->getType( ) << endln ;
  s << endln ;
}

int 
DamageShearInterface::sendSelf(int commitTag, Channel &theChannel) 
{  return -1;  }

int 
DamageShearInterface::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{  return -1;  }


// methods to overwrite the standard implementation in NDMaterial for asking outputs
//----------------------------------------------------------------------------------

Response*
DamageShearInterface::setResponse (const char **argv, int argc, OPS_Stream &output) {

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
DamageShearInterface::getResponse (int responseID, Information &matInfo)
{
  Vector ep_dot(2);

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