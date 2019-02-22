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
                                                                        
// $Revision: 1.1 $
// $Date: 2008/12/09 20:00:16 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/DamageShearInterface.h,v $
                                                                        
#ifndef DamageShearInterface_h
#define DamageShearInterface_h

// Written: Francesco Vanin
//
// Description: This file contains the class definition for 
// AxialShearCombinedMaterial. AxialShearCombinedMaterial provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) DamageShearInterface.h, revA"

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include "CohesiveSurface.h"

//additional includes that I wrote here
//...

class Information;
class Response;
class Channel;
class NDMaterial;
//class CohesiveSurface;


class DamageShearInterface : public NDMaterial
{
  public:
	// full constructor
    DamageShearInterface(int tag, double _E, double _G, double c, double mu, double muR, double beta, double dropDrift, bool elasticSolution);

	// constructor for copies
    DamageShearInterface(int tag, CohesiveSurface* _copyCohesiveSurface);

	//null constructor
    DamageShearInterface();    

	//destructor
    ~DamageShearInterface();

	//methods
    int setTrialStrain (const Vector &strain);
    int setTrialStrain (const Vector &strain, const Vector &rate);
    int setTrialStrainIncr (const Vector &strain);
    int setTrialStrainIncr (const Vector &strain, const Vector &rate);

	const Matrix &getTangent( ) ;
    const Matrix &getInitialTangent (void);
    const Vector &getStress (void);
    const Vector &getStrain (void);
	const Matrix& getCommittedTangent(void);

	int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
	NDMaterial *getCopy (const char *type);
    const char *getType (void) const;
    int getOrder (void) const;
	double getMode();
	

	void Print( OPS_Stream &s, int flag );
    int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	Response* setResponse (const char **argv, int argc, OPS_Stream &output);
	int getResponse (int responseID, Information &matInfo);

  
//--------------------------------------------------------------------------------
  private:	

	  CohesiveSurface* cohesiveSurface;
	  Vector CommittedStrain;   // committed strain, to avoid useless iterations and wrong "elastic" updates
	  Vector epsTrial;

	  Vector sigma;
	  Matrix K;
	  bool elasticSolution;
	
};


#endif



