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
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/WrappedMaterial.h,v $
                                                                        
#ifndef WrappedMaterial_h
#define WrappedMaterial_h

// Written: Francesco Vanin

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

//additional includes that I wrote here
//...

class Information;
class Response;
class Channel;
class NDMaterial;
class UniaxialMaterial;


class WrappedMaterial : public NDMaterial
{
  public:
	// full constructor
    WrappedMaterial(int tag, double E,  UniaxialMaterial* shearMat, double alpha=0.0);

	//null constructor
    WrappedMaterial();    

	//destructor
    ~WrappedMaterial();

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

	void Print( OPS_Stream &s, int flag );
    int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	Response* setResponse (const char **argv, int argc, OPS_Stream &output);
	int getResponse (int responseID, Information &matInfo);
  
//--------------------------------------------------------------------------------
  private:	

    Matrix Kpen;

	double E;

	UniaxialMaterial*  theShearMat;

	Vector stress, stressCommitted;
	Vector u, uCommitted;
	Matrix K, KCommitted;		

	double alpha;
};


#endif



