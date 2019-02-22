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
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/GambarottaLagomarsinoModel.h,v $
                                                                        
#ifndef GambarottaLagomarsinoModel_h
#define GambarottaLagomarsinoModel_h

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


class GambarottaLagomarsinoModel : public NDMaterial
{
  public:
	// full constructor
    GambarottaLagomarsinoModel(int tag, double _E, double _G, double _c, double _mu, double _ct, double _beta, double _L, double _t, double _h, bool _elasticSolution=false);

	//null constructor
    GambarottaLagomarsinoModel();    

	//destructor
    ~GambarottaLagomarsinoModel();

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

	// internal methods
	double toughnessFunction(double alpha);
	double derToughnessFunction(double alpha);
	double sign(double V);

  
//--------------------------------------------------------------------------------
  private:	

    Matrix Kpen;

	double L,t,h;

	double mu;
	double c;
	double ct;
	double beta;

	double alpha, alphaCommitted;
	double s, sCommitted;
	double deltaLambda;
	double deltaAlpha;


	Vector stress, stressCommitted;
	Vector u, uCommitted;
	Matrix K, KCommitted;

	bool elasticSolution;

		
};


#endif



