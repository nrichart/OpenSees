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
                                                                        
// $Revision: 6049 $
// $Date: 2015-07-17 06:56:36 +0200 (Fri, 17 Jul 2015) $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/dispBeamColumn/DispBeamColumn3d.h $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn3d.
// The element displacement field gives rise to constant axial strain,
// linear curvature, and constant twist angle.

#ifndef Macroelement3d_h
#define Macroelement3d_h

//#ifndef _bool_h
//#include "bool.h"
//#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <vector>
#include <string>

class Node;
class SectionForceDeformation;
class NDMaterial;
class Response;

class Macroelement3d : public Element
{
  public:
    Macroelement3d(int tag, int nd1, int nd2, int ndE, SectionForceDeformation *sI, SectionForceDeformation *sE, SectionForceDeformation *sJ, NDMaterial* shearModel, NDMaterial* shearModelOOP, 
		double h, double E_, 
		Vector driftF, Vector driftF_ALR, Vector driftS, Vector driftS_ALR, double Ltfc, double alphaNC_SD, double betaShearSpanF, double betaShearSpanS, double failureFactorF, double failureFactorS,
		Vector axis, Vector oop, Vector intLength, Vector intLengthMasses, Vector massDir, 
		int PDelta=0, double rho=0.0, int cm=0.0, double isGable=false);
    Macroelement3d();
    ~Macroelement3d();

    const char *getClassType(void) const {return "Macroelement3d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    int update(void);
	Vector getBasicDisplacement(Vector uI, Vector uJ, Vector uE);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);    

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    //int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

	//drift model
	void driftModel(double currentDriftF, double currentDriftS, double axialLoadRatio, double H0overL=1.0);

	// override standard damping methods
	//const Matrix &getTangentStiffDamping(void);
	//const Matrix &getSecantStiffDamping(void);
    //const Matrix &getInitialStiffDamping(void);
 
	//const Vector& getRayleighDampingForces(void); 
	//const Vector& getDampingForces(void); 
	//const Matrix &getDamp(void);


  protected:
    
  private:
    const Matrix &getInitialBasicStiff(void);
    const Matrix &getIncrementalCompatibilityMatrix(bool flagIncremental=true);
	int trasformMatrixToGlobal(Matrix& A);

    const int numSections;
	const int numShearModels;

    SectionForceDeformation **theSections; // pointer to the ND material objects
	NDMaterial **theShearModel;
	//UniaxialMaterial **theDampingModel;

    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[3];

    static Matrix K;		// Element stiffness, damping, and mass Matrix.
	Matrix GammaC;	        // Compatibility matrix, including or not P-Delta effects
    static Vector P;		// Element resisting force vector
    Matrix Tgl, Tgl6;       // Rotation 3d matrix from global to local coordinates. Full version and 6x6 version to rotate only one node dofs
    double R[3][3];

    Vector Q;                       // Applied nodal loads
    Vector q;                       // Basic force
	Vector uBasic, uBasicCommitted; // Basic displcaments

    double q0[12];  // Fixed end forces. Ordering can be improved. Used only for element loads (self weight, distributed)
    double p0[12];  // Reactions in basic system. Used only for element loads (self weight, distributed)

    double rho;    // Mass density per unit length
    int cMass;     // consistent mass flag
	Vector massglobalDir;
    int PDelta;
    double deltaW1, deltaV1, deltaW3, deltaV3; // relative displacements used to define the compatibility matrix with P-Delta transformation

	int parameterID;

    enum {maxNumSections = 20};

    static double workArea[];  // I think it is not used, not by me

    double* nodeIInitialDisp;
    double* nodeJInitialDisp;

    double* nodeIOffset;
    double* nodeJOffset;

    Vector xAxis, yAxis, zAxis;  // local orientation axes expressed in global coordinates
	Vector intLength, intLengthMasses;  
    double L;
	double E;

	double driftF, driftS;
	double Ltfc, betaShearSpanF, betaShearSpanS;
	Vector limitDriftF, limitDriftS, limitDriftF_ALR, limitDriftS_ALR;
	double alphaNC_SD;

	bool failedF, failedS;
	bool failedFcommitted, failedScommitted;

	double failureFactorF,failureFactorS;

	bool collapsedF, collapsedS;
	bool collapsedFcommitted, collapsedScommitted;

	bool isGable;
	double wx, wy, wz;

	double committedTime;

};

#endif

