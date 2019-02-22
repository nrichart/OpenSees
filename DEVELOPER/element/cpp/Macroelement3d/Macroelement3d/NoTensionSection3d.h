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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:46:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection3d.h,v $

#ifndef NoTensionSection3d_h
#define NoTensionSection3d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;

class NoTensionSection3d : public SectionForceDeformation
{
 public:
  NoTensionSection3d(int tag, double _k, double _kg, double _L, double _t, double _J, double _fc, int _nSections, bool stronger=false, bool elastic=false, bool crushing=true, bool spandrel=false);
  NoTensionSection3d(void);
  ~NoTensionSection3d(void);
  
  const char *getClassType(void) const {return "NoTensionSection3d";};
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  int setTrialSectionDeformation(const Vector&);
  const Vector &getSectionDeformation(void);
  
  const Vector &getStressResultant(void);
  const Matrix &getSectionTangent(void);
  const Matrix &getInitialTangent(void);
  
  SectionForceDeformation *getCopy(void);
  const ID &getType(void);
  int getOrder(void) const;
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag = 0);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector& getStressResultantSensitivity(int gradIndex,
					      bool conditional);
  const Matrix& getInitialTangentSensitivity(int gradIndex);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &info);
  Vector getCrushingVariables(int sliceNum);

 protected:
  
 private:
  
  double k, kg, L, t, J, fc;
  
  int nSections;
  double  OOPfactor;
  Matrix muX, muY;         // ductility demands for all section discretisations 
  Matrix zetaX, zetaY;     // damaged portions for all section discretisations 
  Matrix muXt, muYt;       // ductility demands for all section discretisations (trial)
  Matrix zetaXt, zetaYt;   // damaged portions for all section discretisations  (trial)

  Vector pos;              // positions of the sections
  Vector weight;           // weight of the sections
  
  Vector eCommitted, e, sCommitted, s;     // strain and stress vectors (committed and trial) 
  Matrix D, Dcommitted, Dtrial;   // stiffness matrix

  //Vector e, eTrial;			// section trial deformations
  
  //static Vector s, sTrial;
  //static Matrix ks;
  //static Matrix kCommitted;
  static ID code;

  bool stronger;
  bool elastic;
  bool crushing;
  bool spandrel;
  double factorStronger;

  int parameterID;
  int sliceOutput;
  double torsionalStiffnessFactor;
};

#endif
