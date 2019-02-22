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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection3d.cpp,v $

#include "NoTensionSection3d.h"

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <Parameter.h>
#include <stdlib.h>
#include <string.h>
#include <MaterialResponse.h>

#include <classTags.h>
#include <elementAPI.h>

//Vector NoTensionSection3d::s(4);
//Matrix NoTensionSection3d::ks(4,4);
//Matrix NoTensionSection3d::kCommitted(4,4);
ID NoTensionSection3d::code(4);


NoTensionSection3d::NoTensionSection3d(void) 
	:SectionForceDeformation(0, SEC_TAG_Elastic3d), k(0.0), kg(0.0), L(0.0), t(0.0), J(0.0), fc(0.0), e(4), 
	eCommitted(4), sCommitted(4), s(4), D(4,4), Dtrial(4,4), Dcommitted(4,4), nSections(0), OOPfactor(1.0), stronger(false), elastic(false), crushing(true)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MY;	// OOP moment (moment around x axis for me)
    code(2) = SECTION_RESPONSE_MZ;	// in-plane moment (moment around y axis for me)
    code(3) = SECTION_RESPONSE_T;	// Torsion (moment around z axis for me)
  }
}

NoTensionSection3d::NoTensionSection3d (int tag, double _k, double _kg, double _L, double _t, double _J, double _fc, int _nSections, bool stronger, bool elastic, bool crushing, bool spandrel)
   :SectionForceDeformation(tag, SEC_TAG_Elastic3d),
   k(_k), kg(_kg), L(_L), t(_t), J(_J), fc(_fc), e(4), eCommitted(4), s(4), sCommitted(4), D(4,4), Dtrial(4,4), Dcommitted(4,4), nSections(_nSections), OOPfactor(0.0),
   muX(_nSections,2), muY(_nSections,2), zetaX(_nSections,2), zetaY(_nSections,2), pos(_nSections), weight(_nSections), sliceOutput(0),
   muXt(_nSections,2), muYt(_nSections,2), zetaXt(_nSections,2), zetaYt(_nSections,2), stronger(stronger), elastic(elastic), crushing(crushing), spandrel(spandrel), factorStronger(0.001),
   torsionalStiffnessFactor(1.0)
{  
	//if (J<0.0) {
		//J = 1/3.*L*t*t*t;  // infinitely narrow section
		//J = L*t*t*t * (1/3. - 0.21*t/L*(1.-pow(t,4)/(12*pow(L,4))));
	//}


    if (L>t) {
		J = 1/3.*L*t*t*t;
		J = L*t*t*t * (1/3. - 0.21*t/L*(1.-pow(t,4)/(12*pow(L,4))));
	} else {
		J = 1/3.*L*L*L*t;
		J = t*L*L*L * (1/3. - 0.21*L/t*(1.-pow(L,4)/(12*pow(t,4))));
	}



	D(0,0) = k*t*L;
    D(1,1) = k*L*t*t*t /12.0;
    D(2,2) = k*L*L*L*t /12.0;

	D(3,3) = kg*J;

	D(3,3)*= torsionalStiffnessFactor;

	if (stronger && !elastic) {
		D(0,0) += factorStronger*k*t*L;
        D(1,1) += factorStronger*k*L*t*t*t /12.0;
        D(2,2) += factorStronger*k*L*L*L*t /12.0;
		D(3,3) += factorStronger*kg*J;
	} 

	Dtrial = D;
	Dcommitted= D;

	s.Zero();
	sCommitted.Zero();

  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MY;	// Mz is the second
    code(2) = SECTION_RESPONSE_MZ;	// My is the third 
    code(3) = SECTION_RESPONSE_T;	// T is the fourth
  }

  if (nSections==1) {
	  pos(0) = 0.0;
	  weight(0) = 1;
  } else {
	 double slice = 1./(nSections);
	 for (int i=0; i<nSections; i++) {
		  pos(i) = -0.5 + slice/2. + slice * i;
		  weight(i) = slice;
	 }
  }


  for (int i=0; i<nSections; i++) {
	muX(i,0) = 1.;
	muX(i,1) = 1.;
	muY(i,0) = 1.;
	muY(i,1) = 1.;
  }

}

NoTensionSection3d::~NoTensionSection3d(void)
{
    return;
}

int 
NoTensionSection3d::commitState(void)
{
	eCommitted = e;
	sCommitted = s;

	Dcommitted = Dtrial;

    muX = muXt;
	muY = muYt;
	zetaX = zetaXt;
	zetaY = zetaYt;

  return 0;
}

int 
NoTensionSection3d::revertToLastCommit(void)
{
	e = eCommitted;
	s = sCommitted;
	Dtrial = Dcommitted;

	muXt = muX;
	muYt = muY;
	zetaXt = zetaX;
	zetaYt = zetaY;


  return 0;
}

int 
NoTensionSection3d::revertToStart(void)
{
	muX.Zero();
	muY.Zero();
	zetaX.Zero();
	zetaY.Zero();

	e.Zero();
	s.Zero();
	eCommitted.Zero();
	sCommitted.Zero();
	Dtrial = D;
	Dcommitted = D;

  return 0;
}

int
NoTensionSection3d::setTrialSectionDeformation (const Vector &def)
{
	Vector increment(4);
	increment = e - def;

	if (fabs(increment.Norm()) < DBL_EPSILON) {
      return 0;
	}

	//opserr << "new strain: " << def; 
    e = def;
	increment = e - eCommitted;

	if (elastic) {
		s = D*e;
		Dtrial = D;
		return 0;
	}
	
	muXt = muX;
	muYt = muY;
	zetaXt = zetaX;
	zetaYt = zetaY;

	//opserr << "e: " << e(0) << ", " << e(2) << " ( " << -abs(e(2)) + 2.0*e(0)/t << " )\n";

	//---------------------------------------------------------------------------------------------------------//
	// set stress
	//---------------------------------------------------------------------------------------------------------//

	s.Zero();
	Dtrial.Zero();

	if (abs(e(2)) < DBL_EPSILON) {
		// zero in plane moment
		if (abs(e(1)) + 2.0*e(0)/t > DBL_EPSILON) {
			if  (-abs(e(1)) + 2.0*e(0)/t > DBL_EPSILON) {
				// all tension, all zeroes, do nothing
			} else {
				
					// case 5
					double w = e(0);   // vertical displacement
					double fiX = e(1);  // rotation around x, fi_x
					double fiY = 0.0;  // rotation around y, fi_y

					if (fiX<0.0) {
						// case 4a,c: both out from different sides, fiX>=0
							 s(0) = (pow(fiY,2)*pow(L,3) + 3*L*pow(fiX*t + 2*w,2))/(24.*fiX);
							 s(1) = -(pow(fiY,2)*pow(L,3)*w + L*(-(fiX*t) + w)*pow(fiX*t + 2*w,2))/(24.*pow(fiX,2));
							 s(2) = (fiY*pow(L,3)*(fiX*t + 2*w))/(24.*fiX);

					} else {
						// case 4b,d: both out from different sides, fiX<0
							 s(0) = -(pow(fiY,2)*pow(L,3) + 3*L*pow(fiX*t - 2*w,2))/(24.*fiX);
							 s(1) = (pow(fiY,2)*pow(L,3)*w + L*pow(fiX*t - 2*w,2)*(fiX*t + w))/(24.*pow(fiX,2));
							 s(2) = (fiY*pow(L,3)*(fiX*t - 2*w))/(24.*fiX);
					}
			}

		} else {
					s(0) = t*L *e(0);
                    s(1) = L*t*t*t /12.0 * e(1);
                    s(2) = L*L*L*t /12.0 * e(2);
		}


	} else {
		double w = e(0);   // vertical displacement
		double fiX = e(1);  // rotation around x, fi_x
		double fiY = e(2);  // rotation around y, fi_y

		double t1 = (2.0*w + t*fiX)/(2.0*fiY);   // x coordinate of the neutral axis for y = +t/2
		double t2 = (2.0*w - t*fiX)/(2.0*fiY);   // x coordinate of the neutral axis for y = -t/2

		int inside1 = 0;
		if (t1<-L/2.0)   inside1=-1;
		if (t1> L/2.0)   inside1=+1;

	    int inside2 = 0;
		if (t2<-L/2.0)   inside2=-1;
		if (t2> L/2.0)   inside2=+1;
		
		if (inside1==0 && inside2==0) {
			if (fiY>0.0) {
			   // case 1a,b: both neutral points inside the long side, fi_y>0
				 s(0) = -(t*(pow(fiX,2)*pow(t,2) + 3*pow(fiY*L - 2*w,2)))/(24.*fiY);
				 s(1) = (fiX*pow(t,3)*(fiY*L - 2*w))/(24.*fiY);
				 s(2) = (t*(pow(fiY,3)*pow(L,3) - 3*pow(fiY,2)*pow(L,2)*w + pow(fiX,2)*pow(t,2)*w + 4*pow(w,3)))/(24.*pow(fiY,2));

			} else {
               // case 1c,d: both neutral points inside the long side, fi_y<0
				 s(0) = (t*(pow(fiX,2)*pow(t,2) + 3*pow(fiY*L + 2*w,2)))/(24.*fiY);
				 s(1) = (fiX*pow(t,3)*(fiY*L + 2*w))/(24.*fiY);
				 s(2) = (t*(pow(fiY,3)*pow(L,3) + 3*pow(fiY,2)*pow(L,2)*w - pow(fiX,2)*pow(t,2)*w - 4*pow(w,3)))/(24.*pow(fiY,2));
			}

		} else {
			if (inside1*inside2>0) {
				if (w - L/2.0*fiY + t/2.0*fiX < DBL_EPSILON) {
				    // case 0a: both neutral points ouside the same side, all compressed
					s(0) = t*L *e(0);
                    s(1) = L*t*t*t /12.0 * e(1);
                    s(2) = L*L*L*t /12.0 * e(2);
				} else {
					// case 0b: both neutral points ouside the same side, all tension
				    s(0) = 0.0;
                    s(1) = 0.0;
                    s(2) = 0.0;
					//opserr << "open in tension\n";
				}
			} else {
				if (inside1*inside2==0) {
					// one in and one out. Cases 2-3
					if (inside2>0) {
						if (fiY>0.0) {
							// case 2a
							 s(0) = -pow(-(fiY*L) + fiX*t + 2*w,3)/(48.*fiX*fiY);
							 s(1) = -((fiY*L + 3*fiX*t - 2*w)*pow(-(fiY*L) + fiX*t + 2*w,3))/(384.*pow(fiX,2)*fiY);
							 s(2) = (pow(-(fiY*L) + fiX*t + 2*w,3)*(3*fiY*L + fiX*t + 2*w))/(384.*fiX*pow(fiY,2));

						} else {
							// case 3a
						     s(0) = L*t*w + pow(-(fiY*L) + fiX*t + 2*w,3)/(48.*fiX*fiY);
							 s(1) = (fiX*L*pow(t,3))/12. + ((fiY*L + 3*fiX*t - 2*w)*pow(-(fiY*L) + fiX*t + 2*w,3))/(384.*pow(fiX,2)*fiY);
							 s(2) = (fiY*pow(L,3)*t)/12. - (pow(-(fiY*L) + fiX*t + 2*w,3)*(3*fiY*L + fiX*t + 2*w))/(384.*fiX*pow(fiY,2));
						}

					} else {
						if (inside2<0) {
							if (fiY<0.0) {
								// case 2c
								 s(0) = pow(fiY*L + fiX*t + 2*w,3)/(48.*fiX*fiY);
								 s(1) = -((fiY*L - 3*fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,3))/(384.*pow(fiX,2)*fiY);
								 s(2) = -((-3*fiY*L + fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,3))/(384.*fiX*pow(fiY,2));;

							} else {
								// case 3c
									 s(0) = L*t*w - pow(fiY*L + fiX*t + 2*w,3)/(48.*fiX*fiY);
									 s(1) = (fiX*L*pow(t,3))/12. + ((fiY*L - 3*fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,3))/(384.*pow(fiX,2)*fiY);
									 s(2) = (fiY*pow(L,3)*t)/12. + ((-3*fiY*L + fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,3))/(384.*fiX*pow(fiY,2));

							}

						} else {
							if (inside1>0) {
								if (fiY>0.0) {
									// case 2b
								    s(0) = -pow(fiY*L + fiX*t - 2*w,3)/(48.*fiX*fiY);
									s(1) = (pow(fiY*L + fiX*t - 2*w,3)*(-(fiY*L) + 3*fiX*t + 2*w))/(384.*pow(fiX,2)*fiY);
									s(2) = (pow(fiY*L + fiX*t - 2*w,3)*(3*fiY*L - fiX*t + 2*w))/(384.*fiX*pow(fiY,2));

								} else {
									// case 3b
								     s(0) = pow(fiY*L + fiX*t - 2*w,3)/(48.*fiX*fiY) + L*t*w;
									 s(1) = (fiX*L*pow(t,3))/12. + ((fiY*L - 3*fiX*t - 2*w)*pow(fiY*L + fiX*t - 2*w,3))/(384.*pow(fiX,2)*fiY);
									 s(2) = (fiY*pow(L,3)*t)/12. + ((-3*fiY*L + fiX*t - 2*w)*pow(fiY*L + fiX*t - 2*w,3))/(384.*fiX*pow(fiY,2));
								}

							} else {
								if (fiY<0.0) {
									// case 2d
									 s(0) = pow(-(fiY*L) + fiX*t - 2*w,3)/(48.*fiX*fiY);
									 s(1) = (pow(fiY*L - fiX*t + 2*w,3)*(fiY*L + 3*fiX*t + 2*w))/(384.*pow(fiX,2)*fiY);
									 s(2) = (pow(-(fiY*L) + fiX*t - 2*w,3)*(3*fiY*L + fiX*t - 2*w))/(384.*fiX*pow(fiY,2));
								} else {
									// case 3d
							 s(0) = L*t*w + pow(fiY*L - fiX*t + 2*w,3)/(48.*fiX*fiY);
							 s(1) = (fiX*L*pow(t,3))/12. + (pow(-(fiY*L) + fiX*t - 2*w,3)*(fiY*L + 3*fiX*t + 2*w))/(384.*pow(fiX,2)*fiY);
							 s(2) = (fiY*pow(L,3)*t)/12. + ((3*fiY*L + fiX*t - 2*w)*pow(fiY*L - fiX*t + 2*w,3))/(384.*fiX*pow(fiY,2));

								}
							}
						}
					}

				} else {
								//	opserr << "entrato sbagliato\n";
					if (fiX<0.0) {
						// case 4a,c: both out from different sides, fiX>=0
							 s(0) = (pow(fiY,2)*pow(L,3) + 3*L*pow(fiX*t + 2*w,2))/(24.*fiX);
							 s(1) = -(pow(fiY,2)*pow(L,3)*w + L*(-(fiX*t) + w)*pow(fiX*t + 2*w,2))/(24.*pow(fiX,2));
							 s(2) = (fiY*pow(L,3)*(fiX*t + 2*w))/(24.*fiX);

					} else {
						// case 4b,d: both out from different sides, fiX<0							 
							 s(0) = -(pow(fiY,2)*pow(L,3) + 3*L*pow(fiX*t - 2*w,2))/(24.*fiX);
							 s(1) = (pow(fiY,2)*pow(L,3)*w + L*pow(fiX*t - 2*w,2)*(fiX*t + w))/(24.*pow(fiX,2));
							 s(2) = (fiY*pow(L,3)*(fiX*t - 2*w))/(24.*fiX);
					}
				}
			}
		}

		}

	    s *= k;
	    s(3) = kg*J*e(3);

		s(3)*= torsionalStiffnessFactor;

		// ----------------------------------------------------------------------------------------------------------------------
		// CRUSHING CORRECTION
		// ----------------------------------------------------------------------------------------------------------------------

		// crushing correction (IP slices) 
		//---------------------------------------------------------------------------------------------------------//
		double wLoc;
	    double fiLoc;
		double muNew;
		double zetaNew;
		double corr;
		double IPfactor;
		Vector dmu0_de(3);
		Vector dmu1_de(3);
		Vector dzeta0_de(3);
		Vector dzeta1_de(3);
		Vector dCorr_de(3);


		// enter only if there is some stress (not if it is open in tension, and so all stresses are zero)
		if ((s.Norm() > DBL_EPSILON) && crushing) {

		if ((sqrt(pow(e(1),2) + pow(e(2),2)) ) > DBL_EPSILON )  // if there is some rotation define the proportion of IP/OOP rotation

		// defines how much the section is bent in plane or out of plane
		// remind that in this case the oop and ip direction are switched (L is the shortest length of the section)

			//if (abs(e(1))*t/2. > 1./100. *(fc/k) || abs(e(2))*L/2. > 1./100. *(fc/k) ) 
			//	OOPfactor = abs(e(1)) / sqrt(pow(e(1),2.) + pow(e(2),2.));

		    OOPfactor = 1.0;  // apply only oop correction
            IPfactor = 1.0-OOPfactor;

			muXt = muX;
			zetaXt = zetaX;

			for (int i=0; i<nSections; i++) {
				
				dmu0_de.Zero();
				dmu1_de.Zero();
				dzeta0_de.Zero();
				dzeta1_de.Zero();
				dCorr_de.Zero();

				wLoc  = e(0) - (L*pos(i))*e(2);  // minus sign consistent with our reference system
				fiLoc = e(1);
	
				if (fiLoc > DBL_EPSILON) {  // negative side (more) compressed. displ[z_] := wLoc + fiLoc*z; 
					muNew = -(wLoc + fiLoc*(-t/2.)) / (fc/k);   // it is mu' in Penna (2014): current ductility demand (at the negative edge)
					zetaNew = (muNew-1.) * (fc/k) / (fiLoc*t);	
					        /*
							opserr <<"zetaNew=" << zetaNew << endln;
							opserr << "e=" << e;
							opserr << "s=" << s;
							opserr << "muNew=" << muNew << endln<< endln;
							*/

					if (muNew > muXt(i,0))   {
						// update ductility demand
						muXt(i,0) = muNew;

						// calculate derivative
						dmu0_de(0) = -k/fc;
						dmu0_de(1) = 0.5*k*t/fc;
						dmu0_de(2) = k*L*pos(i)/fc;
					}

					if (zetaNew > zetaXt(i,0)) {
						// update depth of nonlinear area
						if (zetaNew<=1.0) {
							zetaXt(i,0) = zetaNew;

							dzeta0_de(0) = -1.0/(fiLoc*t);  
							dzeta0_de(1) = (fc + e(0)*k - e(2)*k*L*pos(i)) / (e(1)*e(1) * k*t);
							dzeta0_de(2) = (L*pos(i)) / (e(1) *t);


						} else {					
							zetaXt(i,0) = 1.0;
							zetaXt(i,1) = 1.0;
							/*
							zetaXt(i,0) = 0.999;
							opserr <<"zetaNew=" << zetaNew << endln;
							opserr << "e=" << e;
							opserr << "s=" << s;
							opserr << "muNew=" << muNew << endln<< endln;

							while (zetaXt(i,0) < 1.0) 
								zetaXt(i,0) = 0.999;
							*/

							//dzeta0_de(0) = -1.0/(fiLoc*t);  
							//dzeta0_de(1) = (fc + e(0)*k - e(2)*k*L*pos(i)) / (e(1)*e(1) * k*t);
							//dzeta0_de(2) = (L*pos(i)) / (e(1) *t);
							
							
							// update ductility demand at the other side
							double muOtherSide = muNew *(zetaNew-1.0)/zetaNew;

							if (muOtherSide > muXt(i,1)) {
								muXt(i,1) = muOtherSide;

								dmu1_de(0) = ((-1. + zetaNew)*dmu0_de(0))/zetaNew;
								dmu1_de(1) = ((-1. + zetaNew)*dmu0_de(1))/zetaNew;
								dmu1_de(2) = ((-1. + zetaNew)*dmu0_de(2))/zetaNew;
							}	
							
						}
					}

	
				} else if (fiLoc < -DBL_EPSILON) {               // positive side (more) compressed. displ[z_] := wLoc + fiLoc*z; 
					
					muNew = (-(wLoc + t/2.*fiLoc)) / (fc/k);  
					zetaNew = (muNew-1) * (fc/k) / (-fiLoc*t);			

					if (muNew > muXt(i,1))  {
						muXt(i,1) = muNew;
						
						dmu1_de(0) = -k/fc;
						dmu1_de(1) = -0.5*k*t/fc;
						dmu1_de(2) = k*L*pos(i)/fc;
					}

					if (zetaNew > zetaXt(i,1)) {
						// update depth of nonlinear area
						if (zetaNew<=1.0) {
							zetaXt(i,1) = zetaNew;

							dzeta1_de(0) = 1.0/(fiLoc*t);  
							dzeta1_de(1) = -(fc + e(0)*k - e(2)*k*L*pos(i)) / (fiLoc*fiLoc * k*t);
							dzeta1_de(2) = -(L*pos(i)) / (fiLoc *t);

						} else {
							
							zetaXt(i,0) = 1.0;
							zetaXt(i,1) = 1.0;
							
							// update ductility demand at the other side
							double muOtherSide = muNew *(zetaNew-1.0)/zetaNew;

							if (muOtherSide > muXt(i,0)) {
								muXt(i,0) = muOtherSide;

								dmu0_de(0) = ((-1. + zetaNew)*dmu1_de(0))/zetaNew;
								dmu0_de(1) = ((-1. + zetaNew)*dmu1_de(1))/zetaNew;
								dmu0_de(2) = ((-1. + zetaNew)*dmu1_de(2))/zetaNew;
							}				

							//zetaXt(i,1) = 0.999;
							//dzeta1_de(0) = -1.0/(fiLoc*t);  
							//dzeta1_de(1) = (fc + e(0)*k - e(2)*k*L*pos(i)) / (e(1)*e(1) * k*t);
							//dzeta1_de(2) = (L*pos(i)) / (e(1) *t);

						}
					}


				} else {      // the section rotation is basically null

					//opserr << "zero rotation. Strains set: " << wLoc << ", " << fiLoc << " (" << e;

					muNew = (-wLoc) / (fc/k);  

					//opserr << "muNew = " << muNew << endln;
					// zeta is in theory infinite

					if (muNew > muXt(i,0))  {
						muXt(i,0) = muNew;
						
						dmu0_de(0) = -k/fc;
						dmu0_de(1) = -0.5*k*t/fc;
						dmu0_de(2) = k*L*pos(i)/fc;
					}

					if (muNew > muXt(i,1))  {
						muXt(i,1) = muNew;
						
						dmu1_de(0) = -k/fc;
						dmu1_de(1) = -0.5*k*t/fc;
						dmu1_de(2) = k*L*pos(i)/fc;
					}

					if (zetaXt(i,0) != 1.0 && muNew>1.) {  // was not already all crushed and there is some crushing now
						zetaXt(i,0) = 1.0;
					}
			        if (zetaXt(i,1) != 1.0 && muNew>1.) {  // was not already all crushed and there is some crushing now
						zetaXt(i,1) = 1.0;
					}
				}

				// apply correction term

				if (muXt(i,0)>1.0 && (wLoc + (-t/2.)*fiLoc) < 0.0) {  

					// some nonlinearity was recorded at the negative side (in this or a previous step) and the section is not open at that side 		    
					corr = (muXt(i,0)-1.0) / (2.0*muXt(i,0)) * (zetaXt(i,0)*t) *L*weight(i) * (-wLoc+fiLoc*(t/2))*k;  // positive

					s(0) += OOPfactor*corr;
					s(1) -= corr * ( (t/2.) - (zetaXt(i,0)*t)/3. );
					s(2) -= corr * (L*pos(i));

					dCorr_de(0) = (L*t*weight(i)*(zetaXt(i,0)*((0.5-0.5*muXt(i,0))*muXt(i,0) + (-0.5*e(0) + 0.5*e(2)*L*pos(i) + 0.25*e(1)*t)*dmu0_de(0)) + 
						           muXt(i,0)*(0.5*e(0) - 0.5*e(2)*L*pos(i) - 0.25*e(1)*t + (-0.5*e(0) + 0.5*e(2)*L*pos(i) + 0.25*e(1)*t)*muXt(i,0))*dzeta0_de(0))) / pow(muXt(i,0),2) *k;

					dCorr_de(1) = (L*t*weight(i)*(zetaXt(i,0)*(t*(-0.25 + 0.25*muXt(i,0))*muXt(i,0) + (-0.5*e(0) + 0.5*e(2)*L*pos(i) + 0.25*e(1)*t)*dmu0_de(1)) + 
                                   muXt(i,0)*(0.5*e(0) - 0.5*e(2)*L*pos(i) - 0.25*e(1)*t + (-0.5*e(0) + 0.5*e(2)*L*pos(i) + 0.25*e(1)*t)*muXt(i,0))*dzeta0_de(1)))/pow(muXt(i,0),2) *k;

					dCorr_de(2) = (L*t*weight(i)*(zetaXt(i,0)*(L*pos(i)*(-0.5 + 0.5*muXt(i,0))*muXt(i,0) + (-0.5*e(0) + 0.5*e(2)*L*pos(i) + 0.25*e(1)*t)*dmu0_de(2)) + 
                                   muXt(i,0)*(0.5*e(0) - 0.5*e(2)*L*pos(i) - 0.25*e(1)*t + (-0.5*e(0) + 0.5*e(2)*L*pos(i) + 0.25*e(1)*t)*muXt(i,0))*dzeta0_de(2)))/pow(muXt(i,0),2) *k;

					
					Dtrial(0,0) += OOPfactor*dCorr_de(0) /k;
					Dtrial(0,1) += OOPfactor*dCorr_de(1) /k;
					Dtrial(0,2) += OOPfactor*dCorr_de(2) /k;

					Dtrial(1,0) -= ( dCorr_de(0)* (  (t/2.)-(zetaXt(i,0)*t)/3.) - corr*t/3.*dzeta0_de(0) )/k;
					Dtrial(1,1) -= ( dCorr_de(1)* (  (t/2.)-(zetaXt(i,0)*t)/3.) - corr*t/3.*dzeta0_de(1) )/k;
					Dtrial(1,2) -= ( dCorr_de(2)* (  (t/2.)-(zetaXt(i,0)*t)/3.) - corr*t/3.*dzeta0_de(2) )/k;

					Dtrial(2,0) -= dCorr_de(0)* (L*pos(i)) /k;
					Dtrial(2,1) -= dCorr_de(1)* (L*pos(i)) /k;
					Dtrial(2,2) -= dCorr_de(2)* (L*pos(i)) /k;
					
					

					//opserr << "corr=" << corr << "; mu=" << muXt(i,0)<<"; zeta = "<< zetaXt(i,0) << endln;
					
					
										
				}


				if (muXt(i,1)>1.0 && (wLoc + (t/2.)*fiLoc) < 0.0) {

					// some nonlinearity was recorded at the positive side (in this or a previous step) and the section is not open at that side 		    
					corr = (muXt(i,1)-1.0) / (2.0*muXt(i,1)) * (zetaXt(i,1)*t) *L*weight(i) * (-wLoc-fiLoc*(t/2))*k;  // positive

					s(0) += OOPfactor*corr;
					s(1) += corr * ( (t/2.) - (zetaXt(i,1)*t)/3. );
					s(2) -= corr * (L*pos(i));


					dCorr_de(0) = (L*t*weight(i)*(zetaXt(i,1)*((0.5 - 0.5*muXt(i,1))*muXt(i,1) + (-0.5*e(0) + 0.5*e(2)*L*pos(i) - 0.25*e(1)*t)*dmu1_de(0)) + 
                                  muXt(i,1)*(0.5*e(0) - 0.5*e(2)*L*pos(i) + 0.25*e(1)*t + (-0.5*e(0) + 0.5*e(2)*L*pos(i)- 0.25*e(1)*t)*muXt(i,1))*dzeta1_de(0)))/pow(muXt(i,1),2) *k;

					dCorr_de(1) = (L*t*weight(i)*(zetaXt(i,1)*(t*(0.25 - 0.25*muXt(i,1))*muXt(i,1) + (-0.5*e(0) + 0.5*e(2)*L*pos(i) - 0.25*e(1)*t)*dmu1_de(1)) + 
                                   muXt(i,1)*(0.5*e(0) - 0.5*e(2)*L*pos(i) + 0.25*e(1)*t + (-0.5*e(0) + 0.5*e(2)*L*pos(i) - 0.25*e(1)*t)*muXt(i,1))*dzeta1_de(1)))/pow(muXt(i,1),2) *k;

					dCorr_de(2) = (L*t*weight(i)*(zetaXt(i,1)*(L*pos(i)*(-0.5 + 0.5*muXt(i,1))*muXt(i,1) + (-0.5*e(0) + 0.5*e(2)*L*pos(i) - 0.25*e(1)*t)*dmu1_de(2)) + 
                                    muXt(i,1)*(0.5*e(0) - 0.5*e(2)*L*pos(i) + 0.25*e(1)*t + (-0.5*e(0) + 0.5*e(2)*L*pos(i) - 0.25*e(1)*t)*muXt(i,1))*dzeta1_de(2)))/pow(muXt(i,1),2) *k;							

					/*
					Dtrial(0,0) += OOPfactor*dCorr_de(0) /k;
					Dtrial(0,1) += OOPfactor*dCorr_de(1) /k;
					Dtrial(0,2) += OOPfactor*dCorr_de(2) /k;

					Dtrial(1,0) += ( dCorr_de(0)* (  (t/2.)-(zetaXt(i,1)*t)/3.) - corr*t/3.*dzeta1_de(0) )/k;
					Dtrial(1,1) += ( dCorr_de(1)* (  (t/2.)-(zetaXt(i,1)*t)/3.) - corr*t/3.*dzeta1_de(1) )/k;
					Dtrial(1,2) += ( dCorr_de(2)* (  (t/2.)-(zetaXt(i,1)*t)/3.) - corr*t/3.*dzeta1_de(2) )/k;

					Dtrial(2,0) -= dCorr_de(0)* (L*pos(i)) /k;
					Dtrial(2,1) -= dCorr_de(1)* (L*pos(i)) /k;
					Dtrial(2,2) -= dCorr_de(2)* (L*pos(i)) /k;
					*/
					//opserr << "corr2=" << corr << "; mu2=" << muXt(i,1)<<"; zeta2 = "<< zetaXt(i,1) << endln;
	
				}

				//opserr << "sCorr=" << s;
				//opserr << "fc=" << fc << "; k=" << k<< endln;


			}  // end for sections
	} // end if crushing and no all tension


	// end crushing OOP

	if (stronger) {
		
		s(0) += factorStronger* k*t*L           * e(0);
        s(1) += factorStronger* k*L*t*t*t /12.0 * e(1);
        s(2) += factorStronger* k*L*L*L*t /12.0 * e(2);
		s(3) += factorStronger* kg*J            * e(3);
		
	}

	//if (abs(s(0))<DBL_EPSILON && abs(s(1))<DBL_EPSILON) 
	//	s(3) *= factorStronger; //0.001 ;

	//opserr << "new stress: " << s; 

	//---------------------------------------------------------------------------------------------------------//
	// set tangent
	//---------------------------------------------------------------------------------------------------------//


	if (abs(e(2)) < DBL_EPSILON) {
		// zero in plane moment
		if (abs(e(1)) + 2.0*e(0)/t > DBL_EPSILON) {
			if  (-abs(e(1)) + 2.0*e(0)/t > DBL_EPSILON) {
				// all tension, all zeroes, do nothing
				Dtrial.Zero();
			} else {
					// case 5
					double w = e(0);   // vertical displacement
					double fiX = e(1);  // rotation around x, fi_x
					double fiY = 0.0;  // rotation around y, fi_y

				  if (fiX<0.0) {
						// case 4a,c: both out from different sides, fiX>=0
						 Dtrial(0,0) += L*(t/2. + w/fiX); 
						 Dtrial(1,0) += -(L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(2,0) += (fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(0,1) += -(L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(1,1) += (L*(pow(fiX,3)*pow(t,3) + 2*pow(fiY,2)*pow(L,2)*w + 8*pow(w,3)))/(24.*pow(fiX,3)); 
						 Dtrial(2,1) += -(fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(0,2) += (fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(1,2) += -(fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(2,2) += (pow(L,3)*(fiX*t + 2*w))/(24.*fiX); 

					} else {
						// case 4b,d: both out from different sides, fiX<0
						 Dtrial(0,0) += (L*(t - (2*w)/fiX))/2.; 
						 Dtrial(1,0) += (L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(2,0) += -(fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(0,1) += (L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(1,1) += (L*(pow(fiX,3)*pow(t,3) - 2*pow(fiY,2)*pow(L,2)*w - 8*pow(w,3)))/(24.*pow(fiX,3)); 
						 Dtrial(2,1) += (fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(0,2) += -(fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(1,2) += (fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(2,2) += (pow(L,3)*(fiX*t - 2*w))/(24.*fiX); 
				    }
			}

		} else {
			Dtrial(0,0) += t*L;
            Dtrial(1,1) += L*t*t*t /12.0;
            Dtrial(2,2) += L*L*L*t /12.0;
		}
	} else {
		double w = e(0);   // vertical displacement
		double fiX = e(1);  // rotation around x, fi_x
		double fiY = e(2);  // rotation around y, fi_y

		double t1 = (2.0*w + t*fiX)/(2.0*fiY);   // x coordinate of the neutral axis for y = +t/2
		double t2 = (2.0*w - t*fiX)/(2.0*fiY);   // x coordinate of the neutral axis for y = -t/2

		int inside1 = 0;
		if (t1<-L/2.0)   inside1=-1;
		if (t1> L/2.0)   inside1=+1;

	    int inside2 = 0;
		if (t2<-L/2.0)   inside2=-1;
		if (t2> L/2.0)   inside2=+1;
		
		if (inside1==0 && inside2==0) {
			if (fiY>0.0) {
			   // case 1a,b: both neutral points inside the long side, fi_y>0
			     Dtrial(0,0) += (t*(L - (2*w)/fiY))/2.; 
				 Dtrial(1,0) += -(fiX*pow(t,3))/(12.*fiY); 
				 Dtrial(2,0) += (t*(-3*pow(fiY,2)*pow(L,2) + pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiY,2)); 
				 Dtrial(0,1) += -(fiX*pow(t,3))/(12.*fiY); 
				 Dtrial(1,1) += (pow(t,3)*(fiY*L - 2*w))/(24.*fiY); 
				 Dtrial(2,1) += (fiX*pow(t,3)*w)/(12.*pow(fiY,2)); 
				 Dtrial(0,2) += (t*(-3*pow(fiY,2)*pow(L,2) + pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiY,2)); 
				 Dtrial(1,2) += (fiX*pow(t,3)*w)/(12.*pow(fiY,2)); 
				 Dtrial(2,2) += (t*(pow(fiY,3)*pow(L,3) - 2*pow(fiX,2)*pow(t,2)*w - 8*pow(w,3)))/(24.*pow(fiY,3)); 

			} else {
               // case 1c,d: both neutral points inside the long side, fi_y<0
				 Dtrial(0,0) += (L*t)/2. + (t*w)/fiY; 
				 Dtrial(1,0) += (fiX*pow(t,3))/(12.*fiY); 
				 Dtrial(2,0) += -(t*(-3*pow(fiY,2)*pow(L,2) + pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiY,2)); 
				 Dtrial(0,1) += (fiX*pow(t,3))/(12.*fiY); 
				 Dtrial(1,1) += (pow(t,3)*(fiY*L + 2*w))/(24.*fiY); 
				 Dtrial(2,1) += -(fiX*pow(t,3)*w)/(12.*pow(fiY,2)); 
				 Dtrial(0,2) += -(t*(-3*pow(fiY,2)*pow(L,2) + pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiY,2)); 
				 Dtrial(1,2) += -(fiX*pow(t,3)*w)/(12.*pow(fiY,2)); 
				 Dtrial(2,2) += (t*(pow(fiY,3)*pow(L,3) + 2*pow(fiX,2)*pow(t,2)*w + 8*pow(w,3)))/(24.*pow(fiY,3)); 
			}

		} else {
			if (inside1*inside2>0) {
				if (w - L/2.0*fiY + t/2.0*fiX < DBL_EPSILON) {
				    // case 0a: both neutral points ouside the same side, all compressed
					Dtrial(0,0) += t*L;
                    Dtrial(1,1) += L*t*t*t /12.0;
                    Dtrial(2,2) += L*L*L*t /12.0;
				} else {
					// case 0b: both neutral points ouside the same side, all tension
				    // all zeroes, do nothing
				}
			} else {
				if (inside1*inside2==0) {
					// one in and one out. Cases 2-3
					if (inside2>0) {
						if (fiY>0.0) {
							// case 2a
							Dtrial(0,0) += -pow(-(fiY*L) + fiX*t + 2*w,2)/(8.*fiX*fiY); 
							Dtrial(1,0) += -((fiY*L + 2*fiX*t - 2*w)*pow(-(fiY*L) + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
							Dtrial(2,0) += (pow(-(fiY*L) + fiX*t + 2*w,2)*(2*fiY*L + fiX*t + 2*w))/(48.*fiX*pow(fiY,2)); 
							Dtrial(0,1) += -((fiY*L + 2*fiX*t - 2*w)*pow(-(fiY*L) + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
							Dtrial(1,1) += -(pow(-(fiY*L) + fiX*t + 2*w,2)*(pow(fiY,2)*pow(L,2) + 2*fiX*fiY*L*t + 3*pow(fiX,2)*pow(t,2) - 4*(fiY*L + fiX*t)*w + 4*pow(w,2)))/(192.*pow(fiX,3)*fiY); 
							Dtrial(2,1) += (pow(-(fiY*L) + fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) - 4*fiY*L*w + 4*fiX*t*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
							Dtrial(0,2) += (pow(-(fiY*L) + fiX*t + 2*w,2)*(2*fiY*L + fiX*t + 2*w))/(48.*fiX*pow(fiY,2)); 
							Dtrial(1,2) += (pow(-(fiY*L) + fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) - 4*fiY*L*w + 4*fiX*t*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
							Dtrial(2,2) += -(pow(-(fiY*L) + fiX*t + 2*w,2)*(3*pow(fiY,2)*pow(L,2) + 2*fiY*L*(fiX*t + 2*w) + pow(fiX*t + 2*w,2)))/(192.*fiX*pow(fiY,3)); 

						} else {
							// case 3a
							 Dtrial(0,0) += L*t + pow(-(fiY*L) + fiX*t + 2*w,2)/(8.*fiX*fiY); 
							 Dtrial(1,0) += ((fiY*L + 2*fiX*t - 2*w)*pow(-(fiY*L) + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
							 Dtrial(2,0) += -(pow(-(fiY*L) + fiX*t + 2*w,2)*(2*fiY*L + fiX*t + 2*w))/(48.*fiX*pow(fiY,2)); 
							 Dtrial(0,1) += ((fiY*L + 2*fiX*t - 2*w)*pow(-(fiY*L) + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
							 Dtrial(1,1) += (pow(fiY,4)*pow(L,4) + 3*pow(fiX,4)*pow(t,4) - 8*pow(fiY,3)*pow(L,3)*w + 8*pow(fiX,3)*pow(t,3)*w + 24*pow(fiY,2)*pow(L,2)*pow(w,2) + 16*pow(w,4) + 4*fiY*L*(3*pow(fiX,3)*pow(t,3) - 8*pow(w,3)))/(192.*pow(fiX,3)*fiY); 
							 Dtrial(2,1) += -(pow(-(fiY*L) + fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) - 4*fiY*L*w + 4*fiX*t*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
							 Dtrial(0,2) += -(pow(-(fiY*L) + fiX*t + 2*w,2)*(2*fiY*L + fiX*t + 2*w))/(48.*fiX*pow(fiY,2)); 
							 Dtrial(1,2) += -(pow(-(fiY*L) + fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) - 4*fiY*L*w + 4*fiX*t*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
							 Dtrial(2,2) += (3*pow(fiY,4)*pow(L,4) + 4*pow(fiY,3)*pow(L,3)*(3*fiX*t - 2*w) + pow(fiX*t + 2*w,4))/(192.*fiX*pow(fiY,3)); 
						}

					} else {
						if (inside2<0) {
							if (fiY<0.0) {
								// case 2c
								 Dtrial(0,0) += pow(fiY*L + fiX*t + 2*w,2)/(8.*fiX*fiY); 
								 Dtrial(1,0) += -((fiY*L - 2*fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
								 Dtrial(2,0) += -((-2*fiY*L + fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
								 Dtrial(0,1) += -((fiY*L - 2*fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
								 Dtrial(1,1) += (pow(fiY*L + fiX*t + 2*w,2)*(pow(fiY,2)*pow(L,2) - 2*fiX*fiY*L*t + 3*pow(fiX,2)*pow(t,2) + 4*fiY*L*w - 4*fiX*t*w + 4*pow(w,2)))/(192.*pow(fiX,3)*fiY); 
								 Dtrial(2,1) += -(pow(fiY*L + fiX*t + 2*w,2)*(3*pow(fiY*L - fiX*t,2) + 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
								 Dtrial(0,2) += -((-2*fiY*L + fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
								 Dtrial(1,2) += -(pow(fiY*L + fiX*t + 2*w,2)*(3*pow(fiY*L - fiX*t,2) + 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
								 Dtrial(2,2) += (pow(fiY*L + fiX*t + 2*w,2)*(3*pow(fiY,2)*pow(L,2) - 2*fiY*L*(fiX*t + 2*w) + pow(fiX*t + 2*w,2)))/(192.*fiX*pow(fiY,3)); 

							} else {
								// case 3c
									 Dtrial(0,0) += L*t - pow(fiY*L + fiX*t + 2*w,2)/(8.*fiX*fiY); 
									 Dtrial(1,0) += ((fiY*L - 2*fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(2,0) += ((-2*fiY*L + fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(0,1) += ((fiY*L - 2*fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(1,1) += -(pow(fiY,4)*pow(L,4) + 3*pow(fiX,4)*pow(t,4) + 8*pow(fiY,3)*pow(L,3)*w + 8*pow(fiX,3)*pow(t,3)*w + 24*pow(fiY,2)*pow(L,2)*pow(w,2) + 16*pow(w,4) + 4*fiY*L*(-3*pow(fiX,3)*pow(t,3) + 8*pow(w,3)))/(192.*pow(fiX,3)*fiY); 
									 Dtrial(2,1) += (pow(fiY*L + fiX*t + 2*w,2)*(3*pow(fiY*L - fiX*t,2) + 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(0,2) += ((-2*fiY*L + fiX*t + 2*w)*pow(fiY*L + fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(1,2) += (pow(fiY*L + fiX*t + 2*w,2)*(3*pow(fiY*L - fiX*t,2) + 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(2,2) += -(3*pow(fiY,4)*pow(L,4) + 4*pow(fiY,3)*pow(L,3)*(-3*fiX*t + 2*w) + pow(fiX*t + 2*w,4))/(192.*fiX*pow(fiY,3)); 

							}

						} else {
							if (inside1>0) {
								if (fiY>0.0) {
									// case 2b
									Dtrial(0,0) += pow (fiY*L + fiX*t - 2*w,2)/(8.*fiX*fiY); 
									Dtrial(1,0) += (pow(fiY*L + fiX*t - 2*w,2)*(fiY*L - 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									Dtrial(2,0) += ((-2*fiY*L + fiX*t - 2*w)*pow(fiY*L + fiX*t - 2*w,2))/(48.*fiX*pow(fiY,2)); 
									Dtrial(0,1) += (pow(fiY*L + fiX*t - 2*w,2)*(fiY*L - 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									Dtrial(1,1) += (pow(fiY*L + fiX*t - 2*w,2)*(pow(fiY,2)*pow(L,2) + 3*pow(fiX,2)*pow(t,2) + 4*fiX*t*w + 4*pow(w,2) - 2*fiY*L*(fiX*t + 2*w)))/(192.*pow(fiX,3)*fiY); 
									Dtrial(2,1) += -(pow(fiY*L + fiX*t - 2*w,2)*(3*pow(fiY*L - fiX*t,2) - 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									Dtrial(0,2) += ((-2*fiY*L + fiX*t - 2*w)*pow(fiY*L + fiX*t - 2*w,2))/(48.*fiX*pow(fiY,2)); 
									Dtrial(1,2) += -(pow(fiY*L + fiX*t - 2*w,2)*(3*pow(fiY*L - fiX*t,2) - 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									Dtrial(2,2) += (pow(fiY*L + fiX*t - 2*w,2)*(3*pow(fiY,2)*pow(L,2) + pow(fiX*t - 2*w,2) + fiY*(-2*fiX*L*t + 4*L*w)))/(192.*fiX*pow(fiY,3)); 


								} else {
									// case 3b
									 Dtrial(0,0) += L*t - pow(fiY*L + fiX*t - 2*w,2)/(8.*fiX*fiY); 
									 Dtrial(1,0) += (pow(fiY*L + fiX*t - 2*w,2)*(-(fiY*L) + 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(2,0) += (pow(fiY*L + fiX*t - 2*w,2)*(2*fiY*L - fiX*t + 2*w))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(0,1) += (pow(fiY*L + fiX*t - 2*w,2)*(-(fiY*L) + 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(1,1) += -(pow(fiY,4)*pow(L,4) + 3*pow(fiX,4)*pow(t,4) - 8*pow(fiY,3)*pow(L,3)*w - 8*pow(fiX,3)*pow(t,3)*w + 24*pow(fiY,2)*pow(L,2)*pow(w,2) + 16*pow(w,4) - 4*fiY*L*(3*pow(fiX,3)*pow(t,3) + 8*pow(w,3)))/(192.*pow(fiX,3)*fiY); 
									 Dtrial(2,1) += (pow(fiY*L + fiX*t - 2*w,2)*(3*pow(fiY*L - fiX*t,2) - 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(0,2) += (pow(fiY*L + fiX*t - 2*w,2)*(2*fiY*L - fiX*t + 2*w))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(1,2) += (pow(fiY*L + fiX*t - 2*w,2)*(3*pow(fiY*L - fiX*t,2) - 4*(fiY*L + fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(2,2) += -(3*pow(fiY,4)*pow(L,4) + pow(fiX*t - 2*w,4) - 4*pow(fiY,3)*pow(L,3)*(3*fiX*t + 2*w))/(192.*fiX*pow(fiY,3)); 
								}

							} else {
								if (fiY<0.0) {
									// case 2d
									 Dtrial(0,0) += -pow(fiY*L - fiX*t + 2*w,2)/(8.*fiX*fiY); 
									 Dtrial(1,0) += (pow(fiY*L - fiX*t + 2*w,2)*(fiY*L + 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(2,0) += -((2*fiY*L + fiX*t - 2*w)*pow(fiY*L - fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(0,1) += (pow(fiY*L - fiX*t + 2*w,2)*(fiY*L + 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(1,1) += -(pow(fiY*L - fiX*t + 2*w,2)*(pow(fiY,2)*pow(L,2) + 3*pow(fiX,2)*pow(t,2) + 4*fiX*t*w + 4*pow(w,2) + 2*fiY*L*(fiX*t + 2*w)))/(192.*pow(fiX,3)*fiY); 
									 Dtrial(2,1) += (pow(fiY*L - fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) + 4*(fiY*L - fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(0,2) += -((2*fiY*L + fiX*t - 2*w)*pow(fiY*L - fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(1,2) += (pow(fiY*L - fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) + 4*(fiY*L - fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(2,2) += -((3*pow(fiY,2)*pow(L,2) + 2*fiY*L*(fiX*t - 2*w) + pow(fiX*t - 2*w,2))*pow(fiY*L - fiX*t + 2*w,2))/(192.*fiX*pow(fiY,3)); 
								} else {
									// case 3d
									 Dtrial(0,0) += L*t + pow(fiY*L - fiX*t + 2*w,2)/(8.*fiX*fiY); 
									 Dtrial(1,0) += -(pow(fiY*L - fiX*t + 2*w,2)*(fiY*L + 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(2,0) += ((2*fiY*L + fiX*t - 2*w)*pow(fiY*L - fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(0,1) += -(pow(fiY*L - fiX*t + 2*w,2)*(fiY*L + 2*(fiX*t + w)))/(48.*pow(fiX,2)*fiY); 
									 Dtrial(1,1) += (pow(fiY,4)*pow(L,4) + 3*pow(fiX,4)*pow(t,4) + 8*pow(fiY,3)*pow(L,3)*w - 8*pow(fiX,3)*pow(t,3)*w + 24*pow(fiY,2)*pow(L,2)*pow(w,2) + 16*pow(w,4) + 4*fiY*L*(3*pow(fiX,3)*pow(t,3) + 8*pow(w,3)))/(192.*pow(fiX,3)*fiY); 
									 Dtrial(2,1) += -(pow(fiY*L - fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) + 4*(fiY*L - fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(0,2) += ((2*fiY*L + fiX*t - 2*w)*pow(fiY*L - fiX*t + 2*w,2))/(48.*fiX*pow(fiY,2)); 
									 Dtrial(1,2) += -(pow(fiY*L - fiX*t + 2*w,2)*(3*pow(fiY*L + fiX*t,2) + 4*(fiY*L - fiX*t)*w - 4*pow(w,2)))/(384.*pow(fiX,2)*pow(fiY,2)); 
									 Dtrial(2,2) += (3*pow(fiY,4)*pow(L,4) + pow(fiX*t - 2*w,4) + 4*pow(fiY,3)*pow(L,3)*(3*fiX*t + 2*w))/(192.*fiX*pow(fiY,3)); 
								}
							}
						}
					}

				} else {
					if (fiX<0.0) {
						// case 4a,c: both out from different sides, fiX>=0
						 Dtrial(0,0) += L*(t/2. + w/fiX); 
						 Dtrial(1,0) += -(L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(2,0) += (fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(0,1) += -(L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(1,1) += (L*(pow(fiX,3)*pow(t,3) + 2*pow(fiY,2)*pow(L,2)*w + 8*pow(w,3)))/(24.*pow(fiX,3)); 
						 Dtrial(2,1) += -(fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(0,2) += (fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(1,2) += -(fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(2,2) += (pow(L,3)*(fiX*t + 2*w))/(24.*fiX); 

					} else {
						// case 4b,d: both out from different sides, fiX<0
						 Dtrial(0,0) += (L*(t - (2*w)/fiX))/2.; 
						 Dtrial(1,0) += (L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(2,0) += -(fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(0,1) += (L*(pow(fiY,2)*pow(L,2) - 3*pow(fiX,2)*pow(t,2) + 12*pow(w,2)))/(24.*pow(fiX,2)); 
						 Dtrial(1,1) += (L*(pow(fiX,3)*pow(t,3) - 2*pow(fiY,2)*pow(L,2)*w - 8*pow(w,3)))/(24.*pow(fiX,3)); 
						 Dtrial(2,1) += (fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(0,2) += -(fiY*pow(L,3))/(12.*fiX); 
						 Dtrial(1,2) += (fiY*pow(L,3)*w)/(12.*pow(fiX,2)); 
						 Dtrial(2,2) += (pow(L,3)*(fiX*t - 2*w))/(24.*fiX); 

					}
				}
			}
		}
	}


	
    Dtrial *= k;
	Dtrial(3,3) = kg*J;

	Dtrial(3,3)*= torsionalStiffnessFactor;



    if (stronger) {		
		Dtrial(0,0) += factorStronger*k*t*L;
        Dtrial(1,1) += factorStronger*k*L*t*t*t /12.0;
        Dtrial(2,2) += factorStronger*k*L*L*L*t /12.0;
		Dtrial(3,3) += factorStronger*kg*J;
	} else {
		// fake stiffness instead of zero
		
		double factor = 0.000001; // factorStronger;  // 0.0
		if  (sqrt(pow(s(0),2) + pow(s(1),2) + pow(s(1),2))<DBL_EPSILON) {
		//if  (abs(Dtrial(0,0))<DBL_EPSILON && abs(Dtrial(1,1))<DBL_EPSILON) {
			
			//Dtrial.Zero();
			//s.Zero();

			Dtrial(0,0) += factor *k* t*L;
			Dtrial(1,1) += factor *k* L*t*t*t /12.0;
			Dtrial(2,2) += factor *k* L*L*L*t /12.0;
			//Dtrial(3,3) += factorStronger *kg*J;
			

		// Dtrial = Dcommitted;
			
	     }
	}


	if (spandrel) {
		//s(0) -= factorStronger*k*t*L *e(0);
        //s(1) -= factorStronger*k*L*t*t*t /12.0 * e(1);
        //s(2) -= factorStronger*k*L*L*L*t /12.0 * e(2);	

		//Dtrial = D;
	}
	//opserr << "k= " << k << endln;
	//opserr << "fc= " << fc << endln;

	//opserr << "new tangent: " << Dtrial; 


	return 0;
}

const Vector &
NoTensionSection3d::getSectionDeformation (void)
{
    return e;
}

const Vector &
NoTensionSection3d::getStressResultant (void)
{	        if (abs(s(0))>-DBL_EPSILON && 
				abs(s(1))>-DBL_EPSILON &&
				abs(s(2))>-DBL_EPSILON )		
			return s;
		else
			return s*0.0;

  return s;
}

const Matrix &
NoTensionSection3d::getSectionTangent(void)
{
			if (abs(Dtrial(0,0))>-DBL_EPSILON && 
				abs(Dtrial(0,1))>-DBL_EPSILON &&
				abs(Dtrial(0,2))>-DBL_EPSILON &&
				
				abs(Dtrial(1,0))>-DBL_EPSILON &&
				abs(Dtrial(1,1))>-DBL_EPSILON &&
				abs(Dtrial(1,2))>-DBL_EPSILON &&
				
				abs(Dtrial(2,0))>-DBL_EPSILON &&
				abs(Dtrial(2,1))>-DBL_EPSILON &&
				abs(Dtrial(2,2))>-DBL_EPSILON )		
			return Dtrial;
		else
			return D;

}

const Matrix &
NoTensionSection3d::getInitialTangent(void)
{  
  return D;
}


SectionForceDeformation*
NoTensionSection3d::getCopy ()
{
    // Make a copy of the hinge
    NoTensionSection3d *theCopy =
		new NoTensionSection3d (this->getTag(), k, kg, L, t, J, fc, nSections, stronger, elastic, crushing, spandrel);
    theCopy->parameterID = parameterID;

    return theCopy;
}

const ID&
NoTensionSection3d::getType ()
{
    return code;
}

int
NoTensionSection3d::getOrder () const
{
    return 4;
}

int
NoTensionSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(6);

    int dataTag = this->getDbTag();
    
	data(0) = this->getTag();
    data(1) = k;
    data(2) = kg;    
    data(3) = L;
    data(4) = t;
    data(5) = J;
    
    res += theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection3d::sendSelf -- failed to send data\n";
      return res;
    }
    
    return res;
}

int
NoTensionSection3d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
	static Vector data(6);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection3d::recvSelf -- failed to receive data\n";
      return res;
    }

	this->setTag((int)data(0));
    k = data(1);
    kg = data(2);    
    L = data(3);
    t = data(4);
    J = data(5);

    return res;
}
 
void
NoTensionSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

  } else {
    s << "NoTensionSection3d, tag: " << this->getTag() << endln;
    s << "\t k: " << k << endln;
    s << "\t kg: " << kg << endln;
    s << "\t L: " << L << endln;
    s << "\t t: " << t << endln;
    s << "\t J: " << J << endln;
  }
}

int
NoTensionSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"k") == 0) {
    param.setValue(k);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"kg") == 0) {
    param.setValue(kg);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"L") == 0) {
    param.setValue(L);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"t") == 0) {
    param.setValue(t);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"J") == 0) {
    param.setValue(J);
    return param.addObject(5, this);
  }
  return -1;
}

int
NoTensionSection3d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    k = info.theDouble;
  if (paramID == 2)
    kg = info.theDouble;
  if (paramID == 3)
    L = info.theDouble;
  if (paramID == 4)
    t = info.theDouble;
  if (paramID == 5)
    J = info.theDouble;

  return 0;
}

int
NoTensionSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
NoTensionSection3d::getStressResultantSensitivity(int gradIndex,
						bool conditional)
{
  s.Zero();

  if (parameterID == 1) { // k
    s(0) = L*t*e(0);
    s(1) = L*t*t*t/12.*e(1);
    s(2) = L*L*L*t/12.*e(2);
  }
  if (parameterID == 2) // kg
    s(3) = J*e(3);
  if (parameterID == 3) // L
    s(0) = k* t*e(0);
    s(1) = k* t*t*t /12.0 * e(1);
    s(2) = k* 3.*L*L*t /12.0 * e(2);
  if (parameterID == 4) // t
    s(0) = k* L *e(0);
    s(1) = k* L*3.*t*t /12.0 * e(1);
    s(2) = k* L*L*L /12.0 * e(2);
  if (parameterID == 5) // G
    s(3) = J*e(3);

  return s;
}

const Matrix&
NoTensionSection3d::getInitialTangentSensitivity(int gradIndex)
{
  Dtrial.Zero();

  return Dtrial;
}



// recorder methods
Vector
NoTensionSection3d::getCrushingVariables(int sliceNum) {
	Vector crushingVarSet(4);

	crushingVarSet(0) = muX(sliceNum-1, 0);
	crushingVarSet(1) = zetaX(sliceNum-1, 0);
	crushingVarSet(2) = muX(sliceNum-1, 1);
	crushingVarSet(3) = zetaX(sliceNum-1, 1);

	return crushingVarSet;
}



Response*
NoTensionSection3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();

  Response *theResponse = 0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {

         theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    
         theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 

         theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  // state parameters
  } else if (strcmp(argv[0],"state") == 0 || strcmp(argv[0],"muZeta")==0 || strcmp(argv[0],"crushing") == 0) {
	  
	int sliceNum = nSections/2 +1;
	
	if (argc > 1) 
		sliceNum = atoi(argv[1]);
      
	if (sliceNum > 0 && sliceNum <= nSections) {
			   
		output.tag("SliceOutput");
		output.attr("number",sliceNum);

		sliceOutput = sliceNum;
		theResponse = new MaterialResponse(this, 5, this->getCrushingVariables(sliceNum) );
		output.endTag();
      
    } else if (sliceNum == 0) { // argv[1] was not an int, we want the central slice only
	
		sliceNum = nSections/2 +1;

		output.tag("SliceOutput");
		output.attr("number",sliceNum);

		sliceOutput = sliceNum;
		theResponse = new MaterialResponse(this, 5, this->getCrushingVariables(sliceNum) );
		output.endTag();

	}
   

  }

  output.endTag();
  return theResponse;
}


int 
NoTensionSection3d::getResponse(int responseID, Information &sectInfo)
{
  if (responseID == 5) { 
    Vector data(4);
		return sectInfo.setVector(this->getCrushingVariables(sliceOutput));	
  } else {  
		return SectionForceDeformation::getResponse(responseID, sectInfo);
  }
}