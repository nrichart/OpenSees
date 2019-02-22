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
                                                                        
// $Revision: 6421 $
// $Date: 2016-09-10 04:34:49 +0200 (Sat, 10 Sep 2016) $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/dispBeamColumn/DispBeamColumn3d.cpp $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for Macroelement3d.


#include "Macroelement3d.h"
#include "NoTensionSection3d.h"
#include "DamageShearInterface.h"
#include "GambarottaLagomarsinoModel.h"
#include "WrappedMaterial.h"

#include <Node.h>
#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Parameter.h>
#include <math.h>
#include <algorithm>
#include <elementAPI.h>
#include <string>
#include <TransientIntegrator.h>
#include <CompositeResponse.h>


Matrix Macroelement3d::K(18,18);
Vector Macroelement3d::P(18);
double Macroelement3d::workArea[200];

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numMyMacroelement = 0;

OPS_Export 	void * 
OPS_Macroelement3d()
{
	Vector intLength(3);
	intLength(0) = 1.0;
	intLength(1) = 1.0;
	intLength(2) = 1.0;

	Vector intLengthMasses(3);
	/*
	intLengthMasses(0) = 0.1666;
	intLengthMasses(1) = 0.6666;
	intLengthMasses(2) = 0.1666;
	*/

	intLengthMasses(0) = 0.25;
	intLengthMasses(1) = 0.50;
	intLengthMasses(2) = 0.25;

	Vector driftModelF(0);
	Vector driftModelF_ALR(0);
	Vector driftModelS(0);
	Vector driftModelS_ALR(0);

	double Ltfc = 0.;
	double alphaNC_SD = -4.0/3.0;
	double failureFactorF = 0.001;
	double failureFactorS = 0.001;
	double betaShearSpanF = 0.0;
	double betaShearSpanS = 0.0;

	
	  // print out a message about who wrote this element & any copyright info wanted
  if (numMyMacroelement  == 0) {
    opserr << "Macroelement3d - Written by Francesco Vanin, EPFL, 2018\n";
    numMyMacroelement++;
	//opserr << "Remaining arguments to read: " << OPS_GetNumRemainingInputArgs() << endln;
  }
  int remaining=OPS_GetNumRemainingInputArgs();

    if (remaining < 6) {
	opserr<<"insufficient arguments: eleTag, iNode, jNode, eNode, axisX,axisY,axisZ, oopX,oopY,oopZ -standard sectionTagI,sectionTagE,sectionTagJ,shearModelTag,h,E  <-mass mass> <-cmass> <-pDelta>   \n";
	opserr<<"or alternatively:       eleTag, iNode, jNode, eNode, axisX,axisY,axisZ, oopX,oopY,oopZ -tremuri h,b,t, E,G,fc,ft,GfI,tau0,mu,GfII  <-mass mass> <-cmass> <-pDelta>   \n";
	opserr << "Remaining arguments to read: " << remaining << endln;
	return 0;
    }

    // inputs: 
    int iData[4];
    int numData = 4;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING: invalid integer inputs, first part\n";
	return 0;
    }

    // inputs: 
    double dData[6];
    numData = 6;
    if(OPS_GetDoubleInput(&numData,&dData[0]) < 0) {
	opserr<<"WARNING: invalid double inputs, first part\n";
	return 0;
    }

    Vector axis(3); 
    axis(0) = dData[0];
    axis(1) = dData[1];
    axis(2) = dData[2];

	Vector oop(3);
    oop(0) = dData[3]; 
    oop(1) = dData[4]; 
    oop(2) = dData[5]; 

	SectionForceDeformation* theSectionI;
    SectionForceDeformation* theSectionE;
	SectionForceDeformation* theSectionJ;


    NDMaterial* theShearModel;
	NDMaterial* theShearModelOOP;
	//UniaxialMaterial*  theViscousModel = NULL;
    double E_, h; 
	
	bool gable = false;
	double b = 0.;
	double t = 0.;

    // options for input
	const char* inputStructure = OPS_GetString();
	if (strcmp(inputStructure,"-tremuri") == 0)  {
		// standard input, model fully equivalent to Tremuri (in-plane)
		// arguments: h,b,t, E,G, fc, mu,tau0, Gc,beta
		double dData2[10];
	    numData = 10;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-tremuri' input structure incorrect, invalid double input(s). \nRequired structure: h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>\n";
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		h = dData2[0];
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = dData2[8];
		double beta = dData2[9];

		Ltfc = abs(fc*b*t);
			
		theSectionI = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, dData2[5], 5, false, false, true);   // true for stronger, true for elastic, true for crushing
		
		theSectionE = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, dData2[5], 5, true, true, true);

		theSectionJ = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, dData2[5], 5, false, false, true);

		// Gambarotta Lagomarsino model for shear
		theShearModelOOP = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, true);  // true-false for elastic solution or not. The 5/6 factor is dropped to make it equivalent
		
		theShearModel    = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, false);


		//theShearModel = new DamageShearInterface(0, E_*b*t/h, G/(h)*b*t, 
		//	                                        tau0, mu, mu, Gc, 0.008*h, false);		


				
     } 

		else if (strcmp(inputStructure,"-tremuriIP") == 0) {
		// standard input, model fully equivalent to Tremuri (only in-plane response)
		// arguments: h,b,t, E,G, fc, mu,tau0, Gc,beta
		double dData2[10];
	    numData = 10;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-tremuri' input structure incorrect, invalid double input(s). \nRequired structure: h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>\n";
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		intLengthMasses = intLength;

		h = dData2[0];
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = dData2[8];
		double beta = dData2[9];

		Ltfc = abs(fc*b*t);
			
		theSectionI = new NoTensionSection3d(0, 100.*E_, 100.*G,  
			                                0.01*t, b, -1.0, 100.*dData2[5], 5, false, false, true);   // true for stronger, true for elastic, true for crushing
		
		theSectionE = new NoTensionSection3d(0, 100.*E_, 100.*G,  
			                                0.01*t, b, -1.0, 100.*dData2[5], 5, false, true, true);

		theSectionJ = new NoTensionSection3d(0, 100.*E_, 100.*G,  
			                                0.01*t, b, -1.0, 100.*dData2[5], 5, false, false, true);

		// Gambarotta Lagomarsino model for shear
		theShearModelOOP = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, true);  // true-false for elastic solution or not. The 5/6 factor is dropped to make it equivalent
		
		theShearModel    = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, false);

		//theShearModel = new DamageShearInterface(0, E_*b*t/h, G/(h)*b*t, 
		//	                                        tau0, mu, mu, Gc, 0.008*h, false);	
		
     } 

		else if (strcmp(inputStructure,"-fiberSection") == 0) {
		// input for the use of one or more fiber section models for the rocking behaviour, coupled to a standard shear model
	    // two-line input for standard rectangular prismatic elements-> line 1: section definition, line 2: macroelement definition 
	    // arguments: sectionI, sectionJ, sectionK, h,b,t, E,G, fc, mu,tau0, Gc,beta

		int iData2[3];
		int numData = 3;
		if (OPS_GetIntInput(&numData,&iData2[0]) < 0) {
			opserr<<"WARNING: '-fiberSection' input structure incorrect, invalid integer input(s). \nRequired structure: sectionI, sectionE, sectionJ, h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>\n";
			return 0;
		}

		
		// create copies of the sectional models
		int secTag = iData2[0];
	    theSectionI = OPS_GetSectionForceDeformation(secTag);
		if (theSectionI == 0) {
			opserr<<"SectionI: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[1];
		theSectionE = OPS_GetSectionForceDeformation(secTag);
		if (theSectionE == 0) {
			opserr<<"SectionE: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[2];
	    theSectionJ = OPS_GetSectionForceDeformation(secTag);
		if (theSectionJ == 0) {
			opserr<<"SectionJ: section model tagged " << secTag << " not found." << endln;
			return 0;
		}


		double dData2[10];
	    numData = 10;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-fiberSection' input structure incorrect, invalid double input(s). \nRequired structure: sectionI, sectionE, sectionJ, h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>\n";
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		h = dData2[0];
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = dData2[8];
		double beta = dData2[9];

		Ltfc = abs(fc*b*t);

		// Gambarotta Lagomarsino model for shear
		theShearModelOOP = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, true);  // true-false for elastic solution or not. The 5/6 factor is dropped to make it equivalent
		
		theShearModel    = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, false);

     } 



	else if (strcmp(inputStructure,"-fiberSectionShearModel1d") == 0) {
		// input for the use of one or more fiber section models for the rocking behaviour, plus a uniaxial shear model, uncoupled from the axial load
	    // arguments: sectionI, sectionJ, sectionK, ShearModel1dIP, ShearModel1dOOP

		int iData2[5];
		int numData = 5;
		if (OPS_GetIntInput(&numData,&iData2[0]) < 0) {
			opserr<<"WARNING: '-fiberSectionShearModel1d' input structure incorrect, invalid integer input(s). \nRequired structure: sectionI, sectionE, sectionJ, shearModel1dIP, shearModel1dOOP, alpha, h <-flags>\n";
			return 0;
		}

		
		// create copies of the sectional models
		int secTag = iData2[0];
	    theSectionI = OPS_GetSectionForceDeformation(secTag);
		if (theSectionI == 0) {
			opserr<<"SectionI: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[1];
		theSectionE = OPS_GetSectionForceDeformation(secTag);
		if (theSectionE == 0) {
			opserr<<"SectionE: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[2];
	    theSectionJ = OPS_GetSectionForceDeformation(secTag);
		if (theSectionJ == 0) {
			opserr<<"SectionJ: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		E_ = (theSectionE->getInitialTangent())(0,0);
		

		double dData2[2];
	    numData = 2;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-fiberSectionShearModel1d' input structure incorrect, invalid double input(s). \nRequired structure: sectionI, sectionE, sectionJ, shearModel1dIP, shearModel1dOOP, alpha, h <-flags>\n";
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		t = 0.;
		double alpha = dData2[0];
		h = dData2[1];
			

		// create pointers to uniaxial shear model
		UniaxialMaterial* pointerShear1dIP;
		UniaxialMaterial* pointerShear1dOOP;

		int matTag = iData2[3];
	    pointerShear1dIP = OPS_GetUniaxialMaterial(matTag);
		if (pointerShear1dIP== 0) {
			opserr<<"Shear1dIP: material model tagged " << secTag << " not found." << endln;
			return 0;
		}

		matTag = iData2[4];
	    pointerShear1dOOP = OPS_GetUniaxialMaterial(matTag);
		if (pointerShear1dOOP== 0) {
			opserr<<"Shear1dOOP: material model tagged " << secTag << " not found." << endln;
			return 0;
		}

			
		theShearModel     = new WrappedMaterial(0, E_, pointerShear1dIP, alpha);
		theShearModelOOP  = new WrappedMaterial(0, E_, pointerShear1dOOP);

		//delete [] pointerShear1dIP;
		//delete [] pointerShear1dOOP;

     } 

	
	else if (strcmp(inputStructure,"-pier") == 0)  {
		// standard input for piers, with my shear model and 3 inelastic sections
		// arguments: h, b, t, E, G, fc, mu, tau0, Gc, muR,GfIIB 

		double dData2[11];
	    numData = 11;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-pier' input structure incorrect, invalid double input(s). \nRequired structure: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR <-flags>\n";
			return 0;
		}

		intLength(0) = 1./6.;
		intLength(1) = 2./3.;
		intLength(2) = 1./6.;

		h = dData2[0];
		

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = 1.0 + dData2[8];
		double dropDrift = dData2[9];
		double muR = dData2[10];

		Ltfc = abs(fc*b*t);


		theSectionI = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, fc, 11, false, false, true);   // true for stronger, true for elastic, true for crushing
		
		theSectionE = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, fc, 11, false, false, true);

		theSectionJ = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, fc, 11, false, false, true);

		theShearModelOOP = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, 
			                                        tau0*b*t, mu, muR, Gc, dropDrift*h, true);

		theShearModel = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, 
			                                        tau0*b*t, mu, muR, Gc, dropDrift*h, false);									
	
     }
	
	
	else if (strcmp(inputStructure,"-spandrel") == 0) {
        // standard input for spandrels, with my shear model and the middle section linear elastic (end sections slightly stronger)
		// arguments: h, b, t, E, G, fc, mu, tau0, Gc, muR,GfIIB 

		double dData2[13];
	    numData = 13;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-spandrel' input structure incorrect, invalid double input(s). \nRequired structure: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR  <-flags>\n";
			return 0;
		}

		intLength(0) = 1./6.;
		intLength(1) = 2./3.;
		intLength(2) = 1./6.;

		h = dData2[0];
		

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = 1.0 + dData2[8];
		double dropDrift = dData2[9];
		double muR = dData2[10];

		Ltfc = abs(fc*b*t);


		theSectionI = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, fc, 5, true, false, false);   // true for stronger, true for elastic, true for crushing
		
		theSectionE = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, fc, 5, false, true, false);

		theSectionJ = new NoTensionSection3d(0, E_, G,  
			                                t, b, -1.0, fc, 5, true, false, false);

		theShearModelOOP = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, 
			                                        tau0*b*t, mu, muR, Gc, dropDrift*h, true);

		theShearModel = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, 
			                                        tau0*b*t, mu, muR, Gc, dropDrift*h, false);							
		
     } 

	else if (strcmp(inputStructure,"-gable") == 0)  {

        // standard input for gables, with my shear model
		// arguments: h, b, t, E, G, fc, mu, tau0, Gc, muR,GfIIB 
		gable = true;

		double dData2[11];
	    numData = 11;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr<<"WARNING Macroelement3d:'-gable' input structure incorrect, invalid double input(s). \nRequired structure: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR <-flags>\n";
			return 0;
		}

		intLength(0) = 1./6.;
		intLength(1) = 2./3.;
		intLength(2) = 1./6.;

		h = dData2[0];
		

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = 1.0 + dData2[8];
		double dropDrift = dData2[9];
		double muR = dData2[10];

		Ltfc = abs(fc*0.5*b*t);


		theSectionI = new NoTensionSection3d(0, E_, 0.5*G,  
			                                t, b, -1.0, fc, 5, false, false, true);   // true for stronger, true for elastic, true for crushing
		
		theSectionE = new NoTensionSection3d(0, E_, G,  
			                                t, 0.5*b, -1.0, fc, 5, false, false, true);

		theSectionJ = new NoTensionSection3d(0, E_, G,  
			                                t, 0.1*b, -1.0, fc, 5, false, false, true);

		theShearModelOOP = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, 
			                                        tau0*b*t, mu, muR, Gc, dropDrift*h, true);

		theShearModel = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, 
			                                        tau0*b*t, mu, muR, Gc, dropDrift*h, true);	
     }
	
	
	else {
		// standard input, opensees like. Custom sectional models have to be created before, and their tags are passed to the macrolement.
		// arguments: sectionTag, sectionTagE, shearModelTag, h, E 
		// inputs: 
		int iData2[5];
		int numData = 5;
		if (OPS_GetIntInput(&numData,&iData2[0]) < 0) {
			opserr<<"WARNING: invalid integer inputs, sectional models. \nRequired structure: sectionItag, sectionEtag, sectionJtag, shearModelIP, shearModelOOP, h, E <-flags>\n";
			return 0;
		}

		numData = 1;
		if (OPS_GetDoubleInput(&numData,&h) < 0) {
			opserr<<"WARNING: invalid double inputs, sectional models\n";
			return 0;
		}

		numData = 1;
		if (OPS_GetDoubleInput(&numData,&E_) < 0) {
			opserr<<"WARNING: invalid double inputs, sectional models\n";
			return 0;
		}


		// create copies of the sectional models
		int secTag = iData2[0];
	    theSectionI = OPS_GetSectionForceDeformation(secTag);
		if (theSectionI == 0) {
			opserr<<"SectionI: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[1];
		theSectionE = OPS_GetSectionForceDeformation(secTag);
		if (theSectionE == 0) {
			opserr<<"SectionE: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[2];
	    theSectionJ = OPS_GetSectionForceDeformation(secTag);
		if (theSectionJ == 0) {
			opserr<<"SectionJ: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		int matTag = iData2[3];
		theShearModel = OPS_GetNDMaterial(matTag);
		if (theShearModel == 0) {
			opserr<<"Shear model IP: NDMaterial " << matTag << " not found." << endln;
			return 0;
		}

		int matTagOOP = iData2[4];
		theShearModelOOP = OPS_GetNDMaterial(matTagOOP);
		if (theShearModelOOP == 0) {
			opserr<<"Shear model OOP: NDMaterial " << matTag << " not found." << endln;
			return 0;
		}

		t=0.0;

	 }
	 

    double mass = 0.0;
	double weights[3];
	double massGlobal[3];
	int interfaceDampingTag = 0;
    int cmass = 0;
    int PDelta = 0;
	Vector massDir(3);
	massDir(0) = massDir(1) = massDir(2) = 1.;

	bool read_Ltfc=false;
	if (Ltfc==0) { read_Ltfc=true;}
	bool alreadReadType = false;
	const char* type = NULL;

    
    while(OPS_GetNumRemainingInputArgs() > 0) {
	if (!alreadReadType)	type = OPS_GetString();

	if (strcmp(type,"-cMass") == 0) {
	    cmass = 1;
		alreadReadType = false;
	} 
	else if (strcmp(type,"-mass") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&mass) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid mass\n";
				return 0;
			}
			alreadReadType = false;
	    }
	} 
	else if (strcmp(type,"-AxialCollapseRatio") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&alphaNC_SD) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid AxialCollapseRatio\n";
			}
			if (alphaNC_SD < 1) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", AxialCollapseRatio must be bigger than 1.0 (set equal to 1.0)\n";
				alphaNC_SD = 1.0;
			}
			alreadReadType = false;
	    }
	}

	else if (strcmp(type,"-failureFactor") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&failureFactorF) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid failureFactor\n";
			}
			if ((failureFactorF > 1.0) || (failureFactorF<0.0)) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", failureFactor must be between 0.0 and 1.0 (set equal to 0.001)\n";
				failureFactorF = 0.001;
			}
			failureFactorS = failureFactorF;
			alreadReadType = false;
	    }
	}

	else if (strcmp(type,"-failureFactorFlexure") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&failureFactorF) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid failureFactor\n";
			}
			if ((failureFactorF > 1.0) || (failureFactorF<0.0)) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", failureFactorFlexure must be between 0.0 and 1.0 (set equal to 0.001)\n";
				failureFactorF = 0.001;
			}
			alreadReadType = false;
	    }
	}

     else if (strcmp(type,"-failureFactorShear") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&failureFactorS) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid failureFactor\n";
			}
			if ((failureFactorS > 1.0) || (failureFactorS<0.0)) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", failureFactorShear must be between 0.0 and 1.0 (set equal to 0.001)\n";
				failureFactorS = 0.001;
			}
			alreadReadType = false;
	    }
	}


	else if (strcmp(type,"-density") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&mass) < 0) {
				opserr << "WARNING: Macroelement " << iData[0] <<", invalid mass.\n";
				return 0;
			}
			mass *= b*t; 
			if (b*t<DBL_EPSILON)
				opserr << "WARNING: Macroelement " << iData[0] <<" has zero area, use -mass flag instead of -density.\n";
	    }
		alreadReadType = false;
	} 
	else if (strcmp(type,"-PDelta") == 0  || strcmp(type,"-pDelta") == 0) {
             PDelta = 1;
			 //opserr << "Including PDelta effects\n";
			 alreadReadType = false;
	       }
	       
	/*
	else if (strcmp(type,"-dampingModel") == 0  || strcmp(type,"-OOPdamping") == 0) {
	        if (OPS_GetNumRemainingInputArgs() > 0) {
			   numData = 1;
			   if (OPS_GetIntInput(&numData,&interfaceDampingTag) < 0) {
				   opserr << "WARNING: Macroelement " << iData[0] <<", invalid damping model.\n";
				   return 0;			   
			   }
			   opserr << "WARNING: Macroelement " << iData[0] <<", damping model not yet implemented.\n";
			  // theViscousModel = OPS_GetUniaxialMaterial(interfaceDampingTag);
			  // if (theViscousModel == 0) {
			  //	opserr<<"WARNING: Macroelement " << iData[0] <<", invalid damping model.\n" << endln;
			  //		return 0;
			  //  }
			   
			}
	     }
		 */

	else  if (strcmp(type,"-intWeights") == 0) {
              if(OPS_GetNumRemainingInputArgs() > 2) {
				  numData = 3;
				  if (OPS_GetDoubleInput(&numData,&weights[0]) < 0) {
					  opserr << "WARNING: Macroelement " << iData[0] <<", invalid integration weights, standard weights are applied.\n";
				  } else {
					  intLength(0) = weights[0];
					  intLength(1) = weights[1];
					  intLength(2) = weights[2];
					  // maybe check that they sum to 1
					  double sum = 0.;
					  sum = intLength(0) + intLength(1) + intLength(2);
					  if (abs(sum-1.0) < 0.01)
						  opserr << "WARNING: Macroelement " << iData[0] <<", specified integration weights do not sum exactly to 1.\n";
				  }
			  }
			  alreadReadType = false;
	       }
		else  if (strcmp(type,"-massDir") == 0) {
              if (OPS_GetNumRemainingInputArgs() > 2) {
				  numData = 3;
				  if (OPS_GetDoubleInput(&numData,&massGlobal[0]) < 0) {
					  opserr << "WARNING: Macroelement " << iData[0] <<", invalid mass directions, mass applied in all global directions.\n";
				  } else {
					  massDir(0) = massGlobal[0];
					  massDir(1) = massGlobal[1];
					  massDir(2) = massGlobal[2];
					  // check that they are not bigger than 1, but allow the choice anyway
					  if ( massDir(0)>1.0 || massDir(0)<0.0 || massDir(1)>1.0 || massDir(1)<0.0 || massDir(2)>1.0 || massDir(2)<0.0 )
						  opserr << "WARNING: Macroelement " << iData[0] <<", specified mass directions are negative or bigger than one.\n";
				  }
			  }
			  alreadReadType = false;
	       }
		//drift flexure
		else  if (strcmp(type,"-driftFlexure") == 0) {
			//opserr << "remaining args: " << OPS_GetNumRemainingInputArgs() << endln;
			int otherInputs=OPS_GetNumRemainingInputArgs();
			std::vector<double> parsed;

			bool goOn = true;
			int i=0;
			int numData = 1;
			double singleInput;

			const char* stringRead = NULL;
			char* pEnd;
			
			while (goOn && i<otherInputs) {
				stringRead = OPS_GetString();
				if (strstr(stringRead,"-") ) { // means there is a flag
					goOn = false;					
				} else {
					singleInput = strtod(stringRead, &pEnd);
					parsed.push_back(singleInput);
					i += 1;
				}

			}

			if (read_Ltfc) {
				Ltfc = parsed[parsed.size()-1];
				parsed.pop_back();
			}

			betaShearSpanF = parsed[parsed.size()-1];
			parsed.pop_back();

			// if only one value is given, constant drfit model
			if (parsed.size() ==1) {
				double constantDrift  = parsed[0];
				parsed[0] = 0.0;
				parsed.push_back(constantDrift);
				parsed.push_back(1.0);
				parsed.push_back(constantDrift);
			}

			// check that the law is defined between ALR=0 and ALR = 1, if not extend it (constant up to the next defined point, or to the end)
			if (parsed[0]>0.0) {
				parsed.insert(parsed.begin(), parsed[1]);
				parsed.insert(parsed.begin(), 0.0);
			}

			if (parsed[parsed.size()-2]<1.0) {
				parsed.push_back(1.0),
				parsed.push_back(parsed[parsed.size()-2]);
			}

			int nPoints = (int) (parsed.size()/2);
			driftModelF_ALR.resize(nPoints);
			driftModelF.resize(nPoints);

			for (int i=0; i<nPoints; i++) {
				driftModelF_ALR(i) = parsed[2*i];
				driftModelF(i)     = parsed[2*i+1];
			}

			alreadReadType = true;
			type = stringRead;

			//opserr << "driftF_ALR = " << driftModelF_ALR;
			//opserr << "driftF = " << driftModelF;
		}

		//drift shear
		else  if (strcmp(type,"-driftShear") == 0) {
            int otherInputs=OPS_GetNumRemainingInputArgs();
			std::vector<double> parsedS;

			bool goOn = true;
			int i=0;
			int numData = 1;
			double singleInput;

			const char* stringRead = NULL;
			char* pEnd;
			
			while (goOn && i<otherInputs) {
				stringRead = OPS_GetString();
				if (strstr(stringRead,"-") ) { // means there is a flag
					goOn = false;
				} else {
					singleInput = strtod(stringRead, &pEnd);
					parsedS.push_back(singleInput);
					i += 1;
				}

			}

			if (read_Ltfc) {
				Ltfc = parsedS[parsedS.size()-1];
				parsedS.pop_back();
			}

			betaShearSpanS = parsedS[parsedS.size()-1];
			parsedS.pop_back();

			// if only one value is given, constant drfit model
			if (parsedS.size() ==1) {
				double constantDrift  = parsedS[0];
				parsedS[0] = 0.0;
				parsedS.push_back(constantDrift);
				parsedS.push_back(1.0);
				parsedS.push_back(constantDrift);
			}

			// check that the law is defined between ALR=0 and ALR = 1, if not extend it (constant up to the next defined point, or to the end)
			if (parsedS[0]>0.0) {
				parsedS.insert(parsedS.begin(), parsedS[1]);
				parsedS.insert(parsedS.begin(), 0.0);
			}

			if (parsedS[parsedS.size()-2]<1.0) {
				parsedS.push_back(1.0),
				parsedS.push_back(parsedS[parsedS.size()-2]);
			}

			
			int nPoints = (int) (parsedS.size()/2);
			driftModelS_ALR.resize(nPoints);
			driftModelS.resize(nPoints);

			for (int i=0; i<nPoints; i++) {
				driftModelS_ALR(i) = parsedS[2*i];
				driftModelS(i)     = parsedS[2*i+1];
			}

			alreadReadType = true;
			type = stringRead;

			//opserr << "driftS_ALR = " << driftModelS_ALR;
			//opserr << "driftS = " << driftModelS;
	       }
	
	}  // end while other inputs

	//opserr << "alphaNC= " << alphaNC_SD << endln;

	 Element *theEle =  new Macroelement3d(iData[0],iData[1],iData[2],iData[3],theSectionI,theSectionE,theSectionJ,theShearModel,theShearModelOOP,h,E_,
		                                    driftModelF,driftModelF_ALR, driftModelS, driftModelS_ALR, Ltfc, alphaNC_SD, betaShearSpanF, betaShearSpanS, failureFactorF, failureFactorS,
											axis,oop, intLength, intLengthMasses, massDir, PDelta,mass,cmass, gable);

	 if ((strcmp(inputStructure,"-tremuri")==0) || (strcmp(inputStructure,"-tremuriIP")==0)  || (strcmp(inputStructure,"-pier")==0) || (strcmp(inputStructure,"-spandrel")==0) || (strcmp(inputStructure,"-gable")==0) )  {
		// opserr << "input structure " << inputStructure << endln;
	    if (theSectionI!=NULL)  delete [] theSectionI;
	    if (theSectionE!=NULL)  delete [] theSectionE;
	    if (theSectionJ!=NULL)  delete [] theSectionJ;
	    if (theShearModel!=NULL)  delete [] theShearModel;
	    if (theShearModelOOP!=NULL)  delete [] theShearModelOOP;
	 }
	 
     return theEle;
}


Macroelement3d::Macroelement3d(int tag, int nd1, int nd2, int ndE, 
							   SectionForceDeformation *sI, SectionForceDeformation *sE, SectionForceDeformation *sJ, 
							   NDMaterial* shearModel, NDMaterial* shearModelOOP, double h, double E_, 
							   Vector driftF, Vector driftF_ALR, Vector driftS, Vector driftS_ALR, double Ltfc, double alphaNC_SD, double betaShearSpanF, double betaShearSpanS, double failureFactorF, double failureFactorS,
							   Vector axis, Vector oop, Vector _intLength, Vector _intLengthMasses, Vector massDir, 
							   int PDelta, double rho, int cm, double _isGable)
   :Element (tag, ELE_TAG_DispBeamColumn3d),
   numSections(3), numShearModels(2), theSections(0), theShearModel(0), connectedExternalNodes(3), Q(18), q(12), uBasic(12), uBasicCommitted(12), Tgl(18,18), GammaC(12,18), Tgl6(6,6), rho(rho), cMass(cm), 
   parameterID(0), E(E_), intLength(_intLength), intLengthMasses(_intLengthMasses),  massglobalDir(massDir), driftF(0.), driftS(0.), 
   limitDriftF(driftF), limitDriftS(driftS), limitDriftF_ALR(driftF_ALR), limitDriftS_ALR(driftS_ALR), Ltfc(Ltfc), alphaNC_SD(alphaNC_SD), betaShearSpanF(betaShearSpanF), betaShearSpanS(betaShearSpanS), failureFactorF(failureFactorF), failureFactorS(failureFactorS),
   failedF(false), failedS(false), failedFcommitted(false), failedScommitted(false), 
   collapsedF(false), collapsedS(false), collapsedFcommitted(false), collapsedScommitted(false), 
   isGable(_isGable), wx(0.0), wy(0.0), wz(0.0), nodeIInitialDisp(0), nodeJInitialDisp(0), nodeIOffset(0), nodeJOffset(0), 
   xAxis(3), yAxis(3), zAxis(3), L(h/2.0), PDelta(PDelta), deltaW1(0.), deltaV1(0.), deltaW3(0.), deltaV3(0.), committedTime(0.0)
{

  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
  if (theSections == 0) {
    opserr << "Macroelement3d::Macroelement3d - failed to allocate section model pointer\n";
    exit(-1);
  }
  
  theSections[0] = sI->getCopy();       if (theSections[0] == 0) { opserr << "Macroelement3d::Macroelement3d -- failed to get a copy of section model\n";  exit(-1);  }
  theSections[1] = sE->getCopy();       if (theSections[1] == 0) { opserr << "Macroelement3d::Macroelement3d -- failed to get a copy of section model\n";  exit(-1);  }
  theSections[2] = sJ->getCopy();       if (theSections[2] == 0) { opserr << "Macroelement3d::Macroelement3d -- failed to get a copy of section model\n";  exit(-1);  }


  theShearModel = new NDMaterial *[numShearModels];
  theShearModel[0] = shearModel->getCopy("BeamFiber");
  theShearModel[1] = shearModelOOP->getCopy("BeamFiber");
  if (theShearModel == 0) {
    opserr << "Macroelement3d::Macroelement3d - failed to allocate shear model pointer\n";
    exit(-1);
  }
 
  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = ndE;

  theNodes[0] = 0;
  theNodes[1] = 0;
  theNodes[2] = 0;

  for (int i=0; i<12; i++) {
	  q0[i] = 0.0;
      p0[i] = 0.0;
  }


  // define orientation
  Vector vAxis(3);
  for (int i=0; i<3; i++) {
      xAxis(i) = axis(i);
      vAxis(i) = oop(i);
  }

  // make xAxis unitary
  double xnorm = xAxis.Norm();
  if (xnorm == 0) {
        opserr << "Macroelement3d ("<< tag << "). Defined between two coincident nodes, check geometry." << endln;
        exit(-1);
  }
  
  // y = v cross x
  yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
  yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
  yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);

  double ynorm = yAxis.Norm();
  if (ynorm == 0) {
        opserr << "Macroelement3d ("<< tag << ").The vector v that defines plane xz is parallel to x axis." << endln;
        exit(-1);
  }
    
  yAxis /= ynorm;

  // z = x cross y. Already unitary
  zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
  zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
  zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);

  // Fill in transformation matrix
  R[0][0] = xAxis(0);
  R[0][1] = xAxis(1);
  R[0][2] = xAxis(2);
    
  R[1][0] = yAxis(0);
  R[1][1] = yAxis(1);
  R[1][2] = yAxis(2);
    
  R[2][0] = zAxis(0);
  R[2][1] = zAxis(1);
  R[2][2] = zAxis(2);

  // setup transformation matrix from global to local
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = Tgl(12,12) = Tgl(15,15)  = R[0][0];
    Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = Tgl(12,13) = Tgl(15,16)  = R[0][1];
    Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = Tgl(12,14) = Tgl(15,17)  = R[0][2];
    Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = Tgl(13,12) = Tgl(16,15)  = R[1][0];
    Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = Tgl(13,13) = Tgl(16,16)  = R[1][1];
    Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = Tgl(13,14) = Tgl(16,17)  = R[1][2];
    Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = Tgl(14,12) = Tgl(17,15)  = R[2][0];
    Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = Tgl(14,13) = Tgl(17,16)  = R[2][1];
    Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = Tgl(14,14) = Tgl(17,17)  = R[2][2];
    

    Tgl6(0,0) = Tgl6(3,3)  = R[0][0];
    Tgl6(0,1) = Tgl6(3,4)  = R[0][1];
    Tgl6(0,2) = Tgl6(3,5)  = R[0][2];
    Tgl6(1,0) = Tgl6(4,3)  = R[1][0];
    Tgl6(1,1) = Tgl6(4,4)  = R[1][1];
    Tgl6(1,2) = Tgl6(4,5)  = R[1][2];
    Tgl6(2,0) = Tgl6(5,3)  = R[2][0];
    Tgl6(2,1) = Tgl6(5,4)  = R[2][1];
    Tgl6(2,2) = Tgl6(5,5)  = R[2][2];

    
	this->intLength *= 2*L;
	this->intLengthMasses *= 2*L;
    
}

Macroelement3d::Macroelement3d()
:Element (0, ELE_TAG_DispBeamColumn3d),
          numSections(0), numShearModels(0), theSections(0), theShearModel(0),connectedExternalNodes(3), Q(18), q(12), uBasic(12), uBasicCommitted(12), Tgl(18,18), Tgl6(6,6), GammaC(12,18),
		  rho(0.0), cMass(0), parameterID(0), E(0.), intLength(3), intLengthMasses(3), driftF(0.), driftS(0.), limitDriftF(0), limitDriftS(0), limitDriftF_ALR(0), limitDriftS_ALR(0), 
		  alphaNC_SD(0.0), betaShearSpanF(0.0), betaShearSpanS(0.0), failureFactorF(0.0), failureFactorS(0.0), failedF(false), failedS(false), failedFcommitted(false), failedScommitted(false), 
		  collapsedF(false), collapsedS(false), collapsedFcommitted(false), collapsedScommitted(false),
		  nodeIInitialDisp(0), nodeJInitialDisp(0), nodeIOffset(0), nodeJOffset(0), xAxis(3), yAxis(3), zAxis(3), L(0.0), PDelta(0), deltaW1(0.), deltaV1(0.), deltaW3(0.), deltaV3(0.)
{
  for (int i=0; i<12; i++) {
	  q0[i] = 0.0;
      p0[i] = 0.0;
  }

  for (int i=0; i<3; i++) {
      R[i][0] = 0.0;
      R[i][1] = 0.0;
      R[i][2] = 0.0;
   }

  theNodes[0] = 0;
  theNodes[1] = 0;
  theNodes[2] = 0;

}

Macroelement3d::~Macroelement3d()
{    
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }

  for (int i = 0; i < 2; i++) {
    if (theShearModel[i])
      delete theShearModel[i];
  }
  
  // Delete the array of pointers to SectionForceDeformation pointer arrays and to the Shear model
  if (theSections)
    delete [] theSections;

  if (theShearModel)
    delete [] theShearModel;
}

int
Macroelement3d::getNumExternalNodes() const {
    return 3;
}

const ID&
Macroelement3d::getExternalNodes() {
    return connectedExternalNodes;
}

Node **
Macroelement3d::getNodePtrs() {
    return theNodes;
}

int
Macroelement3d::getNumDOF() {
    return 18;
}

void
Macroelement3d::setDomain(Domain *theDomain) {
	// Check Domain is not null - invoked when object removed from a domain

	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
        theNodes[2] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int NdE = connectedExternalNodes(2);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(NdE);

	

    if (theNodes[0] == 0 || theNodes[1] == 0) {
		return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNdE = theNodes[2]->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6 || dofNdE != 6) {
		return;
    }


    // element projection
    static Vector dx(3);
    
    const Vector &ndICoords = theNodes[0]->getCrds();
    const Vector &ndJCoords = theNodes[1]->getCrds();
    
    dx(0) = ndJCoords(0) - ndICoords(0);
    dx(1) = ndJCoords(1) - ndICoords(1);
    dx(2) = ndJCoords(2) - ndICoords(2);



    // see if there is some initial displacements at nodes
    const Vector &nodeIDisp = theNodes[0]->getDisp();
    const Vector &nodeJDisp = theNodes[0]->getDisp();
    bool addInitialDisplacement = false;

    for (int i=0; i<6; i++) {
            if (nodeIDisp(i) != 0.0 || nodeJDisp(i) != 0.0) {
			     addInitialDisplacement = true;
			}
	}

	if (addInitialDisplacement) {
          nodeIInitialDisp = new double [6];
          for (int j=0; j<6; j++) {
               nodeIInitialDisp[j] = nodeIDisp(j);
		  }
            
          nodeJInitialDisp = new double [6];
          for (int j=0; j<6; j++) {
				  nodeJInitialDisp[j] = nodeJDisp(j);
          }
	}

    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
        dx(2) -= nodeIInitialDisp[2];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
        dx(2) += nodeJInitialDisp[2];
    }



    // define position
    Vector inode(3);
    Vector jnode(3);
    const Vector &ndECoords = theNodes[2]->getCrds();

    inode = ndECoords - L*xAxis;
    jnode = ndECoords + L*xAxis;

     Vector offset(3);
     offset =  inode - ndICoords;
	 //opserr << "Macroelement " << this->getTag() << ": offset I = " << offset;
	 if (offset.Norm()>DBL_EPSILON) {
            nodeIOffset = new double[3];
            nodeIOffset[0] = offset(0);
            nodeIOffset[1] = offset(1);
            nodeIOffset[2] = offset(2);
      }


     offset = jnode - ndJCoords;
	 //opserr << "Macroelement " << this->getTag() << ": offset J = " << offset;

	 if (offset.Norm()>DBL_EPSILON) {
            nodeJOffset = new double[3];
            nodeJOffset[0] = offset(0);
            nodeJOffset[1] = offset(1);
            nodeJOffset[2] = offset(2);
      }

    
    if (nodeJOffset != 0) {
        dx(0) += nodeJOffset[0];
        dx(1) += nodeJOffset[1];
        dx(2) += nodeJOffset[2];
		
    }
    
    if (nodeIOffset != 0) {
        dx(0) -= nodeIOffset[0];
        dx(1) -= nodeIOffset[1];
        dx(2) -= nodeIOffset[2];
		
    }

	
    
    // calculate the element length
    L = dx.Norm();

	// we call L the height of each of the two halves
	L /= 2.;
   
    this->DomainComponent::setDomain(theDomain);
	this->getIncrementalCompatibilityMatrix(false);
	this->update();

	
}

int
Macroelement3d::commitState() {

	//if (this->getTag()==101) opserr << "El. " << this->getTag() << ": called commit\n";

    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
		opserr << "Macroelement3d::commitState () - failed in base class";
    }    

    // Loop over the interfaces and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    // commit shear model
	for (int i = 0; i < numShearModels; i++)
		retVal += theShearModel[i]->commitState();


	failedFcommitted = failedF;
	failedScommitted = failedS;

	collapsedFcommitted = collapsedF;
	collapsedScommitted = collapsedS;

	// commit basic displacements
	uBasicCommitted = uBasic;

    // commit transformation. No need to do anything for P-delta and linear transformations.

	// commit time for transient analysis (could not ask for dt directly)
	committedTime = OPS_GetDomain()->getCurrentTime();

	/*
	if (Kc!=0) {
		KsecantCommitted = this->getSecantStiffDamping();
		*Kc = KsecantCommitted ;
		//opserr << "Kc =" << KsecantCommitted;
	} 
	*/

	return retVal;
}

int
Macroelement3d::revertToLastCommit() {

	//if (this->getTag()==101) opserr << "El. " << this->getTag() << ": called revert to last commit\n";

    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    // revert shear model
	for (int i = 0; i < numShearModels; i++)
		retVal += theShearModel[i]->revertToLastCommit();

	failedF = failedFcommitted;
	failedS = failedScommitted;

	collapsedF = collapsedFcommitted;
	collapsedS = collapsedScommitted;

    // transformation. No need to do anything for P-delta and linear transformations.

    return retVal;
}

int
Macroelement3d::revertToStart()
{
	//if (this->getTag()==101) opserr << "El. " << this->getTag() << ": called revert to start\n";

    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();
    
	// revert shear model
	for (int i = 0; i < numShearModels; i++)
		retVal += theShearModel[i]->revertToStart();

	failedS = false;
	failedF = false;
	failedScommitted = false;
	failedFcommitted = false;

	collapsedS = false;
	collapsedF = false;
	collapsedScommitted = false;
	collapsedFcommitted = false;

    return retVal;
}

int
Macroelement3d::update(void)
{

   //if (this->getTag()==101) opserr << "El. " << this->getTag() << ": called update\n";

  int err = 0;

  const Vector &dispI = theNodes[0]->getTrialDisp();
  const Vector &dispJ = theNodes[1]->getTrialDisp();
  const Vector &dispE = theNodes[2]->getTrialDisp();

  uBasic = this->getBasicDisplacement(dispI, dispJ, dispE);

  double oneOverL = 1./L;

    // Set section deformations
	// divide by integration length
    Vector e(4);
    for (int i=0; i<3; i++)
		e(i) = uBasic(i)  / intLength(0);
	e(3) = uBasic(3) /(2.*L);  
	err += theSections[0]->setTrialSectionDeformation(e);

    e.Zero();
    for (int i=0; i<3; i++)
		e(i) = uBasic(i+4)   / intLength(1);
    err += theSections[1]->setTrialSectionDeformation(e);

	Vector N(4);
	N = theSections[1]->getStressResultant();

	Matrix EshearModel = theShearModel[0]->getInitialTangent();
	e(0) = N(0) / EshearModel(0,0);

	Vector s(2);
    s(0) = e(0);   
    s(1) = uBasic(10);
	err += theShearModel[0]->setTrialStrain(s);

	EshearModel = theShearModel[1]->getInitialTangent();
	e(0) = N(0) / EshearModel(0,0);

    s(0) = e(0);   
    s(1) = uBasic(11);
    err += theShearModel[1]->setTrialStrain(s);


    e.Zero();
    for (int i=0; i<3; i++)
		e(i) = uBasic(i+7)   / intLength(2);
    err += theSections[2]->setTrialSectionDeformation(e);
   

	
	// DRIFT CAPACITY MODEL
	// check for in-plane drift failure
	driftS = abs(-uBasic(10)) /2.0 *oneOverL;
	driftF = max(abs(-(uBasic(1)*(2*L) + uBasic(5)*L) /(2.0*L)), 
		         abs(-(uBasic(8)*(2*L) + uBasic(5)*L) /(2.0*L)));


	double axialLoadRatio;
	if (Ltfc>0.0)
		axialLoadRatio = -N(0)/Ltfc;
	else
	    axialLoadRatio = 0.0;

	this->driftModel(driftF, driftS, axialLoadRatio);
	
	/*
	if (limitDriftS>0.0) {
	    if (abs(driftS)>limitDriftS)	    failedS = true;
	} 
	*/
    
	/*
	if (limitDriftF>0.0) {
	    if (abs(driftF)>limitDriftF)	    failedF = true;
	} 
	*/

    if (err != 0) {
        opserr << "Macroelement3d::update() - failed setTrialSectionDeformations()\n";
    return err;
    }

    return 0;
}

void 
Macroelement3d::driftModel(double currentDriftF, double currentDriftS, double axialLoadRatio, double H0overL) {


	if (axialLoadRatio<DBL_EPSILON)  axialLoadRatio = 0.0;
	if (axialLoadRatio>1.0)  axialLoadRatio = 0.999;

	// flexure model
	if (limitDriftF.Size()>0) {

		int index = 0;
		bool ALR_found = false;

		while (index<limitDriftF_ALR.Size() && !ALR_found) {
			if (axialLoadRatio >= limitDriftF_ALR(index) ) {
				index++;
			} else { 
				ALR_found = true;
			}
		}

		double driftF_capacity = limitDriftF(index-1) + (limitDriftF(index)-limitDriftF(index-1)) / (limitDriftF_ALR(index)-limitDriftF_ALR(index-1)) 
		                     * (axialLoadRatio - limitDriftF_ALR(index-1));
	    
		//check for loss of lateral force capacity
		if (abs(driftF)>driftF_capacity)    failedF = true;

		if (alphaNC_SD>0) {  // check for axial collapse
			if (abs(driftF)>driftF_capacity*alphaNC_SD)    collapsedF = true;
		}


	}
	 
	//shear model
	if (limitDriftS.Size()>0) {

		int index = 0;
		bool ALR_found = false;

		while (index<limitDriftS_ALR.Size() && !ALR_found) {
			if (axialLoadRatio >= limitDriftS_ALR(index) ) 
				index++;
			else 
				ALR_found = true;
		}

		double driftS_capacity = limitDriftS(index-1) + (limitDriftS(index)-limitDriftS(index-1)) / (limitDriftS_ALR(index)-limitDriftS_ALR(index-1)) 
		                     * (axialLoadRatio - limitDriftS_ALR(index-1));
	    
		//check for loss of lateral force capacity
		if (abs(driftS)>driftS_capacity)    failedS = true;

		if (alphaNC_SD>0) {  // check for axial collapse
			if (abs(driftS)>driftS_capacity*alphaNC_SD)    collapsedS = true;
		}
	}

	return;
}

Vector 
Macroelement3d::getBasicDisplacement(Vector dispI, Vector dispJ, Vector dispE) {

  //opserr << "El. " << this->getTag() << ": called get basic displacement\n";

  // passes from global displacements/velocities to local
  Vector uGlobal(18);
  Vector uILocal(6);
  Vector uELocal(6);
  Vector uJLocal(6);

  for (int i=0; i<6; i++) {
       uGlobal(i)    = dispI(i);
       uGlobal(i+6)  = dispJ(i);
       uGlobal(i+12) = dispE(i);

       uILocal(i) = dispI(i);
       uJLocal(i) = dispJ(i);
       uELocal(i) = dispE(i);
  }

    
  if (nodeIInitialDisp != 0) {
        for (int j=0; j<6; j++) {
            uGlobal(j) -= nodeIInitialDisp[j];
            uILocal(j) -= nodeIInitialDisp[j];
		}
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<6; j++) {
            uGlobal(j+6) -= nodeJInitialDisp[j];
            uJLocal(j)   -= nodeIInitialDisp[j];
		}
    }

	// rotate to local
    Vector uLocal(18);
	uILocal = Tgl6*uILocal; 
    uJLocal = Tgl6*uJLocal; 
    uELocal = Tgl6*uELocal; 


	// correct if offset
	static double Wu[3];
    
    if (nodeIOffset) {
        Wu[0] =  nodeIOffset[2]*uGlobal(4) - nodeIOffset[1]*uGlobal(5);
        Wu[1] = -nodeIOffset[2]*uGlobal(3) + nodeIOffset[0]*uGlobal(5);
        Wu[2] =  nodeIOffset[1]*uGlobal(3) - nodeIOffset[0]*uGlobal(4);
        
		uILocal(0) += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
        uILocal(1) += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
        uILocal(2) += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
    }
    
    if (nodeJOffset) {
        Wu[0] =  nodeJOffset[2]*uGlobal(10) - nodeJOffset[1]*uGlobal(11);
        Wu[1] = -nodeJOffset[2]*uGlobal(9)  + nodeJOffset[0]*uGlobal(11);
        Wu[2] =  nodeJOffset[1]*uGlobal(9)  - nodeJOffset[0]*uGlobal(10);
        
		uJLocal(0) += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
        uJLocal(1) += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
        uJLocal(2) += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
    }

	
    for (int i=0; i<6; i++) {
       uLocal(i)    = uILocal(i);
       uLocal(i+6)  = uJLocal(i);
       uLocal(i+12) = uELocal(i);
    }

	double dV1, dV2, dV3;
	double dW1, dW2, dW3;

    // P-Delta transformation: update displacements
	if (PDelta) {   
		dV1 = uELocal(1) - uILocal(1);
		dW1 = uELocal(2) - uILocal(2);

		dV2 = uELocal(4) - uELocal(1);
		dW2 = uELocal(5) - uELocal(2);

		dV3 = uJLocal(1) - uELocal(4);
		dW3 = uJLocal(2) - uELocal(5);
	}

    // drop subscript Local for brevity
    Vector &uI = uILocal;
    Vector &uE = uELocal;
    Vector &uJ = uJLocal;

    double oneOverL = 1./L;

	Vector updatedBasic(12);

    // Get basic deformations. Total compatibility relations (linear) can be expressed also in matrix form as uBasic = Gamma_c * uLocal. 
	if (PDelta) {
		// non-linear total compatibility relations: OpenSees does not incluide this
		// first order approximationd of rotations, second order approximations of axial strain
		updatedBasic(0)  =  uE(0) - uI(0)  + 0.5*oneOverL*(pow(dV1,2) + pow(dW1,2));
		updatedBasic(1)  = -uI(5)  + oneOverL*dV1 ;
		updatedBasic(2)  = -uI(4)  - oneOverL*dW1;
		updatedBasic(3)  =  uJ(3) - uI(3);
		updatedBasic(4)  = -uE(0) + uE(3); // + 0.5*oneOverL*(-deltaV2*(uI(1)-uJ(1)) -deltaW2*(uI(2)-uJ(2)));  term to be included if exact formulation is applied
		updatedBasic(5)  =  oneOverL*(-dV1 + dV3);
		updatedBasic(6)  =  oneOverL*( dW1 - dW3);
		updatedBasic(7)  =  uJ(0) - uE(3) + 0.5*oneOverL*(pow(dV3,2) + pow(dW3,2));
		updatedBasic(8)  =  uJ(5) - oneOverL*dV3;
		updatedBasic(9)  =  uJ(4) + oneOverL*dW3;
		updatedBasic(10) =  dV2;  // the terms for the exact formulation are not included as they need an update of deltaU also
		updatedBasic(11) =  dW2;  // same here
	} else {
		// linear total compatibility relations
		updatedBasic(0)  =  uE(0) - uI(0);
		updatedBasic(1)  = -oneOverL*uI(1) - uI(5)  +oneOverL*uE(1);
		updatedBasic(2)  =  oneOverL*uI(2) - uI(4)  -oneOverL*uE(2);
		updatedBasic(3)  =  uJ(3) - uI(3);
		updatedBasic(4)  = -uE(0) + uE(3);
		updatedBasic(5)  =  oneOverL*(uI(1) + uJ(1) -uE(1) - uE(4));
		updatedBasic(6)  = oneOverL*(-uI(2) - uJ(2) +uE(2) + uE(5));
		updatedBasic(7)  = uJ(0) - uE(3);
		updatedBasic(8)  = -oneOverL*uJ(1) + uJ(5) + oneOverL*uE(4);
		updatedBasic(9)  =  oneOverL*uJ(2) + uJ(4) - oneOverL*uE(5);
		updatedBasic(10) = -uE(1) + uE(4);
		updatedBasic(11) = -uE(2) + uE(5);
	}


	if (PDelta) {
		deltaV1 = dV1;
		deltaW1 = dW1;

		//deltaV2 = dV2;
		//deltaW2 = dW2;

		deltaV3 = dV3;
		deltaW3 = dW3;

		getIncrementalCompatibilityMatrix(true);
	} 
	
	return updatedBasic;
}

const Matrix&
 Macroelement3d::getIncrementalCompatibilityMatrix(bool flagIncremental) {
	//if (this->getTag()==101)  opserr << "El. " << this->getTag() << ": called get incr comp matrix\n";

	GammaC.Zero();
    double oneOverL = 1./L;

	    GammaC(0,0) = -1;
		GammaC(0,12) = 1;
		GammaC(1,1) = -oneOverL;
		GammaC(1,5) = -1;
		GammaC(1,13) = oneOverL;
		GammaC(2,2) = oneOverL;
		GammaC(2,4) = -1;
		GammaC(2,14) = -oneOverL;
		GammaC(3,3) = -1;
		GammaC(3,9) = 1;
		GammaC(4,12) = -1;
		GammaC(4,15) = 1;
		GammaC(5,1) = oneOverL;
		GammaC(5,7) = oneOverL;
		GammaC(5,13) = -oneOverL;
		GammaC(5,16) = -oneOverL;
		GammaC(6,2) = -oneOverL;
		GammaC(6,8) = -oneOverL;
		GammaC(6,14) = oneOverL;
		GammaC(6,17) = oneOverL;
		GammaC(7,6) = 1;
		GammaC(7,15) = -1;
		GammaC(8,7) = -oneOverL;
		GammaC(8,11) = 1;
		GammaC(8,16) = oneOverL;
		GammaC(9,8) = oneOverL;
		GammaC(9,10) = 1;
		GammaC(9,17) = -oneOverL;
		GammaC(10,13) = -1;
		GammaC(10,16) = 1;
		GammaC(11,14) = -1;
		GammaC(11,17) = 1;


    if ((PDelta) && flagIncremental) {
		/*
		GammaC(0,1) = -(deltaV1*oneOverL);
		GammaC(0,2) = -(deltaW1*oneOverL);
		GammaC(0,13) = deltaV1*oneOverL;
		GammaC(0,14) = deltaW1*oneOverL;
		GammaC(4,1) = -(deltaV2*oneOverL)/2.;
		GammaC(4,2) = -(deltaW2*oneOverL)/2.;
		GammaC(4,7) = (deltaV2*oneOverL)/2.;
		GammaC(4,8) = (deltaW2*oneOverL)/2.;
		GammaC(4,13) = ((-deltaV1 - deltaV2 - deltaV3)*oneOverL)/2.;
		GammaC(4,14) = ((-deltaW1 - deltaW2 - deltaW3)*oneOverL)/2.;
		GammaC(4,16) = ((deltaV1 + deltaV2 + deltaV3)*oneOverL)/2.;
		GammaC(4,17) = ((deltaW1 + deltaW2 + deltaW3)*oneOverL)/2.;
		GammaC(7,7) = deltaV3*oneOverL;
		GammaC(7,8) = deltaW3*oneOverL;
		GammaC(7,16) = -(deltaV3*oneOverL);
		GammaC(7,17) = -(deltaW3*oneOverL);
		*/
		GammaC(0,1) = -(deltaV1*oneOverL);
		GammaC(0,2) = -(deltaW1*oneOverL);
		GammaC(0,13) = deltaV1*oneOverL;
		GammaC(0,14) = deltaW1*oneOverL;
		GammaC(7,7) = deltaV3*oneOverL;
		GammaC(7,8) = deltaW3*oneOverL;
		GammaC(7,16) = -(deltaV3*oneOverL);
		GammaC(7,17) = -(deltaW3*oneOverL);
    }

    return GammaC;
  };

const Matrix&
Macroelement3d::getTangentStiff()
{
  //if (this->getTag()==101) opserr << "El. " << this->getTag() << ": called get tangent\n";

  static Matrix kb(12,12);
  
  kb.Zero();
  q.Zero();
  
  double oneOverL = 1.0/L;

  // first interface
  int order = 4;
  const Vector &s  = theSections[0]->getStressResultant();
  const Matrix &ks = theSections[0]->getSectionTangent();
 
  for (int i=0; i<order; i++) {
	  q(i) = s(i);
	  for (int j=0; j<4; j++) 
		  kb(i,j) = ks(i,j)    / intLength(0);
  }
  kb(3,3) = ks(3,3) /(2.*L);


  /*
  if (this->getTag()==3 || this->getTag()==2) {
	  opserr << "El. " << this->getTag() << ", Forces 1:  " << s;
	  opserr << "El. " << this->getTag() << ", Section 1: " << ks;
  }
  */
  
  // second and third interface (theSections[1] and theSections[2])
  order = 3;
  for (int sect=1; sect<3; sect++) {
	  const Vector &s  = theSections[sect]->getStressResultant();
	  const Matrix &ks = theSections[sect]->getSectionTangent();
	  for (int i=0; i<order; i++) {
		  q(1+sect*(order)+i) = s(i);
		  for (int j=0; j<order; j++) 
			  kb(1+sect*(order)+i,1+sect*(order)+j) = ks(i,j)      / intLength(sect);
	  }
	/* 
	  if (this->getTag()==3 || this->getTag()==2) {
		opserr << "El. " << this->getTag() << ", Forces " << sect+1 << ":  " << s;
	    opserr << "El. " << this->getTag() << ", Section " << sect+1 << ": " << ks;
     }
	 */
  }

  // shear interface 1 and 2
  for (int sect=0; sect<2; sect++) {

	Matrix EshearModel = theShearModel[sect]->getInitialTangent();
	double N = q(4);
	double Esm = EshearModel(0,0);

    const Matrix &ks = theShearModel[sect]->getTangent();
    const Vector &s  = theShearModel[sect]->getStress();

    kb(10+sect,10+sect) = ks(1,1);
	q(10+sect) = s(1);

	// update non-diagonal terms of the stiffness matrix (dependency on N2)
	for (int j=0; j<3; j++)
		kb(10+sect, 4+j)  += ks(1,0) /Esm * kb(4, 4+j);
   
  }
 
  // kill lateral capacity if failed
  if (collapsedS || collapsedF) {
	  double failureFactor;
	  if (collapsedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: all
	  int dofsToKill[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 
	  int dof;

	  for (int kDof=0; kDof<12; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
		  for (int j=0; j<kb.noCols(); j++)
			  kb(dof,j) *= failureFactor;
	  }
  } else if (failedS || failedF) {

	  double failureFactor;
	  if (failedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
	  int dofsToKill[8] = {1,2, 5,6, 8,9,  10,11}; 
	  int dof;

	  for (int kDof=0; kDof<8; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
		  for (int j=0; j<kb.noCols(); j++)
			  kb(dof,j) *= failureFactor;
	  }
  }
  
   
  // Transform to local stiffness matrix
  
  Matrix Klocal(18,18); 
  //if (PDelta==1)   getIncrementalCompatibilityMatrix(true);
  Klocal.addMatrixTripleProduct(0.0, GammaC, kb, 1.0); 

  //opserr << "local stiffness = " << Klocal;
  
  // add geometric local stiffness matrix  
  if (PDelta==1) {
        Matrix kgeom(18,18);

	    // add fixed end forces. Update
	    for (int i=0; i<12; i++) {
		   q(i) += q0[i];
	    }

	    double qq0 = q(0);
	    double qq7 = q(7);
	    double qq4 = q(4);

	  	kgeom(1,1) = qq0*oneOverL;
		kgeom(1,13) = -(qq0*oneOverL);
		kgeom(2,2) = qq0*oneOverL;
		kgeom(2,14) = -(qq0*oneOverL);
		kgeom(7,7) = qq7*oneOverL;
		kgeom(7,16) = -(qq7*oneOverL);
		kgeom(8,8) = qq7*oneOverL;
		kgeom(8,17) = -(qq7*oneOverL);
		kgeom(13,1) = -(qq0*oneOverL);
		kgeom(13,13) = qq0*oneOverL;
		kgeom(14,2) = -(qq0*oneOverL);
		kgeom(14,14) = qq0*oneOverL;
		kgeom(16,7) = -(qq7*oneOverL);
		kgeom(16,16) = qq7*oneOverL;
		kgeom(17,8) = -(qq7*oneOverL);
		kgeom(17,17) = qq7*oneOverL;

        Klocal += kgeom;
  }
   

   // Transform to global stiffness matrix   Kglobal = Tgl' * Klocal * Tgl
   // done now otherwise the secant stiffness matrix would be rrturned
	trasformMatrixToGlobal(Klocal);

  return K;
}

const Matrix&
Macroelement3d::getInitialBasicStiff()
{

	//opserr << "El. " << this->getTag() << ": called get K0 basic\n";
  //  opserr << "called initial basic stiffness\n";
  static Matrix kb(12,12);  
  kb.Zero();
  
  double oneOverL = 1.0/L;

  // first interface
  int order = 3;
  const Matrix &ks = theSections[0]->getInitialTangent();
 
  for (int i=0; i<order; i++) {
	  for (int j=0; j<3; j++) 
		  kb(i,j) = ks(i,j)    / intLength(0);
  }
  kb(3,3) = ks(3,3) /(2.*L);
  
  // second and third interface
  order = 3;
  for (int sect=1; sect<3; sect++) {
	  const Matrix &ks = theSections[sect]->getInitialTangent();
	  for (int i=0; i<order; i++) {
		  for (int j=0; j<3; j++) 
			  kb(1+sect*(order)+i,1+sect*(order)+j) = ks(i,j)    / intLength(sect);
	  }
  }

  // shear interface 1 and 2
  for (int sect=0; sect<2; sect++) {
    const Matrix &ks = theShearModel[sect]->getInitialTangent();
    kb(10+sect,10+sect) = ks(1,1);
	// no need to update diagonal terms->they are already uncoupled
  }

  //opserr << "Basic stiffness: \n" << kb;

  if (collapsedS || collapsedF) {

	  double failureFactor;
	  if (collapsedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: all
	  int dofsToKill[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 
	  int dof;

	  for (int kDof=0; kDof<12; kDof++) {
		  dof = dofsToKill[kDof];
		  for (int j=0; j<kb.noCols(); j++)
			  kb(dof,j) *= failureFactor;
	  }
  } else if (failedS || failedF) {

	  double failureFactor;
	  if (failedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
	  int dofsToKill[8] = {1,2, 5,6, 8,9,  10,11}; 
	  int dof;

	  for (int kDof=0; kDof<8; kDof++) {
		  dof = dofsToKill[kDof];
		  for (int j=0; j<kb.noCols(); j++)
			  kb(dof,j) *= failureFactor;
	  }
  }

  
  return kb;
}

const Matrix&
Macroelement3d::getInitialStiff()
{

  //if (this->getTag()==101)  opserr << "El. " << this->getTag() << ": called get K0\n";

  const Matrix &kb = this->getInitialBasicStiff();
 
  // Transform to local stiffness matrix
  //if (PDelta==1)  getIncrementalCompatibilityMatrix(false);
  //getIncrementalCompatibilityMatrix();   //  check this 

  Matrix Klocal(18,18);

  Klocal.addMatrixTripleProduct(0.0, GammaC, kb, 1.0); 
  trasformMatrixToGlobal(Klocal);

  // opserr << "Stiffness matrix: \n" << K;
  
  return K;
}

const Matrix&
Macroelement3d::getMass()
{

 // if (this->getTag()==101)  opserr << "El. " << this->getTag() << ": called get mass\n";
  K.Zero();
  if (rho == 0.0)
    return K;
  
  if (cMass == 0)  {
    // lumped mass matrix. 1/4 of the mass to the end nodes, 1/2 to the center node
	  if (isGable) {
		double m = rho*2.*L;
		K(0,0) = 5./12.*m *massglobalDir(0);
		K(1,1) = 5./12.*m *massglobalDir(1);
		K(2,2) = 5./12.*m *massglobalDir(2);
		
		K(6,6) = 1./12.*m *massglobalDir(0);
		K(7,7) = 1./12.*m *massglobalDir(1);
		K(8,8) = 1./12.*m *massglobalDir(2);

		K(12,12) = 1./3.*m *massglobalDir(0);
		K(13,13) = 1./3.*m *massglobalDir(1);
		K(14,14) = 1./3.*m *massglobalDir(2);

		K(15,15) = 1./6.*m *massglobalDir(0);
		K(16,16) = 1./6.*m *massglobalDir(1);
		K(17,17) = 1./6.*m *massglobalDir(2);
		
	  } else {
		//double m = rho*2.*L /4.;
		//K(0,0) = m *massglobalDir(0);
		//K(1,1) = m *massglobalDir(1);
		//K(2,2) = m *massglobalDir(2);
		
		//K(6,6) = m *massglobalDir(0);
		//K(7,7) = m *massglobalDir(1);
		//K(8,8) = m *massglobalDir(2);

		//K(12,12) = m *massglobalDir(0);
		//K(13,13) = m *massglobalDir(1);
		//K(14,14) = m *massglobalDir(2);
		
		//K(15,15) = m *massglobalDir(0);
		//K(16,16) = m *massglobalDir(1);
		//K(17,17) = m *massglobalDir(2);

		double m = rho;
		K(0,0) = m *intLengthMasses(0) *massglobalDir(0);
		K(1,1) = m *intLengthMasses(0) *massglobalDir(1);
		K(2,2) = m *intLengthMasses(0) *massglobalDir(2);
		
		K(6,6) = m *intLengthMasses(2) *massglobalDir(0);
		K(7,7) = m *intLengthMasses(2) *massglobalDir(1);
		K(8,8) = m *intLengthMasses(2) *massglobalDir(2);

		K(12,12) = m *intLengthMasses(1)/2. *massglobalDir(0);
		K(13,13) = m *intLengthMasses(1)/2. *massglobalDir(1);
		K(14,14) = m *intLengthMasses(1)/2. *massglobalDir(2);
		
		K(15,15) = m *intLengthMasses(1)/2. *massglobalDir(0);
		K(16,16) = m *intLengthMasses(1)/2. *massglobalDir(1);
		K(17,17) = m *intLengthMasses(1)/2. *massglobalDir(2);
	  }

  } else  {
    // consistent mass matrix	  
    static Matrix massLocal(18,18);
    double m = rho*L;
    
	  if (isGable) {
		massLocal(0,0) = 5./12.*2.; 
		massLocal(1,1) = 5./12.*2.; 
		massLocal(2,2) = 5./12.*2.;

		massLocal(6,6) = 1./12.*2.; 
		massLocal(7,7) = 1./12.*2.; 
		massLocal(8,8) = 1./12.*2.;

		massLocal(12,12) = 1./3.*2; 
		massLocal(13,13) = 1./3.*2; 
		massLocal(14,14) = 1./3.*2;

		massLocal(15,15) = 1./6.*2.; 
		massLocal(16,16) = 1./6.*2.; 
		massLocal(17,17) = 1./6.*2.;

	  } else {
		  /*
		massLocal(0,0) = 0.3333333333333333;
		massLocal(0,12) = 0.16666666666666666;
		massLocal(1,1) = 0.3333333333333333;
		massLocal(1,13) = 0.16666666666666666;
		massLocal(2,2) = 0.3333333333333333;
		massLocal(2,14) = 0.16666666666666666;
		massLocal(6,6) = 0.3333333333333333;
		massLocal(6,15) = 0.16666666666666666;
		massLocal(7,7) = 0.3333333333333333;
		massLocal(7,16) = 0.16666666666666666;
		massLocal(8,8) = 0.3333333333333333;
		massLocal(8,17) = 0.16666666666666666;
		massLocal(12,0) = 0.16666666666666666;
		massLocal(12,12) = 0.3333333333333333;
		massLocal(13,1) = 0.16666666666666666;
		massLocal(13,13) = 0.3333333333333333;
		massLocal(14,2) = 0.16666666666666666;
		massLocal(14,14) = 0.3333333333333333;
		massLocal(15,6) = 0.16666666666666666;
		massLocal(15,15) = 0.3333333333333333;
		massLocal(16,7) = 0.16666666666666666;
		massLocal(16,16) = 0.3333333333333333;
		massLocal(17,8) = 0.16666666666666666;
		massLocal(17,17) = 0.3333333333333333;
		*/

		massLocal(1,1) = 0.3333333333333333;
		massLocal(1,13) = 0.16666666666666666;
		massLocal(2,2) = 0.3333333333333333;
		massLocal(2,14) = 0.16666666666666666;
		massLocal(7,7) = 0.3333333333333333;
		massLocal(7,16) = 0.16666666666666666;
		massLocal(8,8) = 0.3333333333333333;
		massLocal(8,17) = 0.16666666666666666;
		massLocal(12,12) = 1;
		massLocal(13,1) = 0.16666666666666666;
		massLocal(13,13) = 0.3333333333333333;
		massLocal(14,2) = 0.16666666666666666;
		massLocal(14,14) = 0.3333333333333333;
		massLocal(15,15) = 1;
		massLocal(16,7) = 0.16666666666666666;
		massLocal(16,16) = 0.3333333333333333;
		massLocal(17,8) = 0.16666666666666666;
		massLocal(17,17) = 0.3333333333333333;

	  }
	massLocal *= m;
	trasformMatrixToGlobal(massLocal);


	  // apply mass directions
  
  for (int i=0; i<3; i++) {
	 for (int j=0; j<18; j++) {
		  K(i,   j)  *= massglobalDir(i);
		  K(i+6, j)  *= massglobalDir(i);
		  K(i+12,j)  *= massglobalDir(i);
		  K(i+15,j)  *= massglobalDir(i);

		  K(j, i)  *= massglobalDir(i);
		  K(j, i+6)  *= massglobalDir(i);
		  K(j, i+12)  *= massglobalDir(i);
		  K(j, i+15)  *= massglobalDir(i);
	  }
   }
  

  }

  return K;
}

void
Macroelement3d::zeroLoad(void)
{

	//opserr << "El. " << this->getTag() << ": called zero load\n";

 Q.Zero();

  for (int i=0; i<12; i++) {
	 q0[i] = 0.0;
	 p0[i] = 0.0;
  }
  return;
}

int 
Macroelement3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{

	//LOAD_TAG_Beam3dUniformLoad
	//LOAD_TAG_SelfWeight

	//if (this->getTag()==101)  opserr << "El. " << this->getTag() << ": called add load\n";


   int type;
   const Vector &data = theLoad->getData(type, loadFactor);

   if (type == LOAD_TAG_SelfWeight) {

	   Vector gravity(6);
	   gravity(0) = data(0) * loadFactor;   // gravity direction in global coordinates
	   gravity(1) = data(1) * loadFactor;
	   gravity(2) = data(2) * loadFactor;

	   // rotate to local
	   gravity = Tgl6*gravity; 

	    double wx = gravity(0)*rho;  // Axial (+direction from node I to J)
		double wy = gravity(1)*rho;  // Transverse y
		double wz = gravity(2)*rho;  // Transverse z

		if (isGable) {
			double Vy = wy*L;
			double Vz = wz*L;
			double Fx = wx*L;   

			// Reactions in basic system
			p0[0] -= 5./12.*Fx;
			p0[1] -= 1./12.*Fx;
			p0[2] -= 1./3. *Fx;
			p0[3] -= 1./6. *Fx;

		    p0[4] -= 5./12.*Vy;
			p0[5] -= 1./12.*Vy;
			p0[6] -= 1./3. *Vy;
			p0[7] -= 1./6. *Vy;

			p0[8]  -= 5./12.*Vz;
			p0[9]  -= 1./12.*Vz;
			p0[10] -= 1./3. *Vz;
			p0[11] -= 1./6. *Vz;
		
			// Fixed end forces in basic system
			q0[0] -= 5./12.*Fx;
			q0[7] += 1./12.*Fx;

		} else {

			double Vy = 0.5*wy*L;
			double Vz = 0.5*wz*L;
			double Fx = wx*L*2.;   

			// Reactions in basic system
			p0[0] -= 0.25*Fx;
			p0[1] -= 0.25*Fx;
			p0[2] -= 0.25*Fx;
			p0[3] -= 0.25*Fx;

			p0[4] -= Vy;
			p0[5] -= Vy;
			p0[6] -= Vy;
			p0[7] -= Vy;

			p0[8] -= Vz;
			p0[9] -= Vz;
			p0[10] -= Vz;
			p0[11] -= Vz;

			// Fixed end forces in basic system
			q0[0] -= 0.25*Fx;
			q0[7] += 0.25*Fx;
		}
  }
   else if (type == LOAD_TAG_Beam3dUniformLoad) {
	    double wy = data(0)*loadFactor;  // Transverse y
		double wz = data(1)*loadFactor;  // Transverse z
		double wx = data(2)*loadFactor;  // Axial 

			double Vy = 0.5*wy*L;
			double Vz = 0.5*wz*L;
			double Fx = wx*L*2.;   

			// Reactions in basic system
			p0[0] -= 0.25*Fx;
			p0[1] -= 0.25*Fx;
			p0[2] -= 0.25*Fx;
			p0[3] -= 0.25*Fx;
		
			p0[4] -= Vy;
			p0[5] -= Vy;
			p0[6] -= Vy;
			p0[7] -= Vy;
		
			p0[8] -= Vz;
			p0[9] -= Vz;
			p0[10] -= Vz;
			p0[11] -= Vz;

			// Fixed end forces in basic system
			q0[0] -= 0.25*Fx;
			q0[7] += 0.25*Fx;
   }
   else {
		opserr << "Macroelement3d::addLoad() -- load type unknown for element with tag: " << this->getTag() << endln;
		opserr << "Available elemental load types: -selfWeight $gX $gY $gZ (global system) and -beamUniform $Wy $Wz $Wx (local system)\n";
		return -1;
  }

	/*
  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1, V2;
    V1 = Py*(1.0-aOverL);
    V2 = Py*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz*(1.0-aOverL);
    V2 = Pz*aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1, M2;
    M1 = -a * b2 * Py * L2;
    M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  }
  */

   return 0;
 }

int 
Macroelement3d::addInertiaLoadToUnbalance(const Vector &accel)
{
//	opserr << "El. " << this->getTag() << ": called add inertia\n";
  // Check for a quick return
  if (rho == 0.0) {

    return 0;
  }
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &RaccelE = theNodes[2]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "Macroelement3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
	if (isGable) {
		double m = rho*2.*L;
		Q(0)  -= 5./12.*m *Raccel1(0) *massglobalDir(0);
		Q(1)  -= 5./12.*m *Raccel1(1) *massglobalDir(1);
		Q(2)  -= 5./12.*m *Raccel1(2) *massglobalDir(2);
		
		Q(6)  -= 1./12.*m *Raccel2(0) *massglobalDir(0);
		Q(7)  -= 1./12.*m *Raccel2(1) *massglobalDir(1);
		Q(8)  -= 1./12.*m *Raccel2(2) *massglobalDir(2);

		Q(12) -= 1./3.*m *RaccelE(0) *massglobalDir(0);
		Q(13) -= 1./3.*m *RaccelE(1) *massglobalDir(1);
		Q(14) -= 1./3.*m *RaccelE(2) *massglobalDir(2);

		Q(15) -= 1./6.*m *RaccelE(3) *massglobalDir(0);
		Q(16) -= 1./6.*m *RaccelE(4) *massglobalDir(1);
		Q(17) -= 1./6.*m *RaccelE(5) *massglobalDir(2);
		
	  } else {
		double m = rho;
		Q(0)  -= m *intLengthMasses(0) *Raccel1(0) *massglobalDir(0);
		Q(1)  -= m *intLengthMasses(0) *Raccel1(1) *massglobalDir(1);
		Q(2)  -= m *intLengthMasses(0) *Raccel1(2) *massglobalDir(2);
		
		Q(6)  -= m *intLengthMasses(2) *Raccel2(0) *massglobalDir(0);
		Q(7)  -= m *intLengthMasses(2) *Raccel2(1) *massglobalDir(1);
		Q(8)  -= m *intLengthMasses(2) *Raccel2(2) *massglobalDir(2);

		Q(12) -= m *intLengthMasses(1)/2. *RaccelE(0) *massglobalDir(0);
		Q(13) -= m *intLengthMasses(1)/2. *RaccelE(1) *massglobalDir(1);
		Q(14) -= m *intLengthMasses(1)/2. *RaccelE(2) *massglobalDir(2);

		Q(15) -= m *intLengthMasses(1)/2. *RaccelE(3) *massglobalDir(0);
		Q(16) -= m *intLengthMasses(1)/2. *RaccelE(4) *massglobalDir(1);
		Q(17) -= m *intLengthMasses(1)/2. *RaccelE(5) *massglobalDir(2);


	  }


  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(18);
    for (int i=0; i<6; i++)  {
      Raccel(i)    = Raccel1(i);
      Raccel(i+6)  = Raccel2(i);
	  Raccel(i+12) = RaccelE(i);
    }

    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector&
Macroelement3d::getResistingForce()
{

 // if (this->getTag()==101) 	opserr << "El. " << this->getTag() << ": called get force\n";

  // Get transpose of the equilibrium matrix
  //if (PDelta==1)  getIncrementalCompatibilityMatrix(false);
  
  P.Zero();
  q.Zero();
  
  double oneOverL = 1.0/L;

  // first interface
  int order = 4;
  const Vector &s  = theSections[0]->getStressResultant();

 
  for (int i=0; i<order; i++) {
	  q(i) = s(i);
  }

  
  if (s.Norm()<=DBL_EPSILON) {
		if ((theSections[0]->getSectionDeformation()).Norm() >= DBL_EPSILON) {
		 //opserr << "Element " << this->getTag() << ": failed in tension at section 1 " << endln;	  
		}
  }
 
  // second and third interface
  order = 3;
  for (int sect=1; sect<3; sect++) {
	  const Vector &s  = theSections[sect]->getStressResultant();
	  for (int i=0; i<order; i++) {
		  q(1+sect*(order)+i) = s(i);
	  }
	  
  	  if (s.Norm()<=DBL_EPSILON) {
		  if ((theSections[sect]->getSectionDeformation()).Norm() >= DBL_EPSILON) {
			 // opserr << "Element " << this->getTag() << ": failed in tension at section " << sect+1 << endln;
		  }
      }
  }

  // shear interface 1 and 2
  for (int sect=0; sect<2; sect++) {
    const Vector &s  = theShearModel[sect]->getStress();
    q(10+sect) = s(1);
  //  opserr<< "shear 12 = " << s;
  }

  /*
  	if (this->getTag()==3 || this->getTag()==2) {
		opserr << "element " << this->getTag() << "  FORCES ------------------------------" << endln;
		opserr << "1st section: " << q(0) << ", " << q(1) << ", " << q(2) << endln;
		opserr << "2nd section: " << q(4) << ", " << q(5) << ", " << q(6) << endln;
		opserr << "3rd section: " << q(7) << ", " << q(8) << ", " << q(9) << endln;
		opserr << "shear spring: " << q(10) << ", " << q(11) << endln;

	}
	*/

    if (collapsedS || collapsedF) {

		double failureFactor;
	  if (collapsedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: all
	  int dofsToKill[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 
	  int dof;

	  for (int kDof=0; kDof<12; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
	  }
   } else if (failedS || failedF) {

	   double failureFactor;
	  if (failedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
	  int dofsToKill[8] = {1,2, 5,6, 8,9,  10,11}; 
	  int dof;

	  for (int kDof=0; kDof<8; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
	  }
   }


  // add fixed end forces . Update
  for (int i=0; i<12; i++) {
	  q(i) += q0[i];
  }
     
  // Transform forces to local reference system
  Vector pLocal(18);

  pLocal = GammaC ^ q;  // pLocal = GammaE * q = GammaC'*q;

  // add fixed reactions
  pLocal(0)  += p0[0];
  pLocal(6)  += p0[1];
  pLocal(12) += p0[2];
  pLocal(15) += p0[3];

  pLocal(1)  += p0[4];
  pLocal(7)  += p0[5];
  pLocal(13) += p0[6];
  pLocal(16) += p0[7];

  pLocal(2)  += p0[8];
  pLocal(8)  += p0[9];
  pLocal(14) += p0[10];
  pLocal(17) += p0[11];
  

  // Transform forces to global reference system
  P = Tgl ^ pLocal;

  // add offset effect
    if (nodeIOffset) {
        P(3) += -nodeIOffset[2]*P(1) + nodeIOffset[1]*P(2);
        P(4) +=  nodeIOffset[2]*P(0) - nodeIOffset[0]*P(2);
        P(5) += -nodeIOffset[1]*P(0) + nodeIOffset[0]*P(1);
    }
    
    if (nodeJOffset) {
        P(9)  += -nodeJOffset[2]*P(7) + nodeJOffset[1]*P(8);
        P(10) +=  nodeJOffset[2]*P(6) - nodeJOffset[0]*P(8);
        P(11) += -nodeJOffset[1]*P(6) + nodeJOffset[0]*P(7);
    }

  return P;
}

const Vector&
Macroelement3d::getResistingForceIncInertia()
{

 // if (this->getTag()==101)  opserr << "El. " << this->getTag() << ": called get force incl inertia\n";

  P = this->getResistingForce();

  //opserr << "Resisting force: " << P << endln;
    
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  P.addVector(1.0, Q, -1.0);


  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accelE = theNodes[2]->getTrialAccel();
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
	if (isGable) {
		double m = rho*2.*L;
		P(0)  += 5./12.*m* accel1(0) *massglobalDir(0);
		P(1)  += 5./12.*m *accel1(1) *massglobalDir(1);
		P(2)  += 5./12.*m *accel1(2) *massglobalDir(2);
		
		P(6)  += 1./12.*m *accel2(0) *massglobalDir(0);
		P(7)  += 1./12.*m *accel2(1) *massglobalDir(1);
		P(8)  += 1./12.*m *accel2(2) *massglobalDir(2);

		P(12) += 1./3.*m *accelE(0) *massglobalDir(0);
		P(13) += 1./3.*m *accelE(1) *massglobalDir(1);
		P(14) += 1./3.*m *accelE(2) *massglobalDir(2);

		P(15) += 1./6.*m *accelE(3) *massglobalDir(0);
		P(16) += 1./6.*m *accelE(4) *massglobalDir(1);
		P(17) += 1./6.*m *accelE(5) *massglobalDir(2);
		
	} else {
		double m = rho;
		P(0)  += m *intLengthMasses(0) *accel1(0) *massglobalDir(0);
		P(1)  += m *intLengthMasses(0) *accel1(1) *massglobalDir(1);
		P(2)  += m *intLengthMasses(0) *accel1(2) *massglobalDir(2);
		
		P(6)  += m *intLengthMasses(2) *accel2(0) *massglobalDir(0);
		P(7)  += m *intLengthMasses(2) *accel2(1) *massglobalDir(1);
		P(8)  += m *intLengthMasses(2) *accel2(2) *massglobalDir(2);

		P(12) += m *intLengthMasses(1)/2. *accelE(0) *massglobalDir(0);
		P(13) += m *intLengthMasses(1)/2. *accelE(1) *massglobalDir(1);
		P(14) += m *intLengthMasses(1)/2. *accelE(2) *massglobalDir(2);

		P(15) += m *intLengthMasses(1)/2. *accelE(3) *massglobalDir(0);
		P(16) += m *intLengthMasses(1)/2. *accelE(4) *massglobalDir(1);
		P(17) += m *intLengthMasses(1)/2. *accelE(5) *massglobalDir(2);
	  }


  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(18);
    for (int i=0; i<6; i++)  {
      accel(i)   = accel1(i);
      accel(i+6) = accel2(i);
	  accel(i+12) = accelE(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }

   
   // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
	     Vector dampingForces = this->getRayleighDampingForces();
         P.addVector(1.0, dampingForces, 1.0);

		 if (this->getTag()==3) {
	       //opserr << dampingForces;
           //opserr << "uBasic=" << uBasic;
	       //opserr << P << endln;
  
		 }
  }

  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0 ) {
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
	  
	}



  }

  //opserr << "Resisting force including inertia: " << P << endln;

  return P;
}

int
Macroelement3d::sendSelf(int commitTag, Channel &theChannel)
{

  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static Vector data(14);

  data(0) = this->getTag();
  data(1) = connectedExternalNodes(0);
  data(2) = connectedExternalNodes(1);
  data(3) = numSections;
  data(4) = 0; // crdTransf->getClassTag();
  int crdTransfDbTag  = 1; //crdTransf->getDbTag();
  /*
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  */
  data(5) = crdTransfDbTag;
  data(6) = 0; //beamInt->getClassTag();
  int beamIntDbTag  = 1; //beamInt->getDbTag();
  /*
  if (beamIntDbTag  == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag  != 0) 
      beamInt->setDbTag(beamIntDbTag);
  }
  */
  data(7) = beamIntDbTag;
  data(8) = rho;
  data(9) = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;
  
  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "DispBeamColumn3d::sendSelf() - failed to send data Vector\n";
     return -1;
  } 
  // check here
  /*
  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumn3d::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumn3d::sendSelf() - failed to send beamInt\n";
    return -1;
  }      
  
  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2*numSections);
  loc = 0;
  for (i = 0; i<numSections; i++) {
    int sectClassTag = theSections[i]->getClassTag();
    int sectDbTag = theSections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      theSections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    
  

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumn3d::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }
  */

  return 0;
}

int
Macroelement3d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
    // still to implement
    /*
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static Vector data(14);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0)  {
    opserr << "DispBeamColumn3d::recvSelf() - failed to recv data Vector\n";
    return -1;
  }
  
  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect = (int)data(3);
  int crdTransfClassTag = (int)data(4);
  int crdTransfDbTag = (int)data(5);
  
  int beamIntClassTag = (int)data(6);
  int beamIntDbTag = (int)data(7);
  
  rho = data(8);
  cMass = (int)data(9);
  
  alphaM = data(10);
  betaK = data(11);
  betaK0 = data(12);
  betaKc = data(13);
  
  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "DispBeamColumn3d::recvSelf() - " <<
	  "failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "DispBeamColumn3d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "DispBeamColumn3d::sendSelf() - failed to recv beam integration\n";
     return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  //
  // now receive the sections
  //
  
  if (numSections != nSect) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (int i=0; i<numSections; i++)
	delete theSections[i];
      delete [] theSections;
    }

    // create a new array to hold pointers
    theSections = new SectionForceDeformation *[nSect];
    if (theSections == 0) {
      opserr << "DispBeamColumn3d::recvSelf() - out of memory creating sections array of size" <<
	nSect << endln;
      exit(-1);
    }    

    // create a section and recvSelf on it
    numSections = nSect;
    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
	opserr << "DispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn3d::recvSelf() - section " <<
	  i << "failed to recv itself\n";
	return -1;
      }     
    }

  } else {

    // 
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //
    
    loc = 0;
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (theSections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete theSections[i];
	theSections[i] = theBroker.getNewSection(sectClassTag);
	if (theSections[i] == 0) {
	  opserr << "DispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn3d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
	return -1;
      }     
    }
  }
  */

  return 0;
}

void
Macroelement3d::Print(OPS_Stream &s, int flag)
{
  s << "\nMacroelement3d, element id:  " << this->getTag() << endln;
  s << "\tConnected external nodes:  " << connectedExternalNodes;
}

// implement display self (only for plotting?)
/*
int
Macroelement3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numModes)
{
  static Vector v1(3);
  static Vector v2(3);

  if (displayMode >= 0) {

      theNodes[0]->getDisplayCrds(v1, fact);
      theNodes[1]->getDisplayCrds(v2, fact);

  } else {

    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);

    // add eigenvector values
    int mode = displayMode * -1;
    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 3; i++) {
	v1(i) += eigen1(i,mode-1)*fact;
	v2(i) += eigen2(i,mode-1)*fact;    
      }    
    }
  }
  return theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag());
}
*/

Response*
Macroelement3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

//	opserr << "El. " << this->getTag() << ": called set response\n";

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","Macroelement3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Fe_11");
      output.tag("ResponseType","Fe_21");
      output.tag("ResponseType","Fe_31");
      output.tag("ResponseType","Fe_12");
      output.tag("ResponseType","Fe_22");
      output.tag("ResponseType","Fe_32");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 1, P);

    // local force -
    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_2");
	  output.tag("ResponseType","N_3");
      output.tag("ResponseType","Mz_3");
      output.tag("ResponseType","My_3");
	  output.tag("ResponseType","V_y");
      output.tag("ResponseType","V_z");

      theResponse = new ElementResponse(this, 2, q);

    // in-plane drifts - shear and flexure
    }  else if (strcmp(argv[0],"drifts") == 0 || strcmp(argv[0],"drift") == 0 
	      || strcmp(argv[0],"driftSF") == 0) {

      output.tag("ResponseType","driftShear");
      output.tag("ResponseType","driftFlexure");

      theResponse = new ElementResponse(this, 3, Vector(2));


	// generic sectional response
	}  else if (strstr(argv[0],"section") != 0) {

    if (argc > 1) {
      
      int sectionNum = atoi(argv[1]);
      
      if (sectionNum > 0 && sectionNum <= 3 && argc > 2) {

		output.tag("SectionOutput");
		output.attr("number",sectionNum);

		theResponse = theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
		output.endTag();
      
      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 
	
		CompositeResponse *theCResponse = new CompositeResponse();
		int numResponse = 0;

		for (int i=0; i<numSections; i++) {
	  
			output.tag("SectionOutput");
			output.attr("number",i+1);
	  
			Response *theSectionResponse = theSections[i]->setResponse(&argv[1], argc-1, output);

			output.endTag();	  

			if (theSectionResponse != 0) {
				numResponse = theCResponse->addResponse(theSectionResponse);
			}
		}
	
		if (numResponse == 0) // no valid responses found
			delete theCResponse;
		else
			theResponse = theCResponse;
	  }
    }
  
	 
	// in-plane shear state  
    }  else if (strcmp(argv[0],"shear") == 0 || strcmp(argv[0],"shearSpring") == 0) {

		output.tag("Shear");
		theResponse = theShearModel[0]->setResponse(&argv[1], argc-1, output);
		output.endTag();



  } else if (strcmp(argv[0],"RayleighForces") == 0 || strcmp(argv[0],"rayleighForces") == 0) {

    theResponse =  new ElementResponse(this, 12, P);

  }   

  output.endTag();
  return theResponse;
}

int 
Macroelement3d::getResponse(int responseID, Information &eleInfo)
{

	//opserr << "El. " << this->getTag() << ": called get response\n";

  double N, V, M1, M2, T;
  double oneOverL = 1.0/L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());
    
  else if (responseID == 2) {
   // P = Tgl ^ (this->getResistingForce());
//	return eleInfo.setVector(P);

	this->getResistingForce();
	for (int i=0; i<12; i++) {
	  q(i) -= q0[i];
    }
	Vector basicForces(12);
	basicForces = q;

    return eleInfo.setVector(basicForces);
  }

  // Drifts
  else if (responseID == 3) {

     Vector driftsVec(2);
	 driftsVec(0) = driftS;
	 driftsVec(1) = driftF;
     
    return eleInfo.setVector(driftsVec);
  }

  // Crushing state variables
  else if (responseID == 4) {
     Vector shearStateVec(2);

	 Vector sectionResponse(4);




     
    return eleInfo.setVector(shearStateVec);
  }

  else
    return -1;
}

int 
Macroelement3d::trasformMatrixToGlobal(Matrix& A) {

	int m,j;
	double tmp[18][18];
	double kl[18][18];

	for (m=0; m<18; m++)
		for (j=0;j<18;j++) 
			kl[m][j] = A(m,j);


    static double RWI[3][3];
        
        if (nodeIOffset) {
            // Compute RWI
            RWI[0][0] = -R[0][1]*nodeIOffset[2] + R[0][2]*nodeIOffset[1];
            RWI[1][0] = -R[1][1]*nodeIOffset[2] + R[1][2]*nodeIOffset[1];
            RWI[2][0] = -R[2][1]*nodeIOffset[2] + R[2][2]*nodeIOffset[1];
            
            RWI[0][1] =  R[0][0]*nodeIOffset[2] - R[0][2]*nodeIOffset[0];
            RWI[1][1] =  R[1][0]*nodeIOffset[2] - R[1][2]*nodeIOffset[0];
            RWI[2][1] =  R[2][0]*nodeIOffset[2] - R[2][2]*nodeIOffset[0];
            
            RWI[0][2] = -R[0][0]*nodeIOffset[1] + R[0][1]*nodeIOffset[0];
            RWI[1][2] = -R[1][0]*nodeIOffset[1] + R[1][1]*nodeIOffset[0];
            RWI[2][2] = -R[2][0]*nodeIOffset[1] + R[2][1]*nodeIOffset[0];
        }
        
        static double RWJ[3][3];
        
        if (nodeJOffset) {
            // Compute RWJ
            RWJ[0][0] = -R[0][1]*nodeJOffset[2] + R[0][2]*nodeJOffset[1];
            RWJ[1][0] = -R[1][1]*nodeJOffset[2] + R[1][2]*nodeJOffset[1];
            RWJ[2][0] = -R[2][1]*nodeJOffset[2] + R[2][2]*nodeJOffset[1];
            
            RWJ[0][1] =  R[0][0]*nodeJOffset[2] - R[0][2]*nodeJOffset[0];
            RWJ[1][1] =  R[1][0]*nodeJOffset[2] - R[1][2]*nodeJOffset[0];
            RWJ[2][1] =  R[2][0]*nodeJOffset[2] - R[2][2]*nodeJOffset[0];
            
            RWJ[0][2] = -R[0][0]*nodeJOffset[1] + R[0][1]*nodeJOffset[0];
            RWJ[1][2] = -R[1][0]*nodeJOffset[1] + R[1][1]*nodeJOffset[0];
            RWJ[2][2] = -R[2][0]*nodeJOffset[1] + R[2][1]*nodeJOffset[0];
        }
        
        // Transform local stiffness to global system
        // First compute kl*T_{lg}

        for (m = 0; m < 18; m++) {
            tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
            tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
            tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];
            
            tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
            tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
            tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];
            
            if (nodeIOffset) {
                tmp[m][3]  += kl[m][0]*RWI[0][0]  + kl[m][1]*RWI[1][0]  + kl[m][2]*RWI[2][0];
                tmp[m][4]  += kl[m][0]*RWI[0][1]  + kl[m][1]*RWI[1][1]  + kl[m][2]*RWI[2][1];
                tmp[m][5]  += kl[m][0]*RWI[0][2]  + kl[m][1]*RWI[1][2]  + kl[m][2]*RWI[2][2];
            }
            
            tmp[m][6] = kl[m][6]*R[0][0] + kl[m][7]*R[1][0]  + kl[m][8]*R[2][0];
            tmp[m][7] = kl[m][6]*R[0][1] + kl[m][7]*R[1][1]  + kl[m][8]*R[2][1];
            tmp[m][8] = kl[m][6]*R[0][2] + kl[m][7]*R[1][2]  + kl[m][8]*R[2][2];
            
            tmp[m][9]  = kl[m][9]*R[0][0] + kl[m][10]*R[1][0] + kl[m][11]*R[2][0];
            tmp[m][10] = kl[m][9]*R[0][1] + kl[m][10]*R[1][1] + kl[m][11]*R[2][1];
            tmp[m][11] = kl[m][9]*R[0][2] + kl[m][10]*R[1][2] + kl[m][11]*R[2][2];
            
            if (nodeJOffset) {
                tmp[m][9]   += kl[m][6]*RWJ[0][0]  + kl[m][7]*RWJ[1][0]  + kl[m][8]*RWJ[2][0];
                tmp[m][10]  += kl[m][6]*RWJ[0][1]  + kl[m][7]*RWJ[1][1]  + kl[m][8]*RWJ[2][1];
                tmp[m][11]  += kl[m][6]*RWJ[0][2]  + kl[m][7]*RWJ[1][2]  + kl[m][8]*RWJ[2][2];
            }

			tmp[m][12] = kl[m][12]*R[0][0] + kl[m][13]*R[1][0]  + kl[m][14]*R[2][0];
            tmp[m][13] = kl[m][12]*R[0][1] + kl[m][13]*R[1][1]  + kl[m][14]*R[2][1];
            tmp[m][14] = kl[m][12]*R[0][2] + kl[m][13]*R[1][2]  + kl[m][14]*R[2][2];
            
            tmp[m][15] = kl[m][15]*R[0][0] + kl[m][16]*R[1][0] + kl[m][17]*R[2][0];
            tmp[m][16] = kl[m][15]*R[0][1] + kl[m][16]*R[1][1] + kl[m][17]*R[2][1];
            tmp[m][17] = kl[m][15]*R[0][2] + kl[m][16]*R[1][2] + kl[m][17]*R[2][2];
            
        }
        
        // Now compute T'_{lg}*(kl*T_{lg})
        for (m = 0; m < 18; m++) {
            K(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
            K(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
            K(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];
            
            K(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
            K(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
            K(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];
            
            if (nodeIOffset) {
                K(3,m) += RWI[0][0]*tmp[0][m]  + RWI[1][0]*tmp[1][m] + RWI[2][0]*tmp[2][m];
                K(4,m) += RWI[0][1]*tmp[0][m]  + RWI[1][1]*tmp[1][m] + RWI[2][1]*tmp[2][m];
                K(5,m) += RWI[0][2]*tmp[0][m]  + RWI[1][2]*tmp[1][m] + RWI[2][2]*tmp[2][m];
            }
            
            K(6,m) = R[0][0]*tmp[6][m] + R[1][0]*tmp[7][m]  + R[2][0]*tmp[8][m];
            K(7,m) = R[0][1]*tmp[6][m] + R[1][1]*tmp[7][m]  + R[2][1]*tmp[8][m];
            K(8,m) = R[0][2]*tmp[6][m] + R[1][2]*tmp[7][m]  + R[2][2]*tmp[8][m];
            
            K(9,m)  = R[0][0]*tmp[9][m] + R[1][0]*tmp[10][m] + R[2][0]*tmp[11][m];
            K(10,m) = R[0][1]*tmp[9][m] + R[1][1]*tmp[10][m] + R[2][1]*tmp[11][m];
            K(11,m) = R[0][2]*tmp[9][m] + R[1][2]*tmp[10][m] + R[2][2]*tmp[11][m];
            
            if (nodeJOffset) {
                K(9,m)  += RWJ[0][0]*tmp[6][m]  + RWJ[1][0]*tmp[7][m] + RWJ[2][0]*tmp[8][m];
                K(10,m) += RWJ[0][1]*tmp[6][m]  + RWJ[1][1]*tmp[7][m] + RWJ[2][1]*tmp[8][m];
                K(11,m) += RWJ[0][2]*tmp[6][m]  + RWJ[1][2]*tmp[7][m] + RWJ[2][2]*tmp[8][m];
            }

			K(12,m) = R[0][0]*tmp[12][m] + R[1][0]*tmp[13][m]  + R[2][0]*tmp[14][m];
            K(13,m) = R[0][1]*tmp[12][m] + R[1][1]*tmp[13][m]  + R[2][1]*tmp[14][m];
            K(14,m) = R[0][2]*tmp[12][m] + R[1][2]*tmp[13][m]  + R[2][2]*tmp[14][m];
            
            K(15,m) = R[0][0]*tmp[15][m] + R[1][0]*tmp[16][m] + R[2][0]*tmp[17][m];
            K(16,m) = R[0][1]*tmp[15][m] + R[1][1]*tmp[16][m] + R[2][1]*tmp[17][m];
            K(17,m) = R[0][2]*tmp[15][m] + R[1][2]*tmp[16][m] + R[2][2]*tmp[17][m];

        }

	return 0;
};

