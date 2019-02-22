#pragma once

#include "Vector.h"
#include "Matrix.h"
#include <string>


class CohesiveSurface
{
public:
	CohesiveSurface(void);
	CohesiveSurface(Matrix& Kpen, double c, double mu, double muR, double beta, double dropDrift, bool elasticSolution);
	// CohesiveSurface(Matrix& Kpen, double mu, double muR, double kStar, double ft, double sN0, double sNT0, double etaN, double etaNT, bool elasticSolution);
	~CohesiveSurface(void);

public:
	double evolveDamage(const Vector& s);
	void updatePlasticStrain(Vector& sigmaTrial);

	Matrix getAlgorithmicTangent();
	Matrix getElasticTangent();
	Matrix getCommittedTangent();
	Vector getSigma();

	int setTrialDisplacement(const Vector& s);

	double signTau(double tau);
	double heaviside(double arg);
	double macauley(double arg);

	double getMu(double kk);
	double getdMu_dk(double kk);

	void commit();
	CohesiveSurface* getCopy();
	double getMode(Vector& epsTrial);
	

protected:
	Matrix Kpen;

	double mu;
	double muR;
	double c;
	double Gc;
	double GfII;  
	/*
	double etaN, etaNT;
	double eta;
	double sN0, sNT0;
	double s02, sc2;
	double s01, sc1;
	*/

	// committed state variables
	Vector s_di, s_di_nplus1;
	double D, Dnplus1;
	double k;
	double redStiffness, redStiffnessComm;
	double x, xCommitted;

	// iteration stress and stiffness
	Vector sigma;
	Matrix Kalg;
	Matrix KCommitted;
	Vector dD_ds;
	double deltaLambda;
	bool flagAlg;
	bool elasticSolution;
};

