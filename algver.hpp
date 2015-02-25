#ifndef _ALGVER_HPP_
#define _ALGVER_HPP_

#include "linalgwrapper.hpp"
#include "rngwrapper.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>


enum CorrType {CorrNone, CorrEqual, CorrAR1, CorrGauss, CorrUndef};
enum KnownParams {KnownLambdaPhi, KnownPhi, KnownNone};
enum PriorCase {PriorNonInform, PriorInform};


class CovMatrix {
    public:
        unsigned int    m_n;        // dimension
        CorrType        m_ct;       // type of correlation
        double          m_lambda;   // correlation length
        double          *m_phi;     // correlation function parameter
        double          *m_chol;    // data dependent on m_ct
    
        void LeftInvMultMat(LAW_Mat *M, LAW_Mat *res);
        void LeftInvMultVec(LAW_Vec *v, LAW_Vec *res);
        void LeftMultMat(LAW_Mat *M, LAW_Mat *res);
        void LeftMultVec(LAW_Vec *v, LAW_Vec *res);
        void Setup(unsigned int, CorrType, double, double*);
        double LogLikelihood(double, double*, LAW_Vec*, LAW_Vec*);
        double LogLikelihood(LAW_Vec*, LAW_Vec*);
        double Likelihood(double, double*, LAW_Vec*, LAW_Vec*);
        double Likelihood(LAW_Vec*, LAW_Vec*);
        void ScaleNormal(double*, double, double*);
        void ScaleNormal(double*);
        void Print(const char *varname, std::ostream &os);
        void InvPrint(const char *varname, std::ostream &os);
        double Det(double *phi);
        double LogDet(double *phi);


        // Constructors / Destructors
        CovMatrix(unsigned int, CorrType, double, double*);
        ~CovMatrix();
 

    private:
        unsigned int    m_len;        // length of data in m_chol
};


class algverOptions {
    public:
        char expname[256];              // string for experiment name
        unsigned int n;             // sample size
        unsigned int d;             // number of covariates
        LAW_Vec *beta;
        double lambda;
        double *phi;
        LAW_Vec *mu;
        CorrType ct;
        KnownParams kp;
        int seed; 

        // Parameters for the covariance for Xn. For simplicity it has to be
        // CorrNone, CorrEqual, or CorrAR1. 
        double x_lambda;
        double *x_phi;
        CorrType x_ct;

        PriorCase pc;
        double q;                   // Sigma_0 parameter for informative prior
        double *r;                  // Sigma parameter for non-informative prior

        // Domain boundaries. For simplicity, just one is used for all beta 
        // parameters
        double betamax;
        double betamin;
        double lambdamax;
        double lambdamin;
        double phimax;
        double phimin;


        void Display();             // Dump the options settings for comparison
        // Constructors / Destructors
        algverOptions(unsigned int n0, unsigned int d0, LAW_Vec *beta0, 
                      double lambda0, double *phi0, LAW_Vec *mu0, 
                      CorrType ct0, KnownParams kp0, int seed0, 
                      double x_lambda0, double *x_phi0, CorrType x_ct0,
                      PriorCase pc0, double q0, double *r0,
                      double betamax0, double betamin0,
                      double lambdamax0, double lambdamin0,
                      double phimax0, double phimin, char *expname0);
        algverOptions(char **argv);
        ~algverOptions();
};


class algver {
    public:
        algverOptions *m_opt;
        LAW_Vec     *m_yn;          // generated data
        LAW_Vec     *m_zn;          // no-noise output
        LAW_Mat     *m_Xn;          // regression matrix
        CovMatrix   *m_Rn;          // measurement error covariance


        void Response(LAW_Vec*);
        void Setup();
        void PrintPrior(char *varname, std::ostream &os);
        void PrintParams(char *betaname, char *lambdaname, char *phiname, std::ostream &os);
        // Constructors / Destructors
        algver(algverOptions*);
        algver(int argc, char **argv);
        ~algver();
};


// Class for computing the parameters for the 
// exact posteriors of a given algver instance 
class algverPosterior {
    public:
        algver          *m_av;
        // M x M where M is the number of beta parameters
        LAW_Mat         *tmpsquaremat;
        // N x M where N is number of samples, M is number of beta parameters
        LAW_Mat         *tmprectmat;
        // M x 1 where M is the number of beta parameters
        LAW_Vec         *tmpvec;


        void betahat(LAW_Vec *bhat); 
        void Xntr_Rni_Xn(LAW_Mat *res);

        // Beta prior info
        void            sigma2(LAW_Mat *s2);
        void            mu2(LAW_Vec *m2, LAW_Mat *s2, LAW_Vec *bh);
        void            PrintBeta(char *varname, std::ostream os);

        // Lambda prior info
        double          a1();
        double          b1();
        void            PrintLambda(char *varname, std::ostream os);

        // Phi prior info
        void            PrintPhi(char *varname, std::ostream os);

        // Generate the densities
        void            MakeDensity();

        // Generate a MATLAB file to generate the distribution parameters

        void            MakeMatlab();

        // Constructors / Destructors
        algverPosterior(algver *av);
        ~algverPosterior();
};

#endif /* _ALGVER_HPP_ */
