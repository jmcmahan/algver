#include "main.hpp"
#include <queso/GammaVectorRV.h>
#include <queso/ConcatenationSubset.h>




/////////// Likelihood Routine Data Structure /////////// 

struct likelihoodRoutine_Data // user defined class
{
    likelihoodRoutine_Data(const QUESO::BaseEnvironment& env, algver *av);
    ~likelihoodRoutine_Data();


    algver *m_av;
    const QUESO::BaseEnvironment* m_env;
};


/////////// Likelihood Routine Data Structure Constructor /////////// 

likelihoodRoutine_Data::likelihoodRoutine_Data(const QUESO::BaseEnvironment& env, algver *av)
    :
    m_env    (&env),
    m_av     (av)
{

}

/////////// Likelihood Routine Data Structure Destructor /////////// 

likelihoodRoutine_Data::~likelihoodRoutine_Data()
{
}


/////////// Likelihood Routine Functions /////////// 
/* Several functions are included for the various cases. This allows the
   conditionsl to be moved outside fo the loop when assigning the callback. */


// Only beta unknown
double likelihoodRoutine_Beta(
    const QUESO::GslVector& paramValues,
    const QUESO::GslVector* paramDirection,
    const void*             functionDataPtr,
    QUESO::GslVector*       gradVector,
    QUESO::GslMatrix*       hessianMatrix,
    QUESO::GslVector*       hessianEffect)
{

    int i;
    const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    algver *av = ((likelihoodRoutine_Data*) functionDataPtr)->m_av;
    LAW_Vec *beta;
    double lambda;
    double *phi;
    unsigned int d;
    double rv;

    d = av->m_opt->d;

    // Compute model response 
    beta = new LAW_Vec(d + 1);
    for (i=0; i<d+1; i++) {
        beta->m_x[i] = paramValues[i];
    }
    av->Response(beta);
    delete beta;

    // Known parameters taken from the algver instance
    lambda = av->m_Rn->m_lambda;
    phi = av->m_Rn->m_phi;

    
    return av->m_Rn->LogLikelihood(lambda, phi, av->m_yn, av->m_zn);
}

// Beta, lambda unknown
double likelihoodRoutine_BetaLambda(
    const QUESO::GslVector& paramValues,
    const QUESO::GslVector* paramDirection,
    const void*             functionDataPtr,
    QUESO::GslVector*       gradVector,
    QUESO::GslMatrix*       hessianMatrix,
    QUESO::GslVector*       hessianEffect)
{

    int i;
    const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    algver *av = ((likelihoodRoutine_Data*) functionDataPtr)->m_av;
    LAW_Vec *beta;
    double lambda;
    double *phi;
    unsigned int d;

    d = av->m_opt->d;

    // Compute model response 
    beta = new LAW_Vec(d + 1);
    for (i=0; i<d+1; i++) {
        beta->m_x[i] = paramValues[i];
    }
    av->Response(beta);
    delete beta;

    // Lambda should be a parameter
    lambda = paramValues[d + 1];

    // Known parameters taken from the algver instance
    phi = av->m_Rn->m_phi;

    return av->m_Rn->LogLikelihood(lambda, phi, av->m_yn, av->m_zn);
}

// Beta, lambda, phi unknown, correlation other than Gaussian
double likelihoodRoutine_BetaLambdaPhiDefault(
    const QUESO::GslVector& paramValues,
    const QUESO::GslVector* paramDirection,
    const void*             functionDataPtr,
    QUESO::GslVector*       gradVector,
    QUESO::GslMatrix*       hessianMatrix,
    QUESO::GslVector*       hessianEffect)
{

    int i;
    const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    algver *av = ((likelihoodRoutine_Data*) functionDataPtr)->m_av;
    LAW_Vec *beta;
    double lambda;
    double phi[1];
    unsigned int d;

    d = av->m_opt->d;

    // Compute model response 
    beta = new LAW_Vec(d + 1);
    for (i=0; i<d+1; i++) {
        beta->m_x[i] = paramValues[i];
    }
    av->Response(beta);
    delete beta;

    // Lambda should be a parameter
    lambda = paramValues[d + 1];

    // Default case has 1-dimensional phi.
    phi[0] =  paramValues[d + 2];

    return av->m_Rn->LogLikelihood(lambda, phi, av->m_yn, av->m_zn);
}


// Beta, lambda, phi unknown, Gaussian correlation
double likelihoodRoutine_BetaLambdaPhiCorrGauss(
    const QUESO::GslVector& paramValues,
    const QUESO::GslVector* paramDirection,
    const void*             functionDataPtr,
    QUESO::GslVector*       gradVector,
    QUESO::GslMatrix*       hessianMatrix,
    QUESO::GslVector*       hessianEffect)
{

    int i;
    const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    algver *av = ((likelihoodRoutine_Data*) functionDataPtr)->m_av;
    LAW_Vec *beta;
    double lambda;
    double *phi;
    double ll;
    unsigned int d;

    d = av->m_opt->d;

    // Compute model response 
    beta = new LAW_Vec(d + 1);
    for (i=0; i<d+1; i++) {
        beta->m_x[i] = paramValues[i];
    }
    av->Response(beta);
    delete beta;

    // Lambda should be a parameter
    lambda = paramValues[d + 1];

    // This case has d-dimensional phi
    phi = new double[d];
    for (i=0; i<d; i++) {
        phi[i] = paramValues[d + 2 + i];
    }
    ll = av->m_Rn->LogLikelihood(lambda, phi, av->m_yn, av->m_zn);
    delete phi;

    return ll;
}



/////////// QUESO Setup and Execution /////////// 

void AlgVerQuesoSolve(const QUESO::FullEnvironment& env, algver *av)
{
    
    int i, j;
    int p, n, numUnknown;
    struct timeval timevalNow;
    LAW_Vec *beta;
    double lambda;
    double *phi;
    double v;
    LAW_Mat *XT;    // transpose of Xn
    LAW_Mat *XTX;   // X^T*X
    LAW_Vec *z;
    LAW_Vec *zp;
    algverPosterior *avp;


    n = av->m_opt->n;
    p = av->m_opt->d + 1;

  
    gettimeofday(&timevalNow, NULL);
    if (env.fullRank() == 0) {
        std::cout << "\nBeginning run of 'Linear Algebraic Verification' at "
            << ctime(&timevalNow.tv_sec)
            << "\n my fullRank = "         << env.fullRank()
            << "\n my subEnvironmentId = " << env.subId()
            << "\n my subRank = "          << env.subRank()
            << "\n my interRank = "        << env.inter0Rank()
            << std::endl << std::endl;
    }

    // Just examples of possible calls
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
        *env.subDisplayFile() 
            << "Beginning run of 'Linear Algebraic Verification' at "
            << ctime(&timevalNow.tv_sec)
            << std::endl;
    }
    env.fullComm().Barrier();
    env.subComm().Barrier();  // Just an example of a possible call
  
    //================================================================
    // Statistical inverse problem (SIP): find posterior PDF for 'g'
    //================================================================
    gettimeofday(&timevalNow, NULL);
    if (env.fullRank() == 0) {
        std::cout 
            << "Beginning 'SIP -> Linear Algebraic Verification' at "
            << ctime(&timevalNow.tv_sec)
            << std::endl;
    }

    //------------------------------------------------------
    // SIP Step 1 of 6: Instantiate the parameter space
    //------------------------------------------------------
    if (av->m_opt->kp == KnownLambdaPhi) {
        numUnknown = p;
    }
    else if (av->m_opt->kp == KnownPhi) {
        numUnknown = p+1;
    }
    else if (av->m_opt->kp == KnownNone && av->m_opt->ct != CorrGauss) {
        numUnknown = p+2;
    }

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        paramSpace(env, "param_", numUnknown, NULL);

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        paramSpaceBetaLambda(env, "beta_lambda_", p+1, NULL);

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        paramSpaceBeta(env, "beta_", p, NULL);
    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        paramSpaceLambda(env, "lambda_", 1, NULL);
    /* This phi is only valid for the correlation cases which are 
       non-Gaussian, since those have a 1-d phi. Need to change when
       that correlation is added. */
    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        paramSpacePhi(env, "phi", 1, NULL);

    //------------------------------------------------------
    // SIP Step 2 of 6: Instantiate the parameter domain
    //------------------------------------------------------
    QUESO::GslVector paramMinValues(paramSpace.zeroVector());
    QUESO::GslVector paramMaxValues(paramSpace.zeroVector());

    QUESO::GslVector paramMinValuesBetaLambda(paramSpaceBetaLambda.zeroVector());
    QUESO::GslVector paramMaxValuesBetaLambda(paramSpaceBetaLambda.zeroVector());

    QUESO::GslVector paramMinValuesBeta(paramSpaceBeta.zeroVector());
    QUESO::GslVector paramMaxValuesBeta(paramSpaceBeta.zeroVector());

    QUESO::GslVector paramMinValuesLambda(paramSpaceLambda.zeroVector());
    QUESO::GslVector paramMaxValuesLambda(paramSpaceLambda.zeroVector());

    QUESO::GslVector paramMinValuesPhi(paramSpacePhi.zeroVector());
    QUESO::GslVector paramMaxValuesPhi(paramSpacePhi.zeroVector());

    // Beta is always unknown, so fill in the domain
    for (i=0; i<p; i++) {
        paramMinValues[i] = av->m_opt->betamin;
        paramMaxValues[i] = av->m_opt->betamax;
        paramMinValuesBeta[i] = av->m_opt->betamin;
        paramMaxValuesBeta[i] = av->m_opt->betamax;
        paramMinValuesBetaLambda[i] = av->m_opt->betamin;
        paramMaxValuesBetaLambda[i] = av->m_opt->betamax;
    }

    if (av->m_opt->kp == KnownPhi || av->m_opt->kp == KnownNone) {
        paramMinValues[p] = av->m_opt->lambdamin;
        paramMaxValues[p] = av->m_opt->lambdamax;
    }
    paramMinValuesLambda[0] = av->m_opt->lambdamin;
    paramMaxValuesLambda[0] = av->m_opt->lambdamax;
    paramMinValuesBetaLambda[p] = av->m_opt->lambdamin;
    paramMaxValuesBetaLambda[p] = av->m_opt->lambdamax;



    /* Needs to be updated for Gaussian correlation case */
    if (av->m_opt->kp == KnownNone) {
        paramMinValues[p+1] = av->m_opt->phimin;
        paramMaxValues[p+1] = av->m_opt->phimax;
    }
    paramMinValuesPhi[0] = av->m_opt->phimin;
    paramMaxValuesPhi[0] = av->m_opt->phimax;
  
  
    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
        paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

 
    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
        paramDomainBeta("beta_", paramSpaceBeta, paramMinValuesBeta, paramMaxValuesBeta);
    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
        paramDomainLambda("lambda_", paramSpaceLambda, paramMinValuesLambda, paramMaxValuesLambda);
    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
        paramDomainPhi("phi_", paramSpacePhi, paramMinValuesPhi, paramMaxValuesPhi);

    //QUESO::ConcatenationSubset<QUESO::GslVector,QUESO::GslMatrix>
    //    paramDomainBetaLambda("beta_lambda_", paramSpaceBetaLambda, 
    //                    paramDomainBeta, paramDomainLambda);
    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
        paramDomainBetaLambda("beta_lambda_", paramSpaceBetaLambda, 
                        paramMinValuesBetaLambda, paramMaxValuesBetaLambda);
    //------------------------------------------------------
    // SIP Step 3 of 6: Instantiate the likelihood function 
    // object to be used by QUESO.
    //------------------------------------------------------

    likelihoodRoutine_Data likelihoodRoutine_Data(env, av);


    double (*likelihoodRoutinePtr)(const QUESO::GslVector& domainVector, 
                                   const QUESO::GslVector* domainDirection, 
                                   const void* routinesDataPtr, 
                                   QUESO::GslVector* gradVector, 
                                   QUESO::GslMatrix* hessianMatrix, 
                                   QUESO::GslVector* hessianEffect);

    if (av->m_opt->kp == KnownLambdaPhi) {
        likelihoodRoutinePtr = &likelihoodRoutine_Beta;
    }
    else if (av->m_opt->kp == KnownPhi) {
        likelihoodRoutinePtr = &likelihoodRoutine_BetaLambda;
    }
    else if (av->m_opt->kp == KnownNone) {
        if (av->m_opt->ct == CorrGauss) {
            likelihoodRoutinePtr = &likelihoodRoutine_BetaLambdaPhiCorrGauss;
        }
        else {
            likelihoodRoutinePtr = &likelihoodRoutine_BetaLambdaPhiDefault;
        }
    }

    QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
        likelihoodFunctionObj("like_", paramDomain, 
            likelihoodRoutinePtr, (void *) &likelihoodRoutine_Data, true); 
    

    //------------------------------------------------------
    // SIP Step 4 of 6: Define the prior RV
    //------------------------------------------------------


    /* Note - the prior cases we are using only affect whether beta has
       an informative or non-informative prior. The others are both 
       non-informative*/

    QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
                *BetaPriorNoninf;//("betaprior_", paramDomain);
    QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix> *BetaPriorInf;
    QUESO::GslVector meanVector( paramSpaceBeta.zeroVector() );
    QUESO::GslMatrix covMatrix = QUESO::GslMatrix(paramSpaceBeta.zeroVector());


    //if (av->m_opt->pc == PriorNonInform) {
        // Create an uniform prior RV for beta
    //    BetaPriorNoninf = 
    //        new QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    //                            ("betaprior_", paramDomain);
    //}
#if 1
    ///// NON-INFORMATIVE PRIOR CASE /////
    if (av->m_opt->pc == PriorNonInform) {
        // Create an uniform prior RV for beta
        BetaPriorNoninf = 
            new QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
                                ("betaprior_",paramDomainBeta);
    }

    ///// INFORMATIVE PRIOR CASE /////
    else {
        covMatrix.cwSet(0.0);
        for (i=0; i<av->m_opt->d+1; i++) {
            meanVector[i] = av->m_opt->mu->m_x[i];
            covMatrix(i, i) = (av->m_opt->q + av->m_opt->r[i])
                                / ( av->m_opt->q * av->m_opt->r[i] 
                                     * av->m_opt->lambda
                                );
        }
        BetaPriorInf = new QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
                ("betaprior_",paramDomainBeta,meanVector,covMatrix);
 
    }

    // Create an Jeffrey's prior RV for lambda
#if 1
    QUESO::JeffreysVectorRV<QUESO::GslVector,QUESO::GslMatrix> *LambdaPrior
         = new QUESO::JeffreysVectorRV<QUESO::GslVector,QUESO::GslMatrix>
                    ("lambdaprior_",paramDomainLambda, av);
#else
    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        parama(env, "whocares", 1, NULL);
    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> 
        paramb(env, "whocares", 1, NULL);
    QUESO::GslVector a(parama.zeroVector());
    QUESO::GslVector b(parama.zeroVector());
    a[0] = 1e-1;
    b[0] = 1e-1;

    std::cout << "a: " << a[0] << " b: " << b[0] << "\n";
    QUESO::GammaVectorRV<QUESO::GslVector,QUESO::GslMatrix> *LambdaPrior
         = new QUESO::GammaVectorRV<QUESO::GslVector,QUESO::GslMatrix>
                    ("lambdaprior_",paramDomainLambda, a, b);
#endif

    // Create an uniform prior RV for phi
    QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix> *PhiPrior
        = new QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>("phiprior_",paramDomainPhi);
 
#endif

    /* Assign the final prior according to the chosen options */

    QUESO::BaseVectorRV<QUESO::GslVector,QUESO::GslMatrix> *priorRv;
    QUESO::BaseVectorRV<QUESO::GslVector,QUESO::GslMatrix> *tempRv;

#if 1
    if (av->m_opt->pc == PriorNonInform) {
        if (av->m_opt->kp == KnownLambdaPhi) {
            priorRv = BetaPriorNoninf;
        }
        else if (av->m_opt->kp == KnownPhi) {
            priorRv = new QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
                      ("beta_lambda_", *BetaPriorNoninf, *LambdaPrior, paramDomain);
        }
        else if (av->m_opt->kp == KnownNone) {
            tempRv = new QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
                      ("beta_lambda_phi", *BetaPriorNoninf, *LambdaPrior, paramDomainBetaLambda);
            priorRv = new QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
                      ("beta_lambda_phi", *tempRv, *PhiPrior, paramDomain);
        }

    }
    /* This should all be cleaned up, but is ok for now since we're in a hurry. */
    else if (av->m_opt->pc == PriorInform) {
        if (av->m_opt->kp == KnownLambdaPhi) {
            priorRv = BetaPriorInf;
        }
        else if (av->m_opt->kp == KnownPhi) {
            //priorRv = new QUESO::InformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
            //                ("beta_lambda_", paramDomain, av);
            priorRv = new QUESO::InformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
                            ("beta_lambda_", paramDomain, av);
        }
        else if (av->m_opt->kp == KnownNone) {
            tempRv = new QUESO::InformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
                            ("beta_lambda_", paramDomainBetaLambda, av);
            //tempRv = new QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
            //          ("beta_lambda_phi", *BetaPriorInf, *LambdaPrior, paramDomainBetaLambda);
            priorRv = new QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix> 
                      ("beta_lambda_phi", *tempRv, *PhiPrior, paramDomain);
        }
    }
#endif

    //------------------------------------------------------
    // SIP Step 5 of 6: Instantiate the inverse problem
    //------------------------------------------------------
    QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
        postRv("post_", paramSpace);
        
    QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
        ip("", NULL, *priorRv, likelihoodFunctionObj, postRv); 

    //------------------------------------------------------
    // SIP Step 6 of 6: Solve the inverse problem, that is,
    // set the 'pdf' and the 'realizer' of the posterior RV
    //------------------------------------------------------
    std::cout << "Solving the SIP with Metropolis Hastings" 
	    << std::endl << std::endl;  

    QUESO::GslVector paramInitials(paramSpace.zeroVector());
    priorRv->realizer().realization(paramInitials);

    QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

    beta = new LAW_Vec(p);
    for (i=0; i<p; i++) {
        beta->m_x[i] = paramInitials[i];
    }
    av->Response(beta);
    delete beta;

    lambda = av->m_Rn->m_lambda;
    phi = av->m_Rn->m_phi;

    

    
    // Initial proposal covariance from least squares if more observations
    // than parameters. Otherwise 
    if (n > p) {
        z = new LAW_Vec(av->m_yn->m_N);
        av->m_yn->subtract(av->m_zn, z);
        v = z->inprod(z) / ((double)(n - p));
        delete z;
    
    
    
        XT = new LAW_Mat(av->m_Xn->m_N, av->m_Xn->m_M);
        for (i=0; i<av->m_Xn->m_M; i++) {
            for (j=0; j<av->m_Xn->m_N; j++) {
                XT->m_a[j*XT->m_N + i] = av->m_Xn->m_a[i*XT->m_M + j];
            }
        }
        XTX = new LAW_Mat(p, p);
        
        // Dividing by the variance is the right thing to do since we invert this 
        // next
        XT->lmult(av->m_Xn, XTX, 1.0/v, 0.0, false, false);
        XTX->cholinv();
    
        z = new LAW_Vec(av->m_yn->m_N, av->m_yn->m_x);
        /* Get least squares initial estimate. This is using the normal equations
           which you're supposed to avoid in practice, but it's probably alright
           since we don't really need a great guess anyway. */
    
        zp = new LAW_Vec(p);
        XT->matvec(z, zp);
        XTX->matvec(zp, zp);
        
    
        for (i=0; i<p; i++) {
            for (j=0; j<p; j++) {
                proposalCovMatrix(i,j) = XTX->m_a[i*p + j];
            }
            paramInitials[i] = zp->m_x[i];
        }
        delete XTX;
        delete XT;
        delete z;
        delete zp;
    }
    // Otherwise, we'll just set it like we do with the other parameters, to 
    // some fraction of the domain.
    else {
        for (i=0; i<p; i++) {
            proposalCovMatrix(i, i) = 0.0625*(av->m_opt->betamax - av->m_opt->betamin);
            paramInitials[i] = av->m_opt->betamin 
                        + 0.5*(av->m_opt->betamax - av->m_opt->betamin);
        }
    }






    if (av->m_opt->kp == KnownPhi || av->m_opt->kp == KnownNone) {
        proposalCovMatrix(p, p) = 0.0625*(av->m_opt->lambdamax - av->m_opt->lambdamin);
        paramInitials[p] = av->m_opt->lambdamin 
                        + 0.5*(av->m_opt->lambdamax - av->m_opt->lambdamin);
    }
    // NOTE: Need to do more than this if Gaussian correlation is used
    if (av->m_opt->kp == KnownNone) {
        proposalCovMatrix(p+1, p+1) = 0.0625*(av->m_opt->phimax - av->m_opt->phimin);
        paramInitials[p+1] = av->m_opt->phimin + 0.5*(av->m_opt->phimax - av->m_opt->phimin);
    }


    ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

    gettimeofday(&timevalNow, NULL);

    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
        *env.subDisplayFile() 
            << "Ending run of 'Linear Algebraic Verification' at "
            << ctime(&timevalNow.tv_sec)
            << std::endl;
    }
  
    if (env.fullRank() == 0) {
        std::cout 
            << "Ending run of 'Linear Algebraic Verification' at "
            << ctime(&timevalNow.tv_sec)
            << std::endl;
    }



    std::cout << "Generating posterior distribution info at "
            << ctime(&timevalNow.tv_sec) << "\n";



    avp = new algverPosterior(av);

    avp->MakeMatlab();
    avp->MakeDensity();

    std::cout << "Finished generating posterior info at "
            << ctime(&timevalNow.tv_sec) << "\n";

    delete avp;

}



int main(int argc, char **argv)
{
    std::string quesoOpt;
    algverOptions *opt; 
    algver *av;


    av = new algver(argc, argv);
#if 0
    av->m_Rn->Print("rn", std::cout);
    av->m_yn->Print("yn", std::cout);
    av->m_Xn->Print("xn", std::cout);
    av->m_opt->beta->Print("beta", std::cout);

    algverPosterior avp(av);

    avp.MakeMatlab();
    //LAW_Mat res(av->m_opt->d+1, av->m_opt->d+1);
    //avp.Xntr_Rni_Xn(&res);
    

    return 0;

#endif

    // The QUESO options file is chosen from the experiment name. A nicer
    // approach would be to set the QUESO options from the options file 
    // for the "algver" program, but this is less work.

    quesoOpt = av->m_opt->expname;
    quesoOpt += ".queso";

    MPI_Init(&argc, &argv);
    QUESO::FullEnvironment* env =
        new QUESO::FullEnvironment(MPI_COMM_WORLD, quesoOpt.c_str(), "", NULL);

    AlgVerQuesoSolve(*env, av);
    delete av;

    MPI_Finalize();


    return 0;
}

