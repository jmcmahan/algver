#include "algver.hpp"



void printUsageAndExit(void)
{
    std::cout << "Usage: algver [-default] options-file\n";
    exit(EXIT_FAILURE);
}

///////////////////////////////////////////////////////
// CovMatrix member functions
///////////////////////////////////////////////////////


// Compute res = inv(R)*M where R is the covariance matrix
void CovMatrix::LeftInvMultMat(LAW_Mat *M, LAW_Mat *res)
{
    int i, j, k;

    double f = *m_phi;
    double a, b;

    if (res->m_M != m_n || M->m_M != m_n || res->m_N != M->m_N ) {
        std::cout << "Wrong size in CovMatrix left matrix multiplication\n";
        printUsageAndExit();
    }

    switch (m_ct) {
        case CorrNone:
            memcpy(res->m_a, M->m_a, (M->m_M)*(M->m_N)*sizeof(double));
            break;

        case CorrEqual:
            // Here a and b are coefficients - a is for the diagonal terms, b is for
            // the rest
            a = 1.0 / (1.0 - f)  - 1.0 / (1.0 - f) * f / (1.0 + (m_n-1.0)*f );
            b = - 1.0 / (1.0 - f) * f / (1.0 + (m_n-1.0)*f );
            for (i=0; i<res->m_M; i++) {
                for (j=0; j<res->m_N; j++) {
                    *(res->ind(i,j)) = 0.0;

                    for (k=0; k<m_n; k++) {
                        *(res->ind(i,j)) +=  ( a*(k == i) +  b*(k != i) )  
                                               *   (*(M->ind(k,j))); 
                    }
                }
            }


            break;

        case CorrAR1:
            a = 1.0 / (1-f*f);
            // This is a tridiagonal matrix, so this can certainly be made more 
            // efficient if needed
            for (i=0; i<res->m_M; i++) {
                for (j=0; j<res->m_N; j++) {
                    *(res->ind(i,j)) = 0.0;
                    for (k=0; k<m_n; k++) {
                        *(res->ind(i,j)) += a * ( - (i-k == 1 || i-k == -1         ) * f 
                                                  + (i==k && k != 0 && k!=m_n-1    ) * (1.0+f*f) 
                                                  + (i==k && (k == 0 || k == m_n-1)) * 1.0 )
                                               *   (*(M->ind(k,j))); 
                    }
                }
            }

            break;

        default:
            break;
    }
}

void CovMatrix::LeftInvMultVec(LAW_Vec *v, LAW_Vec *res)
{
    int i, j;
    double a, b;

    double f = *m_phi;

    if (res->m_N != m_n || v->m_N != m_n) {
        std::cout << "Wrong size in CovMatrix left vector multiplication\n";
        printUsageAndExit();
    }

    switch (m_ct) {
        case CorrNone:
            memcpy(res->m_x, v->m_x, m_n*sizeof(double));
            break;

        case CorrEqual:
            // Here a and b are coefficients - a is for the diagonal terms, b is for
            // the rest
            a = 1.0 / (1.0 - f)  - 1.0 / (1.0 - f) * f / (1.0 + (m_n-1.0)*f );
            b = - 1.0 / (1.0 - f) * f / (1.0 + (m_n-1.0)*f );
            for (i=0; i<m_n; i++) {
                res->m_x[i] = 0.0;
                for (j=0; j<m_n; j++) {
                    res->m_x[i] +=  ( a*(j == i) +  b*(j != i) )  
                                               *   (v->m_x[j]);
                }
            }
            break;

        case CorrAR1:
            // Tridiagonal so can be made more efficient if needed
            a = 1.0 / (1-f*f);
            for (i=0; i<m_n; i++) {
                res->m_x[i] = 0.0;
                for (j=0; j<m_n; j++) {
                    res->m_x[i] += a * ( - (i-j == 1 || i-j == -1         ) * f 
                                     + (i==j && j != 0 && j!=m_n-1    ) * (1.0+f*f) 
                                     + (i==j && (j == 0 || j == m_n-1)) * 1.0 )
                                         *   (v->m_x[j]);
                }
            }
            break;

    }

}



// Compute res = R*M where R is the covariance matrix
void CovMatrix::LeftMultMat(LAW_Mat *M, LAW_Mat *res)
{
    int i, j, k;

    double f = *m_phi;

    if (res->m_M != m_n || M->m_M != m_n || res->m_N != M->m_N ) {
        std::cout << "Wrong size in CovMatrix left matrix multiplication\n";
        printUsageAndExit();
    }

    switch (m_ct) {
        case CorrNone:
            memcpy(res->m_a, M->m_a, (M->m_M)*(M->m_N)*sizeof(double));
            break;

        case CorrEqual:
            for (i=0; i<res->m_M; i++) {
                for (j=0; j<res->m_N; j++) {
                    *(res->ind(i,j)) = 0.0;

                    for (k=0; k<m_n; k++) {
                        *(res->ind(i,j)) +=  ( (k == i) +  (k != i) * f )  
                                               *   (*(M->ind(k,j))); 
                    }
                }
            }


            break;

        case CorrAR1:
            for (i=0; i<res->m_M; i++) {
                for (j=0; j<res->m_N; j++) {
                    *(res->ind(i,j)) = 0.0;
                    for (k=0; k<m_n; k++) {
                        *(res->ind(i,j)) += std::pow(f, std::abs(k-i))
                                               *   (*(M->ind(k,j))); 
                    }
                }
            }

            break;

        default:
            break;
    }
}

void CovMatrix::LeftMultVec(LAW_Vec *v, LAW_Vec *res)
{
    int i, j;

    double f = *m_phi;
    if (res->m_N != m_n || v->m_N != m_n) {
        std::cout << "Wrong size in CovMatrix left vector multiplication\n";
        printUsageAndExit();
    }

    switch (m_ct) {
        case CorrNone:
            memcpy(res->m_x, v->m_x, m_n*sizeof(double));
            break;

        case CorrEqual:
            for (i=0; i<m_n; i++) {
                res->m_x[i] = 0.0;
                for (j=0; j<m_n; j++) {
                    res->m_x[i] += ( (j == i) +  (j != i) * f )  
                                               *   (v->m_x[j]);
                }
            }
            break;

        case CorrAR1:
            for (i=0; i<m_n; i++) {
                res->m_x[i] = 0.0;
                for (j=0; j<m_n; j++) {
                    res->m_x[i] += std::pow(f, std::abs(j-i))
                                               *   (v->m_x[j]);
                }
            }
            break;

    }

}

double CovMatrix::Det(double *phi)
{
    return std::exp(LogDet(phi));
}


double CovMatrix::LogDet(double *phi)
{
    switch (m_ct) {
        case CorrNone:
            return 0.0;
            break;

        case CorrEqual:
            return (m_n - 1.0) * std::log(1.0 - *phi) 
                   + std::log(1.0 + (m_n - 1.0) * (*phi));
            break;

        case CorrAR1:
            return (m_n - 1.0) * std::log(1.0 - (*phi) * (*phi));
            break;

        default:
            break;
    }
}


/* Computes the log-likelihood or least squares weighted by the covariance
   matrix. */
double CovMatrix::LogLikelihood(double lambda, double *phi, LAW_Vec *y, LAW_Vec *yest)
{
    int i;
    double v;           // result
    double retval;
    double p = *phi;    
    double a, b, s;     // helper variables
    double ap, bp;
    double xp, yp;  
    double *zptr;       // pointer to residual array
    LAW_Vec *z;         // stores residual array
    double logdet;


    logdet = LogDet(phi);

    z = new LAW_Vec(y->m_N);
    zptr = z->m_x;


    
    // z = yest - y
    y->subtract(yest, z);



    switch (m_ct) {
        case CorrNone:
            // v = z^T * z
            v = z->inprod(z);
            break;


        case CorrEqual:
            // Initializes the recursion
            a = 1.0;
            b = p;

            // Note : inv(L*L') = inv(L')*inv(L) so this is the right order to do
            // this (although it's worth checking to see if the inversion is being
            // done correctly).
            xp = *zptr / a;
            s = b * xp;
            v = xp * xp;
            
            ap = a;
            bp = b;
            zptr++;

            for (i=1; i<y->m_N; i++) {
                // a = diagonals, b = off-diagonals of current column of
                // Cholesky factor
                a = sqrt(ap*ap - bp*bp);
                b = (ap - bp) * bp / a;

                // Computes current component of inverse Cholesky applied
                // to residual
                xp = (*zptr - s) / a;

                s += b * xp;

                // Sum-of-squares. Note comments below about possible
                // speed-up if needed apply. 
                v += xp * xp;

                ap = a;
                bp = b;
                zptr++;
            }

            break;


        case CorrAR1:
            /* Uses a recursion to apply the inverse of the Cholesky 
               factor and takes sum-of-squares. This avoids having to
               represent the matrix and is O(N) for N = sample size. */
            a = 1.0 / sqrt(1.0 - p * p);
            // First element the same
            xp = *zptr; 
            
            // Note - it's possible it's faster to store the inverted
            // vector as a LAW_Vec and then use the inprod method rather
            // than computing the sum-of-squares manually here. Worth
            // trying if this is slow.
            v = xp * xp;
            zptr++;
            for (i=1; i<y->m_N; i++) {
                // Bad form, but xp is used for two different things here.
                // First, temporary storage to compute the square below.
                yp = *zptr;
                xp = (yp - p*xp) * a; 
                v += xp * xp;

                // Now use to hold previous value.
                xp = yp;
                zptr++;
            }
            break;
        default:
            break;
    }

    delete z;

    retval =  0.5*m_n*std::log(lambda / (2.0*M_PI)) 
             -0.5*logdet
             -0.5*lambda*v;
#if 0 
    retval =  0.5*m_n*std::log(lambda / (2.0*M_PI)) 
             -0.5*m_n*std::log(2.0*M_PI) 
             -0.5*std::log(det)
             -0.5*lambda*v;
#endif

    return retval;
}




/* Computes the log-likelihood or least squares weighted by the covariance
   matrix. This one uses the members as arguments. */
double CovMatrix::LogLikelihood(LAW_Vec *y, LAW_Vec *yest)
{
    return LogLikelihood(m_lambda, m_phi, y, yest);
}


/* These two methods return the likelihood, not its logarithm */
double CovMatrix::Likelihood(double lambda, double *phi, LAW_Vec *y, LAW_Vec *yest)
{
    return exp(LogLikelihood(lambda, phi, y, yest));
}

/* This one does it with the member parameters */
double CovMatrix::Likelihood(LAW_Vec *y, LAW_Vec *yest)
{
    return exp(LogLikelihood(m_lambda, m_phi, y, yest));
}



/* Transform standard normal data according to the set covariance matrix. */
void CovMatrix::ScaleNormal(double *x, double lambda, double *phi)
{
    int i;
    double li = 1.0 / sqrt(lambda);

    switch(m_ct) {
        case CorrNone:
            for (i=0; i<m_n; i++) {
                x[i] *= li;
            }
            break;


        case CorrEqual:
            double d0, d1, o0, o1;   // off-diagonal, diagonal components
            double x0;
            double sum;
            // Initialize the recursion. Comes from the Cholesky 
            // factorization.
            sum = 0.0;
            d0 = 1.0;
            o0 = *phi;

            // First iteration
            d1 = d0;
            o1 = o0;
            x0 = x[0];
            x[0] = (d1*x0 + sum) * li;
            sum += o1*x0;

            // Rest of the iterations
            for (i=1; i<m_n; i++) {
                d1 = sqrt(d0*d0 - o0*o0);
                o1 = (d0 - o0)*o0 / d1; // Note only n-1 of these are used.
                x0 = x[i];
                x[i] = (d1*x0 + sum) * li;
                sum += o1*x0;
                d0 = d1;
                o0 = o1;
            }

             
            break;

        case CorrAR1:
            double a, b;
            a = sqrt( 1.0 - (*phi) * (*phi) );
            b = *phi;
            x0 = x[0];
            x[0] = x[0] * li;
            for (i=1; i<m_n; i++) {
                x0 = b*x0 + a*x[i];
                x[i] = x0 * li;
            }
            

            break;

        default:
            break;
    }

}

/* Scale the normal data using the member parameters. */
void CovMatrix::ScaleNormal(double *x)
{
    ScaleNormal(x, m_lambda, m_phi);
}


/* Initializes members */
void CovMatrix::Setup(unsigned int n, CorrType ct, double lambda, double *phi)
{
    m_n = n;
    m_ct = ct;
    m_lambda = lambda;
    m_phi = phi;

    switch (m_ct) {
        // Don't really use m_chol for these cases for now. If things
        // are slow when computing the likelihoods, this can be used to
        // pre-compute stuff that depends on phi for the cases where phi
        // is known.
        case CorrNone:
        case CorrEqual:
        case CorrAR1:
            m_len = 1;
            break;
        default:
            break;
    }

    m_chol = new double[m_len];
}


/* Output MATLAB code for verifying with varname being the variable name */
void CovMatrix::Print(const char *varname, std::ostream &os)
{

    switch(m_ct) {
        case CorrNone:
            os << varname << "=" << "eye(" << m_n 
               << ");lambda=" << m_lambda << ";\n";
            break;

        case CorrEqual:
            os << varname << "=(" << *m_phi << "*ones(" << m_n << ")+(1-"
               << *m_phi << ")*eye(" << m_n << "));lambda=" << m_lambda << ";\n";
            break;

        case CorrAR1:
            os << varname << "=zeros(" << m_n << ");for i=0:" << m_n << "-1;" 
               << varname << "(i+1,:)=" << *m_phi << ".^abs(-i:" << m_n 
               << "-i-1);lambda=" << m_lambda << ";end;\n";
            break;

        default:
            break;
    }
    
}

void CovMatrix::InvPrint(const char *varname, std::ostream &os)
{

    double x, y;
    switch(m_ct) {
        case CorrNone:
            os << varname << "=" << "eye(" << m_n 
               << ");lambda=" << m_lambda << ";\n";
            break;

        case CorrEqual:
            x = m_phi[0] / (1.0 + (m_n - 1.0) * (m_phi[0]));
            y = 1.0 / (1.0 - m_phi[0]);

            os << varname << "=" << y << "*(eye(" << m_n << ") - " 
               << x <<"*ones(" << m_n <<"));lambda=" 
               << m_lambda << ";\n"; 

            break;

        case CorrAR1:
            y = m_phi[0] * m_phi[0];
            x = 1.0 / (1.0 - y);
            y = 1.0 + y;

            os << varname << "=diag(-" << m_phi[0] << "*ones(" << m_n - 1
               << ",1),1);" << varname << "=" << varname << "+" 
               << varname << "' + diag(" << y << "*ones(" << m_n 
               << ",1));"  << varname << "(1,1)=1;" << varname 
               << "(end,end)=1;" 
               << varname << "=" << varname << "*" << x << ";lambda="
               << m_lambda << ";\n";
            break;

        default:
            break;
    }
    
}

// Constructors / Destructors

CovMatrix::CovMatrix(unsigned int n, CorrType ct, double lambda, double *phi)
{
    Setup(n, ct, lambda, phi);
}


CovMatrix::~CovMatrix()
{
    delete m_chol;
}


///////////////////////////////////////////////////////
// algverOptions member functions
///////////////////////////////////////////////////////


/* Display the set options for debugging */
void algverOptions::Display()
{
    int i, j;
    std::cout << "Experiment Name: " << expname << "\n";
    std::cout << "Sample Size: " << n << "\n";
    std::cout << "Number of Beta Parameters: " << d + 1 << "\n";
    std::cout << "Prior case: ";
    if (pc == PriorNonInform) {
        std::cout << "Non-informative\n";
    }
    else if (pc == PriorInform) {
        std::cout << "Informative Gaussian\n";
        std::cout << "Prior parameter q: " << q << "\n";
        std::cout << "Prior parameters r[i] (up to first 5): ";
        if (d+1 > 5) {
            j = 5;
        }
        else {
            j = d+1;
        }
        for (i=0; i<j; i++) {
            std::cout << r[i] << " ";
        }
        std::cout << "\n";
    }
    else {
        std::cout << "Not Defined (should not happen)\n";
    }
    std::cout << "Lambda: " << lambda << "\n";
    std::cout << "Phi: " << *phi << "\n";
    std::cout << "Observation Noise Correlation Type: ";
    if (ct == CorrNone) {
        std::cout << "None\n";
    }
    else if (ct == CorrEqual) {
        std::cout << "Equal\n";
    }
    else if (ct == CorrAR1) {
        std::cout << "AR(1)\n";
    }
    else {
        std::cout << "Not Defined (should not happen)\n";
    }
    std::cout << "Known Parameters: ";
    if (kp == KnownLambdaPhi) {
        std::cout << "Lambda, Phi\n";
    }
    else if (kp == KnownPhi) {
        std::cout << "Phi\n";
    }
    else if (kp == KnownNone) {
        std::cout << "None\n";
    }
    else {
        std::cout << "Not Defined (should not happen)\n";
    }
    std::cout << "Seed (for observation error and regression matrix): " << seed << "\n";
    std::cout << "Regression Matrix Lambda: " << x_lambda << "\n";
    std::cout << "Regression Matrix Phi: " << *x_phi << "\n";
    std::cout << "Regression Matrix Correlation Type: ";
    if (x_ct == CorrNone) {
        std::cout << "None\n";
    }
    else if (x_ct == CorrEqual) {
        std::cout << "Equal\n";
    }
    else if (x_ct == CorrAR1) {
        std::cout << "AR(1)\n";
    }
    else {
        std::cout << "Not Defined (should not happen)\n";
    }
    std::cout << "Beta Domain: [" << betamin << "," << betamax << "]\n";
    std::cout << "Lambda Domain: [" << lambdamin << "," << lambdamax << "]\n";
    std::cout << "Phi Domain: [" << phimin << "," << phimax << "]\n";
    
}

// Constructors / Destructors
/* Directly set the options */
algverOptions::algverOptions(unsigned int n0, unsigned int d0, LAW_Vec *beta0, 
                             double lambda0, double *phi0, LAW_Vec *mu0, 
                             CorrType ct0, KnownParams kp0, int seed0, 
                             double x_lambda0, double *x_phi0, CorrType x_ct0,
                             PriorCase pc0, double q0, double *r0,
                             double betamax0, double betamin0, double lambdamax0,
                             double lambdamin0, double phimax0, double phimin0,
                             char *expname0)
{
    n = n0;
    d = d0;

    beta = new LAW_Vec(d, beta0->m_x);

    lambda = lambda0;
    phi = phi0;

    if (ct0 != CorrGauss) {
        phi = new double[1];
        *phi = *phi0;
    }

    mu = new LAW_Vec(d, mu0->m_x);

    ct = ct0;
    kp = kp0;
    seed = seed0;

    x_lambda = x_lambda0;

    /* Don't allow Gaussian correlation in this case so just make array of one */
    x_phi = new double[1];
    *x_phi = *x_phi0;

    x_ct = x_ct0;

    pc = pc0;
    q = q0;

    r = new double[d];
    memcpy(r, r0, sizeof(double)*d);

    betamax = betamax0;
    betamin = betamin0;
    lambdamax = lambdamax0;
    lambdamin = lambdamin0;
    phimax = phimax0;
    phimin = phimin0;
    strcpy(expname, expname0);
}

/* Load options from file. argv should be the commandline used to
   run the program with argv[1] specifying which type of options file
   to use (-default is a simplified version) and argv[2] giving the
   name of the options file */ 
algverOptions::algverOptions(char **argv)
{
    int i, j;
    int failed;
    std::string optfile(argv[2]);
    std::string ExpNameString;
    std::string stringOpt;
    std::string stringComp;
    double betaAmp;
    double betaBias;
    double muAmp;
    double muBias;

    stringComp = "-default";

    if (!stringComp.compare(argv[1])) {
        /* Default problem setup fills in the true beta parameters as one cycle
           of a sinusoid */
        std::cout << "\n*****************************************\n";
        std::cout << "Default setup selected.\n";

        std::ifstream optin(argv[2]);
        if (!optin) {
            std::cerr << "Error: can't open options file \"" << optfile <<"\".\n";
            printUsageAndExit();
        }
        
        /* Used as a default file name when other files hold things like variable
           data */
        optin >> stringComp;    // For now just get rid of these
        if (stringComp.compare("ExperimentName:")) {
            std::cerr << "Read in option: " << stringComp << "\n";
            std::cerr << "Error: Can't parse \"ExperimentName:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> ExpNameString;
        strcpy(expname, ExpNameString.c_str());


        /* Sample size */
        optin >> stringComp;
        if (stringComp.compare("SampleSize:")) {
            std::cerr << "Error: Can't parse \"SampleSize:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> n; 

        /* Number of parameters in Beta */
        optin >> stringComp;
        if (stringComp.compare("NumberBetaParams:")) {
            std::cerr << "Error: Can't parse \"NumberBetaParams:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> d;
        /* d is the dimension of the co-variates in the regression matrix, which is one
           less than the number of beta parameters */
        d = d - 1;

        /* Informative / Non-informative case */
        optin >> stringComp;
        if (stringComp.compare("PriorCase:")) {
            std::cerr << "Error: Can't parse \"PriorCase:\" option (check for typos)\n";
            printUsageAndExit();
        }

        failed = 1;
        optin >> stringOpt;
        stringComp = "informative";
        if (!stringComp.compare(stringOpt)) {
            pc = PriorInform;
            r = new double[d];

            /* r and q get filled in below */

            failed = 0;
        }

        stringComp = "noninformative";
        if (!stringComp.compare(stringOpt)) {
            pc = PriorNonInform;

            // This isn't used in this case, but avoids having to special when
            // destructing.
            r = new double[1];
            failed = 0;
        }
        if (failed) {
            std::cerr << "Error: invalid prior case specified.\n";
            printUsageAndExit();
        }

        /* Beta amplitude */
        optin >> stringComp;
        if (stringComp.compare("BetaAmplitude:")) {
            std::cerr << "Error: Can't parse \"BetaAmplitude:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> betaAmp;


        /* Beta bias */
        optin >> stringComp;
        if (stringComp.compare("BetaBias:")) {
            std::cerr << "Error: Can't parse \"BetaBias:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> betaBias;

        
        /* Lambda */
        optin >> stringComp;
        if (stringComp.compare("Lambda:")) {
            std::cerr << "Error: Can't parse \"Lambda:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> lambda;

        /* Phi */
        optin >> stringComp;
        if (stringComp.compare("Phi:")) {
            std::cerr << "Error: Can't parse \"Phi:\" option (check for typos)\n";
            printUsageAndExit();
        }

        // For the default case, no Gaussian correlation so just 1-dim phi
        phi = new double[1];
        optin >> *phi;


        /* Correlation Type */
        optin >> stringComp;
        if (stringComp.compare("CorrelationType:")) {
            std::cerr << "Error: Can't parse \"CorrelationType:\" option (check for typos)\n";
            printUsageAndExit();
        }
        failed = 1;
        optin >> stringComp;
        if (!stringComp.compare("none")) {
            ct = CorrNone;
            failed = 0;
        }
        if (!stringComp.compare("equal")) {
            ct = CorrEqual;
            failed = 0;
        }
        if (!stringComp.compare("ar1")) {
            ct = CorrAR1;
            failed = 0;
        }
        if (failed) {
            std::cerr << "Error: invalid correlation type specified.\n";
            printUsageAndExit();
        }
        

        /* Known Parameters */
        optin >> stringComp;
        if (stringComp.compare("KnownParams:")) {
            std::cerr << "Error: Can't parse \"KnownParams:\" option (check for typos)\n";
            printUsageAndExit();
        }
        failed = 1;
        optin >> stringComp;


        if (!stringComp.compare("none")) {
            kp = KnownNone;
            failed = 0;
        }
        if (!stringComp.compare("lambdaphi")) {
            kp = KnownLambdaPhi;
            failed = 0;
        }
        if (!stringComp.compare("phi")) {
            kp = KnownPhi;
            failed = 0;
        }
        if (failed) {
            std::cerr << "Error: invalid known parameters specified.\n";
            printUsageAndExit();
        }

        /* Seed */
        optin >> stringComp;
        if (stringComp.compare("Seed:")) {
            std::cerr << "Error: Can't parse \"Seed:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> seed;

        /* Regression Matrix Lambda */
        optin >> stringComp;
        if (stringComp.compare("RegressionMatrixLambda:")) {
            std::cerr << "Error: Can't parse \"RegressionMatrixLambda:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> x_lambda;

        /* Regression Matrix Phi */
        optin >> stringComp;
        if (stringComp.compare("RegressionMatrixPhi:")) {
            std::cerr << "Error: Can't parse \"RegressionMatrixPhi:\" option (check for typos)\n";
            printUsageAndExit();
        }

        /* No Gaussian correlation for the regression matrix so just 1-dim phi */
        x_phi = new double[1];
        optin >> *x_phi;

        /* Regression Matrix Correlation Type */
        optin >> stringComp;
        if (stringComp.compare("RegressionMatrixCorrelationType:")) {
            std::cerr << "Error: Can't parse \"RegressionMatrixCorrelationType:\" option (check for typos)\n";
            printUsageAndExit();
        }
        failed = 1;
        optin >> stringComp;
        if (!stringComp.compare("none")) {
            x_ct = CorrNone;
            failed = 0;
        }
        if (!stringComp.compare("equal")) {
            x_ct = CorrEqual;
            failed = 0;
        }
        if (!stringComp.compare("ar1")) {
            x_ct = CorrAR1;
            failed = 0;
        }
        if (failed) {
            std::cerr << "Error: invalid regression matrix correlation type specified.\n";
            printUsageAndExit();
        }

        /* BetaMax */
        optin >> stringComp;
        if (stringComp.compare("BetaMax:")) {
            std::cerr << "Error: Can't parse \"BetaMax:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> betamax;

        /* BetaMin */
        optin >> stringComp;
        if (stringComp.compare("BetaMin:")) {
            std::cerr << "Error: Can't parse \"BetaMin:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> betamin;

        /* LambdaMax */
        optin >> stringComp;
        if (stringComp.compare("LambdaMax:")) {
            std::cerr << "Error: Can't parse \"LambdaMax:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> lambdamax;

        /* LambdaMin */
        optin >> stringComp;
        if (stringComp.compare("LambdaMin:")) {
            std::cerr << "Error: Can't parse \"LambdaMin:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> lambdamin;

        /* PhiMax */
        optin >> stringComp;
        if (stringComp.compare("PhiMax:")) {
            std::cerr << "Error: Can't parse \"PhiMax:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> phimax;

        /* PhiMin */
        optin >> stringComp;
        if (stringComp.compare("PhiMin:")) {
            std::cerr << "Error: Can't parse \"PhiMin:\" option (check for typos)\n";
            printUsageAndExit();
        }
        optin >> phimin;

        /* All options read in successfully so generate the known beta */

        beta = new LAW_Vec(d+1);
        mu = new LAW_Vec(d+1);
        //beta = new double[d+1];
        //mu = new double[d+1];

        for (i=0; i<d + 1; i++) {
            beta->m_x[i] = betaBias + betaAmp * sin(2.0 * M_PI * i / (d + 1.0));
        }

        optin.close();

        /* If the prior type is informative, do the needed setup */
        if (pc == PriorInform) {
            stringComp = ExpNameString;
            stringComp += ".prior";
            std::ifstream optin(stringComp.c_str());
            if (!optin) {
                std::cerr << "Error: can't open prior file \"" << stringComp <<"\".\n";
                printUsageAndExit();
            }
            optin >> muAmp;
            optin >> muBias;
            optin >> q;
            for (i=0; i<d+1; i++) {
                optin >> r[i];
            }

            optin.close();
            for (i=0; i<d + 1; i++) {
                mu->m_x[i] = muBias + muAmp * sin(2.0 * M_PI * i / (d + 1.0));
            }
        }
 

        std::cout << "Successfully read in algver options file!\n";
        std::cout << "Option Settings (check for correctness!)\n";
        std::cout << "*****************************************\n";
        //std::cout << "Default setup read in experiment name: " << expname << " \n";
        std::cout << "Default setup read in beta bias: " << betaBias << " \n";
        std::cout << "Default setup read in beta amplitude: " << betaAmp << " \n";
        if (pc == PriorInform) {
            std::cout << "Default setup read in mu bias: " << muBias << " \n";
            std::cout << "Default setup read in mu amplitude: " << muAmp << " \n";
        }
        Display();
        std::cout << "\n\n";
        

    }

    else {
        std::cerr << "Error: unrecognized options type\n";
        printUsageAndExit();
    }
}


algverOptions::~algverOptions()
{
    delete beta;
    delete mu;
    delete r;
    delete phi;
    delete x_phi;
}

///////////////////////////////////////////////////////
// algver member functions
///////////////////////////////////////////////////////

/* Print MATLAB code for the prior */
void algver::PrintPrior(char *varname, std::ostream &os)
{
    int i;
    int p = m_opt->d+1;

    // Note - this is the prior excluding the factor of lambda
    if (m_opt->pc == PriorInform) {
        os << varname << "= eye(" << p << ")/" << m_opt->q 
           << "+ diag([";
        for (i=0; i<p; i++) {
            os << "1/" << m_opt->r[i] << ";";
        }
        os << "]);\n";
    }

}

void algver::PrintParams(char *betaname, char *lambdaname, char *phiname, std::ostream &os)
{
    m_opt->beta->Print(betaname, os);

    os << "\n" << lambdaname << "=" << m_opt->lambda << ";\n";
    // Note - this needs to be changed when Gaussian correlation is added
    os << phiname << "=" << m_opt->phi[0] << ";\n";
}

/* Generate error-free resposne to parameter input */
void algver::Response(LAW_Vec *beta)
{
    m_Xn->matvec(beta, m_zn);
}


void algver::Setup()
{
    int i;
    unsigned int n, d, p;
    CovMatrix *C;

    n = m_opt->n;
    d = m_opt->d;
    p = d + 1;


    RandString rs(n + n*d, m_opt->seed);

    // Scale n values for the measurement errors
    m_Rn->ScaleNormal(rs.m_v);

    // Create covariance matrix for Xn
    C = new CovMatrix(d, m_opt->x_ct, m_opt->x_lambda, m_opt->x_phi);

    // Scale d values n times for Xn
    if (d) {
        for (i=0; i<n; i++) {
            C->ScaleNormal(rs.m_v + n + i*d);

            // First column is all 1's
            m_Xn->m_a[i*p] = 1.0;

            // Other columns get the random data
            memcpy(m_Xn->m_a + i*p + 1, rs.m_v + n + i*d, sizeof(double)*d);
        }
    }
    else {
        for (i=0; i<n; i++) {
            // Only column is all 1's
            m_Xn->m_a[i] = 1.0;
        }
    }

    // yn = Xn*beta
    m_Xn->matvec(m_opt->beta, m_yn);

    // Add the measurement error to generate the data
    for (i=0; i<n; i++) {
        m_yn->m_x[i] += rs.m_v[i];
    }

    delete C;
    
}


// Constructors / Destructors
algver::algver(algverOptions *opt)
{
    m_opt = opt;
    m_yn = new LAW_Vec(m_opt->n);
    m_zn = new LAW_Vec(m_opt->n);
    m_Xn = new LAW_Mat(m_opt->n, m_opt->d + 1);
    m_Rn = new CovMatrix(m_opt->n, m_opt->ct, m_opt->lambda, m_opt->phi);
    Setup();
}


algver::algver(int argc, char **argv)
{
    if (argc != 3) {
        printUsageAndExit();
    }
    m_opt = new algverOptions(argv);

    m_yn = new LAW_Vec(m_opt->n);
    m_zn = new LAW_Vec(m_opt->n);
    m_Xn = new LAW_Mat(m_opt->n, m_opt->d + 1);
    m_Rn = new CovMatrix(m_opt->n, m_opt->ct, m_opt->lambda, m_opt->phi);
    Setup();
}


algver::~algver()
{
    delete m_yn;
    delete m_zn;
    delete m_Xn;
    delete m_Rn;
    delete m_opt;
}



///////////////////////////////////////////////////////
// algverPosterior member functions
///////////////////////////////////////////////////////


void algverPosterior::betahat(LAW_Vec *bhat)
{
    // Must have already computed Xntr_Rni_Xn so that
    // tmprectmat = inv(R) * Xn
    // and
    // tmpsquaremat = Xn^T * inv(R) * Xn
    // tmpvec = tmprectmat*yn = Xn^T * inv(R) * yn
    tmprectmat->matvec(m_av->m_yn, tmpvec, true);

    // Solve bhat = inv(tmprectmat) * tmpvec
    // which should be the same as
    // bhat = inv(Xn^T * inv(R) * Xn) * Xn^T * inv(R) * yn
    tmpsquaremat->cholsolve(tmpvec, bhat);

}

void algverPosterior::Xntr_Rni_Xn(LAW_Mat *res)
{


    // tmp = inv(R) * Xn
    m_av->m_Rn->LeftInvMultMat(m_av->m_Xn, tmprectmat);
    //tmp->Print("tmp", std::cout);

    //std::cout << "M " << tmp->m_M << " N "  << tmp->m_N << "\n";
    //std::cout << "M " << m_av->m_Xn->m_M << " N "  << m_av->m_Xn->m_N << "\n";
    //std::cout << "M " << res->m_M << " N "  << res->m_N << "\n";

    // res = Xn^T * tmp  = Xn^T * inv(R) * Xn
    m_av->m_Xn->lmult(tmprectmat, res, 1.0, 0.0, true, false);
    //res->Print("res", std::cout);

    // Save a copy in tmpsquaremat, also.
    tmpsquaremat->m_M = res->m_M;
    tmpsquaremat->m_N = res->m_N;
    memcpy(tmpsquaremat->m_a, res->m_a, sizeof(double)*(res->m_M)*(res->m_N));

}

// Beta prior info
void algverPosterior::sigma2(LAW_Mat *s2)
{
    int i;
    int M = m_av->m_opt->d+1;

    // This term is always there
    Xntr_Rni_Xn(s2);

    // If informative, need to add contribution from the prior

    // Prior is I/q + I/r. Adding these is (q+r)/qr*I. What gets
    // added as a term in the posterior covariance is the inverse
    // of this, hence the qr/(q+r).
    if (m_av->m_opt->pc == PriorInform) {
        for (i=0; i<M; i++) {
            s2->m_a[i*M + i] +=  m_av->m_opt->q * m_av->m_opt->r[i] 
                                / ( 
                                  (m_av->m_opt->q + m_av->m_opt->r[i]) 
                                  );

        }

    }

    // Invert to finish
    s2->cholinv();
}

void algverPosterior::mu2(LAW_Vec *m2, LAW_Mat *s2, LAW_Vec *bh)
{
    // This routine expects that
    // tmpsquaremat = Xn^T * inv(R) * Xn
    int i;
    int M = m_av->m_opt->d+1;

    // If non-informative prior, just return betahat
    if (m_av->m_opt->pc != PriorInform) {
        memcpy(m2->m_x, bh->m_x, sizeof(double)*M);
        return;
    }

    // Otherwise, do the computations to adjust for the Gaussian prior
    // info.

    // tmpvec = tmpsquaremat * bhat = Xn^T * inv(R) * Xn * bhat
    tmpsquaremat->matvec(bh, tmpvec, false);
    // Add the prior to the diagonal
    for (i=0; i<M; i++) {
        tmpvec->m_x[i*M + i] +=    m_av->m_opt->mu->m_x[i] 
                                   * m_av->m_opt->q 
                                   * m_av->m_opt->r[i] 
                                   / ( 
                                        (m_av->m_opt->q + m_av->m_opt->r[i]) 
                                     );

    }

    // tmpsquaremat = Xn^T * inv(R) * Xn
    // m2 = s2*tmpvec = s2 * (Xn^T * inv(R) * Xn * bhat + prior_stuff)
    s2->matvec(tmpvec, m2, false);
}

void algverPosterior::PrintBeta(char *varname, std::ostream os)
{
}

    // Lambda prior info
double algverPosterior::a1()
{
}

double algverPosterior::b1()
{
}

void algverPosterior::PrintLambda(char *varname, std::ostream os)
{
}

    // Phi prior info
void algverPosterior::PrintPhi(char *varname, std::ostream os)
{
}


void algverPosterior::MakeDensity()
{
    LAW_Mat *s2; 
    LAW_Vec *bh;
    LAW_Vec *mu;
    std::fstream fs;

    s2 = new LAW_Mat(m_av->m_opt->d + 1, m_av->m_opt->d + 1);
    sigma2(s2);

    bh = new LAW_Vec(m_av->m_opt->d + 1);
    betahat(bh);

    mu = new LAW_Vec(m_av->m_opt->d + 1);
    mu2(mu, s2, bh);

    fs.open("test.m", std::fstream::out);
    fs.precision(9);

    m_av->PrintPrior("sigs", fs);
    m_av->m_opt->mu->Print("mu0", fs);
    m_av->m_yn->Print("yn", fs);
    m_av->m_Xn->Print("xn", fs);
    m_av->m_Rn->Print("rn", fs);

    bh->Print("bhat", fs);
    s2->Print("s2", fs);
    mu->Print("mu2", fs);

    fs.close();

}

void algverPosterior::MakeMatlab()
{
    int i;
    double x, y;
    std::fstream fs;
    std::string outfile(m_av->m_opt->expname);

    outfile = outfile + "/posterior.m";

    // For now, we'll output the data and compute the priors
    // in MATLAB. It'll be better to eventually change this to compute 
    // the prior parameters here, since the data can get large for
    // big sample sizes


    //std::cout << outfile + "\n";
    fs.open((const char *)outfile.c_str(), std::fstream::out);
    fs.precision(9);


    // Here are the variables for the linear model

    if (m_av->m_opt->pc == PriorInform) {
        fs << "priorcase='Informative';\n";
    }
    else if (m_av->m_opt->pc == PriorNonInform){
        fs << "priorcase='Non-Informative';\n";
    }

    if (m_av->m_opt->kp == KnownLambdaPhi) {
        fs << "calvars='Beta';\n";
    }
    else if (m_av->m_opt->kp == KnownPhi) {
        fs << "calvars='Beta,Lambda';\n";
    }
    else if (m_av->m_opt->kp == KnownNone) {
        fs << "calvars='Beta,Lambda,Phi';\n";
    }
    
    if (m_av->m_Rn->m_ct == CorrNone) {
        fs << "corrtype='None';\n";
    }
    if (m_av->m_Rn->m_ct == CorrEqual) {
        fs << "corrtype='Equal';\n";
    }
    if (m_av->m_Rn->m_ct == CorrAR1) {
        fs << "corrtype='AR(1)';\n";
    }

    m_av->m_yn->Print("yn", fs);
    m_av->m_Xn->Print("xn", fs);
    m_av->PrintParams("beta", "l", "phi", fs);
    m_av->m_Rn->Print("rn", fs);
    m_av->m_Rn->InvPrint("rni", fs);
    fs << "xrx=xn'*rni*xn;\n";
    fs << "xrxi=inv(xrx);\n";
    fs << "xr=xn'*rni;\n";
    fs << "betahat=xrxi*xr*yn;\n";
    fs << "N=" << m_av->m_opt->n << ";\n";
    fs << "Nb=" << m_av->m_opt->d+1 << ";\n";

    // This stuff is used in all cases, more-or-less

    if (m_av->m_opt->pc == PriorNonInform) {
        fs << "mu2=betahat; sig2= xrxi;\n";
    }
    else {
        m_av->PrintPrior("sigs", fs);
    
        if (m_av->m_opt->pc == PriorInform) {
            m_av->m_opt->mu->Print("mu0", fs);
        }
        fs << "sig2=inv(inv(sigs) + xrx);\n";
        fs << "mu2=sig2*(xr*yn + inv(sigs)*mu0);\n";
    }

    // THIS is the posterior covariance, remember. Also remember that
    // MATLAB asks for standard deviation not variance for the scalar
    // normal function.
    fs << "betacov=sig2/l;\n";

    fs << "seed=" << m_av->m_opt->seed << ";\n";

    fs << "bmin=" << m_av->m_opt->betamin <<";\n";
    fs << "bmax=" << m_av->m_opt->betamax <<";\n";
    fs << "lmin=" << m_av->m_opt->lambdamin <<";\n";
    fs << "lmax=" << m_av->m_opt->lambdamax <<";\n";
    fs << "fmin=" << m_av->m_opt->phimin <<";\n";
    fs << "fmax=" << m_av->m_opt->phimax <<";\n";

    fs << "ty=linspace(0,1,N);\n";
    fs << "tbeta=(0:Nb-1)/Nb;\n";

    if (m_av->m_opt->kp == KnownPhi) {

            // Lambda Parameters 

            // The prior for lambda in the cases considered is
            // always non-informative. If an informative gamma
            // distribution is used, the parameters can be added. 

            fs << "resid=yn-xn*betahat;\n";
            if (m_av->m_opt->pc == PriorNonInform) {
                fs << "a1=(N-Nb)/2;\n";
                fs << "b1=(resid'*rni*resid)/2;\n";
            }
            else {
                fs << "diffprior=betahat-mu0;\n";
                fs << "sig3=sigs+xrxi;\n";
                fs << "a1=N/2;\n";
                fs << "b1=(resid'*rni*resid + diffprior'*inv(sig3)*diffprior)/2;\n";
            }

            // Use this with Matlab since it seems to use a different convention
            // for Gamma distribution parameters than the specification.
            fs << "ml_a1 = a1; ml_b1 = 1/b1;\n";

            // Beta parameters (multivariate t-distribution)
            fs << "dof=2*a1;\n";
            fs << "loc=mu2;\n";
            fs << "scl=b1*sig2/a1;\n";
    }

    if (m_av->m_opt->kp == KnownNone) {
        if (m_av->m_Rn->m_ct == CorrEqual) {
            fs << "x = @(f) f / (1.0 + (N - 1.0) * (f));\n";
            fs << "y = @(f) 1.0 / (1.0 - f);\n";
            fs << "rni = @(f) y(f) *(eye(N) -  x(f)*ones(N));\n";
            fs << "detrn=@(f) (1-f)^(N-1)*(1+(N-1)*f);\n";
        }
        else if (m_av->m_Rn->m_ct == CorrAR1) {
            fs << "offdvec = @(f) -f*ones(N-1,1);\n";
            fs << "diagvec = @(f) ones(N,1) + [0; ones(N-2,1)*f^2; 0];\n";
            fs << "rni = @(f) (diag(diagvec(f)) + diag(offdvec(f),-1) + diag(offdvec(f),1)) / (1-f^2);\n";
            fs << "detrn=@(f) (1-f^2)^(N-1);\n";
        }

        if (m_av->m_opt->pc == PriorNonInform) {
            fs << "xrx=@(f)xn'*rni(f)*xn;\n";
            fs << "detxrx=@(f)det(xrx(f));\n";
            fs << "xrxi=@(f)inv(xrx(f));\n";
            fs << "xr=@(f)xn'*rni(f);\n";
            fs << "betahat=@(f)xrxi(f)*xr(f)*yn;\n";
            fs << "resid=@(f)yn-xn*betahat(f);\n";
            fs << "a1=(N-Nb)/2;\n";
            fs << "b1=@(f)(resid(f)'*rni(f)*resid(f))/2;\n";
            fs << "ml_a1 = a1; ml_b1 = @(f)1/b1(f);\n";
            fs << "mu2=@(f)betahat(f); sig2= @(f)xrxi(f);\n";
            fs << "dof=2*a1;\n";
            fs << "loc=@(f)mu2(f);\n";
            fs << "scl=@(f)b1(f)*sig2(f)/a1;\n";
            // NOTE: Need to ask Brian about this - does the non-informative 
            // really just remove det(sig3)?
            //fs << "piphi=@(f) 1 / (b1(f)^a1 * sqrt(detrn(f)));\n";
            fs << "piphi=@(f) 1 / (b1(f)^a1 * sqrt(detrn(f)) * sqrt(detxrx(f)));\n";
            // TODO: NEED TO CHANGE THIS SO IT WORKS FOR MULTIVARIABLE T-DISTRIBUTIONS!
            fs << "pibetaphi=@(b,f) (1/sqrt(scl(f)))*tpdf((b - loc(f))/sqrt(scl(f)),dof);\n";
            fs << "pilambdaphi=@(l,f) gampdf(l, ml_a1, ml_b1(f));\n";
            fs << "numq = 100; piphin = zeros(numq,1); fn = linspace(fmin,fmax,numq);\n";
            fs << "for i = 1:numq; piphin(i) = piphi(fn(i)); end;\n";
            fs << "c = trapz(fn, piphin); piphin = piphin / c;\n";
            fs << "if exist('lqmax') & exist('lqmin'); ln = linspace(lqmin, lqmax, numq);\n"; 
            fs << "else; ln = linspace(lmin, lmax, numq); end;\n"; 
            fs << "if exist('bqmax') & exist('bqmin'); bn = linspace(bqmin, bqmax, numq);\n";
            fs << "else; bn = linspace(bmin, bmax, numq); end;\n";
            fs << "tmp1 = zeros(numq,1); tmp2 = tmp1; pibetan = tmp1; pilambdan = tmp1;\n";
            //fs << "%for i=1:numq;\n"; 
            //fs << "%for j=1:numq; \n";
            //fs << "%tmp1(j) = pibetaphi(bn(i), fn(j)); tmp2(j) = pilambdaphi(ln(i), fn(j));\n";
            //fs << "%end; pibetan(i) = trapz(fn, tmp1.*piphin); pilambdan(i) = trapz(fn, tmp2.*piphin); end;\n"; 
            fs << "for i=1:numq;\n"; 
            fs << "for j=1:numq; \n";
            fs << "tmp1(j) = pibetaphi(bn(i), fn(j)); tmp2(j) = pilambdaphi(ln(i), fn(j));\n";
            fs << "end; pibetan(i) = trapz(fn, tmp1.*piphin); pilambdan(i) = trapz(fn, tmp2.*piphin); end;\n"; 
        }
        else {
            fs << "xrx=@(f)xn'*rni(f)*xn;\n";
            fs << "detxrx=@(f)det(xrx(f));\n";
            fs << "xrxi=@(f)inv(xrx(f));\n";
            fs << "sig3=@(f)sigs + xrxi(f);\n";
            fs << "detsig3=@(f) det(sig3(f));\n";
            fs << "xr=@(f)xn'*rni(f);\n";
            fs << "betahat=@(f)xrxi(f)*xr(f)*yn;\n";
            fs << "resid=@(f)yn-xn*betahat(f);\n";
            fs << "diffprior=@(f)betahat(f)-mu0;\n";
            fs << "a1=(N-Nb)/2;\n";
            //fs << "b1=@(f)(resid(f)'*rni(f)*resid(f))/2;\n";
            fs << "b1=@(f)(resid(f)'*rni(f)*resid(f) + diffprior(f)'*inv(sig3(f))*diffprior(f))/2;\n";
            fs << "ml_a1 = a1; ml_b1 = @(f)1/b1(f);\n";
            //fs << "mu2=@(f)betahat(f); sig2= @(f) inv(xrx(f) + inv(sigs));\n";
            fs << "sig2= @(f) inv(xrx(f) + inv(sigs)); mu2=@(f) sig2(f)*(xrx(f)*betahat(f) +inv(sigs)*mu0);\n";
            fs << "dof=2*a1;\n";
            fs << "loc=@(f)mu2(f);\n";
            fs << "scl=@(f)b1(f)*sig2(f)/a1;\n";
            // NOTE: Need to ask Brian about this - does the non-informative 
            // really just remove det(sig3)?
            //fs << "piphi=@(f) 1 / (b1(f)^a1 * sqrt(detrn(f)));\n";
            fs << "piphi=@(f) 1 / (b1(f)^a1 * sqrt(detrn(f)) * sqrt(detxrx(f)) * sqrt(detsig3(f)));\n";
            // TODO: NEED TO CHANGE THIS SO IT WORKS FOR MULTIVARIABLE T-DISTRIBUTIONS!
            fs << "pibetaphi=@(b,f) (1/sqrt(scl(f)))*tpdf((b - loc(f))/sqrt(scl(f)),dof);\n";
            fs << "pilambdaphi=@(l,f) gampdf(l, ml_a1, ml_b1(f));\n";
            fs << "numq = 100; piphin = zeros(numq,1); fn = linspace(fmin,fmax,numq);\n";
            fs << "for i = 1:numq; piphin(i) = piphi(fn(i)); end;\n";
            fs << "c = trapz(fn, piphin); piphin = piphin / c;\n";
            fs << "if exist('lqmax') & exist('lqmin'); ln = linspace(lqmin, lqmax, numq);\n"; 
            fs << "else; ln = linspace(lmin, lmax, numq); end;\n"; 
            fs << "if exist('bqmax') & exist('bqmin'); bn = linspace(bqmin, bqmax, numq);\n";
            fs << "else; bn = linspace(bmin, bmax, numq); end;\n";
            fs << "tmp1 = zeros(numq,1); tmp2 = tmp1; pibetan = tmp1; pilambdan = tmp1;\n";
            fs << "%for i=1:numq;\n"; 
            fs << "%for j=1:numq; \n";
            fs << "%tmp1(j) = pibetaphi(bn(i), fn(j)); tmp2(j) = pilambdaphi(ln(i), fn(j));\n";
            fs << "%end; pibetan(i) = trapz(fn, tmp1.*piphin); pilambdan(i) = trapz(fn, tmp2.*piphin); end;\n"; 
            fs << "for i=1:numq;\n"; 
            fs << "for j=1:numq; \n";
            fs << "tmp1(j) = pibetaphi(bn(i), fn(j)); tmp2(j) = pilambdaphi(ln(i), fn(j));\n";
            fs << "end; pibetan(i) = trapz(fn, tmp1.*piphin); pilambdan(i) = trapz(fn, tmp2.*piphin); end;\n"; 
        }
    }



    // Add in the variables needed 

    fs.close();

}


// Constructors / Destructors
    
algverPosterior::algverPosterior(algver *av)
{
    m_av = av;
    tmprectmat = new LAW_Mat(m_av->m_opt->n, m_av->m_opt->d+1);
    tmpsquaremat = new LAW_Mat(m_av->m_opt->d+1, m_av->m_opt->d+1);
    tmpvec = new LAW_Vec(m_av->m_opt->d+1);
}

algverPosterior:: ~algverPosterior()
{
    delete tmprectmat;
    delete tmpsquaremat;
    delete tmpvec;
}

