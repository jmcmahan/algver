#ifndef _RNGWRAPPER_HPP_
#define _RNGWRAPPER_HPP_


#ifndef ALGVER_USE_STDC_RNG
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#else
#include <random>
#include <chrono>
#endif

#ifndef ALGVER_USE_STDC_RNG

class RandString{
    public:
    unsigned int m_n;                           // Size of the string
    int m_seed;                                 // Random seed
    double *m_v;                                // Values

    gsl_rng *m_rng;         // Random number generator structure



    // Note, if there's no need to add a member function to get
    // new random values (i.e., as-is, the constructor makes
    // the random data and nothing else does) then there's
    // no need to have the implementation-specific member 
    // m_rng here or m_gen, m_dist in the C++11 version. 

    // Constructors / Destructors
    RandString(unsigned int, int);
    ~RandString();
};

#else 

/* This commented out code implements the random number generation 
   without GSL using standard C++ as of C++11. I'm commenting
   out and replacing with a GSL-based routine since the build of gcc
   I'm using at the office doesn't support this. */

#include <random>
#include <chrono>

class RandString {
    public:
    unsigned int m_n;                           // Size of the string
    int m_seed;                                 // Random seed
    double *m_v;                                // Values
    std::default_random_engine m_gen;           // Random generator
    std::normal_distribution<double> m_dist;    // Probability distribution

    // Constructors / Destructors
    RandString(unsigned int, int);
    ~RandString();
};
#endif



#endif /* _RNGWRAPPER_HPP_ */
