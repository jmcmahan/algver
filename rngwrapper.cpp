#include "rngwrapper.hpp"



#ifndef ALGVER_USE_STDC_RNG
/* Create a random string of size n sampled from a standard normal 
   distribution. The seed value initializes the random number
   generator. If seed == 0, the timer is used to pick the seed,
   otherwise it is the seed value used. */
RandString::RandString(unsigned int n, int seed)
{
    int i;
    struct timeval tv;
    
    m_n = n;
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);

    // If seed == 0, the seed is selected according to the timer (i.e., 
    // pseudo-randomly), otherwise the value is used as the seed.
    if (seed == 0) {
        gettimeofday(&tv, 0);
        m_seed = tv.tv_sec + tv.tv_usec;
    }
    else {
        m_seed = seed;
    }
    m_v = new double[m_n];

    gsl_rng_set(m_rng, m_seed);
    
    for (i=0; i<m_n; i++) {
        m_v[i] = gsl_ran_gaussian(m_rng, 1.0);
    }


}


RandString::~RandString()
{
    gsl_rng_free(m_rng);
    delete m_v;
    
}



#else
/* C++11 based code */

/* Create a random string of size n sampled from a standard normal 
   distribution. The seed value initializes the random number
   generator. If seed == 0, the timer is used to pick the seed,
   otherwise it is the seed value used. */
RandString::RandString(unsigned int n, int seed)
{
    int i;
    
    m_n = n;

    // If seed == 0, the seed is selected according to the timer (i.e., 
    // pseudo-randomly), otherwise the value is used as the seed.
    if (seed == 0) {
        m_seed = std::chrono::system_clock::now().time_since_epoch().count();     
    }
    else {
        m_seed = seed;
    }
    m_v = new double[m_n];

    std::default_random_engine m_gen(seed);
    std::normal_distribution<double> m_dist(0.0, 1.0);
    
    for (i=0; i<m_n; i++) {
        m_v[i] = m_dist(m_gen);
    }
    
}
#endif

