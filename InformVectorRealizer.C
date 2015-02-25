// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

//#include <queso/InformVectorRealizer.h>
#include "InformVectorRealizer.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <cmath>

namespace QUESO {
    

// Constructor -------------------------------------
template<class V, class M>
InformVectorRealizer<V,M>::InformVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet,
  algver *av)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering InformVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving InformVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
  m_av = av;
}

// Destructor --------------------------------------
template<class V, class M>
InformVectorRealizer<V,M>::~InformVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
InformVectorRealizer<V,M>::realization(V& nextValues) const
{
  int i;
  double a, b, l;
  double q;
  double *r;
  double *m;

  unsigned int p;
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_unifiedImageSet);
  V iidGaussianVector(m_unifiedImageSet.vectorSpace().zeroVector());

  UQ_FATAL_TEST_MACRO(imageBox == NULL,
                      m_env.worldRank(),
                      "InformVectorRealizer<V,M>::realization()",
                      "only box images are supported right now");
  
  /* This check is commented out, as it's only valid for the lambda
     parameter, not the Gaussian parameters. May want to add something
     back in. */
  /*
  double smallerOfMinValues = imageBox->minValues().getMinValue();
  UQ_FATAL_TEST_MACRO(smallerOfMinValues <= 0.0,
                      m_env.worldRank(), 
                      "InformVectorRealizer<V,M>::realization()\n",
                      "invalid input: Inform distribution is only defined in (0, infinity) and min(m_minValues) < 0. ");

  */


  // Want to do exp(  (log(max) - log(min)) * z + log(min) ) for all components
  // where max/min are the max/min of the parameter domain and z is
  // a uniformly distributed variable on [0,1].
  
  p = nextValues.sizeLocal();



  // Note - the particular density being defined assumes a scalar. The vector stuff
  // is there just because.

  m = m_av->m_opt->mu->m_x;
  q = m_av->m_opt->q;
  r = m_av->m_opt->r;
  a = imageBox->minValues()[p-1];
  b = imageBox->maxValues()[p-1];
  

  bool outOfSupport = true;
  do {
    iidGaussianVector.cwSetGaussian(0.0, 1.0);
    l = std::exp(  (std::log(b) - std::log(a))  * 
                            m_env.rngObject()->uniformSample() + log(a));

    for (i=0; i<p-1; i++) {
        // NEED TO ADD IN THE VARIANCE FROM THE PRIOR HERE, TOO!!!
        // SHOULD JUST BE DIVIDED BY SQUARE ROOT OF THE DIAGONAL ENTRIES,
        // HOWEVER THAT IS STORED. NOTE THIS ISN'T IMPORTANT FOR
        // THE INVERSE PROBLEM SOLVER SINCE IT JUST USES THE LIKELIHOOD
        nextValues[i] = m[i] + iidGaussianVector[i] / 
                                std::sqrt(l * q * r[i] / (q + r[i] ));
    }
    nextValues[p-1] = l;

    outOfSupport = !(this->m_unifiedImageSet.contains(nextValues));
  } while (outOfSupport); // prudenci 2011-Oct-04

  return;
}


}  // End namespace QUESO


template class QUESO::InformVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
