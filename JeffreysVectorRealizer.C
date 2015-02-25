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

//#include <queso/JeffreysVectorRealizer.h>
#include "JeffreysVectorRealizer.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <cmath>

namespace QUESO {
    

// Constructor -------------------------------------
template<class V, class M>
JeffreysVectorRealizer<V,M>::JeffreysVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering JeffreysVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving JeffreysVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

// Destructor --------------------------------------
template<class V, class M>
JeffreysVectorRealizer<V,M>::~JeffreysVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
JeffreysVectorRealizer<V,M>::realization(V& nextValues) const
{
  int i;
  double a, b;

  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_unifiedImageSet);

  UQ_FATAL_TEST_MACRO(imageBox == NULL,
                      m_env.worldRank(),
                      "JeffreysVectorRealizer<V,M>::realization()",
                      "only box images are supported right now");
  
  double smallerOfMinValues = imageBox->minValues().getMinValue();
  UQ_FATAL_TEST_MACRO(smallerOfMinValues <= 0.0,
                      m_env.worldRank(), 
                      "JeffreysVectorRealizer<V,M>::realization()\n",
                      "invalid input: Jeffreys distribution is only defined in (0, infinity) and min(m_minValues) < 0. ");


  // Want to do exp(  (log(max) - log(min)) * z + log(min) ) for all components
  // where max/min are the max/min of the parameter domain and z is
  // a uniformly distributed variable on [0,1].
  

  // For now, I'm just going to put a hack in here. For the final implemntation
  // it might make sense to put the actual implementation in
  // GslVector.C, as that's what all the other random variable classes
  // seem to do.

  //nextValues.cwSetUniform(imageBox->minValues(),imageBox->maxValues());


  // Note - the particular density being defined assumes a scalar. The vector stuff
  // is there just because.
  for (i=0; i<nextValues.sizeLocal(); i++) {
      a = imageBox->minValues()[i];
      b = imageBox->maxValues()[i];
      nextValues[i] = std::exp(  (std::log(b) - std::log(a))  * 
                            m_env.rngObject()->uniformSample() + log(a));
  }
  return;
}


}  // End namespace QUESO


template class QUESO::JeffreysVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
