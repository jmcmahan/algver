//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
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

//#include <queso/InformJointPdf.h>
#include "InformJointPdf.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
InformJointPdf<V,M>::InformJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  algver *av)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"uni").c_str(),
                            domainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering InformJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving InformJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
  m_av = av;
}
// Destructor --------------------------------------
template<class V,class M>
InformJointPdf<V,M>::~InformJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
InformJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "InformJointPdf<V,M>::actualValue()",
                      "invalid input");

  double returnValue;
  //const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_domainSet);
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&(this->m_domainSet.vectorSpace()));

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainDirection) {}; // just to remove compiler warning

 
  if (this->m_domainSet.contains(domainVector) == false) { // prudenci 2011-Oct-04
    returnValue = 0.;
  }
  else {
    returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
  }
  return returnValue;
  
}
//--------------------------------------------------
template<class V, class M>
double
InformJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  unsigned int n;
  int i;
  double returnValue = 0.0;
  double diffValue = 0.0;
  double q;
  double *r;
  double l;
  double d;
  double a, b;


  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainVector[0]) {}; // just to remove compiler warning
  if (domainDirection) {}; // just to remove compiler warning

  n = this->m_domainSet.vectorSpace().dimLocal();
  q = m_av->m_opt->q;
  r = m_av->m_opt->r;
  l = domainVector[n-1];
  d = 1.0;

  a = m_av->m_opt->lambdamin;
  b = m_av->m_opt->lambdamax;
  // Interior of the likelihood for the Gaussian part and computing the determinant
  // of the prior covariance in d
  for (i=0; i<n-1; i++) {
    diffValue = domainVector[i] - m_av->m_opt->mu->m_x[i];
    returnValue += -(l * 0.5 * diffValue*diffValue * q * r[i] / (q + r[i]));
    d *= (q + r[i]) / (q*r[i]);
  }
  // Add the normalization term for the likelihood 
  returnValue += 0.5 * n * std::log(l / (2.0*M_PI)) - 0.5 * log(d);
  
  // Contribution due to the 1 / lambda prior
  returnValue += - std::log(l) - std::log(std::log(b/a));

  return returnValue;
  //return log(volume); // No need to add m_logOfNormalizationFactor [PDF-12]
}
//--------------------------------------------------
template<class V, class M>
double
InformJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering InformJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving InformJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::InformJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
