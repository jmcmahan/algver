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

//#include <queso/JeffreysJointPdf.h>
#include "JeffreysJointPdf.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
JeffreysJointPdf<V,M>::JeffreysJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  algver *av)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"uni").c_str(),
                            domainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering JeffreysJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving JeffreysJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
  m_av = av;
}
// Destructor --------------------------------------
template<class V,class M>
JeffreysJointPdf<V,M>::~JeffreysJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
JeffreysJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "JeffreysJointPdf<V,M>::actualValue()",
                      "invalid input");

  //const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_domainSet);
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&(this->m_domainSet.vectorSpace()));

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainDirection) {}; // just to remove compiler warning

  double a,b;

  if (m_normalizationStyle != 0) {
    return 1.0 / (domainVector[0]);
  }
  else {
    // This is kind of lame, but we're in a hurry
    a = m_av->m_opt->lambdamin;
    b = m_av->m_opt->lambdamax;
    //std::cout << "a: " << a << " b: " << b << "\n";
    return 1.0 / (domainVector[0]*std::log(b/a));
    //return 1.0 / (domainVector[0]);
  }

  //return 1./volume; // No need to multiply by exp(m_logOfNormalizationFactor) [PDF-12]
}
//--------------------------------------------------
template<class V, class M>
double
JeffreysJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainVector[0]) {}; // just to remove compiler warning
  if (domainDirection) {}; // just to remove compiler warning


  return std::log(this->actualValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
  //return log(volume); // No need to add m_logOfNormalizationFactor [PDF-12]
}
//--------------------------------------------------
template<class V, class M>
double
JeffreysJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering JeffreysJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving JeffreysJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::JeffreysJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
