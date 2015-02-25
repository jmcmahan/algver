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

//#include <queso/InformVectorRV.h>
#include "InformVectorRV.h"
#include "InformJointPdf.h"
#include "InformVectorRealizer.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor-------------------------------
template<class V, class M>
InformVectorRV<V,M>::InformVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet,
  algver *av)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"jef").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering InformVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf        = new InformJointPdf<V,M>(m_prefix.c_str(),
                                                 m_imageSet,
                                                 av);
  m_realizer   = new InformVectorRealizer<V,M>(m_prefix.c_str(),
                                                 m_imageSet,
                                                 av);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving InformVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
InformVectorRV<V,M>::~InformVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods --------------------------------------
template <class V, class M>
void
InformVectorRV<V,M>::print(std::ostream& os) const
{
  os << "InformVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::InformVectorRV<QUESO::GslVector,QUESO::GslMatrix>;

