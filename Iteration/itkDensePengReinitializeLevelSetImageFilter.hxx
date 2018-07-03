// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkDensePengReinitializeLevelSetImageFilter_hxx
#define itkDensePengReinitializeLevelSetImageFilter_hxx

#include "itkDensePengReinitializeLevelSetImageFilter.h"

namespace itk {

template< typename TInputImage, typename TOutputImage >
DensePengReinitializeLevelSetImageFilter< TInputImage, TOutputImage >
::DensePengReinitializeLevelSetImageFilter()
{
  // Set difference function pointer
  typename DensePengReinitializeLevelSetFunctionType::RadiusType r;
  r.Fill(1);
  m_ReshapeFunction = DensePengReinitializeLevelSetFunctionType::New();
  m_ReshapeFunction->Initialize(r);
  this->SetDifferenceFunction( static_cast< FiniteDifferenceFunctionType * >(
                                 m_ReshapeFunction.GetPointer() ) );

  // We want to use image spacing because we're dealing with physical units
  this->UseImageSpacingOn();

  // Set reasonable defaults
  this->SetMaximumRMSError(1e-4);
  this->SetNumberOfIterations(1000);
}


} // end namspace itk

#endif /* itkDensePengReinitializeLevelSetImageFilter_hxx */