// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkCurvatureBasedBoneResorptionImageFilter_hxx
#define itkCurvatureBasedBoneResorptionImageFilter_hxx

#include "itkCurvatureBasedBoneResorptionImageFilter.h"
#include "itkImageRegionIterator.h"

namespace itk
{

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::CurvatureBasedBoneResorptionImageFilter()
{
  // Set difference function pointer
  typename CurvatureBasedBoneResorptionFunctionType::RadiusType r;
  r.Fill(1);
  m_ResorptionFunction = CurvatureBasedBoneResorptionFunctionType::New();
  m_ResorptionFunction->Initialize(r);
  this->SetDifferenceFunction( static_cast< FiniteDifferenceFunctionType * >(
                                 m_ResorptionFunction.GetPointer() ) );

  // Setup iterations
  m_TotalTime = NumericTraits< RealType >::ZeroValue();
  m_ElapsedTime = NumericTraits< RealType >::ZeroValue();

  // We want to use image spacing because we're dealing with physical units
  this->UseImageSpacingOn();
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
bool
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::Halt()
{
  /* These are all cherry picked from FiniteDifferenceImageFilter. Note that we
  stop based on time, not iterations. These is an edge case where we somehow get
  zero update (del u = 0) and we have a potential infinite loop */
  if ( m_ElapsedTime != 0 )
  {
    this->UpdateProgress( static_cast< float >( this->m_ElapsedTime )
                          / static_cast< float >( this->m_TotalTime ) );
  }
  
  if ( this->m_ElapsedTime == NumericTraits< RealType >::ZeroValue() )
  {
    return false;
  }
  else if ( Math::FloatAlmostEqual(m_ElapsedTime, m_TotalTime) || (m_ElapsedTime > m_TotalTime) )
  {
    return true;
  }
  else
  {
    return false;
  }
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >::TimeStepType
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::CalculateChange()
{
  TimeStepType dt = Superclass::CalculateChange();

  /* Now we can update. Note that ResolveTimeStep is const */
  if (dt > 0) {
    m_ElapsedTime = m_ElapsedTime + static_cast< RealType >(dt);
  } else {
    itkWarningMacro(<< " time update is zero. Are you sure you set the inputs correctly? Exiting early.")
    m_ElapsedTime = m_TotalTime;
  }

  return dt;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >::TimeStepType
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::ResolveTimeStep(const std::vector< TimeStepType >& timeStepList,
                                       const std::vector< bool >& valid ) const
{
  TimeStepType dt = Superclass::ResolveTimeStep(timeStepList, valid);

  /* If dt will cause us to over shoot our total time, shrink it to the delta.
   * This is OK to do since dt is a lower bound in the CFL condition */
  dt = std::min(dt, m_TotalTime - m_ElapsedTime);

  return dt;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::Initialize()
{
  Superclass::Initialize();

  // Reset elapsed time
  m_ElapsedTime = NumericTraits< RealType >::ZeroValue();

  // Pass image
  m_ResorptionFunction->SetMaskImage(this->GetMaskImage());
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::InitializeIteration()
{
  Superclass::InitializeIteration();

  /* Only need to perform reinitialization if we have curvature */
  /* Perform reinitialization */
  m_ReinitializationFilter = ReinitializationFilterType::New();
  m_ReinitializationFilter->SetInput(this->GetOutput());
  m_ReinitializationFilter->Update();

  this->GetOutput()->Graft(m_ReinitializationFilter->GetOutput());
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
CurvatureBasedBoneResorptionImageFilter< TInputImage, TMaskImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "TotalTime: " << m_TotalTime << std::endl;
  os << indent << "ElapsedTime: " << m_ElapsedTime << std::endl;
  os << indent << "ReinitializationFilter: " << m_ReinitializationFilter << std::endl;
}

} // end namespace itk

#endif /* itkCurvatureBasedBoneResorptionImageFilter_hxx */