// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkDensePengReinitializeLevelSetFunction_hxx
#define itkDensePengReinitializeLevelSetFunction_hxx

#include "itkDensePengReinitializeLevelSetFunction.h"
#include "itkNumericTraits.h"

namespace itk {

template< typename TImageType >
DensePengReinitializeLevelSetFunction< TImageType >
::DensePengReinitializeLevelSetFunction() :
  m_Epsilon(1.0e-9),
  m_CFL(0.5f)
{}

template< typename TImageType >
void
DensePengReinitializeLevelSetFunction< TImageType >
::Initialize(const RadiusType & r)
{
  this->SetRadius(r);

  // Dummy neighborhood.
  NeighborhoodType it;
  it.SetRadius(r);

  // Find the center index of the neighborhood.
  m_Center =  it.Size() / 2;

  // Get the stride length for each axis.
  for ( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_xStride[i] = it.GetStride(i);
  }
}

template< typename TImageType >
typename DensePengReinitializeLevelSetFunction< TImageType >::PixelType
DensePengReinitializeLevelSetFunction< TImageType >
::ComputeUpdate(const NeighborhoodType & it, void *globalData,
                const FloatOffsetType & offset)
{
  // Heavily based on itkLevelSetFunction.hxx
  unsigned int          i;
  const ScalarValueType ZERO = NumericTraits< ScalarValueType >::ZeroValue();
  const ScalarValueType center_value  = it.GetCenterPixel();
  ScalarValueType S;

  const NeighborhoodScalesType neighborhoodScales = this->ComputeNeighborhoodScales();

  // Global data structure
  auto * gd = (GlobalDataStruct *)globalData;

  // Compute the first derivatives
  gd->m_GradMagSqr = m_Epsilon;
  for ( i = 0; i < ImageDimension; i++ )
  {
    const auto positionA = static_cast< unsigned int >( m_Center + m_xStride[i] );
    const auto positionB = static_cast< unsigned int >( m_Center - m_xStride[i] );

    gd->m_dx[i] = 0.5 * ( it.GetPixel(positionA)
                          - it.GetPixel(positionB) ) * neighborhoodScales[i];

    gd->m_dx_forward[i]  = ( it.GetPixel(positionA) - center_value ) * neighborhoodScales[i];

    gd->m_dx_backward[i] = ( center_value - it.GetPixel(positionB) ) * neighborhoodScales[i];

    gd->m_GradMagSqr += gd->m_dx[i] * gd->m_dx[i];
  }

  // Compute minimum spacing
  RealType min_spacing = NumericTraits< RealType >::max();
  for ( i = 0; i < ImageDimension; i++ )
  {
    min_spacing = std::min(1.0/(this->m_ScaleCoefficients[i]), min_spacing);
  }

  // Compute S(phi)
  S = center_value / std::sqrt(Math::sqr(center_value) + gd->m_GradMagSqr * Math::sqr(min_spacing));

  // Compute propagation.
  ScalarValueType propagation_gradient = ZERO;
  if ( S > ZERO ) {
    for ( i = 0; i < ImageDimension; i++ )
    {
      propagation_gradient += itk::Math::sqr( std::max(gd->m_dx_backward[i], ZERO) )
                              + itk::Math::sqr( std::min(gd->m_dx_forward[i],  ZERO) );
    }
  } else {
    for ( i = 0; i < ImageDimension; i++ )
    {
      propagation_gradient += itk::Math::sqr( std::min(gd->m_dx_backward[i], ZERO) )
                              + itk::Math::sqr( std::max(gd->m_dx_forward[i],  ZERO) );
    }
  }

  // Store maximum of absolute value
  gd->m_MaxS = std::max(gd->m_MaxS, std::abs(S));

  return (PixelType)( S*(1.0-propagation_gradient) );
}

template< typename TImageType >
typename DensePengReinitializeLevelSetFunction< TImageType >::TimeStepType
DensePengReinitializeLevelSetFunction< TImageType >
::ComputeGlobalTimeStep(void * GlobalData) const
{
  // For CFL condition, see equation 39 of Peng, Danping, et al. "A PDE-based fast local
  // level set method." Journal of computational physics 155.2 (1999): 410-438.
  TimeStepType dt;
  unsigned int i;

  auto * d = (GlobalDataStruct *)GlobalData;

  // Compute the min spacing
  RealType min_spacing = NumericTraits< RealType >::max();
  for ( i = 0; i < ImageDimension; i++ )
  {
    min_spacing = std::min(1.0/(this->m_ScaleCoefficients[i]), min_spacing);
  }

  // dt (max S / min dx) <= 1/2
  RealType D = 2.0f * std::abs(d->m_MaxS) / min_spacing;

  if ( Math::NotAlmostEquals(D, NumericTraits< RealType >::ZeroValue()) ) {
    dt = (TimeStepType) ( m_CFL / D);
  } else {
    dt = (TimeStepType)0.0f;
  }

  // Reset
  d->m_MaxS = NumericTraits< ScalarValueType >::ZeroValue();

  return dt;
}

template< typename TImageType >
void
DensePengReinitializeLevelSetFunction< TImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Epsilon: " << m_Epsilon << std::endl;
  os << indent << "CFL: " << m_CFL  << std::endl;
}

} // end namspace itk

#endif /* itkDensePengReinitializeLevelSetFunction_hxx */