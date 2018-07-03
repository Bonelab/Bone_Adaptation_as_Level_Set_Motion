// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkCurvatureBasedBoneResorptionFunction_hxx
#define itkCurvatureBasedBoneResorptionFunction_hxx

#include "itkCurvatureBasedBoneResorptionFunction.h"
#include "itkMath.h"

namespace itk
{

template< typename TImageType, typename TMaskImageType >
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >::
CurvatureBasedBoneResorptionFunction() :
  m_PropagationWeight(NumericTraits< RealType >::Zero),
  m_CurvatureWeight(NumericTraits< RealType >::Zero),
  m_Epsilon(1.0e-6),
  m_CFL(0.5f)
{
  m_MaskImage = ITK_NULLPTR;
  m_Interpolator = InterpolatorType::New();
}

template< typename TImageType, typename TMaskImageType >
typename CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >::ScalarValueType
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >
::ComputePropagationTerm( const NeighborhoodType & itkNotUsed(it),
                          GlobalDataStruct *gd,
                          const FloatOffsetType & itkNotUsed(offset))
{
  // Heavily based on itkLevelSetFunction.hxx
  unsigned int i;
  const ScalarValueType ZERO = NumericTraits< ScalarValueType >::ZeroValue();
  ScalarValueType propagation_gradient = ZERO;

  if ( this->GetPropagationWeight() > ZERO ) {
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

  return propagation_gradient;
}

template< typename TImageType, typename TMaskImageType >
typename CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >::ScalarValueType
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >
::ComputeCurvatureTerm( const NeighborhoodType & itkNotUsed(it),
                        GlobalDataStruct *gd,
                        const FloatOffsetType & itkNotUsed(offset))
{
  // Heavily based on itkCurvatureFlowFunction.hxx
  ScalarValueType temp, update = 0.0;
  unsigned int          i, j;

  // accumulate dx^2 * (dyy + dzz) terms
  for ( i = 0; i < ImageDimension; i++ ) {
    temp = 0.0;
    for ( j = 0; j < ImageDimension; j++ ) {
      if ( j == i ) { continue; }
      temp += gd->m_dxy[j][j];
    }

    update += temp * Math::sqr( (double)gd->m_dx[i] );
  }

  // accumlate -2 * dx * dy * dxy terms
  for ( i = 0; i < ImageDimension; i++ ) {
    for ( j = i + 1; j < ImageDimension; j++ ) {
      update -= 2 * gd->m_dx[i] * gd->m_dx[j]
                * gd->m_dxy[i][j];
      }
  }

  update /= gd->m_GradMagSqr;
  return update;
}

template< typename TImageType, typename TMaskImageType >
void
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >
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

template< typename TImageType, typename TMaskImageType >
typename CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >::PixelType
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >
::ComputeUpdate(const NeighborhoodType & it, void *globalData,
                const FloatOffsetType & offset)
{
  // Check if this is a valid pixel. First, we check if an image was actually set.
  const PixelType ZERO_PIXEL = NumericTraits< PixelType >::ZeroValue();
  if (m_MaskImage) {
    // Build continuous index
    IndexType           idx = it.GetIndex();
    ContinuousIndexType cdx;

    for ( unsigned k = 0; k < ImageDimension; ++k )
    {
      cdx[k] = static_cast< double >( idx[k] ) - offset[k];
    }

    // If we are outside buffer, return zero
    if ( !m_Interpolator->IsInsideBuffer(cdx) )
    {
      return ZERO_PIXEL;
    }

    // If the mask value is zero, return zero
    if ( Math::AlmostEquals(m_Interpolator->EvaluateAtContinuousIndex(cdx), NumericTraits< MaskPixelType >::ZeroValue()))
    {
      return ZERO_PIXEL;
    }
  }

  // Heavily based on itkLevelSetFunction.hxx
  unsigned int          i, j;
  const ScalarValueType ZERO = NumericTraits< ScalarValueType >::ZeroValue();
  const ScalarValueType center_value  = it.GetCenterPixel();

  const NeighborhoodScalesType neighborhoodScales = this->ComputeNeighborhoodScales();

  ScalarValueType propagation_term, curvature_term;

  // Global data structure
  auto * gd = (GlobalDataStruct *)globalData;

  // Compute the Hessian matrix and various other derivatives.  Some of these
  // derivatives may be used by overloaded virtual functions.
  gd->m_GradMagSqr = m_Epsilon;
  for ( i = 0; i < ImageDimension; i++ )
  {
    const auto positionA = static_cast< unsigned int >( m_Center + m_xStride[i] );
    const auto positionB = static_cast< unsigned int >( m_Center - m_xStride[i] );

    gd->m_dx[i] = 0.5 * ( it.GetPixel(positionA)
                          - it.GetPixel(positionB) ) * neighborhoodScales[i];

    gd->m_dxy[i][i] = ( it.GetPixel(positionA)
                        + it.GetPixel(positionB) - 2.0 * center_value )
                      * Math::sqr(neighborhoodScales[i]);

    gd->m_dx_forward[i]  = ( it.GetPixel(positionA) - center_value ) * neighborhoodScales[i];

    gd->m_dx_backward[i] = ( center_value - it.GetPixel(positionB) ) * neighborhoodScales[i];

    gd->m_GradMagSqr += gd->m_dx[i] * gd->m_dx[i];

    for ( j = i + 1; j < ImageDimension; j++ )
      {
      const auto positionAa = static_cast< unsigned int >(
        m_Center - m_xStride[i] - m_xStride[j] );
      const auto positionBa = static_cast< unsigned int >(
        m_Center - m_xStride[i] + m_xStride[j] );
      const auto positionCa = static_cast< unsigned int >(
        m_Center + m_xStride[i] - m_xStride[j] );
      const auto positionDa = static_cast< unsigned int >(
        m_Center + m_xStride[i] + m_xStride[j] );

      gd->m_dxy[i][j] = gd->m_dxy[j][i] = 0.25 * ( it.GetPixel(positionAa)
                                                   - it.GetPixel(positionBa)
                                                   - it.GetPixel(positionCa)
                                                   + it.GetPixel(positionDa) )
                                          * neighborhoodScales[i] * neighborhoodScales[j];
      }
  }

  // Compute curvature
  if ( Math::NotAlmostEquals(m_CurvatureWeight, ZERO) ) {
    curvature_term = m_CurvatureWeight * ComputeCurvatureTerm(it, gd, offset);
  } else {
    curvature_term = ZERO;
  }

  // Compute propagation
  if ( Math::NotAlmostEquals(m_PropagationWeight, ZERO) ) {
    propagation_term = m_PropagationWeight * ComputePropagationTerm(it, gd, offset);
  } else {
    propagation_term = ZERO;
  }

  return (PixelType)(curvature_term - propagation_term);
}

template< typename TImageType, typename TMaskImageType >
typename CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >::TimeStepType
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >
::ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const
{
  TimeStepType dt;
  unsigned int i;

  // const NeighborhoodScalesType neighborhoodScales = this->ComputeNeighborhoodScales();

  // // Compute the min spacing
  // RealType min_spacing = NumericTraits< RealType >::max();
  // RealType min_spacing_squared = NumericTraits< RealType >::max();
  // for ( i = 0; i < ImageDimension; i++ ) {
  //   if ( Math::NotAlmostEquals(neighborhoodScales[i], NumericTraits< RealType >::ZeroValue()) ) {
  //     min_spacing = std::min(min_spacing, 1.0/neighborhoodScales[i]);
  //     min_spacing_squared = std::min(min_spacing_squared, 1.0/Math::sqr(neighborhoodScales[i]));
  //   }
  // }

  // Compute the min spacing
  RealType min_spacing = NumericTraits< RealType >::max();
  RealType min_spacing_squared = NumericTraits< RealType >::max();
  for ( i = 0; i < ImageDimension; i++ )
  {
    min_spacing = std::min(1.0/(this->m_ScaleCoefficients[i]), min_spacing);
    min_spacing_squared = std::min(1.0/Math::sqr(this->m_ScaleCoefficients[i]), min_spacing_squared);
  }

  RealType D = std::abs(m_PropagationWeight) / min_spacing + 2 * std::abs(m_CurvatureWeight) / (min_spacing_squared);

  if ( Math::NotAlmostEquals(D, NumericTraits< RealType >::ZeroValue()) ) {
    dt = (TimeStepType) ( m_CFL / D);
  } else {
    dt = (TimeStepType)0.0f;
  }
  return dt;
}

template< typename TImageType, typename TMaskImageType >
void
CurvatureBasedBoneResorptionFunction< TImageType, TMaskImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PropagationWeight: " << m_PropagationWeight << std::endl;
  os << indent << "CurvatureWeight: " << m_CurvatureWeight << std::endl;
  os << indent << "Epsilon: " << m_Epsilon << std::endl;
  os << indent << "CFL: " << m_CFL  << std::endl;
  os << indent << "Interpolator: " << m_Interpolator  << std::endl;
}

} // end namespace itk

#endif /* itkCurvatureBasedBoneResorptionFunction_hxx */