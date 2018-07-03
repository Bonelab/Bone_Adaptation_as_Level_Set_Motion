// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkCurvatureImageFilter_hxx
#define itkCurvatureImageFilter_hxx

#include "itkCurvatureImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkMath.h"
#include "itkImageRegionIterator.h"

namespace itk {

template< typename TInputImage, typename TOutputImage >
CurvatureImageFilter< TInputImage, TOutputImage >
::CurvatureImageFilter() {
}

template< typename TInputImage, typename TOutputImage >
typename CurvatureImageFilter< TInputImage, TOutputImage >::GlobalDataStruct
CurvatureImageFilter< TInputImage, TOutputImage >
::ComputeDerivatives( const NeighborhoodType & it)
{ 
  // Heavily based on itkLevelSetFunction.hxx
  unsigned int          i, j;
  const ScalarValueType ZERO = NumericTraits< ScalarValueType >::ZeroValue();
  const ScalarValueType center_value  = it.GetCenterPixel();
  const NeighborhoodScalesType neighborhoodScales = this->ComputeNeighborhoodScales();

  GlobalDataStruct gd;

  // Compute the Hessian matrix and various other derivatives.  Some of these
  // derivatives may be used by overloaded virtual functions.
  gd.m_GradMagSqr = 1.0e-6;
  for ( i = 0; i < ImageDimension; i++ )
  {
    const auto positionA = static_cast< unsigned int >( m_Center + m_xStride[i] );
    const auto positionB = static_cast< unsigned int >( m_Center - m_xStride[i] );

    gd.m_dx[i] = 0.5 * ( it.GetPixel(positionA)
                          - it.GetPixel(positionB) ) * neighborhoodScales[i];

    gd.m_dxy[i][i] = ( it.GetPixel(positionA)
                        + it.GetPixel(positionB) - 2.0 * center_value )
                      * Math::sqr(neighborhoodScales[i]);

    gd.m_dx_forward[i]  = ( it.GetPixel(positionA) - center_value ) * neighborhoodScales[i];

    gd.m_dx_backward[i] = ( center_value - it.GetPixel(positionB) ) * neighborhoodScales[i];

    gd.m_GradMagSqr += gd.m_dx[i] * gd.m_dx[i];

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

      gd.m_dxy[i][j] = gd.m_dxy[j][i] = 0.25 * ( it.GetPixel(positionAa)
                                                   - it.GetPixel(positionBa)
                                                   - it.GetPixel(positionCa)
                                                   + it.GetPixel(positionDa) )
                                          * neighborhoodScales[i] * neighborhoodScales[j];
      }
  }

  return gd;
}

template< typename TInputImage, typename TOutputImage >
typename CurvatureImageFilter< TInputImage, TOutputImage >::PixelRealType
CurvatureImageFilter< TInputImage, TOutputImage >
::ComputeCurvatureTerm(GlobalDataStruct gd)
{
  // Heavily based on itkCurvatureFlowFunction.hxx
  ScalarValueType temp, laplace = 0.0, cross=0.0;
  unsigned int          i, j;

  // accumulate laplacian
  for ( i = 0; i < ImageDimension; i++ ) {
    laplace += gd.m_dxy[i][i];
  }
  laplace /= std::sqrt((double)gd.m_GradMagSqr);

  // accumulate dx dx dxy terms
  for ( i = 0; i < ImageDimension; i++ ) {
    for ( j = 0; j < ImageDimension; j++ ) {
      cross += gd.m_dx[i] * gd.m_dx[j] * gd.m_dxy[i][j];
    }
  }
  cross /= (std::sqrt((double)gd.m_GradMagSqr)*gd.m_GradMagSqr);

  return (PixelRealType)(laplace+cross);
}

template< typename TInputImage, typename TOutputImage >
void
CurvatureImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  // Setup iterator radius
  m_Radius.Fill(1);

  // Dummy neighborhood.
  NeighborhoodType it;
  it.SetRadius(m_Radius);

  // Find the center index of the neighborhood.
  m_Center =  it.Size() / 2;

  // Get the stride length for each axis.
  for ( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_xStride[i] = it.GetStride(i);
  }

  // Setup scale
  InputImageConstPointer input = this->GetInput();
  if ( input == ITK_NULLPTR )
  {
    itkExceptionMacro("Input image is ITK_NULLPTR");
  }

  typedef typename TInputImage::SpacingType SpacingType;
  const SpacingType & spacing = input->GetSpacing();

  for ( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_ScaleCoefficients[i] = 1.0 / static_cast< double >( spacing[i] );
  }
}

template< typename TInputImage, typename TOutputImage >
void
CurvatureImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType & regionToProcess, ThreadIdType threadId)
{ 
  InputImageConstPointer input = this->GetInput();
  OutputImagePointer output = this->GetOutput();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because input has been copied to output.
  using FaceCalculatorType = typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >;
  using FaceListType = typename FaceCalculatorType::FaceListType;
  FaceCalculatorType faceCalculator;
  using OutputIteratorType = ImageRegionIterator< TOutputImage >;

  FaceListType faceList = faceCalculator(input, regionToProcess, m_Radius);
  typename FaceListType::iterator fIt = faceList.begin();
  typename FaceListType::iterator fEnd = faceList.end();

  // Process the non-boundary region.
  NeighborhoodType    nD(m_Radius, input, *fIt);
  OutputIteratorType  nU(output,  *fIt);
  nD.GoToBegin();
  while ( !nD.IsAtEnd() )
  {
    nU.Value() = static_cast< OutputPixelType > (ComputeCurvatureTerm(ComputeDerivatives(nD)));
    ++nD;
    ++nU;
  }

  // Process each of the boundary faces.
  for ( ++fIt; fIt != fEnd; ++fIt )
  {
    NeighborhoodType    bD(m_Radius, input, *fIt);
    OutputIteratorType  bU(output, *fIt);

    bD.GoToBegin();
    bU.GoToBegin();
    while ( !bD.IsAtEnd() )
    {
      bU.Value() = static_cast< OutputPixelType > (ComputeCurvatureTerm(ComputeDerivatives(bD)));
      ++bD;
      ++bU;
    }
  }
}

template< typename TInputImage, typename TOutputImage >
const typename CurvatureImageFilter< TInputImage, TOutputImage >::NeighborhoodScalesType
CurvatureImageFilter< TInputImage, TOutputImage >
::ComputeNeighborhoodScales() const
{
  NeighborhoodScalesType neighborhoodScales;

  neighborhoodScales.Fill(0.0);
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( this->m_Radius[i] > 0 )
      {
      neighborhoodScales[i] = this->m_ScaleCoefficients[i] / this->m_Radius[i];
      }
    }
  return neighborhoodScales;
}

template< typename TInputImage, typename TOutputImage >
void
CurvatureImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif /* itkCurvatureImageFilter_hxx */