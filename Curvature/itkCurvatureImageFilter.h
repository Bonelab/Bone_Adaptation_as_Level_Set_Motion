// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkCurvatureImageFilter_h
#define itkCurvatureImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkCurvatureBasedBoneResorptionFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNumericTraits.h"
#include "itkVector.h"

namespace itk {

template< typename TInputImage, typename TOutputImage = TInputImage >
class CurvatureImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class type aliases. */
  using Self          = CurvatureImageFilter;
  using Superclass    = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer       = SmartPointer< Self >;
  using ConstPointer  = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CurvatureImageFilter, ImageToImageFilter);

  /** Image types */
  using InputImageType          = TInputImage;
  using InputImageConstPointer  = typename InputImageType::ConstPointer;
  using InputPixelType          = typename InputImageType::PixelType;
  using OutputImageType         = TOutputImage;
  using OutputImagePointer      = typename OutputImageType::Pointer;
  using OutputPixelType         = typename OutputImageType::PixelType;
  using OutputImageRegionType   = typename Superclass::OutputImageRegionType;
  using PixelRealType           = typename NumericTraits< InputPixelType >::RealType;
  using ScalarValueType         = PixelRealType;
  using SizeType                = typename OutputImageType::SizeType;
  itkStaticConstMacro(ImageDimension, unsigned int,  TInputImage::ImageDimension);

  /** Iterator types. */
  using DefaultBoundaryConditionType  = ZeroFluxNeumannBoundaryCondition< InputImageType >;
  using NeighborhoodType              = ConstNeighborhoodIterator< InputImageType, DefaultBoundaryConditionType >;
  using RadiusType                    = typename NeighborhoodType::RadiusType;
  using NeighborhoodScalesType        = Vector< PixelRealType, itkGetStaticConstMacro(ImageDimension) >;

protected:
  CurvatureImageFilter();
  ~CurvatureImageFilter() ITK_OVERRIDE {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** Multithreaded filter */
  virtual void ThreadedGenerateData(const OutputImageRegionType &,
                                    ThreadIdType) ITK_OVERRIDE;
  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

  /** Same data struct for computing derivatives, etc.
   */
  struct GlobalDataStruct {
    /** Hessian matrix */
    vnl_matrix_fixed< ScalarValueType,
                      itkGetStaticConstMacro(ImageDimension),
                      itkGetStaticConstMacro(ImageDimension) > m_dxy;

    /** Array of first derivatives */
    ScalarValueType m_dx[itkGetStaticConstMacro(ImageDimension)];

    ScalarValueType m_dx_forward[itkGetStaticConstMacro(ImageDimension)];
    ScalarValueType m_dx_backward[itkGetStaticConstMacro(ImageDimension)];

    ScalarValueType m_GradMagSqr;
  };

  /** Given the spatial derivatives, compute the curvature term $\kappa$. Note
   *  that this is not the same equation as appears in level sets but instead
   *  just curvature without the gradient magnitude.
   */
  GlobalDataStruct ComputeDerivatives( const NeighborhoodType & neighborhood);

  /** Given the spatial derivatives, compute the curvature term $\kappa$. Note
   *  that this is not the same equation as appears in level sets but instead
   *  just curvature without the gradient magnitude.
   */
  PixelRealType ComputeCurvatureTerm(GlobalDataStruct globalData);

const NeighborhoodScalesType ComputeNeighborhoodScales() const;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(CurvatureImageFilter);

  RadiusType      m_Radius;
  PixelRealType   m_ScaleCoefficients[ImageDimension];
  OffsetValueType m_Center;
  OffsetValueType m_xStride[itkGetStaticConstMacro(ImageDimension)];
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCurvatureImageFilter.hxx"
#endif

#endif /* itkCurvatureImageFilter_h */