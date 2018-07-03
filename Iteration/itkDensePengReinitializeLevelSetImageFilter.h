// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkDensePengReinitializeLevelSetImageFilter_h
#define itkDensePengReinitializeLevelSetImageFilter_h

#include "itkDenseFiniteDifferenceImageFilter.h"
#include "itkDensePengReinitializeLevelSetFunction.h"

namespace itk
{

/** Taken from itkMacro */
#define itkDensePengSetMacro(name, type)                    \
  virtual void Set##name (const type _arg)                  \
  {                                                         \
    itkDebugMacro("setting " #name " to " << _arg);         \
CLANG_PRAGMA_PUSH                                           \
CLANG_SUPPRESS_Wfloat_equal                                 \
    if ( this->Get##name() != _arg )  \
    {                                                       \
      this->m_ReshapeFunction->Set##name( _arg );           \
      this->Modified();                                     \
    }                                                       \
CLANG_PRAGMA_POP                                            \
  }

#define itkDensePengGetConstMacro(name, type)                \
  virtual type Get##name () const                            \
  {                                                          \
    return this->m_ReshapeFunction->Get##name();             \
  }

template< typename TInputImage, typename TOutputImage = TInputImage >
class ITK_TEMPLATE_EXPORT DensePengReinitializeLevelSetImageFilter:
  public DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class type aliases. */
  using Self          = DensePengReinitializeLevelSetImageFilter;
  using Superclass    = DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >;
  using Pointer       = SmartPointer< Self >;
  using ConstPointer  = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DensePengReinitializeLevelSetImageFilter,
               DenseFiniteDifferenceImageFilter);

  /** Image types */
  using InputImageType      = typename Superclass::InputImageType;
  using OutputImageType     = typename Superclass::OutputImageType;
  itkStaticConstMacro(ImageDimension, unsigned int,  TInputImage::ImageDimension);

  /** FiniteDifferenceFunction types */
  using FiniteDifferenceFunctionType              = typename Superclass::FiniteDifferenceFunctionType;
  using DensePengReinitializeLevelSetFunctionType = DensePengReinitializeLevelSetFunction< OutputImageType >;

  /** Other types. */
  using PixelType     = typename TInputImage::PixelType;
  using RealType      = typename NumericTraits< PixelType >::RealType;
  using TimeStepType  = typename Superclass::TimeStepType;
  using ValueType     = typename DensePengReinitializeLevelSetFunctionType::ScalarValueType;


  /** Set/Get epsilon. */
  itkDensePengSetMacro(Epsilon, RealType);
  itkDensePengGetConstMacro(Epsilon, RealType);

  /** Set/Get CFL number. Should be set such that: 0<= CFL <= 1. The smaller the
   *  value, the more iterations. */
  itkDensePengSetMacro(CFL, RealType);
  itkDensePengGetConstMacro(CFL, RealType);

protected:
  DensePengReinitializeLevelSetImageFilter();
  ~DensePengReinitializeLevelSetImageFilter() ITK_OVERRIDE {}

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DensePengReinitializeLevelSetImageFilter);

  /** Function pointer */
  typename DensePengReinitializeLevelSetFunctionType::Pointer m_ReshapeFunction;
};
} // end namspace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDensePengReinitializeLevelSetImageFilter.hxx"
#endif

#endif /* itkDensePengReinitializeLevelSetImageFilter_h */