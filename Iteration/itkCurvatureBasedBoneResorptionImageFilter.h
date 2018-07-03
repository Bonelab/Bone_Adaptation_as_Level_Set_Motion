// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkCurvatureBasedBoneResorptionImageFilter_h
#define itkCurvatureBasedBoneResorptionImageFilter_h

#include "itkDenseFiniteDifferenceImageFilter.h"
#include "itkCurvatureBasedBoneResorptionFunction.h"
#include "itkMath.h"
#include "itkNumericTraits.h"
#include "itkDensePengReinitializeLevelSetImageFilter.h"

namespace itk
{

/** Taken from itkMacro */
#define itkCustomSetMacro(name, type)                       \
  virtual void Set##name (const type _arg)                  \
  {                                                         \
    itkDebugMacro("setting " #name " to " << _arg);         \
CLANG_PRAGMA_PUSH                                           \
CLANG_SUPPRESS_Wfloat_equal                                 \
    if ( this->Get##name() != _arg )  \
    {                                                       \
      this->m_ResorptionFunction->Set##name( _arg );        \
      this->Modified();                                     \
    }                                                       \
CLANG_PRAGMA_POP                                            \
  }

#define itkCustomGetConstMacro(name, type)                   \
  virtual type Get##name () const                            \
  {                                                          \
    return this->m_ResorptionFunction->Get##name();          \
  }

template< typename TInputImage, typename TMaskImage = TInputImage, typename TOutputImage = TInputImage >
class ITK_TEMPLATE_EXPORT CurvatureBasedBoneResorptionImageFilter:
  public DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class type aliases. */
  using Self = CurvatureBasedBoneResorptionImageFilter;
  using Superclass = DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CurvatureBasedBoneResorptionImageFilter,
               DenseFiniteDifferenceImageFilter);

  /** Image types */
  using InputImageType = typename Superclass::InputImageType;
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using MaskImageType = TMaskImage;
  using MaskImageTypePointer = typename MaskImageType::Pointer;
  itkStaticConstMacro(ImageDimension, unsigned int,  TInputImage::ImageDimension);

  /** FiniteDifferenceFunction types */
  using FiniteDifferenceFunctionType = typename Superclass::FiniteDifferenceFunctionType;
  using CurvatureBasedBoneResorptionFunctionType = CurvatureBasedBoneResorptionFunction< OutputImageType, TMaskImage >;

  /** Level set reinitialization types. */
  using ReinitializationFilterType = DensePengReinitializeLevelSetImageFilter< TOutputImage >;

  /** Other types. */
  using PixelType = typename TInputImage::PixelType;
  using RealType = typename NumericTraits< PixelType >::RealType;
  using TimeStepType = typename Superclass::TimeStepType;
  using ValueType = typename CurvatureBasedBoneResorptionFunctionType::ScalarValueType;

  /** Methods to set/get the mask image */
  itkSetInputMacro(MaskImage, MaskImageType);
  itkGetInputMacro(MaskImage, MaskImageType);

  /** Setter and Getter for TotalTime */
  itkSetMacro(TotalTime, RealType);
  itkGetConstMacro(TotalTime, RealType);

  /** Getter for ElapsedTime */
  itkGetConstMacro(ElapsedTime, RealType);

  /** Set/Get the propagation weight. */
  itkCustomSetMacro(PropagationWeight, RealType);
  itkCustomGetConstMacro(PropagationWeight, RealType);

  /** Set/Get the curvature weight. */
  itkCustomSetMacro(CurvatureWeight, RealType);
  itkCustomGetConstMacro(CurvatureWeight, RealType);

  /** Set/Get epsilon. */
  itkCustomSetMacro(Epsilon, RealType);
  itkCustomGetConstMacro(Epsilon, RealType);

  /** Set/Get CFL number. Should be set such that: 0<= CFL <= 1. The smaller the
   *  value, the more iterations. */
  itkCustomSetMacro(CFL, RealType);
  itkCustomGetConstMacro(CFL, RealType);

protected:
  CurvatureBasedBoneResorptionImageFilter();
  ~CurvatureBasedBoneResorptionImageFilter() override {}
  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Explicitly specify the halting time. */
  bool Halt() override;

  /** We need to reinitialize the distance transform. */
  void InitializeIteration() override;

  /** Specific instance of resolving time step */
  TimeStepType ResolveTimeStep(const std::vector< TimeStepType >& timeStepList,
                                       const std::vector< bool >& valid ) const override;

  /** Need to capture and sum up each dt */
  TimeStepType CalculateChange() ITK_OVERRIDE;

  /** Need to setup timing variables */
  void Initialize() override;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(CurvatureBasedBoneResorptionImageFilter);

  /** Define the total time variable */
  RealType  m_TotalTime;
  RealType  m_ElapsedTime;

  /** Reinitialize filter */
  typename ReinitializationFilterType::Pointer m_ReinitializationFilter;

  /** Function pointer */
  typename CurvatureBasedBoneResorptionFunctionType::Pointer m_ResorptionFunction;
};
} // end namspace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCurvatureBasedBoneResorptionImageFilter.hxx"
#endif

#endif /* itkCurvatureBasedBoneResorptionImageFilter_h */