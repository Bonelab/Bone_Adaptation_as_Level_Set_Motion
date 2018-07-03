// B.A. Besler, 2018
// Bone Imaging Laboratory

#ifndef itkCurvatureBasedBoneResorptionFunction_h
#define itkCurvatureBasedBoneResorptionFunction_h

#include "itkFiniteDifferenceFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"

namespace itk
{

template< typename TImageType, typename TMaskImageType = TImageType >
class ITK_TEMPLATE_EXPORT CurvatureBasedBoneResorptionFunction:
  public FiniteDifferenceFunction< TImageType >
{
public:
  /**  Standard class type aliases. */
  using Self          = CurvatureBasedBoneResorptionFunction;
  using Superclass    = FiniteDifferenceFunction< TImageType >;
  using Pointer       = SmartPointer< Self >;
  using ConstPointer  = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(CurvatureBasedBoneResorptionFunction,
               FiniteDifferenceFunction);

  /** Extract superclass dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,  Superclass::ImageDimension);

  /** Convenient typedefs. */
  using ImageType               = typename Superclass::ImageType;
  using PixelType               = typename Superclass::PixelType;
  using PixelRealType           = typename Superclass::PixelRealType;
  using ScalarValueType         = PixelType;
  using RealType                = typename NumericTraits< PixelType >::RealType;
  using RadiusType              = typename Superclass::RadiusType;
  using NeighborhoodType        = typename Superclass::NeighborhoodType;
  using NeighborhoodScalesType  = typename Superclass::NeighborhoodScalesType;
  using FloatOffsetType         = typename Superclass::FloatOffsetType;
  using IndexType               = typename ImageType::IndexType;
  using MaskImageType           = TMaskImageType;
  using MaskPixelType           = typename MaskImageType::PixelType;
  using TimeStepType            = typename Superclass::TimeStepType;

  /** Setup interpolator for mask */
  using InterpolatorType = NearestNeighborInterpolateImageFunction< MaskImageType >;
  using ContinuousIndexType = typename InterpolatorType::ContinuousIndexType;

  /** Returns a pointer to a global data structure that is passed to this
   * object from the solver at each calculation.  The idea is that the
   * solver holds the state of any global values needed to calculate the
   * time step, while the equation object performs the actual calculations.
   *
   * The global data should also be initialized in this method.
   * */
  virtual void * GetGlobalDataPointer() const ITK_OVERRIDE
  {
    GlobalDataStruct *ans = new GlobalDataStruct();
    return ans;
  }

  /** When the finite difference solver filter has finished using a global
   * data pointer, it passes it to this method, which frees the memory.
   * The solver cannot free the memory because it does not know the type
   * to which the pointer points. */
  void ReleaseGlobalDataPointer(void *GlobalData) const ITK_OVERRIDE
  { delete (GlobalDataStruct *)GlobalData; }

  /** Computes the time step for an update given a global data structure.
   * The data used in the computation may take different forms depending on
   * the nature of the equations.  This global data cannot be kept in the
   * instance of the equation object itself since the equation object must
   * remain stateless for thread safety.  The global data is therefore managed
   * for each thread by the finite difference solver filters. */
  virtual TimeStepType ComputeGlobalTimeStep(void *GlobalData) const ITK_OVERRIDE;

  /** Compute the equation value at an index in the image. */
  virtual PixelType ComputeUpdate( const NeighborhoodType & neighborhood,
                                   void *globalData,
                                   const FloatOffsetType & = FloatOffsetType(0.0) ) ITK_OVERRIDE;

  /** This method creates the appropriate member variable operators for the
   * level-set calculations.  The argument to this function is a the radius
   * necessary for performing the level-set calculations. */
  virtual void Initialize(const RadiusType & r);

  /** Setter/Getter for PropagationWeight. This is the absolute weight
   *  of propagation. */
  virtual void SetPropagationWeight(const RealType a)
  { m_PropagationWeight = a; }
  RealType GetPropagationWeight() const
  { return m_PropagationWeight; }

  /** Setter/Getter for CurvatureWeight. This is the absolute weight
   *  of curvature. */
  virtual void SetCurvatureWeight(const RealType a)
  { m_CurvatureWeight = a; }
  RealType GetCurvatureWeight() const
  { return m_CurvatureWeight; }

  /** Setter/Getter for Epsilon. This is a constant added to the 
   *  gradient magnitude squared to prevent divide by zero errors
   *  and keep the problem numerically stable. */
  virtual void SetEpsilon(const RealType a)
  { m_Epsilon = a; }
  RealType GetEpsilon() const
  { return m_Epsilon; }

  /** Setter/Getter for CFL. This is the maximum CFL number for
   *  computing the updating time step. 0.5 is conservative
   *  and 0.9 is a near-optimal. Should always be between 0 and
   *  1. The larger the value, the faster the solution. */
  virtual void SetCFL(const RealType a)
  { m_CFL = a; }
  RealType GetCFL() const
  { return m_CFL; }

  /** Set/Get the mask image for modifying calculations. */
  virtual const MaskImageType * GetMaskImage() const
  {
    return m_MaskImage.GetPointer();
  }
  virtual void SetMaskImage(const MaskImageType *f)
  {
    m_MaskImage = f;
    m_Interpolator->SetInputImage(m_MaskImage);
  }

protected:
  CurvatureBasedBoneResorptionFunction();
  ~CurvatureBasedBoneResorptionFunction() ITK_OVERRIDE {}

  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** A global data type for this class of equations.
   * 
   * See itkLevelSetFunction.h and itkCurvatureFlowImageFilter.h for
   * insight into why this structure.
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

  /** Given the spatial derivatives, compute the propogation term $| \del phi |$ */
  ScalarValueType ComputePropagationTerm( const NeighborhoodType & neighborhood,
                                          GlobalDataStruct *globalData,
                                          const FloatOffsetType &);

  /** Given the spatial derivatives, compute the curvature term $\kappa | \del phi |$ */
  ScalarValueType ComputeCurvatureTerm( const NeighborhoodType & neighborhood,
                                        GlobalDataStruct *globalData,
                                        const FloatOffsetType &);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(CurvatureBasedBoneResorptionFunction);

  /** Members for mask */
  typename MaskImageType::ConstPointer m_MaskImage;
  typename InterpolatorType::Pointer   m_Interpolator;

  /** Parameters */
  RealType  m_PropagationWeight;
  RealType  m_CurvatureWeight;
  RealType  m_Epsilon;
  RealType  m_CFL;

  /** The offset of the center pixel in the neighborhood. */
  OffsetValueType m_Center;

  /** Stride length along the y-dimension. */
  OffsetValueType m_xStride[itkGetStaticConstMacro(ImageDimension)];
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCurvatureBasedBoneResorptionFunction.hxx"
#endif

#endif /* itkCurvatureBasedBoneResorptionFunction_h */