// B.A. Besler, 2018
// Bone Imaging Laboratory

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCurvatureBasedBoneResorptionImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCommand.h"
#include "itkAddImageFilter.h"
#include "itkReinitializeLevelSetImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

const unsigned int DIMENSION = 3;
using LevelSetPixelType = double;
using LevelSetImageType = typename itk::Image< LevelSetPixelType, DIMENSION >;
using LevelSetImageTypePointer = LevelSetImageType::Pointer;
using SegmentationPixelType = unsigned char;
using SegmentationImageType = typename itk::Image< SegmentationPixelType, DIMENSION >;
using SegmentationImageTypePointer = SegmentationImageType::Pointer;
using SegmentationWriterType = typename itk::ImageFileWriter<SegmentationImageType>;
using SegmentationReaderType = typename itk::ImageFileReader<SegmentationImageType>;
using DistanceReaderType = typename itk::ImageFileReader<LevelSetImageType>;
using DistanceWriterType = typename itk::ImageFileWriter<LevelSetImageType>;

using DistanceMapFilterType = typename itk::SignedMaurerDistanceMapImageFilter< LevelSetImageType, LevelSetImageType >;
using AddFilterType = typename itk::AddImageFilter< LevelSetImageType >;
using ResorptionFilterType = typename itk::CurvatureBasedBoneResorptionImageFilter< LevelSetImageType, SegmentationImageType, LevelSetImageType >;
using BinaryThresholdImageFilterType = typename itk::BinaryThresholdImageFilter< LevelSetImageType, SegmentationImageType >;
using InputThresholdImageFilterType = typename itk::BinaryThresholdImageFilter< SegmentationImageType, LevelSetImageType >;
using ReinitializationFilterType = itk::ReinitializeLevelSetImageFilter< LevelSetImageType >;

template<class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  using Self = CommandIterationUpdate;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  void Execute(itk::Object *caller,
               const itk::EventObject & event) override
  {
    Execute( (const itk::Object *) caller, event);
  }
  void Execute(const itk::Object * object,
               const itk::EventObject & event) override
  {
    const auto * filter = static_cast< const TFilter * >( object );
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      { return; }

    std::cout << "\r";
    std::cout << filter->GetElapsedIterations() << ": ";
    std::cout << filter->GetElapsedTime() << " - ";
    std::cout << filter->GetTotalTime();

    // Quick test to see if we are done
    if ( itk::Math::AlmostEquals(filter->GetElapsedTime(), filter->GetTotalTime())) {
      std::cout << std::endl;
    }

    std::cout << std::flush;
  }
};

int main( int argc, char *argv[] )
{
  if( argc < 7 )
  {
		std::cerr << "[USAGE]: " << argv[0] << " ";
		std::cerr << "input_image output_image output_seg ";
		std::cerr << "curvature_weight propogation_weight total_time" << std::endl;
    std::cerr << "Optional:" << std::endl;
    std::cerr << "  input_mask" << std::endl;
    std::cerr << "Notes:" << std::endl;
    std::cerr << "  It is assumed that the image data is isotropic. If not, care must be taken with the initial distance transform." << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_file_name, mask_file_name="", output_file_name, seg_file_name;
  double curvature_weight, propogation_weight, total_time;
  input_file_name = argv[1];
  output_file_name = argv[2];
	seg_file_name = argv[3];
  curvature_weight = atof(argv[4]);
  propogation_weight = atof(argv[5]);
  total_time = atof(argv[6]);
  if (argc > 7) {
		mask_file_name= argv[7];
  }

  std::cout << "Parameters:" << std::endl;
	std::cout << "  input_file_name:		" << input_file_name << std::endl;
	std::cout << "  mask_file_name:			" << mask_file_name << std::endl;
	std::cout << "  output_file_name:		" << output_file_name << std::endl;
	std::cout << "  seg_file_name:			" << seg_file_name << std::endl;
	std::cout << "  curvature_weight:		" << curvature_weight << std::endl;
	std::cout << "  propogation_weight:	" << propogation_weight << std::endl;
	std::cout << "  total_time:					" << total_time << std::endl;
  std::cout << std::endl;

	std::cout << "Setting up reader for " << input_file_name << std::endl;
	DistanceReaderType::Pointer input_reader = DistanceReaderType::New();
	input_reader->SetFileName(input_file_name);
	input_reader->Update();

  std::cout << "Setting up resorption filter with curvature " << curvature_weight << ", propogation " << propogation_weight << ", and total time " << total_time << std::endl;
  ResorptionFilterType::Pointer resorption_filter = ResorptionFilterType::New();
  resorption_filter->SetInput(input_reader->GetOutput());
  resorption_filter->SetPropagationWeight(propogation_weight);
  resorption_filter->SetCurvatureWeight(curvature_weight);
  resorption_filter->SetTotalTime(total_time);

  SegmentationReaderType::Pointer mask_reader = SegmentationReaderType::New();
  if (!mask_file_name.empty()) {
    std::cout << "Setting up mask file " << mask_file_name << std::endl;
    mask_reader->SetFileName(mask_file_name);
    resorption_filter->SetMaskImage(mask_reader->GetOutput());
  } else {
    std::cout << "No mask file provided" << std::endl;
  }

  using CommandType = CommandIterationUpdate<ResorptionFilterType>;
  CommandType::Pointer observer = CommandType::New();
  resorption_filter->AddObserver( itk::IterationEvent(), observer );

  resorption_filter->Update();

  std::cout << "Thresholding to get output" << std::endl;
  BinaryThresholdImageFilterType::Pointer thresh_filter = BinaryThresholdImageFilterType::New();
  thresh_filter->SetInput(resorption_filter->GetOutput());
  thresh_filter->SetInsideValue(1);
  thresh_filter->SetOutsideValue(0);
  thresh_filter->SetUpperThreshold(0);

	std::cout << "Writing seg to " << seg_file_name << std::endl;
	SegmentationWriterType::Pointer writer = SegmentationWriterType::New();
	writer->SetInput(thresh_filter->GetOutput());
	writer->SetFileName(seg_file_name);
	writer->Update();

  std::cout << "Writing output to " << output_file_name << std::endl;
  DistanceWriterType::Pointer dw = DistanceWriterType::New();
  dw->SetFileName(output_file_name);
  dw->SetInput(resorption_filter->GetOutput());
  dw->Update();

  return EXIT_SUCCESS;
}
