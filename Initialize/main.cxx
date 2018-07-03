// B.A. Besler, 2018
// Bone Imaging Laboratory

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDensePengReinitializeLevelSetImageFilter.h"
#include "itkAddImageFilter.h"

const unsigned int DIMENSION = 3;
using LevelSetPixelType = double;
using LevelSetImageType = typename itk::Image< LevelSetPixelType, DIMENSION >;
using LevelSetImageTypePointer = LevelSetImageType::Pointer;
using SegmentationPixelType = unsigned char;
using SegmentationImageType = typename itk::Image< SegmentationPixelType, DIMENSION >;
using SegmentationImageTypePointer = SegmentationImageType::Pointer;
using SegmentationReaderType = typename itk::ImageFileReader<SegmentationImageType>;
using DistanceWriterType = typename itk::ImageFileWriter<LevelSetImageType>;

using DistanceMapFilterType = typename itk::SignedMaurerDistanceMapImageFilter< LevelSetImageType, LevelSetImageType >;
using AddFilterType = typename itk::AddImageFilter< LevelSetImageType >;
using BinaryThresholdImageFilterType = typename itk::BinaryThresholdImageFilter< LevelSetImageType, SegmentationImageType >;
using InputThresholdImageFilterType = typename itk::BinaryThresholdImageFilter< SegmentationImageType, LevelSetImageType >;
using ReinitializationFilterType = itk::DensePengReinitializeLevelSetImageFilter< LevelSetImageType >;

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
		std::cerr << "[USAGE]: " << argv[0] << " ";
		std::cerr << "input_image output_phi" << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_file_name, output_file_name;
  input_file_name = argv[1];
  output_file_name = argv[2];

  std::cout << "Parameters:" << std::endl;
  std::cout << "  input_file_name:    " << input_file_name << std::endl;
  std::cout << "  output_file_name:   " << output_file_name << std::endl;
  std::cout << std::endl;

  std::cout << "Setting up reader for " << input_file_name << std::endl;
  SegmentationReaderType::Pointer input_reader = SegmentationReaderType::New();
  input_reader->SetFileName(input_file_name);
  input_reader->Update();

  std::cout << "Thresholding for distance transform" << std::endl;
  InputThresholdImageFilterType::Pointer input_threshold_filter = InputThresholdImageFilterType::New();
  input_threshold_filter->SetInput(input_reader->GetOutput());
  input_threshold_filter->SetInsideValue(1);
  input_threshold_filter->SetOutsideValue(0);
  input_threshold_filter->SetLowerThreshold(1);

  std::cout << "Setting up distance transform" << std::endl;
  DistanceMapFilterType::Pointer distance_filter = DistanceMapFilterType::New();
  distance_filter->SetInput(input_threshold_filter->GetOutput());
  distance_filter->SquaredDistanceOff();
  distance_filter->UseImageSpacingOn();
  distance_filter->InsideIsPositiveOff();

	// This advances the distance transform to the edge of the image
  AddFilterType::Pointer add_filter = AddFilterType::New();
  add_filter->SetInput1(distance_filter->GetOutput());
  add_filter->SetConstant2(input_reader->GetOutput()->GetSpacing()[0]/-2.0);
  add_filter->Update();

	std::cout << "Reinitializing" << std::endl;
  ReinitializationFilterType::Pointer reinit_filter = ReinitializationFilterType::New();
  reinit_filter->SetInput(add_filter->GetOutput());
	reinit_filter->SetMaximumRMSError(1e-6);
	reinit_filter->SetNumberOfIterations(1000);
  reinit_filter->Update();

	std::cout << "Writing output to " << output_file_name << std::endl;
  DistanceWriterType::Pointer dw = DistanceWriterType::New();
  dw->SetFileName(output_file_name);
  dw->SetInput(reinit_filter->GetOutput());
  //dw->SetInput(add_filter->GetOutput());
  dw->Update();

  return EXIT_SUCCESS;
}
