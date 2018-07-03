// B.A. Besler, 2018
// Bone Imaging Laboratory

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCurvatureImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCommand.h"
#include "itkAddImageFilter.h"
#include "itkReinitializeLevelSetImageFilter.h"
#include "itkDensePengReinitializeLevelSetImageFilter.h"

const unsigned int DIMENSION = 3;
using LevelSetPixelType = double;
using LevelSetImageType = typename itk::Image< LevelSetPixelType, DIMENSION >;
using LevelSetImageTypePointer = LevelSetImageType::Pointer;
using SegmentationPixelType = unsigned char;
using SegmentationImageType = typename itk::Image< SegmentationPixelType, DIMENSION >;
using SegmentationImageTypePointer = SegmentationImageType::Pointer;
using SegmentationReaderType = typename itk::ImageFileReader<SegmentationImageType>;
using DistanceWriterType = typename itk::ImageFileWriter<LevelSetImageType>;

using DistanceMapFilterType = typename itk::SignedMaurerDistanceMapImageFilter< SegmentationImageType, LevelSetImageType >;
using AddFilterType = typename itk::AddImageFilter< LevelSetImageType >;
using InputThresholdImageFilterType = typename itk::BinaryThresholdImageFilter< SegmentationImageType, LevelSetImageType >;
using ReinitializationFilterType = itk::ReinitializeLevelSetImageFilter< LevelSetImageType >;
using DensePengReinitializationFilterType = itk::DensePengReinitializeLevelSetImageFilter< LevelSetImageType >;
using CurvatureImageFilterType = itk::CurvatureImageFilter< LevelSetImageType >;

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    std::cerr << "[USAGE]: " << argv[0] << " ";
    std::cerr << "input_image output_image" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Optional:" << std::endl;
    std::cerr << "  dt_image" << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file_name, output_file_name, dt_file_name="";
  input_file_name = argv[1];
  output_file_name = argv[2];
  if (argc > 3) {
    dt_file_name = argv[3];
  }

  std::cout << "Parameters:" << std::endl;
  std::cout << "  input_file_name:    " << input_file_name << std::endl;
  std::cout << "  output_file_name:   " << output_file_name << std::endl;
  std::cout << "  dt_file_name:       " << dt_file_name << std::endl;
  std::cout << std::endl;

  std::cout << "Setting up reader for " << input_file_name << std::endl;
  SegmentationReaderType::Pointer input_reader = SegmentationReaderType::New();
  input_reader->SetFileName(input_file_name);
  input_reader->Update();

  std::cout << "Thresholding for distance transform" << std::endl;
  InputThresholdImageFilterType::Pointer input_threshold_filter = InputThresholdImageFilterType::New();
  input_threshold_filter->SetInput(input_reader->GetOutput());
  input_threshold_filter->SetInsideValue(0);
  input_threshold_filter->SetOutsideValue(1);
  input_threshold_filter->SetLowerThreshold(1);

  std::cout << "Setting up distance transform" << std::endl;
  // DistanceMapFilterType::Pointer distance_filter = DistanceMapFilterType::New();
  // distance_filter->SetInput(input_threshold_filter->GetOutput());
  // distance_filter->SquaredDistanceOff();
  // distance_filter->UseImageSpacingOn();
  // distance_filter->InsideIsPositiveOff();

  // AddFilterType::Pointer add_filter = AddFilterType::New();
  // add_filter->SetInput1(distance_filter->GetOutput());
  // add_filter->SetConstant2(input_reader->GetOutput()->GetSpacing()[0]/-2.0);
  ReinitializationFilterType::Pointer distance_filter = ReinitializationFilterType::New();
  distance_filter->SetInput(input_threshold_filter->GetOutput());
  distance_filter->NarrowBandingOff();
  distance_filter->SetLevelSetValue(0.5f);
  distance_filter->Update();

  std::cout << "Reinitializing distance transform" << std::endl;
  DensePengReinitializationFilterType::Pointer reinit = DensePengReinitializationFilterType::New();
  // reinit->SetLevelSetValue(0.0f);
  // reinit->NarrowBandingOff();
  // reinit->SetInput(add_filter->GetOutput());
  reinit->SetInput(distance_filter->GetOutput());
  reinit->Update();

  std::cout << "Computing curvature" << std::endl;
  CurvatureImageFilterType::Pointer curvature_filter = CurvatureImageFilterType::New();
  curvature_filter->SetInput(reinit->GetOutput());
  curvature_filter->Update();

  std::cout << "Writing output to " << output_file_name << std::endl;
  DistanceWriterType::Pointer writer = DistanceWriterType::New();
  writer->SetFileName(output_file_name);
  writer->SetInput(curvature_filter->GetOutput());
  writer->Update();

  if (!dt_file_name.empty()) {
    std::cout << "Writing phi " << dt_file_name << std::endl;
    writer->SetFileName(dt_file_name);
    writer->SetInput(reinit->GetOutput());
    writer->Update();
  }

  return EXIT_SUCCESS;
}