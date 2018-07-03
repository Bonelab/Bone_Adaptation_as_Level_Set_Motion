# B.A. Besler 2018
# Bone Imaging Laboratory

import vtk
import os

# Read inputs
if len(os.sys.argv) < 2:
    os.sys.exit('[Usage]: visualize image')
image1_file_name = os.sys.argv[1]

def setup_vis(file_name):
    print('Reading in {}'.format(file_name))

    # Reader
    vtk.vtkImageReader2Factory.RegisterReader(vtk.vtkNIFTIImageReader())
    reader = vtk.vtkImageReader2Factory.CreateImageReader2(file_name)
    if reader is None:
        os.sys.exit('Cannot find reader for {}'.format(image_file_name))
    reader.SetFileName(file_name)
    reader.Update()
    extent = list(reader.GetOutput().GetExtent())
    print(extent)
    for i in range(0, len(extent), 2):
        extent[i] = extent[i] - 1
        extent[i+1] = extent[i+1] + 1
    print(extent)

    # Pad
    padder = vtk.vtkImageConstantPad()
    padder.SetInputConnection(reader.GetOutputPort())
    padder.SetConstant(1)
    padder.SetOutputWholeExtent(extent)

    # Grab surface
    iso = vtk.vtkImageMarchingCubes()
    iso.SetInputConnection(padder.GetOutputPort())
    iso.SetValue(0,0)
    iso.ComputeNormalsOn()
    iso.ComputeGradientsOn()
    iso.ComputeScalarsOn()

    #iso = vtk.vtkFlyingEdges3D()
    #iso.SetInputConnection(padder.GetOutputPort())
    #iso.SetValue(0,0)
    #iso.ComputeNormalsOn()
    #iso.ComputeGradientsOn()
    #iso.ComputeScalarsOn()
    #iso.InterpolateAttributesOff()

    # Smoother
    #smoother = vtk.vtkWindowedSincPolyDataFilter()
    #smoother.SetInputConnection(iso.GetOutputPort())
    #smoother.SetNumberOfIterations(15)
    #smoother.BoundarySmoothingOff()
    #smoother.FeatureEdgeSmoothingOff()
    #smoother.SetFeatureAngle(120.0)
    #smoother.SetPassBand(0.001)
    #smoother.NonManifoldSmoothingOn()
    #smoother.NormalizeCoordinatesOn()

    # Mapper
    isoMapper = vtk.vtkPolyDataMapper()
    #isoMapper.SetInputConnection(smoother.GetOutputPort())
    isoMapper.SetInputConnection(iso.GetOutputPort())
    isoMapper.ScalarVisibilityOff()
    isoActor = vtk.vtkActor()
    isoActor.SetMapper(isoMapper)
    isoActor.GetProperty().SetColor(1,1,1)
    isoActor.GetProperty().SetOpacity(1)

    return isoActor

# Get reader
actor1 = setup_vis(image1_file_name)

# Visualize
renderWindow = vtk.vtkRenderWindow()

renderer = vtk.vtkRenderer()
renderer.AddActor(actor1)
renderer.SetBackground(0, 0, 0)
renderWindow.AddRenderer(renderer)
renderer.ResetCamera()
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
renderWindowInteractor.SetRenderWindow(renderWindow)

renderWindow.Render()
renderWindowInteractor.Start()

