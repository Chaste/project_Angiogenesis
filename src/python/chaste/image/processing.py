import itk

def default_image_type():
    return itk.Image[itk.UC, 3]

def image_type_2d():
    return itk.Image[itk.UC, 2]

def read_image(path, image_type):
    
    """
    Read the image from file
    """
    
    reader = itk.ImageFileReader[image_type].New()
    reader.SetFileName(path)
    reader.Update()
    return reader

def write_image(path, image_type, image_container):
    
    """
    Write the image to file
    """
    
    writer = itk.ImageFileWriter[image_type].New()
    writer.SetInput(image_container.GetOutput())
    writer.SetFileName(path)
    writer.Update()
    
def convert_vtk(image_container, image_type):
    
    itk_vtk_converter = itk.ImageToVTKImageFilter[image_type].New()
    itk_vtk_converter.SetInput(image_container.GetOutput())
    itk_vtk_converter.Update()
    return itk_vtk_converter

def merge_tiles(image_containers, image_type, num_x):
    
    tiler = itk.TileImageFilter[image_type,image_type].New()
    layout = [num_x, num_x, 0]
    tiler.SetLayout(layout)
    for idx, eachContainer in enumerate(image_containers):
        tiler.SetInput(idx, eachContainer.GetOutput())
        
    tiler.Update()
    return tiler

def median_filter(image_container, image_type, radius):
    
    median_filter = itk.MedianImageFilter[image_type, image_type].New()
    median_filter.SetInput(image_container.GetOutput())
    median_filter.SetRadius(radius)
    median_filter.Update()
    
    return median_filter
    
def correct_spacing(image_container, image_type, z_spacing):
    
    correct = itk.ChangeInformationImageFilter[image_type].New()
    correct.SetInput(image_container.GetOutput())
    spacing = image_container.GetOutput().GetSpacing()
    spacing[2] = z_spacing
    correct.SetOutputSpacing(spacing)
    correct.ChangeSpacingOn()
    correct.Update()
    
    return correct

def sum_projection(image_container, image_type):
    output_type = itk.Image[itk.UC, 2]
    
    sum_projection = itk.SumProjectionImageFilter[image_type, output_type].New()
    sum_projection.SetInput(image_container.GetOutput())
    
    sum_projection.Update()
    return sum_projection

def shrink_image(image_type, image_container, factors = [2,2,2]):
    
    """
    Reduce the size of the image in each dimension by the specified factors
    """
    
    shrink = itk.ShrinkImageFilter[image_type, image_type].New()
    shrink.SetInput(image_container.GetOutput())
    shrink.SetShrinkFactor(0, factors[0])
    shrink.SetShrinkFactor(1, factors[1])
    shrink.SetShrinkFactor(2, factors[2])
    shrink.Update()
    return shrink

def get_roi(image_type, image_container, start = [0,0,0], size = [200, 200, 10]):
    
    """
    Select a region of interest in the image
    """
    
    roi = itk.RegionOfInterestImageFilter[image_type, image_type].New()
    roi.SetInput(image_container.GetOutput())

    region = itk.ImageRegion[3]()
    region.SetIndex(start)
    region.SetSize(size)

    roi.SetRegionOfInterest(region)
    roi.Update()
    return roi

def threshold_image(image_type, image_container, lower = 195, upper = 200):
    
    """
    Do a binary thresholding
    """
    
    threshold = itk.BinaryThresholdImageFilter[image_type,image_type].New()
    threshold.SetInput(image_container.GetOutput())
    threshold.SetLowerThreshold(lower)
    threshold.SetUpperThreshold(upper)
    threshold.Update()
    return threshold
    
def get_intensity_histogram(image_type, image_container):
    
    """
    Write an intensity histogram
    """
    
    histogram = itk.ImageToHistogramFilter[image_type].New()
    histogram.SetInput(image_container.GetOutput())
    histogram.Update()
    
    for idx in range(histogram.GetOutput().GetSize()[0]):
        print "Hist Intensity: ", idx, " Frequency: ", histogram.GetOutput().GetFrequency(idx)
    