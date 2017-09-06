# NIRAL Utilities

## What is it?

NIRAL Utilities are open-source applications developped at UNC-Chapel Hill in the Neuro Image Research and Analysis Lab (NIRAL). These utilities are C++ based command line applications that allow image analysis and processing using ITK or VTK libraries.

Specifically the following utilities are contained thus far:

3DTo4DImages - Converts 3D-vector images to 4D images
4DTo3DImages - Converts 4D images to 4D-vector images
CombineDWIs - Combines multiple Diffusion Weighted Image (DWI) files into one file.
convertITKformats - convert 3D images in all ITK formats (NRRD, NIFTI, GIPL, Meta etc)
CropTools - crops 3D and 4D images
DWI_NiftiNrrdConversion - convert DWI and DTI from/to NRRD and NIFTI, works with UNC DTI tools and FSL
ImageMath - the swiss army knife image modification
ImageNoise - Computes Signal to Noise Ratio (SNR) or Contrast to Noise Ratio (CNR)
ImageStat - compute stats on images
IntensityRescaler - rescale/normalize intensities using a prior brain tissue segmentation
ITKTransformTools - Swiss army knife tool to handle ITK transform files
MauerDistanceTransform -  Calculates the Euclidean distance transform of a binary image
MultiAtlasSeg - Segmentation with a Multi-Atlas strategy
NormalizeDWIGradientVectors - Normalizes Diffusion Weighted Image gradient vectors
PolyDataCompression - Reads and writes polydata (*.vtk or *.vtp) to save them as binary data as well as activates the compression of the file.
PolydataMerge - Merges VTK polydata files
PolydataTransform - Transforms polydata files
ReadImageHeader - Prints the header of an ITK image (equivalent of 'unu head' but for any image file format that ITK reads)
TransformDeformationField - concatenates or average deformation fields (H-fields or displacement fields)

Software in Unsupported directory is not build with the main CMakeLists.txt file. It has to be build manually.

## License

See License.txt

## More information

Find the tool on [NITRC](http://www.nitrc.org/projects/niral_utilities/)

