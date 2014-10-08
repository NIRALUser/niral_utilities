#NIRAL Utilities

##What is it?

NIRAL Utilities are open-source applications developped at UNC-Chapel Hill in the Neuro Image Research and Analysis Lab (NIRAL). These utilities are C++ based command line applications that allow image analysis and processing using ITK or VTK libraries.

Specifically the following utilities are contained thus far:
ImageMath - the swiss army knife image modification
ImageStat - compute stats on images
IntensityRescaler - rescale/normalize intensities using a prior brain tissue segmentation
convertITKformats - convert 3D images in all ITK formats (NRRD, NIFTI, GIPL, Meta etc)
DWI_NiftiNrrdConversion - convert DWI and DTI from/to NRRD and NIFTI, works with UNC DTI tools and FSL
CropTools - crops 3D and 4D images
PolydataMerge - Merges VTK polydata files
PolydataTransform - Transforms polydata files
TransformDeformationField - concatenates or average deformation fields (H-fields or displacement fields)
DTIAtlasBuilder - Creates a DTI average from multiple DTI images

Software in Unsupported directory is not build with the main CMakeLists.txt file. It has to be build manually.

##License

See License.txt

##More information

Find the tool on [NITRC](http://www.nitrc.org/projects/niral_utilities/)

