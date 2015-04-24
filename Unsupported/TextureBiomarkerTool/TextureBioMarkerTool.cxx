/*=========================================================================
 *
    DMDData::OrientedImageType::IndexType                   index;                                                     // image index of t2 calc
 * Program :   Insight Segmentation & Registration Toolkit
 * Module  :   DMD MRI Biomarker Tool $ * Language:   C++, ITK, CMake * Date    :   Date: 2010-05-28 18:31:00 $ * Version :   $Revision: 0.20 $ * Authors :   Jiahui Wang, Martin Styner $
 * Update  :   1. Read in and write out images (05-25-2010)
 *             2. Add the series information (not successful) (05-28-2010)
 *             3. Smooth the images by gaussian filter (06-01-10)
 *             4. Assign output images by sequnce description (06-03-2010)
 *             5. Added gaussian smooth for t2 fitting (06-07-2010)
 *             6. Read in gipl images (06-09-2010)
 *             7. Read in manually delineated contiguous fat regions (06-16-2010)
 *             8. Calculate fat only image and fat percentage image (06-22-2010)
 *             9. Created feature calculation class (06-25-2010)
 *            10. Modified bug of image calibration (06-29-2010)
 *            11. Created calibration class (06-29-2010)
 *            12. Added smooth to t2 and t2FS images (06-30-2010)
 *            13. Initialize the output image with image properties of original images
 *            14. Added grayscale erosion to muscle mask with a kernel size of 5x5x5 (07-06-2010)
 *            15. Fixed a bug: correctly convert t2 and t2fs images to nrrd format (07-20-2010)
 *            16. Separate the dicom image conversion and feature calculation (07-21-2010)
 *            17. Temporally remove erosion (07-28-2010)
 *            18. Use bladder to calibrate t2 and t2fs (08-11-2010)
 *            19. Add interpolation of manual segmentation (08-24-2010)
 *            20. When do the interpolation, remove the first few and last few slices (10-20-2010)
 *            21. Call registration module of 3D Slicer (11-10-2010)
 *            22. Added curve fitting to create t2 value images (11-24-2010)
 *            23. Added registration between t2 & t2FS and between t2 & t2Fit (12-03-2010)
 *            24. Added N4 correction (01-14-2011)
 *            25. Arranged the main program (01-20-2011)
 *            26. Added water tube segmentation (02-18-2011)
 *            27. Added correction for T2FS (10-18-2011)
 *            28. Save a copy for the biomarkers. This is avaliable for the biomarkers based on calibrated images (10-21-2011)
 *            29. EPSILON was changed to 0.1 from 1 (11-21-2011)
 *            30. 5 echo full leg T2Fit; Change DMDBiomarker.h ECHO_TIME=20, DATA_PAIRS=5, DMDCurveFit.h MAXPAIRS=4  (12-07-2011)
 *            31. Output VOI  (08-01-2012)
 *            32. Calculate texture features using T2FS  (08-26-2012)
 *            33. Create some synthetic images (09-04-2012)
 *            34. Create T2 value map computation function (02-18-2013)
 * Copyright (c) Neuro Image Research and Analysis Lab.@UNC All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.
 *
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "TextureBioMarkerTool.h"
#include "itkEllipseSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include <iostream>
#include <fstream>

typedef itk::Image< float, 3>  synImageType;
typedef itk::ImageFileWriter< synImageType >                  ImageWriterType;
ImageWriterType::Pointer                    syntheticWriter = ImageWriterType::New();
static void CreateEllipseImage(synImageType::Pointer image);
//static void T2ValueMapCreate(DMDData data, DMDData::OrientedImageType::Pointer t2Calc, int t2mm_slice,DMDData::LfITKType data_array[2][DATA_PAIRS-1], DMDData::StrITKType dataDirec, DMDData::StrITKType caseID, DMDCurveFit curveFit);


inline int ipExistsArgument(char **argv, char *keystr) {
  for (int i = 1; argv[i]; i++)
    if (strstr(argv[i], keystr))
        return 1;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{
    // verify the number of parameters in the command line and print out instructions
    if (argc <= 1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help") || ipExistsArgument(argv, "-h")) {
        std::cerr << "Image Texture Evaluatation Tool 0.1 version (Dec 2013)" << std::endl;
        std::cerr << "USAGE:  " << argv[0] << " [OPTIONS] [INPUT FILE] [ROI MASK] [OUTPUT DIRECTORY]" << "\n\n" << std::endl;
        std::cerr << "EXAMPLE: \n" << argv[0] << " -c example.nrrd example_mask.nrrd /home/tester/work/\n\n" << std::endl;
        std::cerr << "         \n" << argv[0] << " -ch example.nrrd example_mask.nrrd /home/tester/work/ histogram_match_reference_image.nrrd\n\n" << std::endl;
        std::cerr << " -c " << "  calculate the texture features: (i) histogram, (ii) co-occurence, (iii) run length" << std::endl;
        std::cerr << " -r " << "  specify the range of intensity" << std::endl;
        std::cerr << " -h " << "  conduct histogram match to a reference image" << std::endl;
        std::cerr << " -help " << " print out this instruction" << std::endl;
        return EXIT_FAILURE;
    }
    // declare variables
    DMDData::StrITKType                                     dataDirec = argv[2], caseID, inputDirectory, outputFilename, fileExtension, lastCh, commandLine, InputFilename ;
    DMDData::OrientedImageReaderType::Pointer               fatMaskLeft, segMuscleLeft, imagereader, imagereader1, bladderMaskLeft; // reader for muscle mask
    DMDData::segMuscleReaderType::Pointer                   segmusclereader = DMDData::segMuscleReaderType::New();
   // std::vector<DMDData::StrITKType>                        outputFilenameList;   // store the names of files that already processed
    DMDData::PixelType                                      meanT2Left = 0, meanT2FSLeft = 0, meanT2LeftBladder = 0, meanT2FSLeftBladder = 0, fk = 0, fb = 0, integrate = 0, maxMuscleIntT2 = 0, maxMuscleIntT2FS = 0;
    DMDData::IntITKType                                     volLeftFat = 0, volLeftBladder = 0, t2mm_slice = 35;// volLeftFat: volume of subcutaneous fat region; volLeftBladder: volume of bladder region; t2mm_slice: number of slices in each echo of t2mm;
    DMDData::LfITKType                                      data_array[2][DATA_PAIRS-1] = {{0}};
    DMDData::FITKType                                       voxelVol = 0;
    DMDData::OrientedImageType::Pointer                     t2LeftFemur = DMDData::OrientedImageType::New(),
                                                            t2FSLeftFemur = DMDData::OrientedImageType::New(),
                                                            t2LeftFemurSmooth,
                                                            t2FSLeftFemurSmooth,
                                                            fatPercImgLeftFemur = DMDData::OrientedImageType::New(),
                                                            waterPercImgLeftFemur = DMDData::OrientedImageType::New(),
                                                            fatOnlyImgLeftFemur = DMDData::OrientedImageType::New(),
                                                            fatOnlyImgLeftFemurCalib = DMDData::OrientedImageType::New(),
                                                            t2LeftFemurCalib = DMDData::OrientedImageType::New(),
                                                            t2FSLeftFemurCalib = DMDData::OrientedImageType::New(),
                                                            t2FitWaterPerc = DMDData::OrientedImageType::New(),
                                                            t2FitFatPerc = DMDData::OrientedImageType::New(),
                                                            erodeSegMuscleLeft = DMDData::OrientedImageType::New(),
                                                            t2Fit = DMDData::OrientedImageType::New(),
                                                            t2FitSegMuscle = DMDData::OrientedImageType::New(),
                                                            t2Seg = DMDData::OrientedImageType::New(),
                                                            t2Calc = DMDData::OrientedImageType::New(),
                                                            t2CalcSmooth = DMDData::OrientedImageType::New();
    DMDData                                                 t2Data;
 //   DMDCalib                                                t2calib;
    TextureFeatureCal                                       muscleFeatT2, muscleFeatT2FS, muscleFeatFatOnly, muscleFeatFatPerc, muscleFeatT2Fit ;
  //  DMDCurveFit                                             curveFit;
    // end of declaration

    DMDData                                                 data( 3, inputDirectory );

    if ( strchr (argv[1], 'c') != NULL ) { // conduct feature analysis
        DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
        std::stringstream pds, dds, his;
        std::string procDirectory, dataDirectory, hisDirectory;
        std::string maskDirectory(argv[3]);
   //     float                                               mint2val = 0, maxt2val = 0;

        writer->UseCompressionOn();   // turn output image compression on

        dds << argv[2];
        pds << argv[4];
        if ( strchr (argv[1], 'h') != NULL ) { // conduct histogram match
            his << argv[5] ;
        }
        dds >> dataDirectory;
        pds >> procDirectory;
        his >> hisDirectory;

        // check whether the working directory path finished by '/'. If not, add it to the end of the path.
        unsigned found = procDirectory.find_last_of("/\\");
        if (found != procDirectory.length() - 1)
            procDirectory += "/";

        // load segmented muscle regions
        InputFilename = maskDirectory;//dataDirec + "fatReg/" + caseID + "_Seg_5slides_11m.nrrd" ;
        data.dataReader ( InputFilename, segMuscleLeft );

        //make working directory
        commandLine = "mkdir " + procDirectory;
        system(commandLine.c_str());
        commandLine = "rm -r " + procDirectory + "*.*";
  //      std::cout << commandLine << std::endl;
        system(commandLine.c_str());
        std::cout << "0%    20%   40%   60%   80%   100%" << std::endl;
        std::cout << "|-----|-----|-----|-----|-----|" << std::endl;
        std::cout << "*" << std::flush;

        // intensity rescale
        commandLine = "ImageMath " + dataDirectory  + " -rescale 0,2000 " + " -outfile " + procDirectory + "Rescaled.nrrd";
        system(commandLine.c_str());

        // N4 intensity correction
        commandLine = "N4ITKBiasFieldCorrection " + procDirectory + "Rescaled.nrrd" + " " + procDirectory + "N4Correct.nrrd" + " --outputbiasfield " + procDirectory + "BiasField.nii";
       // std::cout << commandLine << std::endl;
      //  system(commandLine.c_str());
    //    std::cout << "N4 inhomogeneous intensity correction finished!" << std::endl;
        // end of N4 correction

        // load intensity corrected images
        // load intensity corrected t2
        //InputFilename = procDirectory + "N4Correct.nrrd" ;

        if ( strchr (argv[1], 'h') != NULL ) { // conduct histogram match
            // histogram matching
            commandLine = "ImageMath " + procDirectory + "Rescaled.nrrd"  + " -matchHistogram " + hisDirectory + " -outfile " + procDirectory + "histogramMatched.nrrd";
            system(commandLine.c_str());
            InputFilename = procDirectory + "histogramMatched.nrrd" ;
        }
        else
            InputFilename = dataDirectory ;
        std::cout << "******" << std::flush;
        data.dataReader ( InputFilename, imagereader );
        t2LeftFemur = imagereader->GetOutput();
     //   std::cout << "all the images were loaded" << std::endl;
        DMDData::OrientedImageType::SpacingType                 spacing = t2LeftFemur->GetSpacing() ;
        voxelVol = spacing[0] * spacing[1] * spacing[2] ;
        // erode muscle mask by 1 voxel -- morphological erosion
        //data.morphMultiGrayErod2DIn3D ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, 0 ) ;
        // smooth t2 weighted left femur image and write it out
        data.smoothGradAnisoDiff( t2LeftFemur, t2LeftFemurSmooth, 5, 0.0625, 1.00 );
        writer->SetInput( t2LeftFemurSmooth );
        outputFilename = procDirectory + "Smooth.nrrd";
        data.dataWriter(writer, outputFilename);
        std::cout << "******" << std::flush;

        #ifdef RUN_LENGTH_MATRIX
        // run length matrix features
        std::string runlengthfeaturefilename = procDirectory + "runlength_features.txt";
        std::ofstream runlengthefile( runlengthfeaturefilename.c_str() , std::ios::app );
    //    runlengthefile << caseID << ":  " << "\n";
        runlengthefile.close();
        if ( strchr (argv[1], 'r') != NULL ) { // specify intesnity range
            muscleFeatT2.runlengthFeat3DROI ( segMuscleLeft->GetOutput(), t2LeftFemurSmooth, runlengthfeaturefilename, voxelVol, atoi(argv[argc - 2]), atoi(argv[argc - 1]) );  // use T2w
        }
        else {
            muscleFeatT2.runlengthFeat3DROI ( segMuscleLeft->GetOutput(), t2LeftFemurSmooth, runlengthfeaturefilename, voxelVol, 0, 0);  // use T2w
        }
        #endif
        std::cout << "******" << std::flush;

        #ifdef COOCURRENCE_MATRIX
        // co-occurrence matrix features
        std::string cooccurrencefeaturefilename = procDirectory + "cooccurrence_features.txt";
        std::ofstream cooccurrenceefile( cooccurrencefeaturefilename.c_str() , std::ios::app );
     //   cooccurrenceefile << caseID << ":  " << "\n";
        cooccurrenceefile.close();
        if ( strchr (argv[1], 'r') != NULL )  // specify intesnity range
            muscleFeatT2.cooccurrenceFeat3DROI ( segMuscleLeft->GetOutput(), t2LeftFemurSmooth, cooccurrencefeaturefilename, voxelVol, atoi(argv[argc - 2]), atoi(argv[argc - 1]));
        else
            muscleFeatT2.cooccurrenceFeat3DROI ( segMuscleLeft->GetOutput(), t2LeftFemurSmooth, cooccurrencefeaturefilename, voxelVol, 0, 0);
        #endif
        std::cout << "******" << std::flush;

        #ifdef HISTOGRAM_TEXTURE
        // calculate muscle features in t2 images
        std::string histogramfeaturefilename = procDirectory + "histogram_features.txt";
        std::ofstream efile( histogramfeaturefilename.c_str() , std::ios::app );
        efile.close();
        if ( strchr (argv[1], 'r') != NULL )  // specify intesnity range
            muscleFeatT2.histogramFeat3DROI ( segMuscleLeft->GetOutput(), segMuscleLeft->GetOutput(), t2LeftFemurSmooth, voxelVol, histogramfeaturefilename, TRUE, atoi(argv[argc - 2]), atoi(argv[argc - 1]));
        else
            muscleFeatT2.histogramFeat3DROI ( segMuscleLeft->GetOutput(), segMuscleLeft->GetOutput(), t2LeftFemurSmooth, voxelVol, histogramfeaturefilename, TRUE, 0, 0);

        #endif
        std::cout << "******" << std::endl;
    }
    return EXIT_SUCCESS;
}

