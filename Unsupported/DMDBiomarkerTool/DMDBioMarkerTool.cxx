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

#include "DMDBioMarkerTool.h"
#include "itkEllipseSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include <iostream>
#include <fstream>

typedef itk::Image< float, 3>  synImageType;
typedef itk::ImageFileWriter< synImageType >                  ImageWriterType;
ImageWriterType::Pointer                    syntheticWriter = ImageWriterType::New();
static void CreateEllipseImage(synImageType::Pointer image);
static void T2ValueMapCreate(DMDData data, DMDData::OrientedImageType::Pointer t2Calc, int t2mm_slice,DMDData::LfITKType data_array[2][DATA_PAIRS-1], DMDData::StrITKType dataDirec, DMDData::StrITKType caseID, DMDCurveFit curveFit);

/////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{
    // verify the number of parameters in the command line and print out instructions
    if ( argc < 3 || !strcmp(argv[1], "-h")) {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " [OPTION] inputFileName " << std::endl;
        std::cerr << " -c convert dicom images to .nrrd format " << std::endl;
        std::cerr << " -f calcuate features in t2w sequence\n" << std::endl;
        std::cerr << " -m number of slices in each echo of t2mm; default = 20\n" << std::endl;
        std::cerr << " -i conduct interpolation\n" << std::endl;
        std::cerr << " -p calculate features in point Dixons sequence\n" << std::endl;
        std::cerr << " -h print out this instruction\n" << std::endl;
        std::cerr << "EXAMPLE: ./DMDBioMarkerTool -f /NDRC/DuchenneDystrophy/data/canine/Data/DD_040_09-03-09/ -m 35\n" << std::endl;
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
    DMDCalib                                                t2calib;
    DMDMuscleFeature                                        muscleFeatT2, muscleFeatT2FS, muscleFeatFatOnly, muscleFeatFatPerc, muscleFeatT2Fit ;
    DMDCurveFit                                             curveFit;
    // end of declaration

    if ( argc > 4 ) {
        if ( strchr (argv[4], 'm') != NULL ) {
            t2mm_slice = atoi(argv[5]);
            std::cout << "number of slices: " << t2mm_slice << std::endl;
        }
    }

    #ifdef POMPE
    caseID = dataDirec.substr(dataDirec.find( "OBS_" ));  // process Pompe
    #else
    caseID = dataDirec.substr(dataDirec.find( "DD_" )); // process DMD
    #endif
    if ( dataDirec[(int)dataDirec.length() - 1] == '/' ) {
        // remove "/" if "/" is exist in the end of input data directory
        caseID.erase(caseID.end() - 1, caseID.end());
    }
    else{
        // add "/" if "/" is not exist in the end of input data directory
        dataDirec = dataDirec + "/";
    }

    inputDirectory = dataDirec + "Dicom";
    DMDData                                                 data( 3, inputDirectory );
    if( strchr (argv[1], 'c') != NULL ) {
        // convert dicom images to nrrd format
        DMDData::IntITKType                                     t2Identifier = 0, t2FSIdentifier = 0, t2CalcIdentifier = 0, pointDixonIdentifier = 0;
        // create Orig under the data directory and save converted data in it
        commandLine = "mkdir " + dataDirec + "Orig";
        system(commandLine.c_str());
        // unzip the dicom images
        commandLine = "unzip " + dataDirec + "*.zip" + " -d " + dataDirec + "Dicom/" ;
        system(commandLine.c_str());
        try {
            // identify the list of DICOM series
            DMDData::NamesGeneratorType::Pointer                nameGenerator = DMDData::NamesGeneratorType::New( );
            nameGenerator->SetUseSeriesDetails( true );
            // transfer input directory to name generator
            nameGenerator->SetDirectory( inputDirectory.c_str() );
            typedef std::vector< std::string >                  SeriesIdContainer;    	
            const                                               SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs( );
            DMDData::StrITKType                                 seriesIdentifier;
	    DMDData::OrientedSeriesReaderType::Pointer          reader;
	    DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
            DMDData::ImageIOType::Pointer                       dicomIO;
            // turn output image compression on
            writer->UseCompressionOn();
            for (SeriesIdContainer::const_iterator it = seriesUID.begin() ; it != seriesUID.end(); it++ ){
                // We instantiate the iterators that will make possible to walk through all the entries of the MetaDataDictionary.
	        std::vector<DMDData::StrITKType>                tagvalue, tagkey;
	        std::vector<DMDData::StrITKType>::iterator      itr;
                // specify the location of sequence name
	        DMDData::StrITKType                             specifiedTag = "0008|103e";
                seriesIdentifier = *it;
                data.seriesDataReader( seriesIdentifier, reader, dicomIO );
                data.readDicomTag(dicomIO, tagvalue, tagkey);
	        itr = find (tagkey.begin(), tagkey.end(), specifiedTag);	
                outputFilename = tagvalue[int(itr - tagkey.begin())];	
                fileExtension = ".nrrd";
                lastCh = outputFilename[outputFilename.size() - 1];
                if ( !strcmp(lastCh.c_str(), " ")) {
                    // if the last charactor is space remove the last space of tagvalue then plus .nrrd
                    outputFilename = outputFilename.erase(outputFilename.size() - 1, 1);
                }
                while ( outputFilename.find(' ') < outputFilename.size() ){
                    // replace space charactors by "_"
	            outputFilename.replace(outputFilename.find(' '), 1, "_");
                }	
                outputFilename = dataDirec + "Orig/" + caseID + "_" + outputFilename;
                // the sequence for t2 fitting?
	        std::size_t found = outputFilename.find("calc");
                if( found != std::string::npos ){
                    std::cout << "this is the sequence for t2 fitting!"  << std::endl;
                    writer->SetInput(reader->GetOutput());
                    if (!t2CalcIdentifier) {
                        outputFilename = dataDirec + "Orig/" + caseID + "_T2calc";
                        outputFilename += fileExtension;
                        data.dataWriter(writer, outputFilename);
                        t2CalcIdentifier++;
                    }
                    else {
                        std::stringstream identifier;
                        identifier << t2CalcIdentifier;
                        outputFilename = dataDirec + "Orig/" + caseID + "_T2calc";
                        outputFilename += identifier.str() + fileExtension;
                        data.dataWriter(writer, outputFilename);
                        t2CalcIdentifier++;
                     }
                     outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
                }
                // the sequence for 3 point dixon?
	        found = outputFilename.find("dixon");
                if( found != std::string::npos ){
                    std::cout << "this is the sequence for 3 point dixon"  << std::endl;
                    writer->SetInput(reader->GetOutput());
                    if (!pointDixonIdentifier) {
                        outputFilename += fileExtension;
                        data.dataWriter(writer, outputFilename);
                        pointDixonIdentifier++;
                    }
                    else {
                        std::stringstream identifier;
                        identifier << pointDixonIdentifier;
                        outputFilename += identifier.str() + fileExtension;
                        data.dataWriter(writer, outputFilename);
                        pointDixonIdentifier++;
                     }
                     outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
                }
                found = outputFilename.find("t2_tse");
                if( found != std::string::npos ){
                    std::cout << "this is the t2 weighted sequence!"  << std::endl;
                    writer->SetInput(reader->GetOutput());
                    if (!t2Identifier) {
                        outputFilename += fileExtension;
                        data.dataWriter(writer, outputFilename);
                        t2Identifier++;
                    }
                    else {
                        std::stringstream identifier;
                        identifier << t2Identifier;
                        outputFilename += identifier.str() + fileExtension;
                        data.dataWriter(writer, outputFilename);
                        t2Identifier++;
                     }
                     outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
                }
                found = outputFilename.find("t2fs_tse");
	        if ( found != std::string::npos ) {
                    std::cout << "this is the t2 fat saturation weighted sequence!"  << std::endl;
                    writer->SetInput(reader->GetOutput());
                    if (!t2FSIdentifier) {
                        outputFilename += fileExtension;
                        data.dataWriter(writer, outputFilename);
                        t2FSIdentifier++;
                    }
                    else {
                        std::stringstream identifier;
                        identifier << t2FSIdentifier;
                        outputFilename += identifier.str() + fileExtension;
                        data.dataWriter(writer, outputFilename);
                        t2FSIdentifier++;
                     }
                     outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
                }
	        found = outputFilename.find("left_femur_T2FS");
	        if ( found != std::string::npos ){
                    std::cout << "this is the left femur aligned t2 weighted fat suppressed sequence!"  << std::endl;
                    DMDData::OrientedImageType::Pointer origImage = reader->GetOutput();
                    outputFilename += fileExtension;
                    // only convert sequence with multiple slices
                    if ( origImage->GetLargestPossibleRegion().GetSize()[2] > 1 ) {
                        writer->SetInput(reader->GetOutput());
                        data.dataWriter(writer, outputFilename);
                    }
                    outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
	        }
                else {
                    // use else, otherwise will produce xxxT2FS.nrrd.nrrd because T2 is a subset of T2FS
	            found = outputFilename.find("left_femur_T2");
	            if ( found != std::string::npos ){
                        std::cout << "this is the left aligned t2 weighted sequence!"  << std::endl;
                        DMDData::OrientedImageType::Pointer origImage = reader->GetOutput();
                        outputFilename += fileExtension;
                        // only convert sequence with multiple slices
                        if ( origImage->GetLargestPossibleRegion().GetSize()[2] > 1 ) {
                            writer->SetInput(reader->GetOutput());
                            data.dataWriter(writer, outputFilename);
                        }
                    //    outputFilename.erase (outputFilename.end() - 5, outputFilename.end());
                        outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
	            }	
                    else {
                        DMDData::OrientedImageType::Pointer origImage = reader->GetOutput();
                        outputFilename += fileExtension;
                        // only convert sequence with multiple slices
                        if ( origImage->GetLargestPossibleRegion().GetSize()[2] > 1 ) {
                            writer->SetInput(reader->GetOutput());
                            data.dataWriter(writer, outputFilename);
                        }
                        outputFilename = outputFilename.erase(outputFilename.size() - 6, 6);
                    }
                }
            }
            commandLine = "rm -r " + dataDirec + "Dicom/" ;
     //       system(commandLine.c_str());
            return EXIT_SUCCESS;
        }
        catch (itk::ExceptionObject &ex){
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
        }
    }

    if ( strchr (argv[1], 'i') != NULL ) {
        // conduct feature analysis
        // load original images continuous fat region and segmented muscle
	DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
        float                                               mint2val = 0, maxt2val = 0;
        writer->UseCompressionOn();   // turn output image compression on
        // load segmented muscle regions
        // conduct interpolation between manually segmented slices
        //InputFilename = dataDirec + "Seg/" + caseID + "_manual_full.nrrd" ;
        InputFilename = dataDirec + "Seg/" + caseID + "_right_femur_seg_every5.nrrd" ;
        data.dataReader ( InputFilename, segMuscleLeft );
        data.interpolate3D(segMuscleLeft->GetOutput(), caseID);
        std::cout << "finish interpolation" << std::endl;
        return EXIT_SUCCESS;
    }

    if ( strchr (argv[1], 'f') != NULL ) { // conduct feature analysis
	DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
        float                                               mint2val = 0, maxt2val = 0;
        writer->UseCompressionOn();   // turn output image compression on
        // load fat mask of left femur
//        InputFilename = dataDirec + "fatReg/" + caseID + "_left_femur_fat.nrrd" ;
        InputFilename = argv[3];
        data.dataReader ( InputFilename, fatMaskLeft );

        // load bladder mask of left femur to correct intensity inconsistance between t2 and t2fs
 //       InputFilename = dataDirec + "fatReg/" + caseID + "_left_femur_fat.nrrd" ;
        InputFilename = argv[3];
        data.dataReader ( InputFilename, bladderMaskLeft );

        // load segmented muscle regions
        #ifdef SEG_INTERPOLATE
        // conduct interpolation between manually segmented slices
   //    InputFilename = dataDirec + "fatReg/" + caseID + "_intervolfull_every5.nrrd" ; // generally used segmentation
      // InputFilename = dataDirec + "fatReg/" + caseID + "_intervolfull_every5_distal.nrrd" ; // generally used segmentation
     //  InputFilename = dataDirec + "fatReg/" + caseID + "_intervolfull_every5_medial.nrrd" ; // generally used segmentation
    //   InputFilename = dataDirec + "fatReg/" + caseID + "_intervolfull_every5_proximal.nrrd" ; // generally used segmentation
   //     InputFilename = dataDirec + "fatReg/" + caseID + "_Seg_5slides_11m.nrrd" ; // 5 central slices
   //     InputFilename = dataDirec + "fatReg/" + caseID + "_left_femur_seg.nrrd" ; // DD_032, DD_033, DD_034, CS, RF, ST,
        InputFilename = argv[3];
        data.dataReader ( InputFilename, segMuscleLeft );
    //    data.interpolate3D(segMuscleLeft->GetOutput(), caseID);
        #endif

        // N4 intensity correction
        commandLine = "N4ITKBiasFieldCorrection --inputimage " + dataDirec + "Orig/" + caseID + "_right_femur_T2.nrrd " + "--outputimage " + dataDirec + "Processed/" + caseID + "_left_femur_T2_N4Correct.nrrd" + " --outputbiasfield " + dataDirec + "Processed/" + caseID + "_T2BiasField.nii";
        system(commandLine.c_str());
        // apply the bias field to t2fs
        //commandLine = "ImageMath " + dataDirec + "Orig/" + caseID + "_left_femur_T2FS.nrrd " + "-div " + dataDirec + "Processed/" + caseID + "_T2BiasField.nii" + " -outfile " + dataDirec + "Processed/" + caseID + "_left_femur_T2FS_N4Correct.nrrd";
        commandLine = "N4ITKBiasFieldCorrection --inputimage " + dataDirec + "Orig/" + caseID + "_right_femur_T2FS.nrrd " + "--outputimage " + dataDirec + "Processed/" + caseID + "_left_femur_T2FS_N4Correct.nrrd";
       // std::cout << commandLine << std::endl;
        system(commandLine.c_str());
        commandLine = "N4ITKBiasFieldCorrection --inputimage " + dataDirec + "Orig/" + caseID + "_T2calc.nrrd " "--outputimage " + dataDirec + "Processed/" + caseID + "_T2calc_N4Correct.nrrd";
        system(commandLine.c_str());
        std::cout << "N4 inhomogeneous intensity correction finished!" << std::endl;
        // end of N4 correction

        // load intensity corrected images
        // load intensity corrected t2
        InputFilename = dataDirec + "Processed/" + caseID + "_left_femur_T2_N4Correct.nrrd" ;
        data.dataReader ( InputFilename, imagereader );
        t2LeftFemur = imagereader->GetOutput();
        // load intensity corrected t2fs
        InputFilename = dataDirec + "Processed/" + caseID + "_left_femur_T2FS_N4Correct.nrrd" ;
        data.dataReader ( InputFilename, imagereader );
        t2FSLeftFemur = imagereader->GetOutput();
        // load intensity corrected t2calc
        InputFilename = dataDirec + "Processed/" + caseID + "_T2calc_N4Correct.nrrd" ;
        InputFilename = dataDirec + "Orig/" + caseID + "_T2calc.nrrd" ;
        data.dataReader ( InputFilename, imagereader );
        t2Calc = imagereader->GetOutput();
        std::cout << "all the images were loaded" << std::endl;

        // initialize image volumes for fat percentage, water percentage, fat only, intensity calibrated t2 t2fs, and eroded muscle mask
        data.imageInitialize ( t2LeftFemur, fatPercImgLeftFemur, voxelVol );
        data.imageInitialize ( t2LeftFemur, waterPercImgLeftFemur, voxelVol );
        data.imageInitialize ( t2LeftFemur, fatOnlyImgLeftFemur );
        data.imageInitialize ( t2LeftFemur, fatOnlyImgLeftFemurCalib );
        data.imageInitialize ( t2LeftFemur, t2LeftFemurCalib );	
        data.imageInitialize ( t2LeftFemur, t2FSLeftFemurCalib );
        data.imageInitialize ( t2LeftFemur, erodeSegMuscleLeft );
        data.imageInitialize ( t2LeftFemur, t2FitWaterPerc );	
        data.imageInitialize ( t2LeftFemur, t2FitFatPerc );	
        std::cout << "image volumes were initialized" << std::endl;
        // erode muscle mask by 1 voxel -- morphological erosion

        data.morphMultiGrayErod2DIn3D ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, 0 ) ;
        // smooth t2 weighted left femur image and write it out
        data.smoothGradAnisoDiff( t2LeftFemur, t2LeftFemurSmooth, 5, 0.0625, 1.00 );
        writer->SetInput( t2LeftFemurSmooth );
        outputFilename = dataDirec + "Processed/" + caseID + "_left_femur_T2_Smooth.nrrd";
        data.dataWriter(writer, outputFilename);
        // smooth t2fs images
        t2Data.smoothGradAnisoDiff( t2FSLeftFemur, t2FSLeftFemurSmooth, 5, 0.0625, 1.00 );
        writer->SetInput( t2FSLeftFemurSmooth );
        outputFilename = dataDirec + "Processed/" + caseID + "_left_femur_T2FS_Smooth.nrrd";
        data.dataWriter(writer, outputFilename);

        #ifdef RUN_LENGTH_MATRIX
        // run length matrix features
        std::string runlengthfeaturefilename = "../data/" + caseID + "_runlength_features7m.txt";
        std::ofstream runlengthefile( runlengthfeaturefilename.c_str() , std::ios::app );
        runlengthefile << caseID << ":  " << "\n";
        runlengthefile.close();
        //muscleFeatT2.runlengthFeat3DROI ( erodeSegMuscleLeft, t2LeftFemurSmooth, runlengthfeaturefilename, voxelVol, caseID );  // use T2w
        muscleFeatT2.runlengthFeat3DROI ( erodeSegMuscleLeft, t2FSLeftFemurSmooth, runlengthfeaturefilename, voxelVol, caseID );  // use T2fs
        return EXIT_SUCCESS;
        #endif

        #ifdef COOCURRENCE_MATRIX
        // co-occurrence matrix features
        std::string cooccurrencefeaturefilename = "../data/" + caseID + "_cooccurrence_features7m.txt";
        std::ofstream cooccurrenceefile( cooccurrencefeaturefilename.c_str() , std::ios::app );
        cooccurrenceefile << caseID << ":  " << "\n";
        cooccurrenceefile.close();
        muscleFeatT2.cooccurrenceFeat3DROI ( erodeSegMuscleLeft, t2LeftFemurSmooth, cooccurrencefeaturefilename, voxelVol, caseID );
        return EXIT_SUCCESS;
        #endif

        // begin create t2 maps
        T2ValueMapCreate(data, t2Calc, t2mm_slice, data_array, dataDirec, caseID, curveFit);
        // end of create t2 value image

        commandLine = "BRAINSFit --fixedVolume " + dataDirec + "Processed/" + caseID + "_left_femur_T2_Smooth.nrrd --movingVolume " + dataDirec + "Processed/" + caseID + "_left_femur_T2FS_Smooth.nrrd --outputVolume " + dataDirec + "Processed/" + caseID + "_T2FS_registered.nrrd --transformType Rigid";   // transform for t2FS
        system(commandLine.c_str());
        commandLine = "BRAINSFit --fixedVolume " + dataDirec + "Processed/" + caseID + "_left_femur_T2_Smooth.nrrd --movingVolume " + dataDirec + "Processed/" + caseID + "_T2FIT.nrrd --outputVolume " + dataDirec + "Processed/" + caseID + "_T2FIT_registered.nrrd --transformType Rigid";   // transform for t2FIT
        system(commandLine.c_str());
        commandLine = "cp " + dataDirec + "Processed/" + caseID + "_left_femur_T2_Smooth.nrrd " + dataDirec + "Processed/" + caseID + "_T2_registered.nrrd";
        system(commandLine.c_str());
        std::cout << "image registration was finished" << std::endl;

        // read tranformed t2fs and t2fit
        // load left femur fs
        InputFilename = dataDirec + "Processed/" + caseID + "_T2FS_registered.nrrd" ;
        std::cout << InputFilename << std::endl;
        data.dataReader ( InputFilename, imagereader );
        t2FSLeftFemur = imagereader->GetOutput();
        // load t2 map
        InputFilename = dataDirec + "Processed/" + caseID + "_T2FIT_registered.nrrd" ;
        data.dataReader ( InputFilename, imagereader );
        t2Fit = imagereader->GetOutput();
        // load t2
        InputFilename = dataDirec + "Processed/" + caseID + "_T2_registered.nrrd" ;
        data.dataReader ( InputFilename, imagereader );
        t2LeftFemur = imagereader->GetOutput();

      //  data.waterTubeSegmentation(t2FSLeftFemur);
     //   std::cout << "thresholding" << std::endl;
        // create iterator of t2 weighted left femur image
       // DMDData::OrientedImageType::Pointer origImage = t2LeftFemurSmooth;
        DMDData::ConstIteratorType constOrigT2Iterator( t2LeftFemurSmooth, t2LeftFemurSmooth->GetRequestedRegion() );		
        // create iterator of subcutaneous fat
        DMDData::OrientedImageType::Pointer leftMaskImage = fatMaskLeft->GetOutput();
        DMDData::ConstIteratorType constLeftMaskIterator( leftMaskImage, leftMaskImage->GetRequestedRegion() );
        // create iterator of bladder region
        DMDData::OrientedImageType::Pointer leftBladderMaskImage = bladderMaskLeft->GetOutput();
        DMDData::ConstIteratorType constLeftBladderMaskIterator( leftBladderMaskImage, leftBladderMaskImage->GetRequestedRegion() );
        DMDData::ConstIteratorType constLeftSegMuscleMaskIterator( erodeSegMuscleLeft, erodeSegMuscleLeft->GetRequestedRegion());
        for ( constOrigT2Iterator.GoToBegin(), constLeftMaskIterator.GoToBegin(), constLeftBladderMaskIterator.GoToBegin(), constLeftSegMuscleMaskIterator.GoToBegin(); !constOrigT2Iterator.IsAtEnd(); ++constOrigT2Iterator, ++constLeftMaskIterator, ++constLeftBladderMaskIterator, ++constLeftSegMuscleMaskIterator ) {
	    if ( constLeftMaskIterator.Get() > 0 ) { // calculate sum of pixel value in t2 in subcutaneous fat region
	        meanT2Left += constOrigT2Iterator.Get();
                volLeftFat++;
            }
            if ( constLeftBladderMaskIterator.Get() > 0 ) { // calculate sum of pixel value in t2FS in subcutaneous fat region
                meanT2LeftBladder += constOrigT2Iterator.Get();
                volLeftBladder++;
            }   	
	    if ( constLeftSegMuscleMaskIterator.Get() > 0 ) { // calculate sum of pixel value in t2 in subcutaneous fat region
                if ( constOrigT2Iterator.Get() > maxMuscleIntT2)
                    maxMuscleIntT2 = constOrigT2Iterator.Get();
            }
        }

        // create iterator of t2 fs left femur image
       // DMDData::OrientedImageType::Pointer origImage = t2FSLeftFemurSmooth;
        DMDData::ConstIteratorType constOrigFSIterator( t2FSLeftFemur, t2FSLeftFemur->GetRequestedRegion() );		   // create iterator of subcutaneous fat
      //  DMDData::OrientedImageType::Pointer leftMaskImage = fatMaskLeft->GetOutput();
       // DMDData::ConstIteratorType constLeftMaskIterator( leftMaskImage, leftMaskImage->GetRequestedRegion() );
        // iterator of bladder region
        //DMDData::OrientedImageType::Pointer leftBladderMaskImage = bladderMaskLeft->GetOutput();
        //DMDData::ConstIteratorType constLeftBladderMaskIterator( leftBladderMaskImage, leftBladderMaskImage->GetRequestedRegion() );
        for ( constOrigFSIterator.GoToBegin(), constLeftMaskIterator.GoToBegin(), constLeftBladderMaskIterator.GoToBegin(), constLeftSegMuscleMaskIterator.GoToBegin(); !constOrigFSIterator.IsAtEnd(); ++constOrigFSIterator, ++constLeftMaskIterator, ++constLeftSegMuscleMaskIterator, ++constLeftBladderMaskIterator ) {	
            if ( constLeftMaskIterator.Get() > 0 ) { // calculate sum of pixel value in t2FS in subcutaneous fat region
	        meanT2FSLeft += constOrigFSIterator.Get();
            }
            if ( constLeftBladderMaskIterator.Get() > 0 ) { // calculate sum of pixel value in t2FS in subcutaneous fat region
                meanT2FSLeftBladder += constOrigFSIterator.Get();
            }   	
	    if ( constLeftSegMuscleMaskIterator.Get() > 0 ) { // calculate sum of pixel value in t2 in subcutaneous fat region
                if ( constOrigFSIterator.Get() > maxMuscleIntT2)
                    maxMuscleIntT2FS = constOrigFSIterator.Get();
            }
        }
        // calculate histogram of t2FS image
        meanT2Left /= volLeftFat;
        meanT2FSLeft /= volLeftFat;
        meanT2LeftBladder /= volLeftBladder;
        meanT2FSLeftBladder /= volLeftBladder;
        // begin calibration
        integrate = meanT2LeftBladder - meanT2FSLeftBladder;
        std::cout << "difference between t2 and t2fs:  " << integrate << "  " <<  meanT2LeftBladder << "  " << meanT2FSLeftBladder << std::endl;
        //getchar();
        fb = ( BRIGHT_CALIB * ( meanT2Left - meanT2FSLeft ) - ( BRIGHT_CALIB - DARK_CALIB ) * meanT2Left ) / ( meanT2Left - meanT2FSLeft );
        fk = ( BRIGHT_CALIB - DARK_CALIB ) / ( meanT2Left - meanT2FSLeft );
        //fb += 75.732;
        std::cout << "the size of region is:  " << volLeftFat << "\n mean of the region in T2 image is:  " << meanT2Left << "\n mean of the region in T2FS image is:  " << meanT2FSLeft << "  " << fk << "  " << fb << std::endl;
        // calibrate T2 images
        if (strcmp(caseID.c_str(), "DD_082_04-30-10") == 0 || strcmp(caseID.c_str(), "DD_083_04-30-10") == 0 || strcmp(caseID.c_str(), "DD_084_04-30-10") == 0) {
            t2calib.linearCalib ( t2LeftFemurSmooth, t2LeftFemurCalib, fatMaskLeft, meanT2Left, meanT2FSLeft, fk, fb, integrate ) ;
            std::cout << "calibrated!!! " << std::endl;
        }
        else {
            t2calib.linearCalib ( t2LeftFemurSmooth, t2LeftFemurCalib, fatMaskLeft, meanT2Left, meanT2FSLeft, fk, fb, 0 ) ;
            std::cout << "don't need calibration!!! " << std::endl;
        }
        writer->SetInput(  t2LeftFemurCalib );
        outputFilename = dataDirec + "Processed/" + caseID + "_left_femur_T2_calib.nrrd";
        data.dataWriter(writer, outputFilename);		
        // calibrate T2FS images
        // DMDData::OrientedImageType::Pointer leftMaskImage = fatMaskLeft->GetOutput();
        //fb -= 75.732;
        t2calib.linearCalib ( t2FSLeftFemur, t2FSLeftFemurCalib, fatMaskLeft, meanT2Left, meanT2FSLeft, fk, fb, 0) ;
        writer->SetInput( t2FSLeftFemurCalib );
        outputFilename =  dataDirec + "Processed/" + caseID + "_left_femur_T2FS_calib.nrrd";
        data.dataWriter(writer, outputFilename);
        // end of calibration
        // iterators for fat only, water only, and fat percentage image
        DMDData::IteratorType fatPercLeftFemurIterator( fatPercImgLeftFemur, fatPercImgLeftFemur->GetRequestedRegion() );
        DMDData::IteratorType waterPercLeftFemurIterator( waterPercImgLeftFemur, waterPercImgLeftFemur->GetRequestedRegion() );
        DMDData::IteratorType fatOnlyImgLeftFemurIterator( fatOnlyImgLeftFemur, fatOnlyImgLeftFemur->GetRequestedRegion() );
        DMDData::IteratorType fatOnlyImgLeftFemurCalibIterator( fatOnlyImgLeftFemurCalib, fatOnlyImgLeftFemurCalib->GetRequestedRegion() );
        DMDData::IteratorType t2ValWaterPercIterator( t2FitWaterPerc, t2FitWaterPerc->GetRequestedRegion() );
        DMDData::IteratorType t2ValFatPercIterator( t2FitFatPerc, t2FitFatPerc->GetRequestedRegion() );
        DMDData::ConstIteratorType constt2FitIterator( t2Fit, t2Fit->GetRequestedRegion() );
        DMDData::ConstIteratorType constLeftT2Iterator( t2LeftFemur, t2LeftFemur->GetRequestedRegion() );
        DMDData::ConstIteratorType constLeftT2FSIterator( t2FSLeftFemur, t2FSLeftFemur->GetRequestedRegion() );
        DMDData::ConstIteratorType constLeftT2CalibIterator( t2LeftFemurCalib, t2LeftFemurCalib->GetRequestedRegion() );
        DMDData::ConstIteratorType constLeftT2FSCalibIterator( t2FSLeftFemurCalib, t2FSLeftFemurCalib->GetRequestedRegion() );
        DMDData::OrientedWriterType::Pointer                 fatonlywriter = DMDData::OrientedWriterType::New();
        DMDData::OrientedWriterType::Pointer                 fatpercwriter = DMDData::OrientedWriterType::New();
        DMDData::OrientedWriterType::Pointer                 waterpercwriter = DMDData::OrientedWriterType::New();
        fatonlywriter->UseCompressionOn();   // turn output image compression on
        fatpercwriter->UseCompressionOn();   // turn output image compression on
        waterpercwriter->UseCompressionOn();   // turn output image compression on
        for ( constLeftT2Iterator.GoToBegin(), constLeftT2FSIterator.GoToBegin(), constLeftT2CalibIterator.GoToBegin(), constLeftT2FSCalibIterator.GoToBegin(), constt2FitIterator.GoToBegin(), fatPercLeftFemurIterator.GoToBegin(), waterPercLeftFemurIterator.GoToBegin(), fatOnlyImgLeftFemurIterator.GoToBegin(), fatOnlyImgLeftFemurCalibIterator.GoToBegin(), t2ValWaterPercIterator, t2ValFatPercIterator; !constLeftT2Iterator.IsAtEnd(); ++constLeftT2Iterator, ++constLeftT2FSIterator, ++constLeftT2CalibIterator, ++constLeftT2FSCalibIterator, ++constt2FitIterator, ++fatPercLeftFemurIterator, ++waterPercLeftFemurIterator, ++fatOnlyImgLeftFemurIterator, ++fatOnlyImgLeftFemurCalibIterator, ++t2ValWaterPercIterator, ++t2ValFatPercIterator) {
            // calculate fat only image
            fatOnlyImgLeftFemurIterator.Set( constLeftT2Iterator.Get() - constLeftT2FSIterator.Get());
            fatOnlyImgLeftFemurCalibIterator.Set( constLeftT2CalibIterator.Get() - constLeftT2FSCalibIterator.Get());
            // calculate fat percentage image
            float fraction = constLeftT2Iterator.Get() + EPSILON;
            if (fraction == 0)
                fraction = EPSILON;
            float tmp = (constLeftT2Iterator.Get() - constLeftT2FSIterator.Get()) / fraction;
            if (tmp >= 0){
                if (tmp > 1)
                    tmp = 1;
                fatPercLeftFemurIterator.Set( tmp * 100 );
                //create t2 value fat image
                t2ValFatPercIterator.Set(tmp * constt2FitIterator.Get());
            }
            else{
                fatPercLeftFemurIterator.Set( 0 );
                t2ValFatPercIterator.Set(0);
            }
            // calculate water percentage image
            tmp = constLeftT2FSIterator.Get() / fraction;
            if (tmp >= 0){
                if (tmp > 1)
                    tmp = 1;
                waterPercLeftFemurIterator.Set( tmp * 100 );
                //create t2 value water image
                t2ValWaterPercIterator.Set(tmp * constt2FitIterator.Get());
            }
            else{
                waterPercLeftFemurIterator.Set( 0 );
                t2ValWaterPercIterator.Set(0);
            }
        }

        // output fat percentage image
        fatpercwriter->SetInput( fatPercImgLeftFemur );
        outputFilename = dataDirec + "Processed/" + caseID + "_fatPercImgLeftFemur.nrrd";
        data.dataWriter(fatpercwriter, outputFilename);
        // output fat only image
	fatonlywriter->SetInput(fatOnlyImgLeftFemur );
        outputFilename = dataDirec + "Processed/" + caseID + "_fatOnlyImgLeftFemur.nrrd";
        data.dataWriter(fatonlywriter, outputFilename);
        // output calibrated fat only image
	fatonlywriter->SetInput(fatOnlyImgLeftFemurCalib );
        outputFilename = dataDirec + "Processed/" + caseID + "_fatOnlyImgLeftFemurCalib.nrrd";
        data.dataWriter(fatonlywriter, outputFilename);
        // output water percentage  image
	waterpercwriter->SetInput(waterPercImgLeftFemur );
        outputFilename = dataDirec + "Processed/" + caseID + "_waterPercImgLeftFemur.nrrd";
        data.dataWriter(waterpercwriter, outputFilename);

        std::string featurefilename = "../data/" + caseID + "_features7m.txt";
        commandLine = "rm " + featurefilename ;
        system(commandLine.c_str()); // remove the old feature file before updating
        std::cout << featurefilename.c_str() << std::endl;
        std::ofstream efile( featurefilename.c_str() , std::ios::app );
        efile << caseID << ":  " << "\n";
        efile.close();
        #ifdef HISTOGRAM_TEXTURE
        // calculate muscle features in t2 images
        std::string histogramfeaturefilename = "../data/" + caseID + "_histogram_features7m.txt";
        std::cout << "features in t2 images:  " << std::endl;
        efile.open( histogramfeaturefilename.c_str() , std::ios::app );
        efile.close();
        //muscleFeatT2.features ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, t2LeftFemurCalib, voxelVol, featurefilename, TRUE ); // from T2w
        muscleFeatT2.features ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, t2FSLeftFemurCalib, voxelVol, histogramfeaturefilename, TRUE );  // form T2fs
        return EXIT_SUCCESS;
        #endif
        #ifdef _HISTOGRAM_TEXTURE
        // calculate muscle features in t2 images
        std::cout << "features in t2 images:  " << std::endl;
        efile.open( featurefilename.c_str() , std::ios::app );
        efile.close();
        muscleFeatT2.features ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, t2LeftFemurCalib, voxelVol, featurefilename, TRUE );
        // calculate muscle features in t2FS images
        std::cout << "features in t2FS images:  " << std::endl;
        efile.open( featurefilename.c_str() , std::ios::app );
        efile << "water_only:     ";
        efile.close();
        muscleFeatT2FS.features ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, t2FSLeftFemurCalib, voxelVol, featurefilename, FALSE );
        // calculate muscle features in calibrated fat map images
        std::cout << "features in calibrated fat only images:  " << std::endl;
        efile.open( featurefilename.c_str(), std::ios::app );
        efile << "calibrated_fat_only:           ";
        efile.close();
        muscleFeatFatOnly.features ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, fatOnlyImgLeftFemurCalib, voxelVol, featurefilename, FALSE );
        // calculate muscle features in t2 fit images
        std::cout << "features in t2 fit images:  " << std::endl;
        efile.open( featurefilename.c_str(), std::ios::app );
        efile << "t2_value:           ";
        efile.close();
        muscleFeatT2Fit.features (  segMuscleLeft->GetOutput(), erodeSegMuscleLeft, t2Fit, voxelVol, featurefilename, FALSE );
        #endif
    }

    if ( strchr (argv[1], 'p') != NULL ) {
        // conduct feature analysis on point Dixon
        // load original images continuous fat region and segmented muscle
	DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
        float                                               mint2val = 0, maxt2val = 0;
        writer->UseCompressionOn();   // turn output image compression on
        // load segmented muscle regions
        #ifdef SEG_INTERPOLATE
        // conduct interpolation between manually segmented slices
        InputFilename = dataDirec + "Mask/" + caseID + "_seg.nii.gz" ; // generally used segmentation
        std::cout << InputFilename << std::endl;
        //getchar();
        data.dataReader ( InputFilename, segMuscleLeft );
        #endif
        // load interpolated volume
        // InputFilename = dataDirec + "fatReg/" + caseID + "_intervolfull_every8.nrrd" ;
        // data.dataReader ( InputFilename, segMuscleLeft );
        // N4 intensity correction
        // begin of correction
        //commandLine = "N4ITKBiasFieldCorrection --inputimage " + dataDirec + "Orig/" + caseID + "_left_femur_T2.nrrd " + "--outputimage " + dataDirec + "Processed/" + caseID + "_left_femur_T2_N4Correct.nrrd" + " --outputbiasfield " + dataDirec + "Processed/" + caseID + "_T2BiasField.nii";
       // system(commandLine.c_str());
        //std::cout << "N4 inhomogeneous intensity correction finished!" << std::endl;
        // end of correction
        // load intensity corrected images
        // load intensity corrected t2
        InputFilename = dataDirec + "Orig/" + caseID + "_3PD.nii.gz" ;
        data.dataReader ( InputFilename, imagereader );
        t2LeftFemur = imagereader->GetOutput();
        std::cout << "all the images were loaded" << std::endl;

        // initialize image volumes for fat percentage, water percentage, fat only, intensity calibrated t2 t2fs, and eroded muscle mask
        data.imageInitialize ( t2LeftFemur, fatPercImgLeftFemur, voxelVol );
        data.imageInitialize ( t2LeftFemur, waterPercImgLeftFemur, voxelVol );
        data.imageInitialize ( t2LeftFemur, fatOnlyImgLeftFemur );
        data.imageInitialize ( t2LeftFemur, fatOnlyImgLeftFemurCalib );
        data.imageInitialize ( t2LeftFemur, t2LeftFemurCalib );	
        data.imageInitialize ( t2LeftFemur, t2FSLeftFemurCalib );
        data.imageInitialize ( t2LeftFemur, erodeSegMuscleLeft );
        data.imageInitialize ( t2LeftFemur, t2FitWaterPerc );	
        data.imageInitialize ( t2LeftFemur, t2FitFatPerc );	
        std::cout << "image volumes were initialized" << std::endl;
        // erode muscle mask by 1 voxel -- morphological erosion
        data.morphMultiGrayErod2DIn3D ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, 1 ) ;
        // smooth t2 weighted left femur image and write it out
        data.smoothGradAnisoDiff( t2LeftFemur, t2LeftFemurSmooth, 5, 0.0625, 1.00 );
        writer->SetInput( t2LeftFemurSmooth );

        DMDData::OrientedImageType::IndexType                   index;                                                     // image index of t2 calc
        DMDData::OrientedImageType::SizeType                    t2Size = t2LeftFemur->GetLargestPossibleRegion().GetSize(); // get the image information

        for(index[2] = 0; index[2] < (int)t2Size[2]; index[2]++) {
            for(index[1] = 0; index[1] < (int)t2Size[1]; index[1]++) {
                for(index[0] = 0; index[0] < (int)t2Size[0]; index[0]++) {
                    t2LeftFemur->SetPixel(index, t2LeftFemur->GetPixel(index) * 1024);
                }
            }
        }

        #ifdef RUN_LENGTH_MATRIX
       // run length matrix features
        std::string runlengthfeaturefilename = "../data/" + caseID + "_runlength_features7m.txt";
        std::ofstream runlengthefile( runlengthfeaturefilename.c_str() , std::ios::app );
        runlengthefile << caseID << ":  " << "\n";
        runlengthefile.close();
        muscleFeatT2.runlengthFeat3DROI ( erodeSegMuscleLeft, t2LeftFemur, runlengthfeaturefilename, voxelVol, caseID );
        //muscleFeatT2.runlengthFeat2DROI ( erodeSegMuscleLeft, t2LeftFemurSmooth, runlengthfeaturefilename, voxelVol );
        //exit(1);
        #endif

        std::string featurefilename = "../data/" + caseID + "_features7m.txt";
        commandLine = "rm " + featurefilename ;
        system(commandLine.c_str()); // remove the old feature file before updating
        std::cout << featurefilename.c_str() << std::endl;
        std::ofstream efile( featurefilename.c_str() , std::ios::app );
        efile << caseID << ":  " << "\n";
        efile.close();
        // calculate muscle features in t2 images
        #ifdef HISTOGRAM_TEXTURE
        std::cout << "features in t2 images:  " << std::endl;
        efile.open( featurefilename.c_str() , std::ios::app );
        efile.close();
        muscleFeatT2.features ( segMuscleLeft->GetOutput(), erodeSegMuscleLeft, t2LeftFemur, voxelVol, featurefilename, TRUE );
        #endif
    }
    return EXIT_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateEllipseImage(synImageType::Pointer image)
{
    typedef itk::EllipseSpatialObject< 3 >   EllipseType;
    typedef itk::SpatialObjectToImageFilter<EllipseType, synImageType >   SpatialObjectToImageFilterType;
    std::string commandLine;
    std::string filename;
    int x, y, z, r, noBig = 20;
    SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
    synImageType::SizeType size;

    size[ 0 ] =  256;
    size[ 1 ] =  256;
    size[ 2 ] =  180;

    imageFilter->SetSize( size );
    synImageType::SpacingType spacing;
    spacing.Fill(1);
    imageFilter->SetSpacing(spacing);

    EllipseType::Pointer ellipse    = EllipseType::New();
    EllipseType::ArrayType radiusArray;

    for (int i = 0; i < 200; i++) {
        // size
        if(i < noBig) { //20~30
            r = rand() % 11 + 20;
       //     y = rand() % 6 + 10;
        //    z = rand() % 6 + 10;
        }
        else {//1~15
            r = rand() % 6 + 5;
          //  y = rand() % 5 + 1;
          //  z = rand() % 5 + 1;
        }
        std::cout << "size:  " << r << "  " << r << "  " << r << std::endl;
        radiusArray[0] = r;
        radiusArray[1] = r;
        radiusArray[2] = r;
        ellipse->SetRadius(radiusArray);

        typedef EllipseType::TransformType                 TransformType;
        TransformType::Pointer transform = TransformType::New();
        transform->SetIdentity();

        TransformType::OutputVectorType  translation;
        TransformType::CenterType        center;
        //location
        if(i < noBig) { // 30-254
            x = rand() % 205 + 50;
            y = rand() % 205 + 50;
            z = rand() % 205 + 50;
        }
        else{        // 15~254
            x = rand() % 225 + 30;
            y = rand() % 225 + 30;
            z = rand() % 225 + 30;
        }
        std::cout << "location:  " << x << "  " << y << "  " << z << std::endl;

        for(int j = z - r; j <= z + r; j++) {
            for(int k = y - r; k <= y + r; k++) {
                for(int l = x - r; l <= x + r; l++) {
                    float d = sqrt((l - x) * (l - x) + (k - y) * (k - y) + (j - z) * (j - z));
                    if (d <= r) {
                        float e = 500 * pow(E_CONSTANT, -d);
                        synImageType::IndexType index;
                        index[0] = l;
                        index[1] = k;
                        index[2] = j;
                        if(e < 10)
                            e = 10;
                        image->SetPixel(index, e);
                    }
                }
            }
        }
    }
    syntheticWriter->SetInput(image);
    std::ostringstream strImageNo;
    strImageNo << 1;
    filename = "synImage" + strImageNo.str() + ".nrrd";
    syntheticWriter->SetFileName(filename.c_str());
    syntheticWriter->UseCompressionOn();
    syntheticWriter->Update();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void T2ValueMapCreate(DMDData data, DMDData::OrientedImageType::Pointer t2Calc, int t2mm_slice,DMDData::LfITKType data_array[2][DATA_PAIRS-1], DMDData::StrITKType dataDirec, DMDData::StrITKType caseID, DMDCurveFit curveFit)
{
    DMDData::OrientedImageType::Pointer                     t2CalcSmooth = DMDData::OrientedImageType::New();
    DMDData::OrientedImageType::Pointer                     t2Fit = DMDData::OrientedImageType::New();
    DMDData::OrientedImageType::SizeType                    t2CalcSize = t2Calc->GetLargestPossibleRegion().GetSize(); // get the image information
    DMDData::OrientedImageType::IndexType                   index;                                                     // image index of t2 calc
    DMDData::IntITKType                                     data_pairs = DATA_PAIRS, no_multi_slice = t2mm_slice;              // number of echoes
    DMDData::LfITKType                                      polyfitcoe = 0, t2calcThresh = 200;                        // the coefficient of the poly
    // load header information of t2 calc image for re-orientation of t2Fit
    DMDData::OrientedImageType::RegionType                  region = t2Calc->GetLargestPossibleRegion() ;
    DMDData::OrientedImageType::PointType                   origin = t2Calc->GetOrigin() ;
    DMDData::OrientedImageType::SpacingType                 spacing = t2Calc->GetSpacing() ;
    DMDData::OrientedImageType::DirectionType               direction = t2Calc->GetDirection() ;
    DMDData::OrientedImageType::SizeType                    t2MapSize = t2Calc->GetLargestPossibleRegion().GetSize();
    float                                                   mint2val = 0, maxt2val = 0;
    DMDData::OrientedWriterType::Pointer                    writer = DMDData::OrientedWriterType::New();
    writer->UseCompressionOn();   // turn output image compression on
    DMDData::StrITKType                                     outputFilename;

    t2MapSize[2] = no_multi_slice; // set the size in depth direction of t2 value image
    t2Fit->SetRegions( t2MapSize );
    t2Fit->SetOrigin( origin );
    t2Fit->SetSpacing( spacing );
    t2Fit->SetDirection( direction );
    t2Fit->Allocate();
    for (int i = 1; i < data_pairs; i++)       // set time for each echo
        data_array[0][i - 1] = 0.001 * ECHO_TIME * (i + 1);
    // smooth the t2 calc with gaussian filter
    data.smoothGaussian( t2Calc, t2CalcSmooth, 1 );
    for(int j = 0 ; j < (int)t2CalcSize[2] / data_pairs; j++) {
        for(index[1] = 0; index[1] < (int)t2CalcSize[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)t2CalcSize[0]; index[0]++) {
                bool marker = 1;      // if one of the pixel values less than t2calcThresh,  mark = 0, otherwise mark = 1
                for (int i = 1; i < data_pairs; i++) {
                    if((i - 1) == 0){
                        index[2] = j + no_multi_slice * (i - 1);
                        if( t2CalcSmooth->GetPixel(index) <= t2calcThresh ) {  // thresholding (200) remove noise
                            marker = 0;
                        }
                    }
                    index[2] = j + no_multi_slice * i;
                    // assign value to ten data points
                    if (marker == 1){
                        data_array[1][i - 1] = t2CalcSmooth->GetPixel(index);
                        data_array[1][i - 1] = log(data_array[1][i - 1]);
                    }
                    else {
                        data_array[1][i - 1] = 0;
                    }
                }
                if(marker == 1){
                    polyfitcoe = curveFit.lslcurvefit( data_array );
                    data_array[0][0] = 0.001 * ECHO_TIME * 2 ;
                    index[2] = j;
                    t2Fit->SetPixel(index, (-1000 / polyfitcoe));
                    if(t2Fit->GetPixel(index) < 0)
                        t2Fit->SetPixel(index, 0);
                    if( index[0] == 0 && index[1] == 0 && index[2] == 0)
                        mint2val = maxt2val = t2Fit->GetPixel(index);
                    else {
                        if (mint2val > t2Fit->GetPixel(index))
                            mint2val = t2Fit->GetPixel(index);
                        if (maxt2val < t2Fit->GetPixel(index))
                            maxt2val = t2Fit->GetPixel(index);
                    }
                }
                else {
                    index[2] = j;
                    t2Fit->SetPixel(index, 0);
                }
            }
        }
    }
    writer->SetInput( t2Fit );
    // write out t2 map
    outputFilename = dataDirec + "Processed/" + caseID + "_T2FIT.nrrd";
    data.dataWriter(writer, outputFilename);
    std::cout << "t2 value image was created" << std::endl;
}
