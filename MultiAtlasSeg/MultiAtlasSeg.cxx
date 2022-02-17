/*=========================================================================

  Version:   $  1.20 $
             1. Compute the performance level for all the labels (v 1.13, 12-13-2011)
             2. Iterative weighted mojority voting '-w' (v 1.13, 12-13-2011)
             3. Added computation of harmonic energy (v 1.13, 01-20-2012)
             4. Energy constant was set to 0.1 (v 1.13, 02-16-2012)
             5. (1) run batFloydSearchAll1 batFloydSearchAll2 batFloydSearchAll3
                (2) remove all folders under ~DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/
                (3) run process at ~DMD/Applications/DMDMultiAtlas/data/finalDeformedImage/
                (4) run batWMajorityVoting1 and batWMajorityVoting2 at NeuroLib/ImageMath/ for weighted majority voting
                (5) run batEvaluationWeightedMV for evaluating
                    (v 1.13, 02-17-2012)
             6. Calculate shape energy (circularity) (v 1.13, 03-02-2012)
                circularity = 4PI * (Area / (Perimeter ^ 2))
             7. Add threshold for the length of path determined by Floyd algorithm (04-27-2012)
             8. Only use the closest cases to target image (v 1.20 05-02-2012)
             9. Use single deformation field file to calculate harmonic energy (v 1.21 08-15-2012)
            10. Normalize image before calculating intensity energy (v 1.21 08-16-2012)
            11. Added distance normalization (v 1.21 08-17-2012)
            12. Add Beysian correction (v 1.3 11-28-2012)

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#define REC_FEM 1
#define SEMIT 4
#define CRAN_SART 5
#define BIC_FEM 2
#define ADDUCTOR 3
#define GRACILIS 6
#define BACKGROUND 0
#define NUMBER_OF_CASE_USED 3
#define NUMBER_OF_MUSCLE 6
#define _INDIVIDUAL_MUSCLE
#define ENERGY_CONST 0.0
#define MAX_NODE -1
#define WEIGHT_THRESHOLD 100.0
#define FLOYD_PATH_LENGTH 1
#define PI 3.1415926

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkOrientedImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>

#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkTranslationTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkAffineTransform.h"
#include "itkTransformFactory.h"
#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkWarpImageFilter.h"
#include "itkCommand.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkWarpHarmonicEnergyCalculator.h"
#include "itkGridForwardWarpImageFilter.h"
#include "itkVectorCentralDifferenceImageFunction.h"

void SortStringList(std::string *strList, int size)
{
    std::vector<std::string> stringarray;

    for (int i = 0; i < size ; i++){
        stringarray.push_back(strList[i]);
    }
    std::sort(stringarray.begin(), stringarray.end());

    int j = 0;
    for (std::vector<std::string>::iterator it=stringarray.begin(); it!=stringarray.end(); ++it) {
        std::string tmp = *it;
        strList[j] = tmp;
        j++;
    }
}


inline int ipExistsArgument(char **argv, const char *keystr) {
  for (int i = 1; argv[i]; i++)
    if (strstr(argv[i], keystr) == argv[i])
        return 1;
  return 0;
}

int main( int argc, char * argv[] )
{
    if (argc <= 1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
        std::cout << "Multi-Atlas Segmentation Tool 1.13 version (Feb 2012)" << std::endl;
        std::cout << "Usage: " << std::endl;
        std::cout << argv[0] << "  fixedImage  movingImage" << std::endl;
        std::cout << " -z" << "  chose the most similar template" << std::endl;
        std::cout << " -b" << "  chose the most similar template" << std::endl;
        std::cout << " -g" << "  establish graph and search optimal route between image pairs through graph" << std::endl;
        std::cout << " -m" << "  run majority voting" << std::endl;
        std::cout << " -v" << "  run weighted majority voting" << std::endl;
        std::cout << " -c" << "  run STAPLE label fusion" << std::endl;
        std::cout << " -p" << "  calculate metrics" << std::endl;
        std::cout << " -u" << "  calculate shape energy (circularity)" << std::endl;
        std::cout << " -e" << "  calculate harmonic energy" << std::endl;
        std::cout << " -n" << "  normalize the weighting factors" << std::endl;
        std::cout << " -N" << "  normalize the weighting factors, supplying 2 files" << std::endl;
        std::cout << " -help " << " print out this instruction" << std::endl;
        std::cout << " fixedImage  movingImage NumberOfAtlases NumberOfCases" << " print out this instruction" << std::endl;
        return EXIT_FAILURE;
    }

    int NUMBER_OF_CASE = atoi(argv[argc - 1]);
    int NUMBER_OF_ATLAS = atoi(argv[argc - 2]);

    bool option_V = ipExistsArgument(argv, "-v");
    bool option_S = ipExistsArgument(argv, "-s");
    bool option_W = ipExistsArgument(argv, "-w");
    bool option_Z = ipExistsArgument(argv, "-z");
    bool option_B = ipExistsArgument(argv, "-b");
    bool option_G = ipExistsArgument(argv, "-g");
    bool option_M = ipExistsArgument(argv, "-m");
    bool option_C = ipExistsArgument(argv, "-c");
    bool option_P = ipExistsArgument(argv, "-p");
    bool option_U = ipExistsArgument(argv, "-u");
    bool option_E = ipExistsArgument(argv, "-e");
    bool option_N = ipExistsArgument(argv, "-n");
    bool option_N2 = ipExistsArgument(argv, "-N");
    const     unsigned int   Dimension = 3;
    typedef   unsigned short  PixelType;
    typedef itk::Image< PixelType, Dimension >   ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    int DATASET_SIZE = 25;
    float p[DATASET_SIZE];  // sensitivity
    float q[DATASET_SIZE];   // specificity
    float alpha = 0, beta = 0, gama = 0;

    alpha = 0.5;
    beta = 0.5;

    ReaderType::Pointer fixedReader  = ReaderType::New();
    ReaderType::Pointer movingReader = ReaderType::New();
    ReaderType::Pointer fixedSegReader  = ReaderType::New();
    ReaderType::Pointer movingSegReader = ReaderType::New();
    ReaderType::Pointer fixedReaderT  = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    std::vector<ReaderType::Pointer> fixedReaders(DATASET_SIZE);
    std::vector<ReaderType::Pointer> moveReaders(DATASET_SIZE);

    for(int i = 0; i < DATASET_SIZE; i++) {
        fixedReaders[i]  = ReaderType::New();
        moveReaders[i]  = ReaderType::New();
    }

    if(  option_C ) {
        //label 1
        float p[DATASET_SIZE];  // sensitivity
        float q[DATASET_SIZE];   // specificity
        float A = 0, B = 0, R = 0, f = 0, fT1 = 0, fT0 = 0, f1 = 1, f0 = 1, senThresh = 0.95;
        typedef itk::Image< float, 3 >                         OrientedImageType;
        char filename[1024];
        int cases[25] = {39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 58, 59, 69, 71, 73, 88, 95, 107}, usedCase[25] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1, 1, 1};
        int segLabel = atoi(argv[2]), fixCase = atoi(argv[3]), numIteration = 1, numCaseUsed = DATASET_SIZE, numCaseWillUse = DATASET_SIZE;
        bool terminateCriteria = 1;
    //    char *quotationMark = "\"";
        std::string commandLine;
        std::ifstream readParameter;
        std::ostringstream strFixCase;
        strFixCase << fixCase;
        // remove performance text file created before
        if(fixCase < 100)
            commandLine = "rm ~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0" + strFixCase.str() + "/performance0" + strFixCase.str() + "_syn_100x50x25_label_AllMuscle_*.txt";
        else
            commandLine = "rm ~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/" + strFixCase.str() + "/performance" + strFixCase.str() + "_syn_100x50x25_label_AllMuscle_*.txt";
        std::cout << commandLine << std::endl;
        system(commandLine.c_str()); //warp image
        while (terminateCriteria != 0) {
            float sen[DATASET_SIZE], averageSen[DATASET_SIZE];  // sensitivity
            float spe[DATASET_SIZE];   // specificity
            std::ostringstream numIter;
            numCaseUsed = numCaseWillUse;
            numCaseWillUse = 0;
            numIter << numIteration;
            //exclude target image
            for (int i = 0; i < DATASET_SIZE; i++){
                if (cases[i] == fixCase)
                    usedCase[i] = 0;
            }
            commandLine = "crlSTAPLE";
            for (int i = 0; i < NUMBER_OF_CASE; i++){
                std::ostringstream strTargetCase;
                if(usedCase[i]) {
                    strTargetCase << cases[i];
                    if(fixCase < 100 && cases[i] < 100)
                        commandLine = commandLine + " ~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0" + strFixCase.str() + "/deformedImage_0" + strTargetCase.str() + "to0" + strFixCase.str() + "_seg.nrrd";
                    else if(fixCase >= 100 && cases[i] < 100)
                        commandLine = commandLine + " ~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/" + strFixCase.str() + "/deformedImage_0" + strTargetCase.str() + "to" + strFixCase.str() + "_seg.nrrd";
                    else if(fixCase < 100 && cases[i] >= 100)
                        commandLine = commandLine + " ~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0" + strFixCase.str() + "/deformedImage_" + strTargetCase.str() + "to0" + strFixCase.str() + "_seg.nrrd";
                    else if(fixCase >= 100 && cases[i] >= 100)
                        commandLine = commandLine + " ~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/" + strFixCase.str() + "/deformedImage_" + strTargetCase.str() + "to" + strFixCase.str() + "_seg.nrrd";
                }
            }

            std::cout << commandLine << std::endl;
            system(commandLine.c_str());//run STAPLE
            exit(1);
            if(fixCase < 100)
                commandLine = "crlIndexOfMaxComponent seg_0" + strFixCase.str() + "_fromAllScans_syn_100x50x25_AllMuscle.nrrd seg_0" + strFixCase.str() + "_fromAllScans_syn_100x50x25_warp_AllMuscle_" + numIter.str() + ".nrrd";
            else
                commandLine = "crlIndexOfMaxComponent seg_" + strFixCase.str() + "_fromAllScans_syn_100x50x25_AllMuscle.nrrd seg_" + strFixCase.str() + "_fromAllScans_syn_100x50x25_warp_AllMuscle_" + numIter.str() + ".nrrd";
            std::cout << commandLine << std::endl;
            system(commandLine.c_str()); //warp image
            // load performance level
            if(fixCase < 100)
                sprintf(filename, "/biomed-resimg/NDRC/DuchenneDystrophy/data/canine/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0%d/performance0%d_syn_100x50x25_label_AllMuscle_%d.txt", fixCase, fixCase, numIteration);
            else
                sprintf(filename, "/biomed-resimg/NDRC/DuchenneDystrophy/data/canine/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/%d/performance%d_syn_100x50x25_label_AllMuscle_%d.txt", fixCase, fixCase, numIteration);

            std::cout << filename << std::endl;
            readParameter.open(filename);
            if (readParameter.is_open()) {
                int j = 0;
                while (!readParameter.eof() && j < (NUMBER_OF_CASE - 1)) {
                    char pa[100];
                    float CNum = 0;
                    averageSen[j] = 0;
                    while (usedCase[j] != 1 && j < (NUMBER_OF_CASE - 1)){
                        j++;
                        averageSen[j] = 0;
                    }
                    CNum  = 0;
                    while (CNum <= 0 || CNum >= 1) {
                        readParameter.getline(pa, 256, ',');
                        CNum = atof(pa);
                    }
                    spe[j] = CNum;
                    std::cout << j << "  " ;
                    for (int k = 0; k < NUMBER_OF_MUSCLE; k++){
                        CNum  = 0;
                        while (CNum <= 0 || CNum >= 1) {
                            readParameter.getline(pa, 256, ',');
                            CNum = atof(pa);
                        }
                        sen[j * NUMBER_OF_MUSCLE + k] = CNum;
                        averageSen[j] += sen[j * NUMBER_OF_MUSCLE + k];
                        std::cout << sen[j * NUMBER_OF_MUSCLE + k] <<  " ";
                    }
                    averageSen[j] /= NUMBER_OF_MUSCLE;
                    std::cout << "  specificity: " << spe[j] << " average sensitivity: " << averageSen[j] << std::endl;
                    j++;
                }
            }
            readParameter.close();
        // check which cases are able to use
            for (int i = 0; i < DATASET_SIZE; i++){
                if((usedCase[i] == 0) || (averageSen[i] < 0.90)){
                    usedCase[i] = 0;
                }
                else {
                    numCaseWillUse++;
                }
            }
            std::cout << numCaseWillUse << "  " << numCaseUsed << std::endl;
           // check termination criteria
            if ((numCaseWillUse == numCaseUsed) || (numCaseWillUse < 6))
                terminateCriteria = 0;
            std::cout << "terminate criteria:  " << terminateCriteria << std::endl;
            numIteration++; //next iteration
        }
        std::cout << "STAPLE finished!" << std::endl;
        exit(1);
    }

    if( option_B ) {
        //label 1
        float A = 0, B = 0, R = 0, f = 0, b = 0, fT1 = 0, fT0 = 0, f1 = 1, f0 = 1, senThresh = 0.01;//senThresh = 0.95
        typedef itk::Image< float, 3 >                         OrientedImageType;
        char filename[1024];
        int cases[25] = {39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 58, 59, 69, 71, 73, 88, 95, 107}, usedCase[25] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        int segLabel = atoi(argv[2]), fixCase = atoi(argv[3]), numIteration = 1, numCaseUsed = DATASET_SIZE, numCaseWillUse = DATASET_SIZE;
        bool terminateCriteria = 1;
     //   char *quotationMark = "\"";
        std::string commandLine;
        std::ifstream readParameter;
        std::ostringstream strFixCase;
        strFixCase << fixCase;

        std::vector<ImageType::Pointer> segmentations(DATASET_SIZE);
        ImageType::Pointer segmentation = ImageType::New();
        ImageType::Pointer T = ImageType::New();
        for (int i = 0; i < DATASET_SIZE; i++)
            segmentations[i] = ImageType::New();
        // create final segmentation based on the Baysian Decision Theory
        // load segmentation 1
        while (R < NUMBER_OF_CASE_USED) {
            R = 0;
            int j = 0;
            if (cases[j] == fixCase) { //used for STAPLE
                j++;
            }
            for (int i = 0; i < DATASET_SIZE - 1; i++){ //used for STAPLE
                #ifdef INDIVIDUAL_MUSCLE
                sprintf(filename, "../../../Multi_Atlas_Build/Data/Prepare_Tree_Atlas_NatHist_Normal/ants/syn_100x50x25_0%dto0%d/deformedImage_0%dto0%d_seg_%d.nrrd", cases[j], fixCase, cases[j], fixCase, segLabel);
                #else
                if(fixCase < 100){
                    if(cases[j] < 100)
                        sprintf(filename, "/home/jiahuiw/DMD/Multi_Atlas_Build/Data/ANTS_Fine_reg/DeformedImage/0%d/deformedImage_0%dto0%d_seg.nrrd", fixCase, cases[j], fixCase);
                    else
                        sprintf(filename, "/home/jiahuiw/DMD/Multi_Atlas_Build/Data/ANTS_Fine_reg/DeformedImage/0%d/deformedImage_%dto0%d_seg.nrrd", fixCase, cases[j], fixCase);
                }
                else{
                    if(cases[j] < 100)
                        sprintf(filename, "/home/jiahui/DMD/Multi_Atlas_Build/Data/ANTS_Fine_reg/DeformedImage/%d/deformedImage_0%dto%d_seg.nrrd", fixCase, cases[j], fixCase);
                    else
                        sprintf(filename, "/home/jiahui/DMD/Multi_Atlas_Build/Data/ANTS_Fine_reg/DeformedImage/%d/deformedImage_%dto%d_seg.nrrd", fixCase, cases[j], fixCase);
                }
                #endif
                std::cout << filename << std::endl;
                fixedReaders[i]->SetFileName(filename);
                fixedReaders[i]->Update();
                segmentations[i] = fixedReaders[i]->GetOutput();
                if (p[i] >= senThresh)
                    R++;
                j++;
                while (cases[j] == fixCase) {
                    j++;
                }
            }
            senThresh -= 0.05;
        }
        senThresh += 0.05;
        std::cout << "number of scan used: " << R << "  sensitivity: " << senThresh << std::endl;
        getchar();
        //used for weighted majority voting
        if(fixCase < 100)
            sprintf(filename, "/biomed-resimg/NDRC/DuchenneDystrophy/data/canine/Research/NeuroLib/niral_utilities/trunk/ImageMath/WeightedMVOutput_seg_0%d.nrrd", fixCase);
        else
            sprintf(filename, "/biomed-resimg/NDRC/DuchenneDystrophy/data/canine/Research/NeuroLib/niral_utilities/trunk/ImageMath/WeightedMVOutput_seg_%d.nrrd", fixCase);
        fixedReaderT->SetFileName(filename);
        fixedReaderT->Update();
        segmentation = fixedReaderT->GetOutput();
        ImageType::SizeType                    size = segmentation->GetLargestPossibleRegion().GetSize();
        ImageType::IndexType                   index;
        for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
            for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
                for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                   // segmentation->SetPixel(index, 0);
                    fT1 = fT0 = 0;
                    f1 = f0 = 1;
                    A = 0; B = 0; f = 0;
                    // calculate f(Ti = 1) and f(Ti = 0)
                    for (int i = 0; i < DATASET_SIZE - 1; i++){ //used for STAPLE
                        if (segmentations[i]->GetPixel(index) == segLabel && p[i] >= senThresh) {
                            fT1++;
                        }
                    }
                    fT1 /= R;
                    fT0 = 1 - fT1;
                    // calculate f(D|T)
                    f0 = f1 = 1;
                    for (int i = 0; i < DATASET_SIZE - 1; i++){ //used for STAPLE
                        if (p[i] >= senThresh) {
                            if (segmentations[i]->GetPixel(index) == segLabel) {
                                f1 = f1 * p[i];
                                f0 = f0 * q[i];
                            }
                            else{
                                f1 = f1 * (1 - p[i]);
                                f0 = f0 * (1 - q[i]);
                            }
                        }
                    }
                    A = f1 * fT1;
                    B = f0 * fT0;
                    if((A + B) > 0){
                        f = A / (A + B) * 1;
                        b = B / (A + B) * 1;
                    }
                    else{
                        f = 0;
                        b = 1;
                    }
                    if (f >= b){
                        f = segLabel; //foreground
                    }
                    else
                        f = 0; //background
                    if(segmentation->GetPixel(index) == segLabel){
                        segmentation->SetPixel(index, f);
                    }
                    else
                        segmentation->SetPixel(index, 0);
                }
            }
        }

        #ifdef INDIVIDUAL_MUSCLE
        sprintf(filename, "segmentation0%d_%d.nrrd", fixCase, segLabel);
        #else
        if (fixCase < 100)
            sprintf(filename, "segmentation0%d_AllMuscle_%d.nrrd", fixCase, segLabel);
        else
            sprintf(filename, "segmentation%d_AllMuscle_%d.nrrd", fixCase, segLabel);
        #endif
        writer->UseCompressionOn();
        writer->SetInput(segmentation);
        writer->SetFileName(filename);
        try {
            writer->Update();
        }
        catch (itk::ExceptionObject &ex){
            std::cout << "write:" << ex << std::endl;
            exit(0);
        }
        writer->Update();
    }

    if( option_W ) {
        // weighted majority voting
        fixedReader->SetFileName(argv[2]);
        fixedReader->Update();
        ImageType::SizeType                  size = fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize();
        ImageType::IndexType                 indexInput;
        float averagePerformance = 0;
        int usedSegmentation = DATASET_SIZE;
        int cases[25] = {39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 58, 59, 69, 71, 73, 88, 95, 107}, usedCase[25] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        float weightParameter[25], sumWeight = 1, sumWeightLastIteration = 0;
        int j = 0, k = 0;
        bool firstIteration = 1;
        std::ostringstream strFixCase;
        strFixCase << argv[argc - 3];
        std::string commandLine;
        std::ostringstream strTargetCase;

        while (usedSegmentation >= 6 && sumWeight != sumWeightLastIteration){
            sumWeightLastIteration = sumWeight;
            ReaderType::Pointer segmentationReader  = ReaderType::New();
            segmentationReader->SetFileName(argv[ 2 ]);
            segmentationReader->Update();
            j = 0; k = 0; usedSegmentation = 0;
            commandLine = "rm caseUsedIterativeWMV" + strFixCase.str() + ".txt";
            system(commandLine.c_str());
            commandLine = "rm performanceIterativeWMV" + strFixCase.str() + ".txt";
            system(commandLine.c_str());

        for (int i = 0; i < DATASET_SIZE; i++){
            float RecFemMove = 0, CraSarMove = 0, AddMove = 0, BicFemMove = 0, GraMove = 0, SemTenMove = 0, BackgroundMove = 0;
            float RecFemInt = 0, CraSarInt = 0, AddInt = 0, BicFemInt = 0, GraInt = 0, SemTenInt = 0, BackgroundInt = 0;
            float RecFemUni = 0, CraSarUni = 0, AddUni = 0, BicFemUni = 0, GraUni = 0, SemTenUni = 0, BackgroundUni = 0;
            float RecFemSum = 0, CraSarSum = 0, AddSum = 0, BicFemSum = 0, GraSum = 0, SemTenSum = 0, BackgroundSum = 0;
            float specificity = 0;

            if(atoi(argv[argc - 3]) != cases[i]){
                moveReaders[k]->SetFileName(argv[k + 3]);
                moveReaders[k]->Update();
            for(indexInput[2] = 0; indexInput[2] < (int)size[2]; indexInput[2]++) {
                for(indexInput[1] = 0; indexInput[1] < (int)size[1]; indexInput[1]++) {
                    for(indexInput[0] = 0; indexInput[0] < (int)size[0]; indexInput[0]++) {
                        float fixPixel = segmentationReader->GetOutput()->GetPixel(indexInput);
                        float  movePixel = moveReaders[k]->GetOutput()->GetPixel(indexInput);
                        if (fixPixel == BACKGROUND && movePixel == BACKGROUND)
                            BackgroundInt++;
                        if (fixPixel == BACKGROUND || movePixel == REC_FEM)
                            BackgroundUni++;
                        if (fixPixel == BACKGROUND){
                            BackgroundSum++;
                            BackgroundMove++;
                        }
                        if (movePixel == BACKGROUND)
                            BackgroundSum++;

                        if (fixPixel == REC_FEM && movePixel == REC_FEM)
                            RecFemInt++;
                        if (fixPixel == REC_FEM || movePixel == REC_FEM)
                            RecFemUni++;
                        if (fixPixel == REC_FEM){
                            RecFemSum++;
                            RecFemMove++;
                        }
                        if (movePixel == REC_FEM)
                            RecFemSum++;

                        if (fixPixel == BIC_FEM && movePixel == BIC_FEM)
                            BicFemInt++;
                        if (fixPixel == BIC_FEM || movePixel == BIC_FEM)
                            BicFemUni++;
                        if (fixPixel == BIC_FEM){
                            BicFemSum++;
                            BicFemMove++;
                        }
                        if (movePixel == BIC_FEM)
                            BicFemSum++;

                        if (fixPixel == ADDUCTOR && movePixel == ADDUCTOR)
                            AddInt++;
                        if (fixPixel == ADDUCTOR || movePixel == ADDUCTOR)
                            AddUni++;
                        if (fixPixel == ADDUCTOR){
                            AddSum++;
                            AddMove++;
                        }
                        if (movePixel == ADDUCTOR)
                            AddSum++;

                        if (fixPixel == SEMIT && movePixel == SEMIT )
                            SemTenInt++;
                        if (fixPixel == SEMIT || movePixel == SEMIT)
                            SemTenUni++;
                        if (fixPixel == SEMIT){
                            SemTenSum++;
                            SemTenMove++;
                        }
                        if (movePixel == SEMIT)
                            SemTenSum++;

                        if (fixPixel == CRAN_SART && movePixel == CRAN_SART)
                            CraSarInt++;
                        if (fixPixel == CRAN_SART || movePixel == CRAN_SART)
                            CraSarUni++;
                        if (fixPixel == CRAN_SART){
                            CraSarSum++;
                            CraSarMove++;
                        }
                        if (movePixel == CRAN_SART)
                            CraSarSum++;

                        if (fixPixel == GRACILIS && movePixel == GRACILIS)
                            GraInt++;
                        if (fixPixel == GRACILIS || movePixel == GRACILIS)
                            GraUni++;
                        if (fixPixel == GRACILIS){
                            GraSum++;
                            GraMove++;
                        }
                        if (movePixel == GRACILIS)
                            GraSum++;
                    }
                }
            }
            k++;
            averagePerformance = (RecFemInt / RecFemMove + BicFemInt / BicFemMove + AddInt / AddMove + SemTenInt / SemTenMove + CraSarInt / CraSarMove + GraInt / GraMove) / 6;
            specificity = BackgroundInt / BackgroundMove;
            }
            if(atoi(argv[argc - 3]) == cases[j]){
                usedCase[j] = 0;
                weightParameter[j] = 0;
            }
            else{
                if(averagePerformance > 0.8 && usedCase[j] == 1) {
                    usedCase[j] = 1;
                    weightParameter[j] = averagePerformance;
                }
                else{
                    usedCase[j] = 0;
                    weightParameter[j] = 0;
                }
            }
            std::string filename;
            filename = "caseUsedIterativeWMV" + strFixCase.str() + ".txt";
            std::ofstream pfile(filename.c_str(), std::ios::app );
            if(usedCase[j] == 1) {
                pfile << cases[j] << "\n";
                usedSegmentation++;
            }
            pfile.close();
            std::cout << "the sensitivety for " << cases[j] << " is " << averagePerformance << "  the specificity is " << specificity << "  " << usedCase[j] << std::endl;
            std::cout << "case used: " << usedSegmentation << std::endl;

            if(firstIteration){
                filename = "performanceIterativeWMV" + strFixCase.str() + "_old.txt";
                std::ofstream efile(filename.c_str() , std::ios::app );
                if(atoi(argv[argc - 3]) == cases[j])
                    efile << specificity << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << "\n";
                else
                    efile << specificity << " " << RecFemInt / RecFemMove << " " << BicFemInt / BicFemMove << " " <<  AddInt / AddMove << " " << SemTenInt / SemTenMove << " " << CraSarInt / CraSarMove << " " << GraInt / GraMove << "\n";
                efile.close();
            }
                filename = "performanceIterativeWMV" + strFixCase.str() + ".txt";
                std::ofstream efile(filename.c_str() , std::ios::app );
                if(usedCase[j] == 1)
                    efile << specificity << " " << RecFemInt / RecFemMove << " " << BicFemInt / BicFemMove << " " <<  AddInt / AddMove << " " << SemTenInt / SemTenMove << " " << CraSarInt / CraSarMove << " " << GraInt / GraMove << "\n";
                else
                    efile << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << "\n";
                efile.close();
            j++;
        }
        firstIteration = 0;
        if(usedSegmentation >= 6){
            commandLine = "~/NeuroLib/ImageMath/ImageMath ";
            if(atoi(argv[argc - 3]) < 100) {
                int i = 0;
                while(usedCase[i] == 0){
                    i++;
                }
                std::ostringstream strTargetCase;
                if(cases[i] < 100){
                    strTargetCase << cases[i];
                    commandLine = commandLine + "~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0" + strFixCase.str() + "/deformedImage_0" + strTargetCase.str() + "to0" + strFixCase.str() + "_seg.nrrd ";
                }
                else{
                    strTargetCase << cases[i];
                    commandLine = commandLine + "~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0" + strFixCase.str() + "/deformedImage_" + strTargetCase.str() + "to0" + strFixCase.str() + "_seg.nrrd ";
                }
            }
            else{
                std::ostringstream strTargetCase;
                strTargetCase << cases[0];
                commandLine = commandLine + "~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/" + strFixCase.str() + "/deformedImage_0" + strTargetCase.str() + "to" + strFixCase.str() + "_seg.nrrd ";
            }
            commandLine = commandLine + "-weightedMajorityVoting " ;
            for (int i = 0; i < DATASET_SIZE; i++){
                std::ostringstream strTargetCase;
                strTargetCase << cases[i];
                if(usedCase[i] == 1){
                    if(atoi(argv[argc - 3]) < 100) {
                        commandLine = commandLine + "~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/0" + strFixCase.str() + "/deformedImage_0" + strTargetCase.str() + "to0" + strFixCase.str() + "_seg.nrrd ";
                    }
                    else{
                        commandLine = commandLine + "~/DMD/Multi_Atlas_Build/Data/NatHist_Atlas/STAPLE/" + strFixCase.str() + "/deformedImage_0" + strTargetCase.str() + "to" + strFixCase.str() + "_seg.nrrd ";
                    }
                }
            }
            //normalize weighting factors
            sumWeight = 0;
            for (int i = 0; i < DATASET_SIZE; i++){
                if(usedCase[i] == 1){
                    sumWeight += weightParameter[i];
                }
            }
            std::cout << sumWeight << std::endl;
            for (int i = 0; i < DATASET_SIZE; i++){
                if(usedCase[i] == 1 && sumWeight != 0){
                    weightParameter[i] /= sumWeight;
                }
            }
            commandLine = commandLine + "-weights " ;
            bool cammerMark = 0;
            for (int i = 0; i < DATASET_SIZE; i++){
                std::ostringstream strWeight;
                strWeight << weightParameter[i];
                if(usedCase[i] == 1){
                    if(cammerMark)
                        commandLine = commandLine + ",";
                    else
                        cammerMark = 1;
                    commandLine = commandLine + strWeight.str();
                }
            }

            if(atoi(argv[argc - 3]) < 100)
                commandLine = commandLine + " -outfile ~/NeuroLib/ImageMath/WeightedMVOutput_seg_0" + strFixCase.str() + ".nrrd ";

            else
                commandLine = commandLine + " -outfile ~/NeuroLib/ImageMath/WeightedMVOutput_seg_" + strFixCase.str() + ".nrrd ";
            std::cout << commandLine << std::endl;
            system(commandLine.c_str());
            commandLine = "cp performanceIterativeWMV" + strFixCase.str() + ".txt " + "performanceIterativeWMV" + strFixCase.str() + "_old.txt" ;
            system(commandLine.c_str());
        }
        else{

            commandLine = "cp performanceIterativeWMV" + strFixCase.str() + "_old.txt " + "performanceIterativeWMV" + strFixCase.str() + ".txt" ;
            system(commandLine.c_str());
            commandLine = "rm performanceIterativeWMV" + strFixCase.str() + "_old.txt";
            system(commandLine.c_str());
        }
        }
    }

    if( option_P ) {
        fixedReader->SetFileName(  argv[ 2 ] );
        movingReader->SetFileName( argv[ 3 ] );
    }

    if( option_S ) {
        fixedSegReader->SetFileName(  argv[ 2 ] );
        movingSegReader->SetFileName( argv[ 3 ] );
    }
    try {
        if( option_P ) {
            fixedReader->Update();
            movingReader->Update();
        }
        if( option_S ) {
           fixedSegReader->Update();
           movingSegReader->Update();
        }
    }
    catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception catched !" << std::endl;
    std::cerr << excep << std::endl;
    }

    if( option_S ) {
        ImageType::SizeType                  size = fixedSegReader->GetOutput()->GetLargestPossibleRegion().GetSize();
        ImageType::IndexType                 indexInput;
        float RecFemInt = 0, CraSarInt = 0, AddInt = 0, BicFemInt = 0, GraInt = 0, SemTenInt = 0;
        float RecFemUni = 0, CraSarUni = 0, AddUni = 0, BicFemUni = 0, GraUni = 0, SemTenUni = 0;
        float RecFemSum = 0, CraSarSum = 0, AddSum = 0, BicFemSum = 0, GraSum = 0, SemTenSum = 0;
        float averagePerformance = 0;
        int evaluateCase = atoi(argv[4]);

        for(indexInput[2] = 0; indexInput[2] < (int)size[2]; indexInput[2]++) {
            for(indexInput[1] = 0; indexInput[1] < (int)size[1]; indexInput[1]++) {
                for(indexInput[0] = 0; indexInput[0] < (int)size[0]; indexInput[0]++) {
                    float fixPixel = fixedSegReader->GetOutput()->GetPixel(indexInput);
                    float  movePixel = movingSegReader->GetOutput()->GetPixel(indexInput);
                    if (fixPixel == REC_FEM && movePixel == REC_FEM)
                        RecFemInt++;
                    if (fixPixel == REC_FEM || movePixel == REC_FEM)
                        RecFemUni++;
                    if (fixPixel == REC_FEM)
                        RecFemSum++;
                    if (movePixel == REC_FEM)
                        RecFemSum++;

                    if (fixPixel == BIC_FEM && movePixel == BIC_FEM)
                        BicFemInt++;
                    if (fixPixel == BIC_FEM || movePixel == BIC_FEM)
                        BicFemUni++;
                    if (fixPixel == BIC_FEM)
                        BicFemSum++;
                    if (movePixel == BIC_FEM)
                        BicFemSum++;

                    if (fixPixel == ADDUCTOR && movePixel == ADDUCTOR)
                        AddInt++;
                    if (fixPixel == ADDUCTOR || movePixel == ADDUCTOR)
                        AddUni++;
                    if (fixPixel == ADDUCTOR)
                        AddSum++;
                    if (movePixel == ADDUCTOR)
                        AddSum++;

                    if (fixPixel == SEMIT && movePixel == SEMIT )
                        SemTenInt++;
                    if (fixPixel == SEMIT || movePixel == SEMIT)
                        SemTenUni++;
                    if (fixPixel == SEMIT)
                        SemTenSum++;
                    if (movePixel == SEMIT)
                        SemTenSum++;

                    if (fixPixel == CRAN_SART && movePixel == CRAN_SART)
                        CraSarInt++;
                    if (fixPixel == CRAN_SART || movePixel == CRAN_SART)
                        CraSarUni++;
                    if (fixPixel == CRAN_SART)
                        CraSarSum++;
                    if (movePixel == CRAN_SART)
                        CraSarSum++;

                    if (fixPixel == GRACILIS && movePixel == GRACILIS)
                        GraInt++;
                    if (fixPixel == GRACILIS || movePixel == GRACILIS)
                        GraUni++;
                    if (fixPixel == GRACILIS)
                        GraSum++;
                    if (movePixel == GRACILIS)
                        GraSum++;
                }
            }
        }

     //   std::ofstream efile( "overlapBayesian.txt" , std::ios::app );
      //  std::ofstream efile( "overlapMajorityVoting.txt" , std::ios::app );
        //std::ofstream efile( "overlapWeightedMV.txt" , std::ios::app );
        std::ofstream efile( "overlapMostSimilar.txt" , std::ios::app );
       // std::ofstream efile( "overlapIterativeStaple.txt" , std::ios::app );
   //     std::ofstream efile( "overlapStaple.txt" , std::ios::app );
        if(evaluateCase == REC_FEM)
            efile <<"\nRec_Fem_dice_overlap: " << RecFemInt/RecFemUni << " " << RecFemInt << " " << RecFemUni << " " << RecFemSum << " dice_overlap: " << 2 * RecFemInt / RecFemSum ;
        if(evaluateCase == BIC_FEM)
            efile <<"\nBic_Fem_dice_overlap: " << BicFemInt/BicFemUni << " " << BicFemInt << " " << BicFemUni << " " << BicFemSum << " dice_overlap: " << 2 * BicFemInt / BicFemSum ;
        if(evaluateCase == ADDUCTOR)
            efile <<"\nAdd_dice_overlap: " << AddInt/AddUni << " " <<  AddInt << " " << AddUni << " " << AddSum << " dice_overlap: " << 2 * AddInt / AddSum;
        if(evaluateCase == SEMIT)
            efile <<"\nSemTen_dice_overlap: " << SemTenInt/SemTenUni << " " << SemTenInt << " " << SemTenUni << " " << SemTenSum << " dice_overlap: " << 2 * SemTenInt / SemTenSum;
        if(evaluateCase == CRAN_SART)
            efile <<"\nCra_Sar_dice_overlap: " << CraSarInt/CraSarUni << " " << CraSarInt << " " << CraSarUni << " " << CraSarSum << " dice_overlap: " << 2 * CraSarInt / CraSarSum ;
        if(evaluateCase == GRACILIS)
            efile <<"\nGra_dice_overlap: " << GraInt/GraUni << " " << GraInt << " " << GraUni << " " << GraSum << " dice_overlap: " << 2 * GraInt / GraSum << "\n\n";
        efile.close();
        averagePerformance = (2 * RecFemInt / RecFemSum + 2 * BicFemInt / BicFemSum + 2 * AddInt / AddSum + 2 * SemTenInt / SemTenSum + 2 * CraSarInt / CraSarSum + 2 * GraInt / GraSum) / 6;
        std::cout << "the sensitivety is " << averagePerformance << std::endl;
    }

    if( option_E ) {
        // compute harmonic energy
      typedef itk::Vector<double,3> DeformationPixelType;

        typedef itk::Image< DeformationPixelType, Dimension >         OrientedImageType;
        typedef itk::ImageFileReader<OrientedImageType> OrientedReaderType;
        OrientedReaderType::Pointer orientedreader = OrientedReaderType::New();
        OrientedImageType::Pointer deformationField;
        OrientedImageType::IndexType index;
        OrientedImageType::SizeType size;
        double HE = 0;
        std::string filename;
        std::ostringstream strFixCase;

        orientedreader->SetFileName(argv[2]);
        try {
            orientedreader->Update();
            deformationField = orientedreader->GetOutput();
            size = deformationField->GetLargestPossibleRegion().GetSize();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
        for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
	  for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
	    for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
	      // calculate harmonic energy
	      DeformationPixelType value = deformationField->GetPixel(index);
	      HE += sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
	    }
	  }
        }
        std::cout << "harmonic energy overall: " << HE << std::endl;
	HE = HE / (size[0] * size[1] * size[2]);

        std::cout << "harmonic energy: " << HE << std::endl;
        filename = argv[3];
        std::ofstream efile( filename.c_str() , std::ios::app );
        efile << HE << "\n";
        efile.close();
        return 0;
    }

    if( option_U ) {
        // compute shape energy
        typedef itk::Image< float, Dimension >         OrientedImageType;
        typedef itk::ImageFileReader<OrientedImageType> OrientedReaderType;
        typedef itk::ImageFileWriter<OrientedImageType> OrientedWriterType;
     //   OrientedReaderType::Pointer binaryMask = OrientedReaderType::New();
        OrientedReaderType::Pointer maskReader = OrientedReaderType::New();

        OrientedReaderType::Pointer orientedreaderx = OrientedReaderType::New();
        OrientedReaderType::Pointer orientedreadery = OrientedReaderType::New();
        OrientedReaderType::Pointer orientedreaderz = OrientedReaderType::New();
        OrientedWriterType::Pointer orientedwriter = OrientedWriterType::New();
        OrientedImageType::Pointer deformationFieldx, deformationFieldy, deformationFieldz;
        OrientedImageType::Pointer binaryMask;
        OrientedImageType::IndexType                   index, indexForSearch;
        OrientedImageType::SizeType                    size;
        float SE = 0, perimeter = 0, area = 0;

        maskReader->SetFileName( argv[ 2 ] );
        maskReader->Update();
        binaryMask = maskReader->GetOutput();
        size = binaryMask->GetLargestPossibleRegion().GetSize();
        index[2] = size[2] / 2;
        //for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                if (binaryMask->GetPixel(index) != 0) {
                    area++;
                    if (index[0] != 0 && index[1] != 0 && index[0] != (size[0] - 1) && index[1] != (size[1] - 1)) {
                        float a, b, c, d;
                        indexForSearch[0] = index[0] - 1; indexForSearch[1] = index[1]; indexForSearch[2] = index[2];
                        a = binaryMask->GetPixel(indexForSearch);
                        indexForSearch[0] = index[0] + 1; indexForSearch[1] = index[1]; indexForSearch[2] = index[2];
                        b = binaryMask->GetPixel(indexForSearch);
                        indexForSearch[0] = index[0]; indexForSearch[1] = index[1] - 1; indexForSearch[2] = index[2];
                        c = binaryMask->GetPixel(indexForSearch);
                        indexForSearch[0] = index[0]; indexForSearch[1] = index[1] + 1; indexForSearch[2] = index[2];
                        d = binaryMask->GetPixel(indexForSearch);
                        if (a == 0 || b == 0 || c == 0 || d == 0)
                            perimeter++;
                    }
                }
            }
        }
        //}

        SE = 4 * PI * (area / (perimeter * perimeter));
        std::cout << "shape energy: " << SE << "  " << area << "  " << perimeter << std::endl;
        std::ofstream efile( "shapeEnergy.txt" , std::ios::app );
        efile << SE << "\n";
        efile.close();
        return 0;
    }

    if( option_Z ) {
        float graph[NUMBER_OF_CASE][NUMBER_OF_CASE][2], circularity[NUMBER_OF_CASE]; //include deformation from two directions
        float intensityE, harmonicE, shapeE, tmpHE = 0;
        int startNode = atoi(argv[2]), endNode = atoi(argv[3]), floydRoute[NUMBER_OF_CASE];
        std::ostringstream strFixCase;
        strFixCase << argv[argc - 3];
        std::ifstream efile( "harmonicEnergyNormalizedIdentical.txt");
        std::ifstream invefile( "harmonicEnergyInverse.txt");
        std::ifstream iefile( "intensityEnergyNormalizedIdentical.txt");
        std::ifstream sefile( "shapeEnergyNormalizedIdentical.txt");
        std::string filename = "templateForSeg_" + strFixCase.str() + ".txt";
        std::string commandLine;
        std::ofstream templatefile(filename.c_str(), std::ios::app );
        int cases[25] = {39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 58, 59, 69, 71, 73, 88, 95, 107}, caseUsed[NUMBER_OF_CASE];
        int graphHistogram[120] = {0};


        for (int i = 0; i < NUMBER_OF_CASE; i++){
            floydRoute[i] = MAX_NODE;
            caseUsed[i] = 0;
        }

        int m = 0;
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            sefile >> shapeE;
            circularity[i] = shapeE;
        }
        sefile.close();
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            for (int j = 0; j < NUMBER_OF_CASE; j++){
                if(i == j){
                    graph[i][j][0] = -1;
                    graph[i][j][1] = -1;
                }
                else {
                    m++;
                    efile >> harmonicE;
                    iefile >> intensityE;
                    graph[i][j][0] = alpha * harmonicE;
                    graph[i][j][0] += beta * intensityE;
                    graph[i][j][0] += gama * fabs(circularity[i] - circularity[j]) + ENERGY_CONST;
                    int tmp = round(graph[i][j][0] * 100);
                    if((tmp - 1) < 0)
                        graphHistogram[0]++;
                    else if ((tmp - 1) > 120)
                        graphHistogram[119]++;
                    else
                        graphHistogram[tmp - 1]++;
                    invefile >> tmpHE;
                    graph[i][j][1] = tmpHE;
                }
            }
        }
        efile.close();
        invefile.close();
        iefile.close();
        //search shortest route from startNode to endNode using Dijkstra
        int k = 0;
        float *Dijkstra = new float[NUMBER_OF_CASE];
        for (int i = 0; i < NUMBER_OF_CASE; i++)
            Dijkstra[i] = 0;

        while (startNode != endNode){
            float minDist = 0;
            int minNode = 0, DijkstraNode = 0;
            k++;
            for (int i = 0; i < NUMBER_OF_CASE; i++){
                if (i == 0){
                    if(graph[startNode][i][0] != -1) {
                        while (Dijkstra[i] == 1)
                            i++;
                        minDist = graph[startNode][i][0];
                        Dijkstra[startNode] = 1;
                        DijkstraNode = i;
                    }
                    else {
                        while (Dijkstra[i + 1] == 1)
                            i++;
                        minDist = graph[startNode][i + 1][0];
                        Dijkstra[startNode] = 1;
                        DijkstraNode = i + 1;
                    }
                }
                else{
                    if(Dijkstra[i] != 1){
                        if (graph[startNode][i][0] < minDist && graph[startNode][i][0] != -1){
                            minDist = graph[startNode][i][0];
                            DijkstraNode = i;
                        }
                    }
                }
            }
            Dijkstra[DijkstraNode] = 1;
            startNode = DijkstraNode;
        }
        //search shortest route from startNode to endNode using Floyd algorithm
        startNode = atoi(argv[2]); endNode = atoi(argv[3]);
        int floydNode = startNode, lengthPathFloyd = 0;
        bool changes = 1;
        float directDistance = graph[startNode][endNode][0], currentDistance = graph[startNode][endNode][0], cumDistance = 0;
        k = 0;
        floydRoute[k] = startNode;
        k++;
        while (changes == 1) {
            changes = 0;
            float tmpCumDistance =  directDistance;
            for (int i = 0; i < NUMBER_OF_CASE; i++){
                currentDistance = graph[startNode][i][0];
                if (tmpCumDistance > currentDistance && startNode != i && endNode != i && caseUsed[i] == 0) {
                    tmpCumDistance = currentDistance;
                    floydNode = i;
                    floydRoute[k] = i;
                }
            }
        }
        k++;
        floydRoute[k] = endNode;
        //end of Floyd algorithm
        k = 1;
        std::cout << argv[2] << "  " << argv[3] << ": ";
        while(floydRoute[k + 1] != MAX_NODE){
            std::cout << floydRoute[k] << "  ";
            k++;
        }
        std::cout << "\n";
      //  templatefile << floydRoute[k - 2];
        templatefile << "\n";
        templatefile.close();
    }

    if( option_N ) {
        std::string filename;
        std::ostringstream strFixCase;
        //strFixCase << argv[argc - 1];

        filename = argv[2];
        std::ifstream infile1( filename.c_str() );
        std::ifstream infile2( filename.c_str() );
        float min = 0, max = 0;
        int sizeOfDataset = 0;

        while (!infile1.eof()) {
           float tmp;
           infile1 >> tmp;
           sizeOfDataset++;
        }
        infile1.close();
        sizeOfDataset--;
        std::cout << "size of data: " << sizeOfDataset << std::endl;

        float *weightValue = new float[sizeOfDataset];
        sizeOfDataset = 0;
        while (!infile2.eof()) {
           infile2 >> weightValue[sizeOfDataset];
           sizeOfDataset++;
        }
        infile2.close();
        sizeOfDataset--;

        for (int i = 0; i < sizeOfDataset; i++) {
            if (i == 0) {
                min = weightValue[i];
                max = min;
            }
            else {
              //  infile >> weightValue[i];
                if (weightValue[i] > max)
                    max = weightValue[i];
                if (weightValue[i] < min)
                    min = weightValue[i];
            }
        }
        for (int i = 0; i < sizeOfDataset; i++){
           if(max - min != 0)
               weightValue[i] = (weightValue[i] - min) / (max - min);
           else
               weightValue[i] = 0;
        }
        filename = argv[3];
        std::ofstream outfile( filename.c_str() , std::ios::app );
        for (int i = 0; i < sizeOfDataset; i++){
            outfile << weightValue[i] << "\n";
        }
        outfile.close();
    }

    if( option_N2 ) {
      std::string filename1, filename2;
        std::ostringstream strFixCase;
        float min = 0, max = 0;
        int sizeOfDataset = 0;

        filename1 = argv[2];
        filename2 = argv[3];
        std::ifstream infile1( filename1.c_str() );
        std::ifstream infile2( filename1.c_str() );
        std::ifstream infile3( filename2.c_str() );
        std::ifstream infile4( filename2.c_str() );

        while (!infile1.eof()) {
           float tmp;
           infile1 >> tmp;
           sizeOfDataset++;
        }
        infile1.close();
        while (!infile3.eof()) {
           float tmp;
           infile3 >> tmp;
           sizeOfDataset++;
        }
        infile3.close();
        sizeOfDataset--;
        std::cout << "size of data: " << sizeOfDataset << std::endl;

        float *weightValue = new float[sizeOfDataset];
        sizeOfDataset = 0;
        while (!infile2.eof()) {
           infile2 >> weightValue[sizeOfDataset];
           sizeOfDataset++;
        }
        infile2.close();
        sizeOfDataset = 0;
        while (!infile4.eof()) {
           infile4 >> weightValue[sizeOfDataset];
           sizeOfDataset++;
        }
        infile4.close();
        sizeOfDataset--;

        for (int i = 0; i < sizeOfDataset; i++) {
            if (i == 0) {
                min = weightValue[i];
                max = min;
            }
            else {
              //  infile >> weightValue[i];
                if (weightValue[i] > max)
                    max = weightValue[i];
                if (weightValue[i] < min)
                    min = weightValue[i];
            }
        }
        for (int i = 0; i < sizeOfDataset; i++){
           if(max - min != 0)
               weightValue[i] = (weightValue[i] - min) / (max - min);
           else
               weightValue[i] = 0;
        }

        filename1 = argv[4];
        std::ofstream outfile( filename1.c_str() , std::ios::app );
        for (int i = 0; i < sizeOfDataset; i++){
            outfile << weightValue[i] << "\n";
        }
        outfile.close();
    }

    if( option_G ) {
        float graph[NUMBER_OF_CASE][NUMBER_OF_CASE][2], circularity[NUMBER_OF_CASE]; //include deformation from two directions
        float intensityE, harmonicE, shapeE, tmpHE = 0, maxDistance = 0, minDistance = 100000;
        int startNode = 0, endNode = 0;
        int floydRoute[NUMBER_OF_CASE];
        std::ostringstream strFixCase;
        std::string filename;

        filename = argv[2];
        std::ifstream iefile( filename.c_str() ); // intensity energy

        filename = argv[3];
        std::ifstream efile( filename.c_str() ); // harmonic energy

        filename = argv[4];

        std::string commandLine;
        commandLine = "rm " + filename;
        system(commandLine.c_str());
        std::ofstream templatefile(filename.c_str(), std::ios::app );
        int cases[25] = {39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 58, 59, 69, 71, 73, 88, 95, 107};
	int caseUsed[NUMBER_OF_CASE];
        int graphHistogram[120] = {0};
    //    gama = 0;

        beta = atof(argv[5]);
        alpha = atof(argv[6]);
        gama = atof(argv[7]);

        for (int i = 0; i < NUMBER_OF_CASE; i++){
            floydRoute[i] = MAX_NODE;
            caseUsed[i] = 0;
        }

        int m = 0;
       /* for (int i = 0; i < NUMBER_OF_CASE; i++){
            sefile >> shapeE;
            circularity[i] = shapeE;
        }
        sefile.close();*/
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            for (int j = 0; j < NUMBER_OF_CASE; j++){
                    graph[i][j][0] = -1;
                    graph[i][j][1] = -1;
            }
        }

        for (int i = 0; i < NUMBER_OF_ATLAS; i++){ //start on the route
            for (int j = 0; j < NUMBER_OF_ATLAS; j++){ // end on the route
                if(i == j){
                    graph[i][j][0] = -1;
                    graph[i][j][1] = -1;
                }
                else {
                    m++;
                    efile >> harmonicE;
                    iefile >> intensityE;
                    graph[j][i][0] = alpha * harmonicE;
                    graph[j][i][0] += beta * intensityE;
               //     std::cout << graph[i][j][0] << std::endl;
               //     graph[j][i][0] += gama * fabs(circularity[i] - circularity[j]) + ENERGY_CONST;
                    if(graph[j][i][0] < minDistance)
                        minDistance = graph[j][i][0];
                    if(graph[j][i][0] > maxDistance)
                        maxDistance = graph[j][i][0];
                    int tmp = round(graph[j][i][0] * 100);
                    if((tmp - 1) < 0)
                        graphHistogram[0]++;
                    else if ((tmp - 1) >= 120)
                        graphHistogram[119]++;
                    else
                        graphHistogram[tmp - 1]++;
             //       invefile >> tmpHE;
               //     graph[i][j][1] = tmpHE;
                }
            }
        }
	//std::cout << m << std::endl;

        for (int i = NUMBER_OF_ATLAS; i < NUMBER_OF_CASE; i++){
            for (int j = 0; j < NUMBER_OF_ATLAS; j++){
                    m++;
                    efile >> harmonicE;
                    iefile >> intensityE;
                    graph[j][i][0] = alpha * harmonicE;
                    graph[j][i][0] += beta * intensityE;
                //    graph[j][i][0] += gama * fabs(circularity[i] - circularity[j]) + ENERGY_CONST;
                    if(graph[j][i][0] < minDistance)
                        minDistance = graph[j][i][0];
                    if(graph[j][i][0] > maxDistance)
                        maxDistance = graph[j][i][0];
                    int tmp = round(graph[j][i][0] * 100);
                    if((tmp - 1) < 0)
                        graphHistogram[0]++;
                    else if ((tmp - 1) >= 120)
                        graphHistogram[119]++;
                    else
                        graphHistogram[tmp - 1]++;
            }
        }
	//std::cout << m << std::endl;

	/*	std::cout << "Graph: " << std::endl;
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            for (int j = 0; j < NUMBER_OF_CASE; j++){
	      std::cout << graph[j][i][0] << "," ;
            }
	    std::cout << std::endl;
	    }*/

        //normalize the distance
	
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            minDistance = 1000000;
            maxDistance = 0;
            for (int j = 0; j < NUMBER_OF_ATLAS; j++){
                if(i == j){
                }
                else {
                    if(graph[j][i][0] < minDistance)
                        minDistance = graph[j][i][0];
                    if(graph[j][i][0] > maxDistance)
                        maxDistance = graph[j][i][0];
                }
            }
            for (int j = 0; j < NUMBER_OF_ATLAS; j++){
                if(i == j){
                }
                else {
                    graph[j][i][0] = (graph[j][i][0] - minDistance) / (maxDistance - minDistance);
                }
            }
	    }
        /*for (int i = NUMBER_OF_ATLAS; i < NUMBER_OF_CASE; i++){
            minDistance = 1000000;
            maxDistance = 0;
            for (int j = 0; j < NUMBER_OF_ATLAS; j++){
                if(graph[j][i][0] < minDistance)
                    minDistance = graph[j][i][0];
                if(graph[j][i][0] > maxDistance)
                    maxDistance = graph[j][i][0];
            }
            for (int j = 0; j < NUMBER_OF_ATLAS; j++){
                graph[j][i][0] = (graph[j][i][0] - minDistance) / (maxDistance - minDistance);
            }
        }*/
        //end of normalization
        efile.close();
//        invefile.close();
        iefile.close();


        //search shortest route from startNode to endNode using Floyd algorithm
        //startNode = atoi(argv[5]);
        //endNode = atoi(argv[6]);
        endNode = NUMBER_OF_ATLAS;
        for (int i = 0; i < NUMBER_OF_ATLAS; i++){
             startNode = i;
             int floydNode = startNode, lengthPathFloyd = 0;
             bool changes = 1;
             float directDistance = graph[startNode][endNode][0], currentDistance = graph[startNode][endNode][0], cumDistance = 0;
             int k = 0;
             for (int i = 0; i < NUMBER_OF_CASE; i++){
                 floydRoute[i] = MAX_NODE;
                 caseUsed[i] = 0;
             }

             floydRoute[k] = startNode;
             while (changes == 1) {
                 changes = 0;
                 float tmpCumDistance =  directDistance;
                 for (int i = 0; i < NUMBER_OF_CASE; i++){
                     currentDistance = cumDistance + graph[startNode][i][0] + graph[i][endNode][0];
                     if (tmpCumDistance > currentDistance && startNode != i && endNode != i && caseUsed[i] == 0) {
                         changes = 1;
                         floydNode = i;
                         tmpCumDistance = currentDistance - graph[i][endNode][0];
                     }
                 }
                 if(changes == 1) {
                     cumDistance = tmpCumDistance;
                     startNode = floydNode;
                     k++;
                     floydRoute[k] = startNode;
                     caseUsed[floydNode] = 1;
                     lengthPathFloyd++;
                 }
             }
             k++;
             floydRoute[k] = endNode;
             //end of Floyd algorithm
             k = 0;
             while(floydRoute[k] != MAX_NODE){
                 std::cout << floydRoute[k] << "  ";
                 k++;
             }
             std::cout << "\n";
             templatefile << floydRoute[k - 2];
             templatefile << "\n";
        }
        templatefile.close();
    }

    if( option_M ) {
        std::cout << "running majority voting" << std::endl;
        std::ostringstream strFixCase;

        std::string filename;
        filename = argv[4];
        std::ifstream templatefile( filename.c_str());
	std::string outfolder (argv[5]);

        std::string commandLine;
        int cases[NUMBER_OF_CASE];
        int onlyOneAtlas = 0;
        gama = 0;

        beta = atof(argv[argc - 5]);
        alpha = atof(argv[argc - 4]);
        gama = atof(argv[argc - 3]);

        DIR *dir;
        struct dirent *ent;
        int sizeLabelList = 0;
        char is_a_label[5];

        if ((dir = opendir (outfolder.c_str())) != NULL) {
            std::string tmpFilename;
            while ((ent = readdir (dir)) != NULL) {
                tmpFilename = ent->d_name;
                tmpFilename.copy(is_a_label, 5, 0);
		is_a_label[5] = '\0';
		//std::cout << tmpFilename << "," << is_a_label << "," <<  argv[8] << std::endl;
                if(tmpFilename.at(0) == '.')    // skip . and ..
                    continue;
                else{
		  if (!strcmp(is_a_label, argv[8])) {
		    //std::cout << tmpFilename << "," << is_a_label << std::endl;
		    sizeLabelList++;
		  }
                }
            }
        }
        closedir (dir);
	if (sizeLabelList == 0) {
	  std::cout << "no label/parcellation files found, exiting" << std::endl;
	  exit(-1);
	}
        std::string *labelList = new std::string[sizeLabelList];

        if ((dir = opendir (outfolder.c_str())) != NULL) {
            std::string tmpFilename;
            int i = 0;
            while ((ent = readdir (dir)) != NULL) {
                tmpFilename = ent->d_name;
                tmpFilename.copy(is_a_label, 5, 0);
		is_a_label[5] = '\0';
                if(tmpFilename.at(0) == '.')    // skip . and ..
                    continue;
                else{
                    if (!strcmp(is_a_label, argv[8])){
                        labelList[i] = tmpFilename;
                        i++;
                    }
                }
            }
            SortStringList(labelList, sizeLabelList);
        }
        closedir (dir);

        bool *caseFlag = new bool[NUMBER_OF_CASE];
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            caseFlag[i] = 1;
        }

        commandLine = "ImageMath " ;
        commandLine += outfolder;
        commandLine += labelList[0].c_str() ;
        commandLine += " -majorityVoting ";
        int k = 0;
        for (int i = 0; i < NUMBER_OF_ATLAS; i++){
          commandLine += outfolder + labelList[i] + " ";
        }
        commandLine += "-outfile " ;
        commandLine += argv[6] ;
        commandLine += argv[7];
        std::cout << commandLine.c_str() << std::endl;
        system(commandLine.c_str());
    }

    if( option_V ) {
        std::cout << "running weighted majority voting" << std::endl;
        //run weighted majority voting
        std::ostringstream strFixCase;
      //  strFixCase << argv[argc - 1];

        std::string filename;
        filename = argv[2];
        std::ifstream iefile( filename.c_str() );
        filename = argv[3];
        std::ifstream efile( filename.c_str() );
        filename = argv[4];
        std::ifstream templatefile( filename.c_str());
	std::string labelfolder (argv[5]);

        std::string commandLine;
        float harmonicE, intensityE, shapeE, weightFactor[NUMBER_OF_CASE * (NUMBER_OF_CASE - 1)], circularity[NUMBER_OF_CASE];
        int cases[NUMBER_OF_CASE];// = {1000, 1001, 1002, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1017, 1036, 1003};
        int onlyOneAtlas = 0;
        gama = 0;

        beta = atof(argv[argc - 5]);
        alpha = atof(argv[argc - 4]);
        gama = atof(argv[argc - 3]);

        std::cout << "weights: " << beta << "  " << alpha << "  " << gama << std::endl;
        std::cout << "labelfolder: " << labelfolder << std::endl;

	//std::cout << argv[8] << std::endl;

        DIR *dir;
        struct dirent *ent;
        int sizeLabelList = 0;
        char is_a_label[6];

        if ((dir = opendir (labelfolder.c_str())) != NULL) {
            std::string tmpFilename;
            while ((ent = readdir (dir)) != NULL) {
                tmpFilename = ent->d_name;
                tmpFilename.copy(is_a_label, 5, 0);
		is_a_label[5] = '\0';
		//		std::cout << tmpFilename << "," << is_a_label << std::endl;
                if(tmpFilename.at(0) == '.')    // skip . and ..
                    continue;
                else{
		  if (!strcmp(is_a_label, argv[8])) {
		    //std::cout << tmpFilename << "," << is_a_label << std::endl;
		    sizeLabelList++;
		  }
                }
            }
        }
        closedir (dir);

	if (sizeLabelList == 0) {
	  std::cout << "no label/parcellation files found, exiting" << std::endl;
	  exit(-1);
	}
        std::cout << "size of label list: " << sizeLabelList << std::endl;
        std::string *labelList = new std::string[sizeLabelList];

        if ((dir = opendir (labelfolder.c_str())) != NULL) {
            std::string tmpFilename;
            int i = 0;
            while ((ent = readdir (dir)) != NULL) {
                tmpFilename = ent->d_name;
                tmpFilename.copy(is_a_label, 5, 0);
		is_a_label[5] = '\0';
                if(tmpFilename.at(0) == '.')    // skip . and ..
                    continue;
                else{
                    if (!strcmp(is_a_label, argv[8])){
                        labelList[i] = tmpFilename;
                        i++;
                    }
                }
            }
            SortStringList(labelList, sizeLabelList);
        }
        closedir (dir);

        bool *caseFlag = new bool[NUMBER_OF_CASE];
        for (int i = 0; i < NUMBER_OF_CASE; i++)
            caseFlag[i] = 0;


        for (int i = 0; i < NUMBER_OF_CASE; i++){
            int temp;
            templatefile >> temp;
            caseFlag[temp] = 1;
        }
        templatefile.close();
        for (int i = 0; i < NUMBER_OF_CASE; i++){
            cases[i] = 1;
            onlyOneAtlas++;
            if (i != (NUMBER_OF_CASE - 1) && caseFlag[i] == 0) {
                cases[i] = 0;
                onlyOneAtlas--;
            }
            if (i == (NUMBER_OF_CASE - 1) )
                onlyOneAtlas--;
        }
        int m = 0;
        float minDistance = 1000000, maxDistance = 0;
        for ( int i = 0; i < (NUMBER_OF_ATLAS * (NUMBER_OF_ATLAS - 1) + NUMBER_OF_ATLAS); i++){
            efile >> harmonicE;
            iefile >> intensityE;
            if( i >= (NUMBER_OF_ATLAS * (NUMBER_OF_ATLAS - 1)) ){
                weightFactor[m] = alpha * harmonicE + beta * intensityE;
                if(weightFactor[m] < minDistance)
                    minDistance = weightFactor[m];
                if(weightFactor[m] > maxDistance)
                    maxDistance = weightFactor[m];
                m++;
            }
        }
        for (int i = 0; i < NUMBER_OF_ATLAS; i++){
            weightFactor[i] = (weightFactor[i] - minDistance) / (maxDistance - minDistance);
        }
        efile.close();
        iefile.close();

        commandLine = "ImageMath " ;
        commandLine += labelfolder;
        commandLine += labelList[0] ;
        commandLine += " -weightedMajorityVoting ";
        int k = 0;
        for (int i = 0; i < NUMBER_OF_ATLAS; i++){
            if (cases[i] != 0) {
                commandLine += labelfolder + labelList[i] + " ";
            }
        }
        commandLine += "-weights ";
        bool firstWeightFlag = 0;
        k = 0;
        for (int i = 0; i < NUMBER_OF_ATLAS; i++){
            std::ostringstream strWFactor;
            if(cases[i] != 0) {
                if (firstWeightFlag == 1 && k > 0){
                    commandLine += ",";
                }
                if (onlyOneAtlas == 1)
                    strWFactor << 1;
                else
                    strWFactor << 1 - weightFactor[k];
                commandLine += strWFactor.str();
                firstWeightFlag = 1;
                k++;
            }
        }
        commandLine += " -outfile " ;
        commandLine += argv[6] ;
        commandLine += argv[7];
        std::cout << commandLine.c_str() << std::endl;
        system(commandLine.c_str());

    }

    if( option_P ) {
        typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType >  MSMMetricType;
        typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType >  NCCMetricType;
        typedef itk::MutualInformationImageToImageMetric< ImageType, ImageType >  MIMetricType;
        typedef itk::NormalizedMutualInformationHistogramImageToImageMetric< ImageType, ImageType > NMIMetricType;
        typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;

        MSMMetricType::Pointer MSMmetric = MSMMetricType::New();
        NCCMetricType::Pointer NCCmetric = NCCMetricType::New();
        MIMetricType::Pointer MImetric = MIMetricType::New();
        NMIMetricType::Pointer NMImetric = NMIMetricType::New();

        RescaleFilterType::Pointer rescaleFixFilter = RescaleFilterType::New();
        RescaleFilterType::Pointer rescaleMoveFilter = RescaleFilterType::New();

        rescaleFixFilter->SetInput( fixedReader->GetOutput() );
        rescaleFixFilter->SetOutputMinimum( 0 );
        rescaleFixFilter->SetOutputMaximum( 1023 );
        rescaleFixFilter->Update();

        rescaleMoveFilter->SetInput( movingReader->GetOutput() );
        rescaleMoveFilter->SetOutputMinimum( 0 );
        rescaleMoveFilter->SetOutputMaximum( 1023 );
        rescaleMoveFilter->Update();

        ImageType::ConstPointer fixedImage  = rescaleFixFilter->GetOutput();
        ImageType::ConstPointer movingImage = rescaleMoveFilter->GetOutput();

        typedef itk::AffineTransform<double,3> AffineTransformType;
        AffineTransformType::Pointer transform = AffineTransformType::New();

        typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >  NNInterpolatorType;
        typedef itk::LinearInterpolateImageFunction< ImageType > LinearInterpolatorType;
        LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
        transform->SetIdentity();

        // mean square metric
        MSMmetric->SetTransform( transform );
        MSMmetric->SetInterpolator( interpolator );
        MSMmetric->SetFixedImage(  fixedImage  );
        MSMmetric->SetMovingImage( movingImage );
        // normalized cross correlation
        NCCmetric->SetTransform( transform );
        NCCmetric->SetInterpolator( interpolator );
        NCCmetric->SetFixedImage(  fixedImage  );
        NCCmetric->SetMovingImage( movingImage );
        // mutual information
        MImetric->SetTransform( transform );
        MImetric->SetInterpolator( interpolator );
        MImetric->SetFixedImage(  fixedImage  );
        MImetric->SetMovingImage( movingImage );
        // normalized mutual information
        NMImetric->SetTransform( transform );
        NMImetric->SetInterpolator( interpolator );
        NMImetric->SetFixedImage(  fixedImage  );
        NMImetric->SetMovingImage( movingImage );

        MSMmetric->SetFixedImageRegion( fixedImage->GetLargestPossibleRegion( ) );
        try {
            MSMmetric->Initialize();
        }
        catch( itk::ExceptionObject & excep ){
            std::cerr << "Exception catched !" << std::endl;
            std::cerr << excep << std::endl;
            return EXIT_FAILURE;
        }
        float msm;
	msm = MSMmetric->GetValue( transform->GetParameters( ) );
	//msm = NCCmetric->GetValue( transform->GetParameters( ) );
        std::string filename;
        //std::ostringstream strFixCase;
        filename = argv[4];
        std::ofstream efile( filename.c_str() , std::ios::app );
        efile << msm << "\n";
        efile.close();
        std::cout << "mean square metric value:            " << msm << std::endl;
    }

    return EXIT_SUCCESS;
}


