/*=========================================================================
 *
 * Program :   $Insight Segmentation & Registration Toolkit $
 * Module  :   $DMDBiomarkerTool: DMDMuscleFeature.h $
 * Purpose :   $The class of calculating features in segmented muscle regions $
 * Language:   $C++ $
 * Date    :   $Date: 2010-06-24 16:35:04 $
 * Version :   $Revision: 0.10 $
 * Authors :   $Jiahui Wang, Martin Styner $
 * Update  :   1. Created the class (06-24-10)
 * Copyright (c) Neuro Image Research and Analysis Lab.@UNC All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.
 *
=========================================================================*/

#ifndef DMDMUSCLEFEAT
#define DMDMUSCLEFEAT    // prevent for redefine
#include "DMDData.h"
#include <iostream>
#include <fstream>
//#include "itkScalarImageToCooccurrenceMatrixFilter.h"

#define HISTOGRAM_BIN_SIZE 1000
#define RUN_INTENSITY_LEVEL 64
#define RUN_LENGTH_LEVEL 64
#define INTENSITY_RESCALE_FACTOR 64
#define SEARCHED_PIXEL_VALUE -1
#define BACKGROUND 0
#define SIZE_LABELING_BIN 10000
#define HISTOGRAM_LIMIT 10000
#define INTENSITY_CUT_RATE 0.05
#define NUMBER_MUSCLE 11
#define ROI_SIZE 32
#define CO_OCCURRENCE_LEVEL 64

char *muscleName[] = { "MUSCLES", "CRAN_SART", "REC_FEM", "BIC_FEM", "GRACILIS", "SEMI_MEM", "ADDUCTOR", "SEMIT", "VAS_LAT", "VAS_INT", "VAS_MED", "CAUD_SART"};
//////////////////////////////////////////////////////////////////////////
class DMDMuscleFeature : public virtual DMDData
{
    private:
//        float meanCran, meanRect, meanSemit, meanAdductor, meanGracilis, meanBic, meanSemimem, meanVaslat, meanVasmed, meanVasint, meanCaud, sdCran, sdRect, sdSemit, sdAdductor, sdGracilis, sdBic, sdSemimem, sdVaslat, sdVasmed, sdVasint, sdCaud, volCran, volRect, volSemit, volAdductor, volGracilis, volBic, volSemimem, volVaslat, volVasmed, volVasint, volCaud ;
        float meanCran, meanRect, meanSemit, meanAdductor, meanGracilis, meanBic, meanSemimem, meanVaslat, meanVasmed, meanVasint, meanCaud, sdCran, sdRect, sdSemit, sdAdductor, sdGracilis, sdBic, sdSemimem, sdVaslat, sdVasmed, sdVasint, sdCaud, skewCran, skewRect, skewSemit, skewAdductor, skewGracilis, skewBic, skewSemimem, skewVaslat, skewVasmed, skewVasint, skewCaud, kortCran, kortRect, kortSemit, kortAdductor, kortGracilis, kortBic, kortSemimem, kortVaslat, kortVasmed, kortVasint, kortCaud, volCran, volRect, volSemit, volAdductor, volGracilis, volBic, volSemimem, volVaslat, volVasmed, volVasint, volCaud, entropyCran, entropyRect, entropySemit, entropyAdductor, entropyGracilis, entropyBic, entropySemimem, entropyVaslat, entropyVasmed, entropyVasint, entropyCaud;
        float muscleVolume[NUMBER_OF_MUSCLE];
    public:
       // DMDMuscleFeature():meanCran(0), meanRect(0), meanSemit(0), meanAdductor(0), meanGracilis(0), meanBic(0), meanSemimem(0), meanVaslat(0), meanVasmed(0), meanVasint(0), meanCaud(0), sdCran(0), sdRect(0), sdSemit(0), sdAdductor(0), sdGracilis(0), sdBic(0), sdSemimem(0), sdVaslat(0), sdVasmed(0), sdVasint(0), sdCaud(0), volCran(0), volRect(0), volSemit(0), volAdductor(0), volGracilis(0), volBic(0), volSemimem(0), volVaslat(0), volVasmed(0), volVasint(0), volCaud(0) {
        DMDMuscleFeature():meanCran(0), meanRect(0), meanSemit(0), meanAdductor(0), meanGracilis(0), meanBic(0), meanSemimem(0), meanVaslat(0), meanVasmed(0), meanVasint(0), meanCaud(0), sdCran(0), sdRect(0), sdSemit(0), sdAdductor(0), sdGracilis(0), sdBic(0), sdSemimem(0), sdVaslat(0), sdVasmed(0), sdVasint(0), sdCaud(0), skewCran(0), skewRect(0), skewSemit(0), skewAdductor(0), skewGracilis(0), skewBic(0), skewSemimem(0), skewVaslat(0), skewVasmed(0), skewVasint(0), skewCaud(0), kortCran(0), kortRect(0), kortSemit(0), kortAdductor(0), kortGracilis(0), kortBic(0), kortSemimem(0), kortVaslat(0), kortVasmed(0), kortVasint(0), kortCaud(0), volCran(0), volRect(0), volSemit(0), volAdductor(0), volGracilis(0), volBic(0), volSemimem(0), volVaslat(0), volVasmed(0), volVasint(0), volCaud(0), entropyCran(0), entropyRect(0), entropySemit(0), entropyAdductor(0), entropyGracilis(0), entropyBic(0), entropySemimem(0), entropyVaslat(0), entropyVasmed(0), entropyVasint(0), entropyCaud(0) {

	}
        typedef itk::ImageRegionConstIterator< DMDData::OrientedImageType > ConstIteratorType;
        void features( DMDData::OrientedImageType::Pointer mask, DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, float voxelSize, std::string featurefilename, bool );
        void runlengthFeat( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename);
        void runlengthFeat3D( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType);
        void runlengthFeat2DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType);
        void runlengthFeat3DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType, StrITKType);
        void calRunLengthFeatures(DMDData::FITKType &SRE, DMDData::FITKType &LRE, DMDData::FITKType &GLN, DMDData::FITKType &RLN, DMDData::FITKType &RP, FITKType np, int runlengthmatrix[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], StrITKType);
        void cooccurrenceFeat3DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType, StrITKType);
        void calCooccurrenceFeatures( DMDData::FITKType &EntropyTmp, DMDData::FITKType &EnergyTmp, DMDData::FITKType &ContrastTmp, DMDData::FITKType &HomoGeneityTmp, float  cooccurrencematrix[CO_OCCURRENCE_LEVEL][CO_OCCURRENCE_LEVEL], std::string);

};
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::calRunLengthFeatures(DMDData::FITKType &SRE, DMDData::FITKType &LRE, DMDData::FITKType &GLN, DMDData::FITKType &RLN, DMDData::FITKType &RP, FITKType np, int runlengthmatrix[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], StrITKType caseID)
{
    // calculate run length matrix features
    float nr = 0, pr = 0, pg = 0;  // pr: number of runs with same length but different intensity; pg: number of runs with same intensity but different length; nr: total number of runs.
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrix[j][i];
             nr += runlengthmatrix[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrix[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (np != 0)
        RP = nr / np;
    else
        RP = 0;
   // std::cout << "number of run: " << nr << std::endl;
    std::string featurefilename = "../data/voi" + caseID + ".txt" ;
    std::ofstream efile( featurefilename.c_str() , std::ios::app );

    efile << SRE << "  " << LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";

    efile.close();
}
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::calCooccurrenceFeatures( DMDData::FITKType &EntropyTmp, DMDData::FITKType &EnergyTmp, DMDData::FITKType &ContrastTmp, DMDData::FITKType &HomoGeneityTmp, float  cooccurrencematrix[CO_OCCURRENCE_LEVEL][CO_OCCURRENCE_LEVEL], std::string caseID )
{
    // calculate run length matrix features
    EntropyTmp = 0;
    EnergyTmp = 0;
    ContrastTmp = 0;
    HomoGeneityTmp = 0;
    for (int j = 0; j < CO_OCCURRENCE_LEVEL; j++) {
        for (int i = 0; i < CO_OCCURRENCE_LEVEL; i++) {
           if(cooccurrencematrix[j][i]){
               cooccurrencematrix[j][i] /= (CO_OCCURRENCE_LEVEL * CO_OCCURRENCE_LEVEL);
               EntropyTmp += (-cooccurrencematrix[j][i] * log (cooccurrencematrix[j][i]));
            //   std::cout << EnergyTmp << std::endl;
               EnergyTmp += (cooccurrencematrix[j][i] * cooccurrencematrix[j][i]);
           //    if(EnergyTmp < 0) {
               //    std::cout << EnergyTmp << "  " << cooccurrencematrix[j][i] << std::endl;
                //   getchar();
          //     }

               ContrastTmp += ((i - j) * (i - j) * cooccurrencematrix[j][i]);
               HomoGeneityTmp += (cooccurrencematrix[j][i] / (1 + abs(i - j)));
           }
        }
    }
    // std::cout << "number of run: " << nr << std::endl;
    std::string featurefilename = "../data/voiCOM_" + caseID + ".txt" ;
    std::ofstream efile( featurefilename.c_str() , std::ios::app );
    efile << EntropyTmp << "  " << EnergyTmp << "  " << ContrastTmp << "  " << HomoGeneityTmp << "\n";
    efile.close();
}

//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::runlengthFeat( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename)
{
    int runlengthmatrixRecFem[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixSemit[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixCranSart[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixBicFem[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixGracilis[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixAdductor[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixSemiMem[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixVasLat[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixVasMed[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixVasInt[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixCaudSart[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL];
    int tempIntensity = 0, tempLength = 0;
    DMDData ddata;
    DMDData::OrientedImageType::Pointer rescaledData = DMDData::OrientedImageType::New(), tempMaskX = DMDData::OrientedImageType::New(), tempMaskY = DMDData::OrientedImageType::New(), tempMaskZ = DMDData::OrientedImageType::New();
    float SRE = 0, LRE = 0, GLN = 0, RLN = 0, RP = 0;  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float npRecFem = 0, npSemit = 0, npCranSart = 0, npBicFem = 0, npGracilis = 0, npAdductor = 0, npSemiMem = 0, npVasLat = 0, npVasMed = 0, npVasInt = 0, npCaudSart = 0;

    for (int j = 0; j < RUN_LENGTH_LEVEL; j++) {
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
            runlengthmatrixRecFem[j][i] = 0;
            runlengthmatrixSemit[j][i] = 0;
            runlengthmatrixCranSart[j][i] = 0;
            runlengthmatrixBicFem[j][i] = 0;
            runlengthmatrixGracilis[j][i] = 0;
            runlengthmatrixAdductor[j][i] = 0;
            runlengthmatrixSemiMem[j][i] = 0;
            runlengthmatrixVasLat[j][i] = 0;
            runlengthmatrixVasMed[j][i] = 0;
            runlengthmatrixVasInt[j][i] = 0;
            runlengthmatrixCaudSart[j][i] = 0;
        }
    }

    ddata.imageInitialize ( data, rescaledData );
    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;
    IteratorType ReScaledDataIterator( rescaledData, rescaledData->GetRequestedRegion() ) ;

    //rescale images
    for ( constDataIterator.GoToBegin(), ReScaledDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator,++ReScaledDataIterator ) {
        tempIntensity = constDataIterator.Get() / INTENSITY_RESCALE_FACTOR;
        ReScaledDataIterator.Set(tempIntensity);
    }

    //establish run length matrix
    DMDData::OrientedImageType::SizeType   size = data->GetLargestPossibleRegion().GetSize();
    DMDData::OrientedImageType::IndexType  index;
    int startRunRecFem = -1, runLengthRecFem = 0, startRunSemit = -1, runLengthSemit = 0, startRunCranSart = -1, runLengthCranSart = 0, startRunBicFem = -1, runLengthBicFem = 0, startRunGracilis = -1, runLengthGracilis = 0, startRunAdductor = -1, runLengthAdductor = 0, startRunSemiMem = -1, runLengthSemiMem = 0, startRunVasLat = -1, runLengthVasLat = 0, startRunVasMed = -1, runLengthVasMed = 0, startRunVasInt = -1, runLengthVasInt = 0, startRunCaudSart = -1, runLengthCaudSart = 0;
    //x direction
    for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
        for(index[1] = 0; index[1] < size[1]; index[1]++) {
            for(index[0] = 0; index[0] < size[0]; index[0]++) {
                ddata.connectedComponentRegionGrowing (rescaledData->GetPixel(index),rescaledData->GetPixel(index), rescaledData, index);
                if ( erodemask->GetPixel(index) == REC_FEM ) {
                    npRecFem++;
                    if (startRunRecFem == -1){
                        startRunRecFem = rescaledData->GetPixel(index);
                        runLengthRecFem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunRecFem){
                            runLengthRecFem++;
                        }
                        else {
                            if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                                runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                            if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                                startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                            startRunRecFem = -1;
                            runLengthRecFem = 0;
                        }
                    }
                }
                else {
                    if (startRunRecFem != -1){
                        if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                            runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                        if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                            startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                        startRunRecFem = -1;
                        runLengthRecFem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == SEMIT ) {
                    npSemit++;
                    if (startRunSemit == -1){
                        startRunSemit = rescaledData->GetPixel(index);
                        runLengthSemit++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunSemit){
                            runLengthSemit++;
                        }
                        else {
                            if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                                runLengthSemit = RUN_LENGTH_LEVEL - 1;
                            if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                                startRunSemit = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                            startRunSemit = -1;
                            runLengthSemit = 0;
                        }
                    }
                }
                else {
                    if (startRunSemit != -1){
                        if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                            runLengthSemit = RUN_LENGTH_LEVEL - 1;
                        if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                            startRunSemit = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                        startRunSemit = -1;
                        runLengthSemit = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == CRAN_SART ) {
                    npCranSart++;
                    if (startRunCranSart == -1){
                        startRunCranSart = rescaledData->GetPixel(index);
                        runLengthCranSart++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunCranSart){
                            runLengthCranSart++;
                        }
                        else {
                            if(runLengthCranSart > RUN_LENGTH_LEVEL - 1)
                                runLengthCranSart = RUN_LENGTH_LEVEL - 1;
                            if(startRunCranSart > RUN_INTENSITY_LEVEL - 1)
                                startRunCranSart = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixCranSart[runLengthCranSart][startRunCranSart]++;
                            startRunCranSart = -1;
                            runLengthCranSart = 0;
                        }
                    }
                }
                else {
                    if (startRunCaudSart != -1){
                        if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                            runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                        if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                            startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                        startRunCaudSart = -1;
                        runLengthCaudSart = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == BIC_FEM ) {
                    npBicFem++;
                    if (startRunBicFem == -1){
                        startRunBicFem = rescaledData->GetPixel(index);
                        runLengthBicFem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunBicFem){
                            runLengthBicFem++;
                        }
                        else {
                            if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                                runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                            if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                                startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;
                            startRunBicFem = -1;
                            runLengthBicFem = 0;
                        }
                    }
                }
                else {
                    if (startRunBicFem != -1){
                        if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                            runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                        if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                            startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;
                        startRunBicFem = -1;
                        runLengthBicFem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == GRACILIS ) {
                    npGracilis++;
                    if (startRunGracilis == -1){
                        startRunGracilis = rescaledData->GetPixel(index);
                        runLengthGracilis++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunGracilis){
                            runLengthGracilis++;
                        }
                        else {
                            if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                                runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                            if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                                startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;
                            startRunGracilis = -1;
                            runLengthGracilis = 0;
                        }
                    }
                }
                else {
                    if (startRunGracilis != -1){
                        if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                            runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                        if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                            startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;
                        startRunGracilis = -1;
                        runLengthGracilis = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == ADDUCTOR ) {
                    npAdductor++;
                    if (startRunAdductor == -1){
                        startRunAdductor = rescaledData->GetPixel(index);
                        runLengthAdductor++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunAdductor){
                            runLengthAdductor++;
                        }
                        else {
                            if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                                runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                            if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                                startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                            startRunAdductor = -1;
                            runLengthAdductor = 0;
                        }
                    }
                }
                else {
                    if (startRunAdductor != -1){
                        if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                            runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                        if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                            startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                        startRunAdductor = -1;
                        runLengthAdductor = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == SEMI_MEM ) {
                    npSemiMem++;
                    if (startRunSemiMem == -1){
                        startRunSemiMem = rescaledData->GetPixel(index);
                        runLengthSemiMem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunSemiMem){
                            runLengthSemiMem++;
                        }
                        else {
                            if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                                runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                            if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                                startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                            startRunSemiMem = -1;
                            runLengthSemiMem = 0;
                        }
                    }
                }
                else {
                    if (startRunSemiMem != -1){
                        if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                            runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                        if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                            startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                        startRunSemiMem = -1;
                        runLengthSemiMem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_LAT ) {
                    npVasLat++;
                    if (startRunVasLat == -1){
                        startRunVasLat = rescaledData->GetPixel(index);
                        runLengthVasLat++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasLat){
                            runLengthVasLat++;
                        }
                        else {
                            if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                                runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                                startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;
                            startRunVasLat = -1;
                            runLengthVasLat = 0;
                        }
                    }
                }
                else {
                    if (startRunVasLat != -1){
                        if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                            runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                            startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;
                        startRunVasLat = -1;
                        runLengthVasLat = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_MED ) {
                    npVasMed++;
                    if (startRunVasMed == -1){
                        startRunVasMed = rescaledData->GetPixel(index);
                        runLengthVasMed++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasMed){
                            runLengthVasMed++;
                        }
                        else {
                            if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                                runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                                startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                            startRunVasMed = -1;
                            runLengthVasMed = 0;
                        }
                    }
                }
                else {
                    if (startRunVasMed != -1){
                        if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                            runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                            startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                        startRunVasMed = -1;
                        runLengthVasMed = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_INT ) {
                    npVasInt++;
                    if (startRunVasInt == -1){
                        startRunVasInt = rescaledData->GetPixel(index);
                        runLengthVasInt++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasInt){
                            runLengthVasInt++;
                        }
                        else {
                            if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                                runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                                startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                            startRunVasInt = -1;
                            runLengthVasInt = 0;
                        }
                    }
                }
                else {
                    if (startRunVasInt != -1){
                        if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                            runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                            startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                        startRunVasInt = -1;
                        runLengthVasInt = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == CAUD_SART ) {
                    npCaudSart++;
                    if (startRunCaudSart == -1){
                        startRunCaudSart = rescaledData->GetPixel(index);
                        runLengthCaudSart++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunCaudSart){
                            runLengthCaudSart++;
                        }
                        else {
                            if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                                runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                            if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                                startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                            startRunCaudSart = -1;
                            runLengthCaudSart = 0;
                        }
                    }
                }
                else {
                    if (startRunCaudSart != -1){
                        if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                            runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                        if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                            startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                        startRunCaudSart = -1;
                        runLengthCaudSart = 0;
                    }
                }
            }
        }
    }
    //y direction
    for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
        for(index[0] = 0; index[0] < size[0]; index[0]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                if ( erodemask->GetPixel(index) == REC_FEM ) {
                    if (startRunRecFem == -1){
                        startRunRecFem = rescaledData->GetPixel(index);
                        runLengthRecFem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunRecFem){
                            runLengthRecFem++;
                        }
                        else {
                            if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                                runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                            if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                                startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                            startRunRecFem = -1;
                            runLengthRecFem = 0;
                        }
                    }
                }
                else {
                    if (startRunRecFem != -1){
                        if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                            runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                        if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                            startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                        startRunRecFem = -1;
                        runLengthRecFem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == SEMIT ) {
                    if (startRunSemit == -1){
                        startRunSemit = rescaledData->GetPixel(index);
                        runLengthSemit++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunSemit){
                            runLengthSemit++;
                        }
                        else {
                            if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                                runLengthSemit = RUN_LENGTH_LEVEL - 1;
                            if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                                startRunSemit = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                            startRunSemit = -1;
                            runLengthSemit = 0;
                        }
                    }
                }
                else {
                    if (startRunSemit != -1){
                        if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                            runLengthSemit = RUN_LENGTH_LEVEL - 1;
                        if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                            startRunSemit = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                        startRunSemit = -1;
                        runLengthSemit = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == CRAN_SART ) {
                    if (startRunCranSart == -1){
                        startRunCranSart = rescaledData->GetPixel(index);
                        runLengthCranSart++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunCranSart){
                            runLengthCranSart++;
                        }
                        else {
                            if(runLengthCranSart > RUN_LENGTH_LEVEL - 1)
                                runLengthCranSart = RUN_LENGTH_LEVEL - 1;
                            if(startRunCranSart > RUN_INTENSITY_LEVEL - 1)
                                startRunCranSart = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixCranSart[runLengthCranSart][startRunCranSart]++;
                            startRunCranSart = -1;
                            runLengthCranSart = 0;
                        }
                    }
                }
                else {
                    if (startRunCaudSart != -1){
                        if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                            runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                        if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                            startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                        startRunCaudSart = -1;
                        runLengthCaudSart = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == BIC_FEM ) {
                    if (startRunBicFem == -1){
                        startRunBicFem = rescaledData->GetPixel(index);
                        runLengthBicFem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunBicFem){
                            runLengthBicFem++;
                        }
                        else {
                            if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                                runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                            if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                                startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;
                            startRunBicFem = -1;
                            runLengthBicFem = 0;
                        }
                    }
                }
                else {
                    if (startRunBicFem != -1){
                        if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                            runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                        if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                            startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;
                        startRunBicFem = -1;
                        runLengthBicFem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == GRACILIS ) {
                    if (startRunGracilis == -1){
                        startRunGracilis = rescaledData->GetPixel(index);
                        runLengthGracilis++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunGracilis){
                            runLengthGracilis++;
                        }
                        else {
                            if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                                runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                            if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                                startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;
                            startRunGracilis = -1;
                            runLengthGracilis = 0;
                        }
                    }
                }
                else {
                    if (startRunGracilis != -1){
                        if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                            runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                        if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                            startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;
                        startRunGracilis = -1;
                        runLengthGracilis = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == ADDUCTOR ) {
                    if (startRunAdductor == -1){
                        startRunAdductor = rescaledData->GetPixel(index);
                        runLengthAdductor++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunAdductor){
                            runLengthAdductor++;
                        }
                        else {
                            if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                                runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                            if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                                startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                            startRunAdductor = -1;
                            runLengthAdductor = 0;
                        }
                    }
                }
                else {
                    if (startRunAdductor != -1){
                        if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                            runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                        if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                            startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                        startRunAdductor = -1;
                        runLengthAdductor = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == SEMI_MEM ) {
                    if (startRunSemiMem == -1){
                        startRunSemiMem = rescaledData->GetPixel(index);
                        runLengthSemiMem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunSemiMem){
                            runLengthSemiMem++;
                        }
                        else {
                            if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                                runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                            if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                                startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                            startRunSemiMem = -1;
                            runLengthSemiMem = 0;
                        }
                    }
                }
                else {
                    if (startRunSemiMem != -1){
                        if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                            runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                        if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                            startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                        startRunSemiMem = -1;
                        runLengthSemiMem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_LAT ) {
                    if (startRunVasLat == -1){
                        startRunVasLat = rescaledData->GetPixel(index);
                        runLengthVasLat++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasLat){
                            runLengthVasLat++;
                        }
                        else {
                            if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                                runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                                startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;
                            startRunVasLat = -1;
                            runLengthVasLat = 0;
                        }
                    }
                }
                else {
                    if (startRunVasLat != -1){
                        if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                            runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                            startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;
                        startRunVasLat = -1;
                        runLengthVasLat = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_MED ) {
                    if (startRunVasMed == -1){
                        startRunVasMed = rescaledData->GetPixel(index);
                        runLengthVasMed++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasMed){
                            runLengthVasMed++;
                        }
                        else {
                            if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                                runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                                startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                            startRunVasMed = -1;
                            runLengthVasMed = 0;
                        }
                    }
                }
                else {
                    if (startRunVasMed != -1){
                        if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                            runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                            startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                        startRunVasMed = -1;
                        runLengthVasMed = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_INT ) {
                    if (startRunVasInt == -1){
                        startRunVasInt = rescaledData->GetPixel(index);
                        runLengthVasInt++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasInt){
                            runLengthVasInt++;
                        }
                        else {
                            if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                                runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                                startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                            startRunVasInt = -1;
                            runLengthVasInt = 0;
                        }
                    }
                }
                else {
                    if (startRunVasInt != -1){
                        if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                            runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                            startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                        startRunVasInt = -1;
                        runLengthVasInt = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == CAUD_SART ) {
                    if (startRunCaudSart == -1){
                        startRunCaudSart = rescaledData->GetPixel(index);
                        runLengthCaudSart++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunCaudSart){
                            runLengthCaudSart++;
                        }
                        else {
                            if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                                runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                            if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                                startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                            startRunCaudSart = -1;
                            runLengthCaudSart = 0;
                        }
                    }
                }
                else {
                    if (startRunCaudSart != -1){
                        if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                            runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                        if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                            startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                        startRunCaudSart = -1;
                        runLengthCaudSart = 0;
                    }
                }
            }
        }
    }
    //z direction
    for(index[1] = 0; index[1] < size[1]; index[1]++) {
        for(index[0] = 0; index[0] < size[0]; index[0]++) {
            for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
                if ( erodemask->GetPixel(index) == REC_FEM ) {
                    if (startRunRecFem == -1){
                        startRunRecFem = rescaledData->GetPixel(index);
                        runLengthRecFem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunRecFem){
                            runLengthRecFem++;
                        }
                        else {
                            if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                                runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                            if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                                startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                            startRunRecFem = -1;
                            runLengthRecFem = 0;
                        }
                    }
                }
                else {
                    if (startRunRecFem != -1){
                        if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                            runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                        if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                            startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                        startRunRecFem = -1;
                        runLengthRecFem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == SEMIT ) {
                    if (startRunSemit == -1){
                        startRunSemit = rescaledData->GetPixel(index);
                        runLengthSemit++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunSemit){
                            runLengthSemit++;
                        }
                        else {
                            if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                                runLengthSemit = RUN_LENGTH_LEVEL - 1;
                            if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                                startRunSemit = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                            startRunSemit = -1;
                            runLengthSemit = 0;
                        }
                    }
                }
                else {
                    if (startRunSemit != -1){
                        if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                            runLengthSemit = RUN_LENGTH_LEVEL - 1;
                        if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                            startRunSemit = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                        startRunSemit = -1;
                        runLengthSemit = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == CRAN_SART ) {
                    if (startRunCranSart == -1){
                        startRunCranSart = rescaledData->GetPixel(index);
                        runLengthCranSart++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunCranSart){
                            runLengthCranSart++;
                        }
                        else {
                            if(runLengthCranSart > RUN_LENGTH_LEVEL - 1)
                                runLengthCranSart = RUN_LENGTH_LEVEL - 1;
                            if(startRunCranSart > RUN_INTENSITY_LEVEL - 1)
                                startRunCranSart = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixCranSart[runLengthCranSart][startRunCranSart]++;
                            startRunCranSart = -1;
                            runLengthCranSart = 0;
                        }
                    }
                }
                else {
                    if (startRunCaudSart != -1){
                        if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                            runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                        if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                            startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                        startRunCaudSart = -1;
                        runLengthCaudSart = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == BIC_FEM ) {
                    if (startRunBicFem == -1){
                        startRunBicFem = rescaledData->GetPixel(index);
                        runLengthBicFem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunBicFem){
                            runLengthBicFem++;
                        }
                        else {
                            if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                                runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                            if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                                startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;
                            startRunBicFem = -1;
                            runLengthBicFem = 0;
                        }
                    }
                }
                else {
                    if (startRunBicFem != -1){
                        if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                            runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                        if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                            startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;
                        startRunBicFem = -1;
                        runLengthBicFem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == GRACILIS ) {
                    if (startRunGracilis == -1){
                        startRunGracilis = rescaledData->GetPixel(index);
                        runLengthGracilis++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunGracilis){
                            runLengthGracilis++;
                        }
                        else {
                            if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                                runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                            if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                                startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;
                            startRunGracilis = -1;
                            runLengthGracilis = 0;
                        }
                    }
                }
                else {
                    if (startRunGracilis != -1){
                        if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                            runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                        if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                            startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;
                        startRunGracilis = -1;
                        runLengthGracilis = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == ADDUCTOR ) {
                    if (startRunAdductor == -1){
                        startRunAdductor = rescaledData->GetPixel(index);
                        runLengthAdductor++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunAdductor){
                            runLengthAdductor++;
                        }
                        else {
                            if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                                runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                            if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                                startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                            startRunAdductor = -1;
                            runLengthAdductor = 0;
                        }
                    }
                }
                else {
                    if (startRunAdductor != -1){
                        if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                            runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                        if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                            startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                        startRunAdductor = -1;
                        runLengthAdductor = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == SEMI_MEM ) {
                    if (startRunSemiMem == -1){
                        startRunSemiMem = rescaledData->GetPixel(index);
                        runLengthSemiMem++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunSemiMem){
                            runLengthSemiMem++;
                        }
                        else {
                            if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                                runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                            if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                                startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                            startRunSemiMem = -1;
                            runLengthSemiMem = 0;
                        }
                    }
                }
                else {
                    if (startRunSemiMem != -1){
                        if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                            runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                        if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                            startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                        startRunSemiMem = -1;
                        runLengthSemiMem = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_LAT ) {
                    if (startRunVasLat == -1){
                        startRunVasLat = rescaledData->GetPixel(index);
                        runLengthVasLat++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasLat){
                            runLengthVasLat++;
                        }
                        else {
                            if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                                runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                                startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;
                            startRunVasLat = -1;
                            runLengthVasLat = 0;
                        }
                    }
                }
                else {
                    if (startRunVasLat != -1){
                        if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                            runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                            startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;
                        startRunVasLat = -1;
                        runLengthVasLat = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_MED ) {
                    if (startRunVasMed == -1){
                        startRunVasMed = rescaledData->GetPixel(index);
                        runLengthVasMed++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasMed){
                            runLengthVasMed++;
                        }
                        else {
                            if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                                runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                                startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                            startRunVasMed = -1;
                            runLengthVasMed = 0;
                        }
                    }
                }
                else {
                    if (startRunVasMed != -1){
                        if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                            runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                            startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                        startRunVasMed = -1;
                        runLengthVasMed = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == VAS_INT ) {
                    if (startRunVasInt == -1){
                        startRunVasInt = rescaledData->GetPixel(index);
                        runLengthVasInt++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunVasInt){
                            runLengthVasInt++;
                        }
                        else {
                            if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                                runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                            if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                                startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                            startRunVasInt = -1;
                            runLengthVasInt = 0;
                        }
                    }
                }
                else {
                    if (startRunVasInt != -1){
                        if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                            runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                        if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                            startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                        runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                        startRunVasInt = -1;
                        runLengthVasInt = 0;
                    }
                }
                if ( erodemask->GetPixel(index) == CAUD_SART ) {
                    if (startRunCaudSart == -1){
                        startRunCaudSart = rescaledData->GetPixel(index);
                        runLengthCaudSart++;
                    }
                    else {
                        if(rescaledData->GetPixel(index) == startRunCaudSart){
                            runLengthCaudSart++;
                        }
                        else {
                            if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                                runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                            if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                                startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                            runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                            startRunCaudSart = -1;
                            runLengthCaudSart = 0;
                        }
                    }
                }
                else {
                    if (startRunCaudSart != -1){
                        if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                            runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                        if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                            startRunCaudSart = RUN_INTENSITY_LEVEL -1;
                        runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                        startRunCaudSart = -1;
                        runLengthCaudSart = 0;
                    }
                }
            }
        }
    }
    // calculate run length matrix features
    std::ofstream efile( featurefilename.c_str() , std::ios::app );
    float nr = 0, pr = 0, pg = 0;
    // RefFem
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixRecFem[j][i];
             nr += runlengthmatrixRecFem[j][i];
 //            std::cout <<  runlengthmatrixRecFem[j][i] << "  " ;
        }
 //       std::cout << "\n";
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
//    getchar();
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixRecFem[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npRecFem != 0)
        RP = nr / npRecFem;
    else
        RP = 0;
    efile << "RecFem  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // Semit
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixSemit[j][i];
             nr += runlengthmatrixSemit[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixSemit[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npSemit != 0)
        RP = nr / npSemit;
    else
        RP = 0;
    efile << "Semit  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // CranSart
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixCranSart[j][i];
             nr += runlengthmatrixCranSart[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixCranSart[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npCranSart != 0)
        RP = nr / npCranSart;
    else
        RP = 0;
    efile << "CranSart  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // BicFem
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixBicFem[j][i];
             nr += runlengthmatrixBicFem[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixBicFem[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npBicFem != 0)
        RP = nr / npBicFem;
    else
        RP = 0;
    efile << "BicFem  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // Gracilis
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixGracilis[j][i];
             nr += runlengthmatrixGracilis[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixGracilis[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npGracilis != 0)
        RP = nr / npGracilis;
    else
        RP = 0;
    efile << "Gracilis  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // Adductor
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixAdductor[j][i];
             nr += runlengthmatrixAdductor[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixAdductor[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npAdductor != 0)
        RP = nr / npAdductor;
    else
        RP = 0;
    efile << "Adductor  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // SemiMem
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixSemiMem[j][i];
             nr += runlengthmatrixSemiMem[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixSemiMem[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npSemiMem != 0)
        RP = nr / npSemiMem;
    else
        RP = 0;
    efile << "SemiMem  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // VasLat
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixVasLat[j][i];
             nr += runlengthmatrixVasLat[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixVasLat[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npVasLat != 0)
        RP = nr / npVasLat;
    else
        RP = 0;
    efile << "VasLat  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // VasMed
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixVasMed[j][i];
             nr += runlengthmatrixVasMed[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixVasMed[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npVasMed != 0)
        RP = nr / npVasMed;
    else
        RP = 0;
    efile << "VasMed  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // VasInt
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixVasInt[j][i];
             nr += runlengthmatrixVasInt[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixVasInt[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npVasInt != 0)
        RP = nr / npVasInt;
    else
        RP = 0;
    efile << "VasInt  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // CaudSart
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixCaudSart[j][i];
             nr += runlengthmatrixCaudSart[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixCaudSart[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npCaudSart != 0)
        RP = nr / npCaudSart;
    else
        RP = 0;
    efile << "CaudSart  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    efile.close();
}
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::runlengthFeat3D( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType voxelSize)
{
    //muscle based run length features
    int runlengthmatrixRecFem[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixSemit[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixCranSart[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixBicFem[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixGracilis[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixAdductor[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixSemiMem[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixVasLat[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixVasMed[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixVasInt[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL], runlengthmatrixCaudSart[RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL];
    int tempIntensity = 0, tempLength = 0, intensityMax = 0, intensityMin = 0, intensityMaxCut = 0, intensityMinCut = 0, histogram[HISTOGRAM_LIMIT] = {0};
    float intensityRescaleFactor = 0;
    DMDData ddata;
    DMDData::OrientedImageType::Pointer rescaledData = DMDData::OrientedImageType::New(), rescaledData_REC_FEM = DMDData::OrientedImageType::New(), rescaledData_SEMIT = DMDData::OrientedImageType::New(), rescaledData_CRAN_SART = DMDData::OrientedImageType::New(), rescaledData_BIC_FEM = DMDData::OrientedImageType::New(), rescaledData_ADDUCTOR = DMDData::OrientedImageType::New(), rescaledData_GRACILIS = DMDData::OrientedImageType::New(), rescaledData_SEMI_MEM = DMDData::OrientedImageType::New(), rescaledData_VAS_LAT = DMDData::OrientedImageType::New(), rescaledData_VAS_MED = DMDData::OrientedImageType::New(), rescaledData_VAS_INT = DMDData::OrientedImageType::New(), rescaledData_CAUD_SART = DMDData::OrientedImageType::New(), tempMaskX = DMDData::OrientedImageType::New(), tempMaskY = DMDData::OrientedImageType::New(), tempMaskZ = DMDData::OrientedImageType::New();
    float SRE = 0, LRE = 0, GLN = 0, RLN = 0, RP = 0;  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float npRecFem = 0, npSemit = 0, npCranSart = 0, npBicFem = 0, npGracilis = 0, npAdductor = 0, npSemiMem = 0, npVasLat = 0, npVasMed = 0, npVasInt = 0, npCaudSart = 0;
    double imageVolume = 0;

    for (int j = 0; j < RUN_LENGTH_LEVEL; j++) {
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
            runlengthmatrixRecFem[j][i] = 0;
            runlengthmatrixSemit[j][i] = 0;
            runlengthmatrixCranSart[j][i] = 0;
            runlengthmatrixBicFem[j][i] = 0;
            runlengthmatrixGracilis[j][i] = 0;
            runlengthmatrixAdductor[j][i] = 0;
            runlengthmatrixSemiMem[j][i] = 0;
            runlengthmatrixVasLat[j][i] = 0;
            runlengthmatrixVasMed[j][i] = 0;
            runlengthmatrixVasInt[j][i] = 0;
            runlengthmatrixCaudSart[j][i] = 0;
        }
    }

    ddata.imageInitialize ( data, rescaledData );
    ddata.imageInitialize ( data, rescaledData_REC_FEM );
    ddata.imageInitialize ( data, rescaledData_SEMIT );
    ddata.imageInitialize ( data, rescaledData_CRAN_SART );
    ddata.imageInitialize ( data, rescaledData_BIC_FEM );
    ddata.imageInitialize ( data, rescaledData_ADDUCTOR );
    ddata.imageInitialize ( data, rescaledData_GRACILIS );
    ddata.imageInitialize ( data, rescaledData_SEMI_MEM );
    ddata.imageInitialize ( data, rescaledData_VAS_LAT );
    ddata.imageInitialize ( data, rescaledData_VAS_MED );
    ddata.imageInitialize ( data, rescaledData_VAS_INT );
    ddata.imageInitialize ( data, rescaledData_CAUD_SART );

    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;
    IteratorType ReScaledDataIterator( rescaledData, rescaledData->GetRequestedRegion() ) ;

    // identify minum and maxmum intensity of the image
    constDataIterator.GoToBegin();
    intensityMin = intensityMax = constDataIterator.Get();
    for ( constDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator ) {
        int hisTmp = constDataIterator.Get();
        if (hisTmp >= HISTOGRAM_LIMIT)
            hisTmp = HISTOGRAM_LIMIT - 1;
        histogram[hisTmp]++;
        if(hisTmp > 0)
            imageVolume++;
        if (constDataIterator.Get() < intensityMin)
            intensityMin = constDataIterator.Get();
        if (constDataIterator.Get() > intensityMax)
            intensityMax = constDataIterator.Get();
    }
    // identify 5% minum and maxmum intensity of the image
    float cumImageVolume = 0;
    for (int i = 1; i <= intensityMax; i++){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMinCut = i;
            break;
        }
    }
    std::cout << imageVolume << "  " << cumImageVolume << "  ";
    cumImageVolume = 0;
    for (int i = intensityMax; i >= 1; i--){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMaxCut = i;
            break;
        }
    }
    std::cout << cumImageVolume << std::endl;
    // cut the high and low %5 intensity
 //   intensityMinCut = intensityMin + (intensityMax - intensityMin) * INTENSITY_CUT_RATE;
//    intensityMaxCut = intensityMax - (intensityMax - intensityMin) * INTENSITY_CUT_RATE;
    //rescale images
    intensityRescaleFactor = (float)(intensityMaxCut - intensityMinCut) / RUN_INTENSITY_LEVEL;
    std::cout << intensityMin << "  " << intensityMax << "  " << intensityMinCut << "  " << intensityMaxCut << "  " << intensityRescaleFactor << std::endl;
    for ( constDataIterator.GoToBegin(), ReScaledDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator, ++ReScaledDataIterator ) {
        tempIntensity = constDataIterator.Get();
        if (tempIntensity < intensityMinCut)
            tempIntensity = intensityMinCut;
        if (tempIntensity > intensityMaxCut)
            tempIntensity = intensityMaxCut;
       // std::cout << tempIntensity << "  " << tempIntensity/ intensityRescaleFactor << "  " << round(tempIntensity/ intensityRescaleFactor) << std::endl;
        //getchar();
        tempIntensity = round(tempIntensity / intensityRescaleFactor);
        if(tempIntensity < 0)
            tempIntensity = 0;
        if(tempIntensity >= RUN_INTENSITY_LEVEL)
            tempIntensity = RUN_INTENSITY_LEVEL - 1;
        ReScaledDataIterator.Set(tempIntensity);
    }
    std::cout << "rescale finished! " << std::endl;
    //establish run length matrix
    DMDData::OrientedImageType::SizeType   size = data->GetLargestPossibleRegion().GetSize();
    DMDData::OrientedImageType::IndexType  index;
    int upperBound_REC_FEM = 0, lowerBound_REC_FEM = HISTOGRAM_LIMIT, upperBound_SEMIT = 0, lowerBound_SEMIT = HISTOGRAM_LIMIT, upperBound_CRAN_SART = 0, lowerBound_CRAN_SART = HISTOGRAM_LIMIT, upperBound_BIC_FEM = 0, lowerBound_BIC_FEM = HISTOGRAM_LIMIT, upperBound_GRACILIS = 0, lowerBound_GRACILIS = HISTOGRAM_LIMIT, upperBound_ADDUCTOR = 0, lowerBound_ADDUCTOR = HISTOGRAM_LIMIT, upperBound_SEMI_MEM = 0, lowerBound_SEMI_MEM = HISTOGRAM_LIMIT, upperBound_VAS_LAT = 0, lowerBound_VAS_LAT = HISTOGRAM_LIMIT, upperBound_VAS_MED = 0, lowerBound_VAS_MED = HISTOGRAM_LIMIT, upperBound_VAS_INT = 0, lowerBound_VAS_INT = HISTOGRAM_LIMIT, upperBound_CAUD_SART = 0, lowerBound_CAUD_SART = HISTOGRAM_LIMIT;
    int startRunRecFem = -1, runLengthRecFem = 0, startRunSemit = -1, runLengthSemit = 0, startRunCranSart = -1, runLengthCranSart = 0, startRunBicFem = -1, runLengthBicFem = 0, startRunGracilis = -1, runLengthGracilis = 0, startRunAdductor = -1, runLengthAdductor = 0, startRunSemiMem = -1, runLengthSemiMem = 0, startRunVasLat = -1, runLengthVasLat = 0, startRunVasMed = -1, runLengthVasMed = 0, startRunVasInt = -1, runLengthVasInt = 0, startRunCaudSart = -1, runLengthCaudSart = 0;
    // establish rescaled image for each muscle
    for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
        for(index[1] = 0; index[1] < size[1]; index[1]++) {
            for(index[0] = 0; index[0] < size[0]; index[0]++) {
                // REC_FEM
                if ( erodemask->GetPixel(index) == REC_FEM ) {
                    rescaledData_REC_FEM->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_REC_FEM) {
                        lowerBound_REC_FEM = rescaledData->GetPixel(index);
     //                   std::cout << "lower bound:  " << lowerBound_REC_FEM << std::endl;
                    }
                    if(rescaledData->GetPixel(index) > upperBound_REC_FEM){
                        upperBound_REC_FEM = rescaledData->GetPixel(index);
       //                 std::cout << "upper bound:  " << upperBound_REC_FEM << std::endl;
                    }
                    npRecFem++;
                }
                else
                    rescaledData_REC_FEM->SetPixel(index, BACKGROUND);
                // SEMIT
                if ( erodemask->GetPixel(index) == SEMIT ) {
                    rescaledData_SEMIT->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_SEMIT)
                        lowerBound_SEMIT = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_SEMIT)
                        upperBound_SEMIT = rescaledData->GetPixel(index);
                    npSemit++;
                }
                else
                    rescaledData_SEMIT->SetPixel(index, BACKGROUND);
                // CRAN_SART
                if ( erodemask->GetPixel(index) == CRAN_SART ) {
                    rescaledData_CRAN_SART->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_CRAN_SART)
                        lowerBound_CRAN_SART = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_CRAN_SART)
                        upperBound_CRAN_SART = rescaledData->GetPixel(index);
                    npCranSart++;
                }
                else
                    rescaledData_CRAN_SART->SetPixel(index, BACKGROUND);
                // BIC_FEM
                if ( erodemask->GetPixel(index) == BIC_FEM ) {
                    rescaledData_BIC_FEM->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_BIC_FEM)
                        lowerBound_BIC_FEM = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_BIC_FEM)
                        upperBound_BIC_FEM = rescaledData->GetPixel(index);
                    npBicFem++;
                }
                else
                    rescaledData_BIC_FEM->SetPixel(index, BACKGROUND);
                // GRACILIS
                if ( erodemask->GetPixel(index) == GRACILIS ) {
                    rescaledData_GRACILIS->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_GRACILIS)
                        lowerBound_GRACILIS = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_GRACILIS)
                        upperBound_GRACILIS = rescaledData->GetPixel(index);
                    npGracilis++;
                }
                else
                    rescaledData_GRACILIS->SetPixel(index, BACKGROUND);
                // ADDUCTOR
                if ( erodemask->GetPixel(index) == ADDUCTOR ) {
                    rescaledData_ADDUCTOR->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_ADDUCTOR)
                        lowerBound_ADDUCTOR = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_ADDUCTOR)
                        upperBound_ADDUCTOR = rescaledData->GetPixel(index);
                    npAdductor++;
                }
                else
                    rescaledData_ADDUCTOR->SetPixel(index, BACKGROUND);
                // SEMI_MEM
                if ( erodemask->GetPixel(index) == SEMI_MEM ) {
                    rescaledData_SEMI_MEM->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_SEMI_MEM)
                        lowerBound_SEMI_MEM = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_SEMI_MEM)
                        upperBound_SEMI_MEM = rescaledData->GetPixel(index);
                    npSemiMem++;
                }
                else
                    rescaledData_SEMI_MEM->SetPixel(index, BACKGROUND);
                // VAS_LAT
                if ( erodemask->GetPixel(index) == VAS_LAT ) {
                    rescaledData_VAS_LAT->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_VAS_LAT)
                        lowerBound_VAS_LAT = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_VAS_LAT)
                        upperBound_VAS_LAT = rescaledData->GetPixel(index);
                    npVasLat++;
                }
                else
                    rescaledData_VAS_LAT->SetPixel(index, BACKGROUND);
                // VAS_MED
                if ( erodemask->GetPixel(index) == VAS_MED ) {
                    rescaledData_VAS_MED->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_VAS_MED)
                        lowerBound_VAS_MED = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_VAS_MED)
                        upperBound_VAS_MED = rescaledData->GetPixel(index);
                    npVasMed++;
                }
                else
                    rescaledData_VAS_MED->SetPixel(index, BACKGROUND);
                // VAS_INT
                if ( erodemask->GetPixel(index) == VAS_INT ) {
                    rescaledData_VAS_INT->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_VAS_INT)
                        lowerBound_VAS_INT = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_VAS_INT)
                        upperBound_VAS_INT = rescaledData->GetPixel(index);
                    npVasInt++;
                }
                else
                    rescaledData_VAS_INT->SetPixel(index, BACKGROUND);
                // CAUD_SART
                if ( erodemask->GetPixel(index) == CAUD_SART ) {
                    rescaledData_CAUD_SART->SetPixel(index, rescaledData->GetPixel(index));
                    if(rescaledData->GetPixel(index) < lowerBound_CAUD_SART)
                        lowerBound_CAUD_SART = rescaledData->GetPixel(index);
                    if(rescaledData->GetPixel(index) > upperBound_CAUD_SART)
                        upperBound_CAUD_SART = rescaledData->GetPixel(index);
                    npCaudSart++;
                }
                else
                    rescaledData_CAUD_SART->SetPixel(index, BACKGROUND);
            }
        }
    }
    std::cout << "upper bound: " << upperBound_REC_FEM << "  lower bound: " << lowerBound_REC_FEM << std::endl;
    // establish run length matrix
    for (int i = lowerBound_REC_FEM; i <= upperBound_REC_FEM; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // REC_FEM
                    if ( erodemask->GetPixel(index) == REC_FEM && rescaledData->GetPixel(index) == i) {
                        rescaledData_REC_FEM->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_REC_FEM->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_REC_FEM);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_REC_FEM->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixRecFem[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixRecFem[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Rec Fem finished" << std::endl;
/*    getchar();
    for (int j = 0; j < RUN_LENGTH_LEVEL; j++)  {
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++)  {
            std::cout << runlengthmatrixRecFem[j][i] << " ";
        }
        std::cout << "\n" ;
    }
*/
    for (int i = lowerBound_SEMIT; i <= upperBound_SEMIT; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // SEMIT
                    if ( erodemask->GetPixel(index) == SEMIT && rescaledData->GetPixel(index) == i) {
                        rescaledData_SEMIT->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_SEMIT->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_SEMIT);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_SEMIT->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixSemit[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixSemit[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Semit finished" << std::endl;

    for (int i = lowerBound_BIC_FEM; i <= upperBound_BIC_FEM; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // BIC_FEM
                    if ( erodemask->GetPixel(index) == BIC_FEM && rescaledData->GetPixel(index) == i) {
                        rescaledData_BIC_FEM->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_BIC_FEM->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_BIC_FEM);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_BIC_FEM->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixBicFem[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixBicFem[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Bic Fem finished" << std::endl;

    for (int i = lowerBound_CRAN_SART; i <= upperBound_CRAN_SART; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // CRAN_SART
                    if ( erodemask->GetPixel(index) == CRAN_SART && rescaledData->GetPixel(index) == i) {
                        rescaledData_CRAN_SART->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_CRAN_SART->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_CRAN_SART);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_CRAN_SART->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixCranSart[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixCranSart[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Cran Sart finished" << std::endl;

    for (int i = lowerBound_GRACILIS; i <= upperBound_GRACILIS; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // GRACILIS
                    if ( erodemask->GetPixel(index) == GRACILIS && rescaledData->GetPixel(index) == i) {
                        rescaledData_GRACILIS->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_GRACILIS->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_GRACILIS);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_GRACILIS->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixGracilis[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixGracilis[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Gra finished" << std::endl;

    for (int i = lowerBound_ADDUCTOR; i <= upperBound_ADDUCTOR; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // ADDUCTOR
                    if ( erodemask->GetPixel(index) == ADDUCTOR && rescaledData->GetPixel(index) == i) {
                        rescaledData_ADDUCTOR->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_ADDUCTOR->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_ADDUCTOR);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_ADDUCTOR->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixAdductor[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixAdductor[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Add finished" << std::endl;
/*
    for (int i = lowerBound_SEMI_MEM; i <= upperBound_SEMI_MEM; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // SEMI_MEM
                    if ( erodemask->GetPixel(index) == SEMI_MEM && rescaledData->GetPixel(index) == i) {
                        rescaledData_SEMI_MEM->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_SEMI_MEM->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_SEMI_MEM);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_SEMI_MEM->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixSemiMem[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixSemiMem[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Semi Mem finished" << std::endl;

    for (int i = lowerBound_VAS_LAT; i <= upperBound_VAS_LAT; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // VAS_LAT
                    if ( erodemask->GetPixel(index) == VAS_LAT && rescaledData->GetPixel(index) == i) {
                        rescaledData_VAS_LAT->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_VAS_LAT->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_VAS_LAT);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_VAS_LAT->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixVasLat[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixVasLat[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Vas Lat finished" << std::endl;

    for (int i = lowerBound_VAS_MED; i <= upperBound_VAS_MED; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // VAS_MED
                    if ( erodemask->GetPixel(index) == VAS_MED && rescaledData->GetPixel(index) == i) {
                        rescaledData_VAS_MED->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_VAS_MED->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_VAS_MED);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_VAS_MED->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixVasMed[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixVasMed[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Vas Med finished" << std::endl;

    for (int i = lowerBound_VAS_INT; i <= upperBound_VAS_INT; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // VAS_INT
                    if ( erodemask->GetPixel(index) == VAS_INT && rescaledData->GetPixel(index) == i) {
                        rescaledData_VAS_INT->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_VAS_INT->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_VAS_INT);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_VAS_INT->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixVasInt[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixVasInt[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Vas Int finished" << std::endl;

    for (int i = lowerBound_CAUD_SART; i <= upperBound_CAUD_SART; i++) {
        int objCount = 0;
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    // CAUD_SART
                    if ( erodemask->GetPixel(index) == CAUD_SART && rescaledData->GetPixel(index) == i) {
                        rescaledData_CAUD_SART->SetPixel(index, rescaledData->GetPixel(index));
                    }
                    else
                        rescaledData_CAUD_SART->SetPixel(index, BACKGROUND);
                }
            }
        }
        objCount = connectedComponentLabeling (rescaledData_CAUD_SART);
        int objVolume[SIZE_LABELING_BIN] = {0};
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    int k = rescaledData_CAUD_SART->GetPixel(index);
                    if(k){
                        objVolume[k]++;
                    }
                }
            }
        }
        for (int j = 0; j < SIZE_LABELING_BIN; j++) {
            objVolume[j] = round(objVolume[j] * voxelSize);
            if (objVolume[j]) {
                if(objVolume[j] >= RUN_LENGTH_LEVEL){
                    runlengthmatrixCaudSart[RUN_LENGTH_LEVEL - 1][i]++;
                }
                else{
                    if(objVolume[j]){
                        runlengthmatrixCaudSart[objVolume[j]][i]++;
                    }
                }
            }
        }
    }
    std::cout << "Caud Sart finished" << std::endl;
*/
    /*
    for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
        for(index[1] = 0; index[1] < size[1]; index[1]++) {
            for(index[0] = 0; index[0] < size[0]; index[0]++) {
                if ( erodemask->GetPixel(index) == REC_FEM && rescaledData_REC_FEM->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    //runLengthRecFem = ddata.connectedComponentRegionGrowing (rescaledData_REC_FEM->GetPixel(index),rescaledData_REC_FEM->GetPixel(index), rescaledData_REC_FEM, index);
                    connectedComponentLabeling (rescaledData_REC_FEM);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthRecFem  << std::endl;
                    std::cout << "labelling finished" << std::endl;
                    getchar();
                    startRunRecFem = rescaledData_REC_FEM->GetPixel(index);
                    npRecFem += runLengthRecFem;
                    if(runLengthRecFem > RUN_LENGTH_LEVEL - 1)
                        runLengthRecFem = RUN_LENGTH_LEVEL - 1;
                    if(startRunRecFem > RUN_INTENSITY_LEVEL - 1)
                        startRunRecFem = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixRecFem[runLengthRecFem][startRunRecFem]++;
                }
                if ( erodemask->GetPixel(index) == SEMIT && rescaledData_SEMIT->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthSemit = ddata.connectedComponentRegionGrowing (rescaledData_SEMIT->GetPixel(index),rescaledData_SEMIT->GetPixel(index), rescaledData_SEMIT, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthSemit << std::endl;
                    startRunSemit = rescaledData_SEMIT->GetPixel(index);
                    npSemit += runLengthSemit;
                    if(runLengthSemit > RUN_LENGTH_LEVEL - 1)
                        runLengthSemit = RUN_LENGTH_LEVEL - 1;
                    if(startRunSemit > RUN_INTENSITY_LEVEL - 1)
                        startRunSemit = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixSemit[runLengthSemit][startRunSemit]++;
                }
                if ( erodemask->GetPixel(index) == CRAN_SART && rescaledData_CRAN_SART->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthCranSart = ddata.connectedComponentRegionGrowing (rescaledData_CRAN_SART->GetPixel(index),rescaledData_CRAN_SART->GetPixel(index), rescaledData_CRAN_SART, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthCranSart << std::endl;
                    startRunCranSart = rescaledData_CRAN_SART->GetPixel(index);
                    npCranSart += runLengthCranSart;
                    if(runLengthCranSart > RUN_LENGTH_LEVEL - 1)
                        runLengthCranSart = RUN_LENGTH_LEVEL - 1;
                    if(startRunCranSart > RUN_INTENSITY_LEVEL - 1)
                        startRunCranSart = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixCranSart[runLengthCranSart][startRunCranSart]++;

                }
                if ( erodemask->GetPixel(index) == BIC_FEM && rescaledData_BIC_FEM->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    //runLengthBicFem = ddata.connectedComponentRegionGrowing (rescaledData_BIC_FEM->GetPixel(index),rescaledData_BIC_FEM->GetPixel(index), rescaledData_BIC_FEM, index);
                    connectedComponentLabeling (rescaledData_BIC_FEM);
                    std::cout << "labelling finished" << std::endl;
                    getchar();
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthBicFem << std::endl;
                    startRunBicFem = rescaledData_BIC_FEM->GetPixel(index);
                    npBicFem += runLengthBicFem;
                    if(runLengthBicFem > RUN_LENGTH_LEVEL - 1)
                        runLengthBicFem = RUN_LENGTH_LEVEL - 1;
                    if(startRunBicFem > RUN_INTENSITY_LEVEL - 1)
                        startRunBicFem = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixBicFem[runLengthBicFem][startRunBicFem]++;

                }
                if ( erodemask->GetPixel(index) == GRACILIS && rescaledData_GRACILIS->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthGracilis = ddata.connectedComponentRegionGrowing (rescaledData_GRACILIS->GetPixel(index),rescaledData_GRACILIS->GetPixel(index), rescaledData_GRACILIS, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthGracilis << std::endl;
                    startRunGracilis = rescaledData_GRACILIS->GetPixel(index);
                    npGracilis += runLengthGracilis;
                    if(runLengthGracilis > RUN_LENGTH_LEVEL - 1)
                        runLengthGracilis = RUN_LENGTH_LEVEL - 1;
                    if(startRunGracilis > RUN_INTENSITY_LEVEL - 1)
                        startRunGracilis = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixGracilis[runLengthGracilis][startRunGracilis]++;

                }
                if ( erodemask->GetPixel(index) == ADDUCTOR && rescaledData_ADDUCTOR->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthAdductor = ddata.connectedComponentRegionGrowing (rescaledData_ADDUCTOR->GetPixel(index),rescaledData_ADDUCTOR->GetPixel(index), rescaledData_ADDUCTOR, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthAdductor << std::endl;
                    startRunAdductor = rescaledData_ADDUCTOR->GetPixel(index);
                    npAdductor += runLengthAdductor;
                    if(runLengthAdductor > RUN_LENGTH_LEVEL - 1)
                        runLengthAdductor = RUN_LENGTH_LEVEL - 1;
                    if(startRunAdductor > RUN_INTENSITY_LEVEL - 1)
                        startRunAdductor = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixAdductor[runLengthAdductor][startRunAdductor]++;
                }
                if ( erodemask->GetPixel(index) == SEMI_MEM && rescaledData_SEMI_MEM->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthSemiMem = ddata.connectedComponentRegionGrowing (rescaledData_SEMI_MEM->GetPixel(index),rescaledData_SEMI_MEM->GetPixel(index), rescaledData_SEMI_MEM, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthSemiMem << std::endl;
                    startRunSemiMem = rescaledData_SEMI_MEM->GetPixel(index);
                    npSemiMem += runLengthSemiMem;
                    if(runLengthSemiMem > RUN_LENGTH_LEVEL - 1)
                        runLengthSemiMem = RUN_LENGTH_LEVEL - 1;
                    if(startRunSemiMem > RUN_INTENSITY_LEVEL - 1)
                        startRunSemiMem = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixSemiMem[runLengthSemiMem][startRunSemiMem]++;
                }
                if ( erodemask->GetPixel(index) == VAS_LAT && rescaledData_VAS_LAT->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthVasLat = ddata.connectedComponentRegionGrowing (rescaledData_VAS_LAT->GetPixel(index),rescaledData_VAS_LAT->GetPixel(index), rescaledData_VAS_LAT, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthVasLat << std::endl;
                    startRunVasLat = rescaledData_VAS_LAT->GetPixel(index);
                    npVasLat += runLengthVasLat;
                    if(runLengthVasLat > RUN_LENGTH_LEVEL - 1)
                        runLengthVasLat = RUN_LENGTH_LEVEL - 1;
                    if(startRunVasLat > RUN_INTENSITY_LEVEL - 1)
                        startRunVasLat = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixVasLat[runLengthVasLat][startRunVasLat]++;

                }
                if ( erodemask->GetPixel(index) == VAS_MED && rescaledData_VAS_MED->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthVasMed = ddata.connectedComponentRegionGrowing (rescaledData_VAS_MED->GetPixel(index),rescaledData_VAS_MED->GetPixel(index), rescaledData_VAS_MED, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthVasMed << std::endl;
                    startRunVasMed = rescaledData_VAS_MED->GetPixel(index);
                    npVasMed += runLengthVasMed;
                    if(runLengthVasMed > RUN_LENGTH_LEVEL - 1)
                        runLengthVasMed = RUN_LENGTH_LEVEL - 1;
                    if(startRunVasMed > RUN_INTENSITY_LEVEL - 1)
                        startRunVasMed = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixVasMed[runLengthVasMed][startRunVasMed]++;
                }
                if ( erodemask->GetPixel(index) == VAS_INT && rescaledData_VAS_INT->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthVasInt = ddata.connectedComponentRegionGrowing (rescaledData_VAS_INT->GetPixel(index),rescaledData_VAS_INT->GetPixel(index), rescaledData_VAS_INT, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthVasInt << std::endl;
                    startRunVasInt = rescaledData_VAS_INT->GetPixel(index);
                    npVasInt += runLengthVasInt;
                    if(runLengthVasInt > RUN_LENGTH_LEVEL - 1)
                        runLengthVasInt = RUN_LENGTH_LEVEL - 1;
                    if(startRunVasInt > RUN_INTENSITY_LEVEL - 1)
                        startRunVasInt = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixVasInt[runLengthVasInt][startRunVasInt]++;
                }
                if ( erodemask->GetPixel(index) == CAUD_SART && rescaledData_CAUD_SART->GetPixel(index) != SEARCHED_PIXEL_VALUE ) {
                    runLengthCaudSart = ddata.connectedComponentRegionGrowing (rescaledData_CAUD_SART->GetPixel(index),rescaledData_CAUD_SART->GetPixel(index), rescaledData_CAUD_SART, index);
                    std::cout << index[0] << "  " << index[1] << "  " << index[2] << "  " << runLengthCaudSart << std::endl;
                    startRunCaudSart = rescaledData_CAUD_SART->GetPixel(index);
                    npCaudSart += runLengthCaudSart;
                    if(runLengthCaudSart > RUN_LENGTH_LEVEL - 1)
                        runLengthCaudSart = RUN_LENGTH_LEVEL - 1;
                    if(startRunCaudSart > RUN_INTENSITY_LEVEL - 1)
                        startRunCaudSart = RUN_INTENSITY_LEVEL - 1;
                    runlengthmatrixCaudSart[runLengthCaudSart][startRunCaudSart]++;
                }
            }
        }
    }
*/
    // calculate run length matrix features
    std::ofstream efile( featurefilename.c_str() , std::ios::app );
    float nr = 0, pr = 0, pg = 0;
    // RefFem
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixRecFem[j][i];
             nr += runlengthmatrixRecFem[j][i];
 //            std::cout <<  runlengthmatrixRecFem[j][i] << "  " ;
        }
 //       std::cout << "\n";
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
//    getchar();
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixRecFem[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npRecFem != 0)
        RP = nr / npRecFem;
    else
        RP = 0;
    efile << "RecFem  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // Semit
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixSemit[j][i];
             nr += runlengthmatrixSemit[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixSemit[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npSemit != 0)
        RP = nr / npSemit;
    else
        RP = 0;
    efile << "Semit  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // CranSart
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixCranSart[j][i];
             nr += runlengthmatrixCranSart[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixCranSart[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npCranSart != 0)
        RP = nr / npCranSart;
    else
        RP = 0;
    efile << "CranSart  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // BicFem
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixBicFem[j][i];
             nr += runlengthmatrixBicFem[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixBicFem[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npBicFem != 0)
        RP = nr / npBicFem;
    else
        RP = 0;
    efile << "BicFem  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // Gracilis
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixGracilis[j][i];
             nr += runlengthmatrixGracilis[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixGracilis[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npGracilis != 0)
        RP = nr / npGracilis;
    else
        RP = 0;
    efile << "Gracilis  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // Adductor
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixAdductor[j][i];
             nr += runlengthmatrixAdductor[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixAdductor[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npAdductor != 0)
        RP = nr / npAdductor;
    else
        RP = 0;
    efile << "Adductor  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // SemiMem
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixSemiMem[j][i];
             nr += runlengthmatrixSemiMem[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixSemiMem[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npSemiMem != 0)
        RP = nr / npSemiMem;
    else
        RP = 0;
    efile << "SemiMem  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // VasLat
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixVasLat[j][i];
             nr += runlengthmatrixVasLat[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixVasLat[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npVasLat != 0)
        RP = nr / npVasLat;
    else
        RP = 0;
    efile << "VasLat  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // VasMed
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixVasMed[j][i];
             nr += runlengthmatrixVasMed[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixVasMed[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npVasMed != 0)
        RP = nr / npVasMed;
    else
        RP = 0;
    efile << "VasMed  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // VasInt
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixVasInt[j][i];
             nr += runlengthmatrixVasInt[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixVasInt[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npVasInt != 0)
        RP = nr / npVasInt;
    else
        RP = 0;
    efile << "VasInt  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    // CaudSart
    nr = 0;
    for (int j = 1; j < RUN_LENGTH_LEVEL; j++) {
        pr = 0;
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
             pr += runlengthmatrixCaudSart[j][i];
             nr += runlengthmatrixCaudSart[j][i];
        }
        SRE += pr / (j * j);
        LRE += pr * j * j;
        RLN += pr * pr;
    }
    for (int j = 0; j < RUN_INTENSITY_LEVEL; j++) {
        pg = 0;
        for (int i = 1; i < RUN_LENGTH_LEVEL; i++) {
             pg += runlengthmatrixCaudSart[i][j];
        }
        GLN += pg * pg;
    }
    if(nr != 0){
        SRE = SRE / nr;
        LRE = LRE / nr;
        GLN = GLN / nr;
        RLN = RLN / nr;
    }
    else {
        SRE = 0;
        LRE = 0;
        GLN = 0;
        RLN = 0;
    }
    if (npCaudSart != 0)
        RP = nr / npCaudSart;
    else
        RP = 0;
    efile << "CaudSart  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
    efile.close();
}

//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::runlengthFeat2DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType voxelSize)
{
    //ROI based run length features
    int runlengthmatrix[NUMBER_MUSCLE + 1][RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL];
    int tempIntensity = 0, tempLength = 0, intensityMax = 0, intensityMin = 0, intensityMaxCut = 0, intensityMinCut = 0, histogram[HISTOGRAM_LIMIT] = {0};
    float intensityRescaleFactor = 0;
    DMDData ddata;
    DMDData::OrientedImageType::Pointer rescaledData = DMDData::OrientedImageType::New(), rescaledData_REC_FEM = DMDData::OrientedImageType::New(), rescaledData_SEMIT = DMDData::OrientedImageType::New(), rescaledData_CRAN_SART = DMDData::OrientedImageType::New(), rescaledData_BIC_FEM = DMDData::OrientedImageType::New(), rescaledData_ADDUCTOR = DMDData::OrientedImageType::New(), rescaledData_GRACILIS = DMDData::OrientedImageType::New(), rescaledData_SEMI_MEM = DMDData::OrientedImageType::New(), rescaledData_VAS_LAT = DMDData::OrientedImageType::New(), rescaledData_VAS_MED = DMDData::OrientedImageType::New(), rescaledData_VAS_INT = DMDData::OrientedImageType::New(), rescaledData_CAUD_SART = DMDData::OrientedImageType::New(), tempMaskX = DMDData::OrientedImageType::New(), tempMaskY = DMDData::OrientedImageType::New(), tempMaskZ = DMDData::OrientedImageType::New();
    float SRE = 0, LRE = 0, GLN = 0, RLN = 0, RP = 0;  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float npRecFem = 0, npSemit = 0, npCranSart = 0, npBicFem = 0, npGracilis = 0, npAdductor = 0, npSemiMem = 0, npVasLat = 0, npVasMed = 0, npVasInt = 0, npCaudSart = 0, np = 0;
    double imageVolume = 0;
    DMDData::OrientedImageType::SpacingType    spacing = data->GetSpacing() ;

    for (int j = 0; j < RUN_LENGTH_LEVEL; j++) {
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
            runlengthmatrix[REC_FEM][j][i] = 0;
            runlengthmatrix[SEMIT][j][i] = 0;
            runlengthmatrix[CRAN_SART][j][i] = 0;
            runlengthmatrix[BIC_FEM][j][i] = 0;
            runlengthmatrix[GRACILIS][j][i] = 0;
            runlengthmatrix[ADDUCTOR][j][i] = 0;
            runlengthmatrix[SEMI_MEM][j][i] = 0;
            runlengthmatrix[VAS_LAT][j][i] = 0;
            runlengthmatrix[VAS_MED][j][i] = 0;
            runlengthmatrix[VAS_INT][j][i] = 0;
            runlengthmatrix[CAUD_SART][j][i] = 0;
        }
    }
    ddata.imageInitialize ( data, rescaledData );
    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;
    IteratorType ReScaledDataIterator( rescaledData, rescaledData->GetRequestedRegion() ) ;
    // identify minum and maxmum intensity of the image
    constDataIterator.GoToBegin();
    intensityMin = intensityMax = constDataIterator.Get();
    for ( constDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator ) {
        int hisTmp = constDataIterator.Get();
        if (hisTmp >= HISTOGRAM_LIMIT)
            hisTmp = HISTOGRAM_LIMIT - 1;
        histogram[hisTmp]++;
        if(hisTmp > 0)
            imageVolume++;
        if (constDataIterator.Get() < intensityMin)
            intensityMin = constDataIterator.Get();
        if (constDataIterator.Get() > intensityMax)
            intensityMax = constDataIterator.Get();
    }
    // identify 5% minum and maxmum intensity of the image
    float cumImageVolume = 0;
    for (int i = 1; i <= intensityMax; i++){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMinCut = i;
            break;
        }
    }
    std::cout << imageVolume << "  " << cumImageVolume << "  ";
    cumImageVolume = 0;
    for (int i = intensityMax; i >= 1; i--){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMaxCut = i;
            break;
        }
    }
    //rescale images
    intensityRescaleFactor = (float)(intensityMaxCut - intensityMinCut) / RUN_INTENSITY_LEVEL;
    for ( constDataIterator.GoToBegin(), ReScaledDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator, ++ReScaledDataIterator ) {
        tempIntensity = constDataIterator.Get();
        if (tempIntensity < intensityMinCut)
            tempIntensity = intensityMinCut;
        if (tempIntensity > intensityMaxCut)
            tempIntensity = intensityMaxCut;
        tempIntensity = round(tempIntensity / intensityRescaleFactor);
        if(tempIntensity < 0)
            tempIntensity = 0;
        if(tempIntensity >= RUN_INTENSITY_LEVEL)
            tempIntensity = RUN_INTENSITY_LEVEL - 1;
        ReScaledDataIterator.Set(tempIntensity);
    }
    std::cout << "rescale finished! " << std::endl;
    //establish run length matrix
    DMDData::OrientedImageType::SizeType   size = data->GetLargestPossibleRegion().GetSize();
    DMDData::OrientedImageType::IndexType  index;
    int upperBound[NUMBER_MUSCLE + 1] = {0}, lowerBound[NUMBER_MUSCLE + 1] = {0};
    // establish rescaled image for each muscle
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    if ( erodemask->GetPixel(index) == p ) {
                        if(rescaledData->GetPixel(index) < lowerBound[p]) {
                            lowerBound[p] = rescaledData->GetPixel(index);
                        }
                        if(rescaledData->GetPixel(index) > upperBound[p]){
                            upperBound[p] = rescaledData->GetPixel(index);
                        }
                    }
                }
            }
        }
    }
    DMDData::OrientedImageType::Pointer roi = DMDData::OrientedImageType::New();
    int nroi = 0; //the number of rois used for calculating features for each muscle.
    OrientedImageType::RegionType region;
    OrientedImageType::IndexType start;
    OrientedImageType::SizeType roiSize;
    OrientedImageType::IndexType roiIndex;
    roiSize[0] = roiSize[1] = roiSize[2] = RUN_INTENSITY_LEVEL;
    start[0] = start[1] = start[2] = 0;
    region.SetSize(roiSize);
    region.SetIndex(start);
    roi->SetRegions(region);
    roi->Allocate();
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        std::cout << "lower bound: " << lowerBound[p] << "  upper bound: " << upperBound[p] << std::endl;
        if( lowerBound[p] == 0 && upperBound[p] == 0)
            continue;
        nroi = 0;
        SRE = LRE = GLN = RLN = RP = 0;
        std::ofstream efile( featurefilename.c_str() , std::ios::app );
        for (int m = 0; m < size[2]; m += RUN_LENGTH_LEVEL) {
            for (int n = 0; n < size[1]; n += RUN_LENGTH_LEVEL) {
                for (int r = 0; r < size[0]; r += RUN_LENGTH_LEVEL) {
                    float SRETmp = 0, LRETmp = 0, GLNTmp = 0, RLNTmp = 0, RPTmp = 0;
                    np = 0;
                  //  for (int i = lowerBound[p]; i <= upperBound[p]; i++) {
                    int objCount = 0;
                    bool inMuscleFlag = 0;
                    roi->FillBuffer(itk::NumericTraits<DMDData::OrientedImageType::PixelType>::Zero);
                    inMuscleFlag = 0;
                    for(index[2] = m , roiIndex[2] = 0; index[2] < m + RUN_LENGTH_LEVEL; index[2]++, roiIndex[2]++) {
                        for(index[1] = n, roiIndex[1] = 0; index[1] < n + RUN_LENGTH_LEVEL; index[1]++, roiIndex[1]++) {
                            for(index[0] = r, roiIndex[0] = 0; index[0] < r + RUN_LENGTH_LEVEL; index[0]++, roiIndex[0]++) {
                                if(roiIndex[0] < size[0] && roiIndex[1] < size[1] && roiIndex[2] < size[2]) {
                                    if ( erodemask->GetPixel(index) == p) {
                                        roi->SetPixel(roiIndex, rescaledData->GetPixel(index));
                                        inMuscleFlag = 1;
                                        np++;
                                    }
                                    else
                                        roi->SetPixel(roiIndex, BACKGROUND);
                                }
                            }
                        }
                    }
                    if(inMuscleFlag == 1) {
                        int startRun = -1, runLength = 0;
                        // x direction
                        for(index[2] = m , roiIndex[2] = 0; index[2] < m + RUN_LENGTH_LEVEL; index[2]++, roiIndex[2]++) {
                            for(index[1] = n, roiIndex[1] = 0; index[1] < n + RUN_LENGTH_LEVEL; index[1]++, roiIndex[1]++) {
                                for(index[0] = r, roiIndex[0] = 0; index[0] < r + RUN_LENGTH_LEVEL; index[0]++, roiIndex[0]++) {
                                    if ( erodemask->GetPixel(index) == p ) {
                                        if (startRun == -1){
                                            startRun = roi->GetPixel(roiIndex);
                                            runLength++;
                                        }
                                        else {
                                            if(roi->GetPixel(roiIndex) == startRun){
                                                runLength++;
                                            }
                                            else {
                                                if(runLength > RUN_LENGTH_LEVEL - 1)
                                                    runLength = RUN_LENGTH_LEVEL - 1;
                                                if(startRun > RUN_INTENSITY_LEVEL - 1)
                                                    startRun = RUN_INTENSITY_LEVEL - 1;
                                                runlengthmatrix[p][runLength][startRun]++;
                                                startRun = -1;
                                                runLength = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (startRun != -1){
                                            if(runLength > RUN_LENGTH_LEVEL - 1)
                                                runLength = RUN_LENGTH_LEVEL - 1;
                                            if(startRun > RUN_INTENSITY_LEVEL - 1)
                                                startRun = RUN_INTENSITY_LEVEL - 1;
                                            runlengthmatrix[p][runLength][startRun]++;
                                            startRun = -1;
                                            runLength = 0;
                                        }
                                    }
                                }
                            }
                        }
                        // y direction
                        startRun = -1;
                        runLength = 0;
                        for(index[2] = m , roiIndex[2] = 0; index[2] < m + RUN_LENGTH_LEVEL; index[2]++, roiIndex[2]++) {
                            for(index[0] = r, roiIndex[0] = 0; index[0] < r + RUN_LENGTH_LEVEL; index[0]++, roiIndex[0]++) {
                                for(index[1] = n, roiIndex[1] = 0; index[1] < n + RUN_LENGTH_LEVEL; index[1]++, roiIndex[1]++) {
                                    if ( erodemask->GetPixel(index) == p ) {
                                        if (startRun == -1){
                                            startRun = roi->GetPixel(roiIndex);
                                            runLength++;
                                        }
                                        else {
                                            if(roi->GetPixel(roiIndex) == startRun){
                                                runLength++;
                                            }
                                            else {
                                                 if(runLength > RUN_LENGTH_LEVEL - 1)
                                                     runLength = RUN_LENGTH_LEVEL - 1;
                                                 if(startRun > RUN_INTENSITY_LEVEL - 1)
                                                     startRun = RUN_INTENSITY_LEVEL - 1;
                                                 runlengthmatrix[p][runLength][startRun]++;
                                                 startRun = -1;
                                                 runLength = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (startRun != -1){
                                            if(runLength > RUN_LENGTH_LEVEL - 1)
                                                runLength = RUN_LENGTH_LEVEL - 1;
                                            if(startRun > RUN_INTENSITY_LEVEL - 1)
                                                startRun = RUN_INTENSITY_LEVEL - 1;
                                            runlengthmatrix[p][runLength][startRun]++;
                                            startRun = -1;
                                            runLength = 0;
                                        }
                                    }
                                }
                            }
                        }
                        // z direction
                        startRun = -1;
                        runLength = 0;
                        for(index[1] = n, roiIndex[1] = 0; index[1] < n + RUN_LENGTH_LEVEL; index[1]++, roiIndex[1]++) {
                            for(index[0] = r, roiIndex[0] = 0; index[0] < r + RUN_LENGTH_LEVEL; index[0]++, roiIndex[0]++) {
                                for(index[2] = m , roiIndex[2] = 0; index[2] < m + RUN_LENGTH_LEVEL; index[2]++, roiIndex[2]++) {
                                    if ( erodemask->GetPixel(index) == p ) {
                                        if (startRun == -1){
                                            startRun = roi->GetPixel(roiIndex);
                                            runLength++;
                                        }
                                        else {
                                            if(roi->GetPixel(roiIndex) == startRun){
                                                runLength++;
                                            }
                                            else {
                                                 if(runLength > RUN_LENGTH_LEVEL - 1)
                                                     runLength = RUN_LENGTH_LEVEL - 1;
                                                 if(startRun > RUN_INTENSITY_LEVEL - 1)
                                                     startRun = RUN_INTENSITY_LEVEL - 1;
                                                 runlengthmatrix[p][runLength][startRun]++;
                                                 startRun = -1;
                                                 runLength = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (startRun != -1){
                                            if(runLength > RUN_LENGTH_LEVEL - 1)
                                                runLength = RUN_LENGTH_LEVEL - 1;
                                            if(startRun > RUN_INTENSITY_LEVEL - 1)
                                                startRun = RUN_INTENSITY_LEVEL - 1;
                                            runlengthmatrix[p][runLength][startRun]++;
                                            startRun = -1;
                                            runLength = 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //}
                    np = np * 3 ;
                    if(np != 0) {
                        //calRunLengthFeatures( SRETmp, LRETmp, GLNTmp, RLNTmp, RPTmp, np, runlengthmatrix[p] );
                        nroi++;
                        SRE += SRETmp;
                        LRE += LRETmp;
                        GLN += GLNTmp;
                        RLN += RLNTmp;
                        RP  += RPTmp;
                    }
                    for(int i = 0; i < RUN_LENGTH_LEVEL; i++){
                        for(int j = 0; j < RUN_INTENSITY_LEVEL; j++){
                            runlengthmatrix[p][i][j] = 0;
                        }
                    }
                }
            }
        }
        SRE = SRE / nroi;
        LRE = LRE / nroi;
        GLN = GLN / nroi;
        RLN = RLN / nroi;
        RP = RP / nroi;
        efile << muscleName[p] << "  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
        efile.close();
        std::cout << muscleName[p] << " finished" << std::endl;
    }
}
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::cooccurrenceFeat3DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType voxelSize, StrITKType caseID)
{
    //ROI based run length features
    float  cooccurrencematrix[NUMBER_MUSCLE + 1][CO_OCCURRENCE_LEVEL][CO_OCCURRENCE_LEVEL];
    int tempIntensity = 0, tempLength = 0, intensityMax = 0, intensityMin = 0, intensityMaxCut = 0, intensityMinCut = 0, histogram[HISTOGRAM_LIMIT] = {0};
    float intensityRescaleFactor = 0;
    DMDData ddata;
    DMDData::OrientedImageType::Pointer rescaledData = DMDData::OrientedImageType::New(), rescaledData_REC_FEM = DMDData::OrientedImageType::New(), rescaledData_SEMIT = DMDData::OrientedImageType::New(), rescaledData_CRAN_SART = DMDData::OrientedImageType::New(), rescaledData_BIC_FEM = DMDData::OrientedImageType::New(), rescaledData_ADDUCTOR = DMDData::OrientedImageType::New(), rescaledData_GRACILIS = DMDData::OrientedImageType::New(), rescaledData_SEMI_MEM = DMDData::OrientedImageType::New(), rescaledData_VAS_LAT = DMDData::OrientedImageType::New(), rescaledData_VAS_MED = DMDData::OrientedImageType::New(), rescaledData_VAS_INT = DMDData::OrientedImageType::New(), rescaledData_CAUD_SART = DMDData::OrientedImageType::New(), tempMaskX = DMDData::OrientedImageType::New(), tempMaskY = DMDData::OrientedImageType::New(), tempMaskZ = DMDData::OrientedImageType::New();
    float Entropy = 0, Energy = 0, Contrast = 0, HomoGeneity = 0;  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float EntropyData[NUMBER_OF_MUSCLE] = {0}, EnergyData[NUMBER_OF_MUSCLE] = {0}, ContrastData[NUMBER_OF_MUSCLE] = {0}, HomoGeneityData[NUMBER_OF_MUSCLE] = {0};  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float npRecFem = 0, npSemit = 0, npCranSart = 0, npBicFem = 0, npGracilis = 0, npAdductor = 0, npSemiMem = 0, npVasLat = 0, npVasMed = 0, npVasInt = 0, npCaudSart = 0, np = 0;
    double imageVolume = 0;
    DMDData::OrientedImageType::SpacingType    spacing = data->GetSpacing() ;
    float intervalX = round(ROI_SIZE / spacing[0]), intervalY = round(ROI_SIZE / spacing[1]), intervalZ = round(ROI_SIZE / spacing[2]);

    std::cout << "size of interval: " << intervalX << "  " << intervalY << "  " << intervalZ << std::endl;
    std::cout << "size of interval: " << spacing[0] << "  " << spacing[1] << "  " << spacing[2] << std::endl;
    //getchar();


    for (int j = 0; j < CO_OCCURRENCE_LEVEL; j++) {
        for (int i = 0 ; i < CO_OCCURRENCE_LEVEL ; i++) {
            cooccurrencematrix[REC_FEM][j][i] = 0;
            cooccurrencematrix[SEMIT][j][i] = 0;
            cooccurrencematrix[CRAN_SART][j][i] = 0;
            cooccurrencematrix[BIC_FEM][j][i] = 0;
            cooccurrencematrix[GRACILIS][j][i] = 0;
            cooccurrencematrix[ADDUCTOR][j][i] = 0;
            cooccurrencematrix[SEMI_MEM][j][i] = 0;
            cooccurrencematrix[VAS_LAT][j][i] = 0;
            cooccurrencematrix[VAS_MED][j][i] = 0;
            cooccurrencematrix[VAS_INT][j][i] = 0;
            cooccurrencematrix[CAUD_SART][j][i] = 0;
        }
    }
    ddata.imageInitialize ( data, rescaledData );
    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;
    IteratorType ReScaledDataIterator( rescaledData, rescaledData->GetRequestedRegion() ) ;
    // identify minum and maxmum intensity of the image
    constDataIterator.GoToBegin();
    intensityMin = intensityMax = constDataIterator.Get();
    for ( constDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator ) {
        int hisTmp = constDataIterator.Get();
        if (hisTmp >= HISTOGRAM_LIMIT)
            hisTmp = HISTOGRAM_LIMIT - 1;
        histogram[hisTmp]++;
        if(hisTmp > 0)
            imageVolume++;
        if (constDataIterator.Get() < intensityMin)
            intensityMin = constDataIterator.Get();
        if (constDataIterator.Get() > intensityMax)
            intensityMax = constDataIterator.Get();
    }
    // identify 5% minum and maxmum intensity of the image
    float cumImageVolume = 0;
    for (int i = 1; i <= intensityMax; i++){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMinCut = i;
            break;
        }
    }
    std::cout << imageVolume << "  " << cumImageVolume << "  ";
    cumImageVolume = 0;
    for (int i = intensityMax; i >= 1; i--){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMaxCut = i;
            break;
        }
    }
    //rescale images
    intensityRescaleFactor = (float)(intensityMaxCut - intensityMinCut) / CO_OCCURRENCE_LEVEL;
    for ( constDataIterator.GoToBegin(), ReScaledDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator, ++ReScaledDataIterator ) {
        tempIntensity = constDataIterator.Get();
        if (tempIntensity < intensityMinCut)
            tempIntensity = intensityMinCut;
        if (tempIntensity > intensityMaxCut)
            tempIntensity = intensityMaxCut;
        tempIntensity = round(tempIntensity / intensityRescaleFactor);
        if(tempIntensity < 0)
            tempIntensity = 0;
        if(tempIntensity >= CO_OCCURRENCE_LEVEL)
            tempIntensity = CO_OCCURRENCE_LEVEL - 1;
        ReScaledDataIterator.Set(tempIntensity);
    }
    std::cout << "rescale finished! " << std::endl;

    //establish run length matrix
    DMDData::OrientedImageType::SizeType   size = rescaledData->GetLargestPossibleRegion().GetSize();
    DMDData::OrientedImageType::IndexType  index;
    int upperBound[NUMBER_MUSCLE + 1] = {0}, lowerBound[NUMBER_MUSCLE + 1] = {0};
    // establish rescaled image for each muscle
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
   // for ( int p = CRAN_SART ; p <= SEMI_MEM; p++ ) {
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    if ( erodemask->GetPixel(index) == p ) {
                        if(rescaledData->GetPixel(index) < lowerBound[p]) {
                            lowerBound[p] = rescaledData->GetPixel(index);
                        }
                        if(rescaledData->GetPixel(index) > upperBound[p]){
                            upperBound[p] = rescaledData->GetPixel(index);
                        }
                    }
                }
            }
        }
    }
    DMDData::OrientedImageType::Pointer roi = DMDData::OrientedImageType::New();
    int nroi = 0; //the number of rois used for calculating features for each muscle.
    OrientedImageType::RegionType region;
    OrientedImageType::IndexType start;
    OrientedImageType::SizeType roiSize;
    OrientedImageType::IndexType roiIndex;
    roiSize[0] = intervalX, roiSize[1] = intervalY, roiSize[2] = intervalZ;
    start[0] = start[1] = start[2] = 0;
    region.SetSize(roiSize);
    region.SetIndex(start);
    roi->SetRegions(region);
    roi->Allocate();
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
    //for ( int p = CRAN_SART ; p <= SEMI_MEM; p++ ) {
        std::cout << "lower bound: " << lowerBound[p] << "  upper bound: " << upperBound[p] << std::endl;
        if( lowerBound[p] == 0 && upperBound[p] == 0)
            continue;
        nroi = 0;
        Entropy = Energy = Contrast = HomoGeneity = 0;
      //  std::ofstream efile( featurefilename.c_str() , std::ios::app );
        for (int m = 0; m < size[2]; m += round(intervalZ * 0.5)) {
            for (int n = 0; n < size[1]; n += round(intervalY * 0.5)) {
                for (int r = 0; r < size[0]; r += round(intervalX * 0.5)) {
                    float EntropyTmp = 0, EnergyTmp = 0, ContrastTmp = 0, HomoGeneityTmp = 0;
                    np = 0;
                    bool haveMuscleFlag = 0;
                    for (int i = lowerBound[p]; i <= upperBound[p]; i++) {
                        int objCount = 0;
                        bool inMuscleFlag = 0;
                        roi->FillBuffer(itk::NumericTraits<DMDData::OrientedImageType::PixelType>::Zero);
                        for(index[2] = m , roiIndex[2] = 0; index[2] < m + intervalZ; index[2]++, roiIndex[2]++) {
                            for(index[1] = n, roiIndex[1] = 0; index[1] < n + intervalY; index[1]++, roiIndex[1]++) {
                                for(index[0] = r, roiIndex[0] = 0; index[0] < r + intervalX; index[0]++, roiIndex[0]++) {
                                    if ((index[0] >= 0) && (index[0] < size[0]) && (index[1] >= 0) && (index[1] < size[1]) && (index[2] >= 0) && (index[2] < size[2])) {
                                        if(roiIndex[0] < size[0] && roiIndex[1] < size[1] && roiIndex[2] < size[2]) {
                                            if ( erodemask->GetPixel(index) == p && rescaledData->GetPixel(index) == i) {
                                                roi->SetPixel(roiIndex, rescaledData->GetPixel(index));
                                                inMuscleFlag = 1;
                                            }
                                            else
                                                roi->SetPixel(roiIndex, BACKGROUND);
                                            if(erodemask->GetPixel(index) == p)
                                                haveMuscleFlag = 1;
                                        }
                                        else
                                            roi->SetPixel(roiIndex, BACKGROUND);
                                    }
                                }
                            }
                        }
                        if(haveMuscleFlag == 0 && inMuscleFlag == 0)
                            break;
                        if(haveMuscleFlag == 1 && inMuscleFlag == 0)
                            continue;
                        if(inMuscleFlag == 1) {
                            for(index[2] = m , roiIndex[2] = 0; index[2] < m + intervalZ; index[2]++, roiIndex[2]++) {
                                for(index[1] = n, roiIndex[1] = 0; index[1] < n + intervalY; index[1]++, roiIndex[1]++) {
                                    for(index[0] = r, roiIndex[0] = 0; index[0] < r + intervalX; index[0]++, roiIndex[0]++) {
                                        int vali = 0, valj = 0;
                                        // x direction
                                        if ( (roiIndex[0] + 2) < intervalX ) {
                                            vali = roi->GetPixel(roiIndex);
                                            roiIndex[0] += 2;
                                            valj = roi->GetPixel(roiIndex);
                                            roiIndex[0] -= 2;
                                            if (vali >= 0 && vali < CO_OCCURRENCE_LEVEL && valj >= 0 && valj < CO_OCCURRENCE_LEVEL) {
                                                cooccurrencematrix[p][valj][vali]++;
                                                cooccurrencematrix[p][vali][valj]++;
                                            }
                                        }
                                        // y direction
                                        if ( (roiIndex[1] + 2) < intervalY ) {
                                            vali = roi->GetPixel(roiIndex);
                                            roiIndex[1] += 2;
                                            valj = roi->GetPixel(roiIndex);
                                            roiIndex[1] -= 2;
                                            if (vali >= 0 && vali < CO_OCCURRENCE_LEVEL && valj >= 0 && valj < CO_OCCURRENCE_LEVEL) {
                                                cooccurrencematrix[p][valj][vali]++;
                                                cooccurrencematrix[p][vali][valj]++;
                                            }
                                        }
                                        // z direction
                                        if ( (roiIndex[2] + 2) < intervalZ ) {
                                            vali = roi->GetPixel(roiIndex);
                                            roiIndex[2] += 2;
                                            valj = roi->GetPixel(roiIndex);
                                            roiIndex[2] -= 2;
                                            if (vali >= 0 && vali < CO_OCCURRENCE_LEVEL && valj >= 0 && valj < CO_OCCURRENCE_LEVEL) {
                                                cooccurrencematrix[p][valj][vali]++;
                                                cooccurrencematrix[p][vali][valj]++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(haveMuscleFlag) {
                        calCooccurrenceFeatures( EntropyTmp, EnergyTmp, ContrastTmp, HomoGeneityTmp, cooccurrencematrix[p], caseID );
                        nroi++;
                        //std::cout << "mark1" << SRETmp << "  " << LRETmp << "  " << GLNTmp << "  " << RLNTmp << "  " << RPTmp << std::endl;
                        Entropy += EntropyTmp;
                        Energy += EnergyTmp;
                        Contrast += ContrastTmp;
                        HomoGeneity += HomoGeneityTmp;

                        std::string matrixfilename = "../data/" + caseID + "_matrix.txt";
                        std::ofstream efilex( matrixfilename.c_str() , std::ios::app );
                        for(int i = 0; i < CO_OCCURRENCE_LEVEL; i++){
                            for(int j = 0; j < CO_OCCURRENCE_LEVEL; j++){
                                efilex << cooccurrencematrix[p][i][j] << "  ";
                            }
                            efilex << "\n";
                        }
                        efilex.close();
                        for(int i = 0; i < CO_OCCURRENCE_LEVEL; i++){
                            for(int j = 0; j < CO_OCCURRENCE_LEVEL; j++){
                                cooccurrencematrix[p][i][j] = 0;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "size of " << muscleName[p] << ": " << nroi << std::endl;
        if(nroi) {
            Entropy  /=  nroi;
            Energy /= nroi;
            Contrast /= nroi;
            HomoGeneity /= nroi;
        }
       // efile << muscleName[p] << "  " << Entropy << "  " <<  Energy << "  " << Contrast << "  " << HomoGeneity << "\n";
      //  efile.close();
        EntropyData[p] = Entropy; EnergyData[p] = Energy; ContrastData[p] = Contrast; HomoGeneityData[p] = HomoGeneity;
        std::cout << muscleName[p] << " finished" << std::endl;
    }

    std::ofstream efile( featurefilename.c_str() , std::ios::app );
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << muscleName[p] << "  " ;
    }
    efile << "\nEntropy";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << EntropyData[p];
    }
    efile << "\nEnergy";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << EnergyData[p];
    }
    efile << "\nContrast";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << ContrastData[p];
    }
    efile << "\nHomoGeneity";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << HomoGeneityData[p];
    }
    efile.close();

}
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::runlengthFeat3DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType voxelSize, StrITKType caseID)
{
    //ROI based run length features
    int runlengthmatrix[NUMBER_MUSCLE + 1][RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL];
    int tempIntensity = 0, tempLength = 0, intensityMax = 0, intensityMin = 0, intensityMaxCut = 0, intensityMinCut = 0, histogram[HISTOGRAM_LIMIT] = {0};
    float intensityRescaleFactor = 0;
    DMDData ddata;
    DMDData::OrientedImageType::Pointer rescaledData = DMDData::OrientedImageType::New(), rescaledData_REC_FEM = DMDData::OrientedImageType::New(), rescaledData_SEMIT = DMDData::OrientedImageType::New(), rescaledData_CRAN_SART = DMDData::OrientedImageType::New(), rescaledData_BIC_FEM = DMDData::OrientedImageType::New(), rescaledData_ADDUCTOR = DMDData::OrientedImageType::New(), rescaledData_GRACILIS = DMDData::OrientedImageType::New(), rescaledData_SEMI_MEM = DMDData::OrientedImageType::New(), rescaledData_VAS_LAT = DMDData::OrientedImageType::New(), rescaledData_VAS_MED = DMDData::OrientedImageType::New(), rescaledData_VAS_INT = DMDData::OrientedImageType::New(), rescaledData_CAUD_SART = DMDData::OrientedImageType::New(), tempMaskX = DMDData::OrientedImageType::New(), tempMaskY = DMDData::OrientedImageType::New(), tempMaskZ = DMDData::OrientedImageType::New();
    float SRE = 0, LRE = 0, GLN = 0, RLN = 0, RP = 0;  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float SREData[NUMBER_OF_MUSCLE] = {0}, LREData[NUMBER_OF_MUSCLE] = {0}, GLNData[NUMBER_OF_MUSCLE] = {0}, RLNData[NUMBER_OF_MUSCLE] = {0}, RPData[NUMBER_OF_MUSCLE] = {0};  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float npRecFem = 0, npSemit = 0, npCranSart = 0, npBicFem = 0, npGracilis = 0, npAdductor = 0, npSemiMem = 0, npVasLat = 0, npVasMed = 0, npVasInt = 0, npCaudSart = 0, np = 0;
    double imageVolume = 0;
    DMDData::OrientedImageType::SpacingType    spacing = data->GetSpacing() ;
    float intervalX = round(ROI_SIZE / spacing[0]), intervalY = round(ROI_SIZE / spacing[1]), intervalZ = round(ROI_SIZE / spacing[2]);

    std::cout << "size of interval: " << intervalX << "  " << intervalY << "  " << intervalZ << std::endl;
    std::cout << "size of interval: " << spacing[0] << "  " << spacing[1] << "  " << spacing[2] << std::endl;
    //getchar();


    for (int j = 0; j < RUN_LENGTH_LEVEL; j++) {
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
            runlengthmatrix[REC_FEM][j][i] = 0;
            runlengthmatrix[SEMIT][j][i] = 0;
            runlengthmatrix[CRAN_SART][j][i] = 0;
            runlengthmatrix[BIC_FEM][j][i] = 0;
            runlengthmatrix[GRACILIS][j][i] = 0;
            runlengthmatrix[ADDUCTOR][j][i] = 0;
            runlengthmatrix[SEMI_MEM][j][i] = 0;
            runlengthmatrix[VAS_LAT][j][i] = 0;
            runlengthmatrix[VAS_MED][j][i] = 0;
            runlengthmatrix[VAS_INT][j][i] = 0;
            runlengthmatrix[CAUD_SART][j][i] = 0;
        }
    }
    ddata.imageInitialize ( data, rescaledData );
    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;
    IteratorType ReScaledDataIterator( rescaledData, rescaledData->GetRequestedRegion() ) ;
    // identify minum and maxmum intensity of the image
    constDataIterator.GoToBegin();
    intensityMin = intensityMax = constDataIterator.Get();
    for ( constDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator ) {
        int hisTmp = constDataIterator.Get();
        if (hisTmp >= HISTOGRAM_LIMIT)
            hisTmp = HISTOGRAM_LIMIT - 1;
        histogram[hisTmp]++;
        if(hisTmp > 0)
            imageVolume++;
        if (constDataIterator.Get() < intensityMin)
            intensityMin = constDataIterator.Get();
        if (constDataIterator.Get() > intensityMax)
            intensityMax = constDataIterator.Get();
    }
    // identify 5% minum and maxmum intensity of the image
    float cumImageVolume = 0;
    for (int i = 1; i <= intensityMax; i++){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMinCut = i;
            break;
        }
    }
    std::cout << imageVolume << "  " << cumImageVolume << "  ";
    cumImageVolume = 0;
    for (int i = intensityMax; i >= 1; i--){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMaxCut = i;
            break;
        }
    }
    //rescale images
    intensityRescaleFactor = (float)(intensityMaxCut - intensityMinCut) / RUN_INTENSITY_LEVEL;
    for ( constDataIterator.GoToBegin(), ReScaledDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator, ++ReScaledDataIterator ) {
        tempIntensity = constDataIterator.Get();
        if (tempIntensity < intensityMinCut)
            tempIntensity = intensityMinCut;
        if (tempIntensity > intensityMaxCut)
            tempIntensity = intensityMaxCut;
        tempIntensity = round(tempIntensity / intensityRescaleFactor);
        if(tempIntensity < 0)
            tempIntensity = 0;
        if(tempIntensity >= RUN_INTENSITY_LEVEL)
            tempIntensity = RUN_INTENSITY_LEVEL - 1;
        ReScaledDataIterator.Set(tempIntensity);
    }
    std::cout << "rescale finished! " << std::endl;

    //establish run length matrix
    DMDData::OrientedImageType::SizeType   size = rescaledData->GetLargestPossibleRegion().GetSize();
    DMDData::OrientedImageType::IndexType  index;
    int upperBound[NUMBER_MUSCLE + 1] = {0}, lowerBound[NUMBER_MUSCLE + 1] = {0};
    // establish rescaled image for each muscle
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
   // for ( int p = CRAN_SART ; p <= SEMI_MEM; p++ ) {
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    if ( erodemask->GetPixel(index) == p ) {
                        if(rescaledData->GetPixel(index) < lowerBound[p]) {
                            lowerBound[p] = rescaledData->GetPixel(index);
                        }
                        if(rescaledData->GetPixel(index) > upperBound[p]){
                            upperBound[p] = rescaledData->GetPixel(index);
                        }
                    }
                }
            }
        }
    }
    DMDData::OrientedImageType::Pointer roi = DMDData::OrientedImageType::New();
    int nroi = 0; //the number of rois used for calculating features for each muscle.
    OrientedImageType::RegionType region;
    OrientedImageType::IndexType start;
    OrientedImageType::SizeType roiSize;
    OrientedImageType::IndexType roiIndex;
    roiSize[0] = intervalX, roiSize[1] = intervalY, roiSize[2] = intervalZ;
    start[0] = start[1] = start[2] = 0;
    region.SetSize(roiSize);
    region.SetIndex(start);
    roi->SetRegions(region);
    roi->Allocate();
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
   // for ( int p = CRAN_SART ; p <= SEMI_MEM; p++ ) {
        std::cout << "lower bound: " << lowerBound[p] << "  upper bound: " << upperBound[p] << std::endl;
        if( lowerBound[p] == 0 && upperBound[p] == 0)
            continue;
        nroi = 0;
        SRE = LRE = GLN = RLN = RP = 0;
       // std::ofstream efile( featurefilename.c_str() , std::ios::app );
        for (int m = 0; m < size[2]; m += round(intervalZ * 0.5)) {
            for (int n = 0; n < size[1]; n += round(intervalY * 0.5)) {
                for (int r = 0; r < size[0]; r += round(intervalX * 0.5)) {
                    float SRETmp = 0, LRETmp = 0, GLNTmp = 0, RLNTmp = 0, RPTmp = 0;
                    np = 0;
                    bool haveMuscleFlag = 0;
                    for (int i = lowerBound[p]; i <= upperBound[p]; i++) {
                        int objCount = 0;
                        bool inMuscleFlag = 0;
                        roi->FillBuffer(itk::NumericTraits<DMDData::OrientedImageType::PixelType>::Zero);
                        for(index[2] = m , roiIndex[2] = 0; index[2] < m + intervalZ; index[2]++, roiIndex[2]++) {
                            for(index[1] = n, roiIndex[1] = 0; index[1] < n + intervalY; index[1]++, roiIndex[1]++) {
                                for(index[0] = r, roiIndex[0] = 0; index[0] < r + intervalX; index[0]++, roiIndex[0]++) {
                                    if ((index[0] >= 0) && (index[0] < size[0]) && (index[1] >= 0) && (index[1] < size[1]) && (index[2] >= 0) && (index[2] < size[2])) {
                                        if(roiIndex[0] < size[0] && roiIndex[1] < size[1] && roiIndex[2] < size[2]) {
                                            if ( erodemask->GetPixel(index) == p && rescaledData->GetPixel(index) == i) {
                                                roi->SetPixel(roiIndex, rescaledData->GetPixel(index));
                                                inMuscleFlag = 1;
                                            }
                                            else
                                                roi->SetPixel(roiIndex, BACKGROUND);
                                            if(erodemask->GetPixel(index) == p)
                                                haveMuscleFlag = 1;
                                        }
                                        else
                                            roi->SetPixel(roiIndex, BACKGROUND);
                                    }
                                }
                            }
                        }
                        if(haveMuscleFlag == 0 && inMuscleFlag == 0)
                            break;
                        if(haveMuscleFlag == 1 && inMuscleFlag == 0)
                            continue;
                        if(inMuscleFlag == 1) {
                            objCount = connectedComponentLabeling (roi);
                            int objVolume[SIZE_LABELING_BIN] = {0};
                            for(roiIndex[2] = 0 ; roiIndex[2] < intervalZ; roiIndex[2]++) {
                                for(roiIndex[1] = 0; roiIndex[1] < intervalY; roiIndex[1]++) {
                                    for(roiIndex[0] = 0; roiIndex[0] < intervalX; roiIndex[0]++) {
                                        int k = roi->GetPixel(roiIndex);
                                        if(k){
                                            objVolume[k]++;
                                        }
                                    }
                                }
                            }
                            for (int j = 0; j < SIZE_LABELING_BIN; j++) {
                                if (objVolume[j]) {
                                    objVolume[j] = round(objVolume[j] * voxelSize);
                                }
                                if (objVolume[j]) {
                                    if(objVolume[j] >= RUN_LENGTH_LEVEL){
                                        runlengthmatrix[p][RUN_LENGTH_LEVEL - 1][i]++;
                                        np += objVolume[j];
                                    }
                                    else{
                                        if(objVolume[j]){
                                            runlengthmatrix[p][objVolume[j]][i]++;
                                            np += objVolume[j];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(haveMuscleFlag) {
                        calRunLengthFeatures( SRETmp, LRETmp, GLNTmp, RLNTmp, RPTmp, np, runlengthmatrix[p], caseID );
                        nroi++;
                        //std::cout << "mark1" << SRETmp << "  " << LRETmp << "  " << GLNTmp << "  " << RLNTmp << "  " << RPTmp << std::endl;
                        SRE += SRETmp;
                        LRE += LRETmp;
                        GLN += GLNTmp;
                        RLN += RLNTmp;
                        RP  += RPTmp;
                        std::string matrixfilename = "../data/" + caseID + "_matrix.txt";
                        std::ofstream efilex( matrixfilename.c_str() , std::ios::app );
                        for(int i = 0; i < RUN_LENGTH_LEVEL; i++){
                            for(int j = 0; j < RUN_INTENSITY_LEVEL; j++){
                                efilex << runlengthmatrix[p][i][j] << "  ";
                            }
                            efilex << "\n";
                        }
                        efilex.close();
                        for(int i = 0; i < RUN_LENGTH_LEVEL; i++){
                            for(int j = 0; j < RUN_INTENSITY_LEVEL; j++){
                                runlengthmatrix[p][i][j] = 0;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "size of " << muscleName[p] << ": " << nroi << std::endl;
        if(nroi) {
            SRE = SRE / nroi;
            LRE = LRE / nroi;
            GLN = GLN / nroi;
            RLN = RLN / nroi;
            RP = RP / nroi;
        }
        SREData[p] = SRE; LREData[p] = LRE; GLNData[p] = GLN; RLNData[p] = RLN; RPData[p] = RP;
    //    efile << muscleName[p] << "  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
     //   efile.close();
        std::cout << muscleName[p] << " finished" << std::endl;
    }

    std::ofstream efile( featurefilename.c_str() , std::ios::app );
    efile << caseID << ":\n" ;
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << muscleName[p] << "  " ;
    }
    efile << "\nSRE";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << SREData[p];
    }
    efile << "\nLRE";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << LREData[p];
    }
    efile << "\nGLN";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << GLNData[p];
    }
    efile << "\nRLN";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << RLNData[p];
    }
    efile << "\nRP";
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        efile << "   " << RPData[p];
    }
    efile.close();
}
/*
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::runlengthFeat3DROI( DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, std::string featurefilename, FITKType voxelSize, StrITKType caseID)
{
    //ROI based run length features
    int runlengthmatrix[NUMBER_MUSCLE + 1][RUN_LENGTH_LEVEL][RUN_INTENSITY_LEVEL];
    int tempIntensity = 0, tempLength = 0, intensityMax = 0, intensityMin = 0, intensityMaxCut = 0, intensityMinCut = 0, histogram[HISTOGRAM_LIMIT] = {0};
    float intensityRescaleFactor = 0;
    DMDData ddata;
    DMDData::OrientedImageType::Pointer rescaledData = DMDData::OrientedImageType::New(), rescaledData_REC_FEM = DMDData::OrientedImageType::New(), rescaledData_SEMIT = DMDData::OrientedImageType::New(), rescaledData_CRAN_SART = DMDData::OrientedImageType::New(), rescaledData_BIC_FEM = DMDData::OrientedImageType::New(), rescaledData_ADDUCTOR = DMDData::OrientedImageType::New(), rescaledData_GRACILIS = DMDData::OrientedImageType::New(), rescaledData_SEMI_MEM = DMDData::OrientedImageType::New(), rescaledData_VAS_LAT = DMDData::OrientedImageType::New(), rescaledData_VAS_MED = DMDData::OrientedImageType::New(), rescaledData_VAS_INT = DMDData::OrientedImageType::New(), rescaledData_CAUD_SART = DMDData::OrientedImageType::New(), tempMaskX = DMDData::OrientedImageType::New(), tempMaskY = DMDData::OrientedImageType::New(), tempMaskZ = DMDData::OrientedImageType::New();
    float SRE = 0, LRE = 0, GLN = 0, RLN = 0, RP = 0;  // SRE = Short Run Emphasis, LRE = Long Run Emphasis, GLN = Gray Level Nonuniform, RLN = Run Length Nonuniform, RP = Run Percentage
    float npRecFem = 0, npSemit = 0, npCranSart = 0, npBicFem = 0, npGracilis = 0, npAdductor = 0, npSemiMem = 0, npVasLat = 0, npVasMed = 0, npVasInt = 0, npCaudSart = 0, np = 0;
    double imageVolume = 0;
    DMDData::OrientedImageType::SpacingType    spacing = data->GetSpacing() ;
    float intervalX = round(ROI_SIZE / spacing[0]), intervalY = round(ROI_SIZE / spacing[1]), intervalZ = round(ROI_SIZE / spacing[2]);

   // std::cout << intervalX << "  " << intervalY << "  " << intervalZ << std::endl;
  //  std::cout << spacing[0] << "  " << spacing[1] << "  " << spacing[2] << std::endl;
  //  getchar();

    for (int j = 0; j < RUN_LENGTH_LEVEL; j++) {
        for (int i = 0; i < RUN_INTENSITY_LEVEL; i++) {
            runlengthmatrix[REC_FEM][j][i] = 0;
            runlengthmatrix[SEMIT][j][i] = 0;
            runlengthmatrix[CRAN_SART][j][i] = 0;
            runlengthmatrix[BIC_FEM][j][i] = 0;
            runlengthmatrix[GRACILIS][j][i] = 0;
            runlengthmatrix[ADDUCTOR][j][i] = 0;
            runlengthmatrix[SEMI_MEM][j][i] = 0;
            runlengthmatrix[VAS_LAT][j][i] = 0;
            runlengthmatrix[VAS_MED][j][i] = 0;
            runlengthmatrix[VAS_INT][j][i] = 0;
            runlengthmatrix[CAUD_SART][j][i] = 0;
        }
    }
    ddata.imageInitialize ( data, rescaledData );
    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;
    IteratorType ReScaledDataIterator( rescaledData, rescaledData->GetRequestedRegion() ) ;
    // identify minum and maxmum intensity of the image
    constDataIterator.GoToBegin();
    intensityMin = intensityMax = constDataIterator.Get();
    for ( constDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator ) {
        int hisTmp = constDataIterator.Get();
        if (hisTmp >= HISTOGRAM_LIMIT)
            hisTmp = HISTOGRAM_LIMIT - 1;
        histogram[hisTmp]++;
        if(hisTmp > 0)
            imageVolume++;
        if (constDataIterator.Get() < intensityMin)
            intensityMin = constDataIterator.Get();
        if (constDataIterator.Get() > intensityMax)
            intensityMax = constDataIterator.Get();
    }
    // identify 5% minum and maxmum intensity of the image
    float cumImageVolume = 0;
    for (int i = 1; i <= intensityMax; i++){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMinCut = i;
            break;
        }
    }
    std::cout << imageVolume << "  " << cumImageVolume << "  ";
    cumImageVolume = 0;
    for (int i = intensityMax; i >= 1; i--){
        cumImageVolume += histogram[i];
        if((cumImageVolume / imageVolume) >= INTENSITY_CUT_RATE){
            intensityMaxCut = i;
            break;
        }
    }
    //rescale images
    intensityRescaleFactor = (float)(intensityMaxCut - intensityMinCut) / RUN_INTENSITY_LEVEL;
    for ( constDataIterator.GoToBegin(), ReScaledDataIterator.GoToBegin(); !constDataIterator.IsAtEnd(); ++constDataIterator, ++ReScaledDataIterator ) {
        tempIntensity = constDataIterator.Get();
        if (tempIntensity < intensityMinCut)
            tempIntensity = intensityMinCut;
        if (tempIntensity > intensityMaxCut)
            tempIntensity = intensityMaxCut;
        tempIntensity = round(tempIntensity / intensityRescaleFactor);
        if(tempIntensity < 0)
            tempIntensity = 0;
        if(tempIntensity >= RUN_INTENSITY_LEVEL)
            tempIntensity = RUN_INTENSITY_LEVEL - 1;
        ReScaledDataIterator.Set(tempIntensity);
    }
    std::cout << "rescale finished! " << std::endl;
    //establish run length matrix
    DMDData::OrientedImageType::SizeType   size = data->GetLargestPossibleRegion().GetSize();
    DMDData::OrientedImageType::IndexType  index;
    int upperBound[NUMBER_MUSCLE + 1] = {0}, lowerBound[NUMBER_MUSCLE + 1] = {0};
    // establish rescaled image for each muscle
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        for(index[2] = 0 ; index[2] < size[2]; index[2]++) {
            for(index[1] = 0; index[1] < size[1]; index[1]++) {
                for(index[0] = 0; index[0] < size[0]; index[0]++) {
                    if ( erodemask->GetPixel(index) == p ) {
                        if(rescaledData->GetPixel(index) < lowerBound[p]) {
                            lowerBound[p] = rescaledData->GetPixel(index);
                        }
                        if(rescaledData->GetPixel(index) > upperBound[p]){
                            upperBound[p] = rescaledData->GetPixel(index);
                        }

                    }
                }
            }
        }
    }
    DMDData::OrientedImageType::Pointer roi = DMDData::OrientedImageType::New();
    int nroi = 0; //the number of rois used for calculating features for each muscle.
    OrientedImageType::RegionType region;
    OrientedImageType::IndexType start;
    OrientedImageType::SizeType roiSize;
    OrientedImageType::IndexType roiIndex;
    roiSize[0] = intervalX, roiSize[1] = intervalY, roiSize[2] = intervalZ;
    start[0] = start[1] = start[2] = 0;
    region.SetSize(roiSize);
    region.SetIndex(start);
    roi->SetRegions(region);
    roi->Allocate();
    for ( int p = CRAN_SART ; p <= CAUD_SART; p++ ) {
        std::cout << "lower bound: " << lowerBound[p] << "  upper bound: " << upperBound[p] << std::endl;
        if( lowerBound[p] == 0 && upperBound[p] == 0)
            continue;
        nroi = 0;
        SRE = LRE = GLN = RLN = RP = 0;
        std::ofstream efile( featurefilename.c_str() , std::ios::app );
        for (int m = 0; m < size[2]; m += round(intervalZ * 0.5)) {
            for (int n = 0; n < size[1]; n += round(intervalY * 0.5)) {
                for (int r = 0; r < size[0]; r += round(intervalX * 0.5)) {
                    float SRETmp = 0, LRETmp = 0, GLNTmp = 0, RLNTmp = 0, RPTmp = 0;
                    np = 0;
                    bool haveMuscleFlag = 0;
                    for (int i = lowerBound[p]; i <= upperBound[p]; i++) {
                        int objCount = 0;
                        bool inMuscleFlag = 0;
                    //    haveMuscleFlag = 0;
                        roi->FillBuffer(itk::NumericTraits<DMDData::OrientedImageType::PixelType>::Zero);
                        for(index[2] = m , roiIndex[2] = 0; index[2] < m + intervalZ; index[2]++, roiIndex[2]++) {
                            for(index[1] = n, roiIndex[1] = 0; index[1] < n + intervalY; index[1]++, roiIndex[1]++) {
                                for(index[0] = r, roiIndex[0] = 0; index[0] < r + intervalX; index[0]++, roiIndex[0]++) {
                                    if(roiIndex[0] < size[0] && roiIndex[1] < size[1] && roiIndex[2] < size[2]) {
                                        if ( erodemask->GetPixel(index) == p && rescaledData->GetPixel(index) == i) {
                                            roi->SetPixel(roiIndex, rescaledData->GetPixel(index));
                                            inMuscleFlag = 1;
                                            //np++;
                                        }
                                        else
                                            roi->SetPixel(roiIndex, BACKGROUND);
                                        if(erodemask->GetPixel(index) == p)
                                            haveMuscleFlag = 1;
                                    }
                                }
                            }
                        }
                  //      std::cout << "flags: " << i << "-  " << haveMuscleFlag << "  " << inMuscleFlag << std::endl;
                        if(haveMuscleFlag == 0 && inMuscleFlag == 0)
                            break;
                        if(haveMuscleFlag == 1 && inMuscleFlag == 0)
                            continue;

                        if(inMuscleFlag == 1) {
                            objCount = connectedComponentLabeling (roi);
                            int objVolume[SIZE_LABELING_BIN] = {0};
                            for(roiIndex[2] = 0 ; roiIndex[2] < intervalZ; roiIndex[2]++) {
                                for(roiIndex[1] = 0; roiIndex[1] < intervalY; roiIndex[1]++) {
                                    for(roiIndex[0] = 0; roiIndex[0] < intervalX; roiIndex[0]++) {
                                        int k = roi->GetPixel(roiIndex);
                                        if(k){
                                            objVolume[k]++;
                                        }
                                    }
                                }
                            }
                            for (int j = 0; j < SIZE_LABELING_BIN; j++) {
                                if (objVolume[j]) {
                                    //std::cout << j << ":  " << objVolume[j] << "  "  ;
                                    objVolume[j] = round(objVolume[j] * voxelSize);
                                   // std::cout << "  " << objVolume[j] << "  " << voxelSize  << std::endl;
                                }
                                if (objVolume[j]) {
                                    if(objVolume[j] >= RUN_LENGTH_LEVEL){
                                        runlengthmatrix[p][RUN_LENGTH_LEVEL - 1][i]++;
                                        np += objVolume[j];
                                    }
                                    else{
                                        if(objVolume[j]){
                                            runlengthmatrix[p][objVolume[j]][i]++;
                                            np += objVolume[j];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(haveMuscleFlag) {
                   // np = np * voxelSize;
                   // if(np != 0) {
         //               std::cout << "number of pixel: " << np << std::endl;
                        calRunLengthFeatures( SRETmp, LRETmp, GLNTmp, RLNTmp, RPTmp, np, runlengthmatrix[p], caseID );
                        nroi++;
                        SRE += SRETmp;
                        LRE += LRETmp;
                        GLN += GLNTmp;
                        RLN += RLNTmp;
                        RP  += RPTmp;
                        std::string matrixfilename = caseID + "_matrix.txt";
                        std::ofstream efilex( matrixfilename.c_str() , std::ios::app );
                        for(int i = 0; i < RUN_LENGTH_LEVEL; i++){
                            for(int j = 0; j < RUN_INTENSITY_LEVEL; j++){
                         //       if(runlengthmatrix[p][i][j] != 0)
                        //            std::cout << "rlm element: " << runlengthmatrix[p][i][j] << std::endl;
                                efilex << runlengthmatrix[p][i][j] << "  ";
                            }
                            efilex << "\n";
                        }
                        efilex.close();
                        //std::cout << "flag: " << haveMuscleFlag << std::endl;
                        //std::cout << SRETmp << "  " << LRETmp << "  " << GLN << "  " << RLN << "  " << RP << std::endl;
                        //getchar();
                    //}
                        for(int i = 0; i < RUN_LENGTH_LEVEL; i++){
                            for(int j = 0; j < RUN_INTENSITY_LEVEL; j++){
                                runlengthmatrix[p][i][j] = 0;
                            }
                        }
                    }
                }
            }
        }
        SRE = SRE / nroi;
        LRE = LRE / nroi;
        GLN = GLN / nroi;
        RLN = RLN / nroi;
        RP = RP / nroi;
        efile << muscleName[p] << "  " << SRE << "  " <<  LRE << "  " << GLN << "  " << RLN << "  " << RP << "\n";
        efile.close();
        std::cout << muscleName[p] << " finished" << std::endl;
    }
}
*/
//////////////////////////////////////////////////////////////////////////
void DMDMuscleFeature::features( DMDData::OrientedImageType::Pointer mask, DMDData::OrientedImageType::Pointer erodemask, DMDData::OrientedImageType::Pointer data, float voxelSize, std::string featurefilename, bool volOutputMark )
{
    // if volOutputMark = TRUE, then output the volume of muscles
    int volSemitReal = 0, volRectReal = 0, volCranReal = 0, volAdductorReal = 0, volGracilisReal = 0, volBicReal = 0, volSemimemReal = 0, volVaslatReal = 0, volVasmedReal = 0, volVasintReal = 0, volCaudReal = 0;
    int volSemitNonErod = 0, volRectNonErod = 0, volCranNonErod = 0, volAdductorNonErod = 0, volGracilisNonErod = 0, volBicNonErod = 0, volSemimemNonErod = 0, volVaslatNonErod = 0, volVasmedNonErod = 0, volVasintNonErod = 0, volCaudNonErod = 0;
    float quanlitySemit = 0, quanlityRect = 0, quanlityCran = 0, quanlityAdductor = 0, quanlityGracilis = 0, quanlityBic = 0, quanlitySemimem = 0, quanlityVaslat = 0, quanlityVasmed = 0, quanlityVasint = 0, quanlityCaud = 0;
    float histogramSemit[HISTOGRAM_BIN_SIZE] = {0}, histogramRect[HISTOGRAM_BIN_SIZE] = {0},histogramCran[HISTOGRAM_BIN_SIZE] = {0}, histogramAdductor[HISTOGRAM_BIN_SIZE] = {0}, histogramGracilis[HISTOGRAM_BIN_SIZE] = {0}, histogramBic[HISTOGRAM_BIN_SIZE] = {0}, histogramSemimem[HISTOGRAM_BIN_SIZE] = {0}, histogramVaslat[HISTOGRAM_BIN_SIZE] = {0}, histogramVasmed[HISTOGRAM_BIN_SIZE] = {0}, histogramVasint[HISTOGRAM_BIN_SIZE] = {0}, histogramCaud[HISTOGRAM_BIN_SIZE] = {0} ;
    std::ofstream efile( featurefilename.c_str() , std::ios::app );
    ConstIteratorType constMaskIterator( mask, mask->GetRequestedRegion() ) ;
    ConstIteratorType constErodeMaskIterator( erodemask, erodemask->GetRequestedRegion() ) ;
    ConstIteratorType constDataIterator( data, data->GetRequestedRegion() ) ;

    volSemit = volRect = volCran = volAdductor = volGracilis = volBic = volSemimem = volVaslat = volVasmed = volVasint = volCaud = 0;

    // calculate the mean and volume of muscles using eroded muscle mask
    for ( constMaskIterator.GoToBegin(),constDataIterator.GoToBegin(); !constMaskIterator.IsAtEnd(); ++constMaskIterator, ++constDataIterator ) {
        if ( constMaskIterator.Get() == SEMIT ) {
            if (constDataIterator.Get() > 0) {
                volSemitNonErod++;
            }
        }
        if ( constMaskIterator.Get() == REC_FEM ) {
            if (constDataIterator.Get() > 0) {
                volRectNonErod++;
            }
        }
        if ( constMaskIterator.Get() == CRAN_SART ) {
            if (constDataIterator.Get() > 0) {
                volCranNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == ADDUCTOR ) {
            if (constDataIterator.Get() > 0) {
                volAdductorNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == GRACILIS ) {
            if (constDataIterator.Get() > 0) {
                volGracilisNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == BIC_FEM ) {
            if (constDataIterator.Get() > 0) {
                volBicNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == SEMI_MEM ) {
            if (constDataIterator.Get() > 0) {
                volSemimemNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == VAS_LAT ) {
            if (constDataIterator.Get() > 0) {
                volVaslatNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == VAS_MED ) {
            if (constDataIterator.Get() > 0) {
                volVasmedNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == VAS_INT ) {
            if (constDataIterator.Get() > 0) {
                volVasintNonErod++;
            }
        }   	       	
        if ( constMaskIterator.Get() == CAUD_SART ) {
        //    if (constDataIterator.Get() > 0) {
                volCaudNonErod++;
         //   }
        }   	       	
    }
    // calculate the mean and volume of muscles using eroded muscle mask
    // if the T2 fit only covers the mid section, only caluculate the biomarkers in covered part
    int tmp = 0;
    for ( constErodeMaskIterator.GoToBegin(),constDataIterator.GoToBegin(); !constErodeMaskIterator.IsAtEnd(); ++constErodeMaskIterator, ++constDataIterator ) {
        if ( constErodeMaskIterator.Get() == SEMIT ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanSemit += constDataIterator.Get();
                volSemitReal++;
                if((tmp) < 0)
                    histogramSemit[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramSemit[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramSemit[tmp]++;
            }
            volSemit++;
        }
        if ( constErodeMaskIterator.Get() == REC_FEM ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanRect += constDataIterator.Get();
                volRectReal++;
                if((tmp) < 0)
                    histogramRect[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramRect[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramRect[tmp]++;

            }
            volRect++;
        }
        if ( constErodeMaskIterator.Get() == CRAN_SART ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanCran += constDataIterator.Get();
                volCranReal++;
                if((tmp) < 0)
                    histogramCran[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramCran[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramCran[tmp]++;
            }
            volCran++;
        }   	       	
        if ( constErodeMaskIterator.Get() == ADDUCTOR ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanAdductor += constDataIterator.Get();
                volAdductorReal++;
                if((tmp) < 0)
                    histogramAdductor[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramAdductor[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramAdductor[tmp]++;
            }
            volAdductor++;
        }   	       	
        if ( constErodeMaskIterator.Get() == GRACILIS ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanGracilis += constDataIterator.Get();
                volGracilisReal++;
                if((tmp) < 0)
                    histogramGracilis[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramGracilis[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramGracilis[tmp]++;
            }
            volGracilis++;
        }   	       	
        if ( constErodeMaskIterator.Get() == BIC_FEM ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanBic += constDataIterator.Get();
                volBicReal++;
                if((tmp) < 0)
                    histogramBic[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramBic[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramBic[tmp]++;
            }
            volBic++;
        }   	       	
        if ( constErodeMaskIterator.Get() == SEMI_MEM ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanSemimem += constDataIterator.Get();
                volSemimemReal++;
                if((tmp) < 0)
                    histogramSemimem[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramSemimem[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramSemimem[tmp]++;
            }
            volSemimem++;
        }   	       	
        if ( constErodeMaskIterator.Get() == VAS_LAT ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanVaslat += constDataIterator.Get();
                volVaslatReal++;
                if((tmp) < 0)
                    histogramVaslat[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramVaslat[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramVaslat[tmp]++;
            }
            volVaslat++;
        }   	       	
        if ( constErodeMaskIterator.Get() == VAS_MED ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanVasmed += constDataIterator.Get();
                volVasmedReal++;
                if((tmp) < 0)
                    histogramVasmed[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramVasmed[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramVasmed[tmp]++;
            }
            volVasmed++;
        }   	       	
        if ( constErodeMaskIterator.Get() == VAS_INT ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanVasint += constDataIterator.Get();
                volVasintReal++;
                if((tmp) < 0)
                    histogramVasint[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramVasint[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramVasint[tmp]++;
            }
            volVasint++;
        }   	       	
        if ( constErodeMaskIterator.Get() == CAUD_SART ) {
            if (constDataIterator.Get() > 0) {
                tmp = constDataIterator.Get();
	        meanCaud += constDataIterator.Get();
                volCaudReal++;
                if((tmp) < 0)
                    histogramCaud[0]++;
                else if (tmp >= HISTOGRAM_BIN_SIZE)
                    histogramCaud[HISTOGRAM_BIN_SIZE - 1]++;
                else
                    histogramCaud[tmp]++;
            }
            volCaud++;
        }   	       	
    }
    //calculate histogram entropy
    for (int i = 0; i < HISTOGRAM_BIN_SIZE; i++){
        if(histogramSemit[i] > 0) {
            entropySemit += -1 * histogramSemit[i] / volSemitReal * log (histogramSemit[i] / volSemitReal);
        }
        if(histogramRect[i] > 0) {
            entropyRect += -1 * histogramRect[i] / volRectReal * log (histogramRect[i] / volRectReal);
        }
        if(histogramCran[i] > 0) {
            entropyCran += -1 * histogramCran[i] / volCranReal * log (histogramCran[i] / volCranReal);
          //  std::cout << i << ":  " << histogramCran[i] << "  " << volCranReal << "  " << entropyCran << std::endl;
           // getchar();
        }
        if(histogramAdductor[i] > 0) {
            entropyAdductor += -1 * histogramAdductor[i] / volAdductorReal * log (histogramAdductor[i] / volAdductorReal);
        }
        if(histogramGracilis[i] > 0) {
            entropyGracilis += -1 * histogramGracilis[i] / volGracilisReal * log (histogramGracilis[i] / volGracilisReal);
        }
        if(histogramBic[i] > 0) {
            entropyBic += -1 * histogramBic[i] / volBicReal * log (histogramBic[i] / volBicReal);
        }
        if(histogramSemimem[i] > 0) {
            entropySemimem += -1 * histogramSemimem[i] / volSemimemReal * log (histogramSemimem[i] / volSemimemReal);
        }
        if(histogramVaslat[i] > 0) {
            entropyVaslat += -1 * histogramVaslat[i] / volVaslatReal * log (histogramVaslat[i] / volVaslatReal);
        }
        if(histogramVasmed[i] > 0) {
            entropyVasmed += -1 * histogramVasmed[i] / volVasmedReal * log (histogramVasmed[i] / volVasmedReal);
        }
        if(histogramVasint[i] > 0) {
            entropyVasint += -1 * histogramVasint[i] / volVasintReal * log (histogramVasint[i] / volVasintReal);
        }
        if(histogramCaud[i] > 0) {
            entropyCaud += -1 * histogramCaud[i] / volCaudReal * log (histogramCaud[i] / volCaudReal);
        }
    }
    std::cout << "the volum of the semit region is:  " << volSemit << " original size:  " << volSemit * voxelSize << std::endl;
    std::cout << "the volum of the rect region is:  " << volRect << " original size:  " << volRect * voxelSize << std::endl;
    std::cout << "the volum of the cran region is:  " << volCran << " original size:  " << volCran * voxelSize << std::endl;
    std::cout << "the volum of the Adductor is:  " << volAdductor << " original size:  " << volAdductor * voxelSize << std::endl;
    std::cout << "the volum of the Gracilis is:  " << volGracilis << " original size:  " << volGracilis * voxelSize << std::endl;
    std::cout << "the volum of the Biceps femoris is:  " << volBic << " original size:  " << volBic * voxelSize << std::endl;
    std::cout << "the volum of the Semimem is:  " << volSemimem << " original size:  " << volSemimem * voxelSize << std::endl;
    std::cout << "the volum of the Vaslat is:  " << volVaslat << " original size:  " << volVaslat * voxelSize << std::endl;
    std::cout << "the volum of the Vasmed is:  " << volVasmed << " original size:  " << volVasmed * voxelSize << std::endl;
    std::cout << "the volum of the Vasint is:  " << volVasint << " original size:  " << volVasint * voxelSize << std::endl;
    std::cout << "the volum of the Caud is:  " << volCaud << " original size:  " << volCaud * voxelSize << std::endl;
    quanlitySemit = meanSemit;
    quanlityRect = meanRect;
    quanlityCran = meanCran;
    quanlityAdductor = meanAdductor;
    quanlityGracilis = meanGracilis;
    quanlityBic = meanBic;
    quanlitySemimem = meanSemimem;
    quanlityVaslat = meanVaslat;
    quanlityVasmed = meanVasmed;
    quanlityVasint = meanVasint;
    quanlityCaud = meanCaud;
    if (volSemitReal)
        meanSemit /= volSemitReal ;
    if (volRectReal)
        meanRect /= volRectReal ;
    if (volCranReal)
        meanCran /= volCranReal ;
    if (volAdductorReal)
        meanAdductor /= volAdductorReal ;
    if (volGracilisReal)
        meanGracilis /= volGracilisReal ;
    if (volBicReal)
        meanBic /= volBicReal ;
    if (volSemimemReal)
        meanSemimem /= volSemimemReal ;
    if (volVaslatReal)
        meanVaslat /= volVaslatReal ;
    if (volVasmedReal)
        meanVasmed /= volVasmedReal ;
    if (volVasintReal)
        meanVasint /= volVasintReal ;
    if (volCaudReal)
        meanCaud /= volCaudReal ;
    // calculate the standard deviation of muscles
    for ( constErodeMaskIterator.GoToBegin(),constDataIterator.GoToBegin(); !constErodeMaskIterator.IsAtEnd(); ++constErodeMaskIterator, ++constDataIterator ) {
        if ( constErodeMaskIterator.Get() == SEMIT ) {
            if (constDataIterator.Get() > 0) {
	        sdSemit += (constDataIterator.Get() - meanSemit) * (constDataIterator.Get() - meanSemit);
            }
        }
        if ( constErodeMaskIterator.Get() == REC_FEM ) {
            if (constDataIterator.Get() > 0) {
	        sdRect += (constDataIterator.Get() - meanRect) * (constDataIterator.Get() - meanRect);
            }
        }
        if ( constErodeMaskIterator.Get() == CRAN_SART ) {
            if (constDataIterator.Get() > 0) {
	        sdCran += (constDataIterator.Get() - meanCran) * (constDataIterator.Get() - meanCran);
            }
        }
        if ( constErodeMaskIterator.Get() == ADDUCTOR ) {
            if (constDataIterator.Get() > 0) {
	        sdAdductor += (constDataIterator.Get() - meanAdductor) * (constDataIterator.Get() - meanAdductor);
            }
        }
        if ( constErodeMaskIterator.Get() == GRACILIS ) {
            if (constDataIterator.Get() > 0) {
	        sdGracilis += (constDataIterator.Get() - meanGracilis) * (constDataIterator.Get() - meanGracilis);
            }
        }
        if ( constErodeMaskIterator.Get() == BIC_FEM ) {
            if (constDataIterator.Get() > 0) {
	        sdBic += (constDataIterator.Get() - meanBic) * (constDataIterator.Get() - meanBic);
            }
        }
        if ( constErodeMaskIterator.Get() == SEMI_MEM ) {
            if (constDataIterator.Get() > 0) {
	        sdSemimem += (constDataIterator.Get() - meanSemimem) * (constDataIterator.Get() - meanSemimem);
            }
        }
        if ( constErodeMaskIterator.Get() == VAS_LAT ) {
            if (constDataIterator.Get() > 0) {
	        sdVaslat += (constDataIterator.Get() - meanVaslat) * (constDataIterator.Get() - meanVaslat);
            }
        }
        if ( constErodeMaskIterator.Get() == VAS_MED ) {
            if (constDataIterator.Get() > 0) {
	        sdVasmed += (constDataIterator.Get() - meanVasmed) * (constDataIterator.Get() - meanVasmed);
            }
        }
        if ( constErodeMaskIterator.Get() == VAS_INT ) {
            if (constDataIterator.Get() > 0) {
	        sdVasint += (constDataIterator.Get() - meanVasint) * (constDataIterator.Get() - meanVasint);
            }
        }
        if ( constErodeMaskIterator.Get() == CAUD_SART ) {
            if (constDataIterator.Get() > 0) {
	        sdCaud += (constDataIterator.Get() - meanCaud) * (constDataIterator.Get() - meanCaud);
            }
        }
    }
    sdSemit = sqrt(sdSemit / volSemitReal);
    sdRect = sqrt(sdRect / volRectReal);
    sdCran = sqrt(sdCran / volCranReal);
    sdAdductor = sqrt(sdAdductor / volAdductorReal);
    sdGracilis = sqrt(sdGracilis / volGracilisReal);
    sdBic = sqrt(sdBic / volBicReal);
    sdSemimem = sqrt(sdSemimem / volSemimemReal);
    sdVaslat = sqrt(sdVaslat / volVaslatReal);
    sdVasmed = sqrt(sdVasmed / volVasmedReal);
    sdVasint = sqrt(sdVasint / volVasintReal);
    sdCaud = sqrt(sdCaud / volCaudReal );
    // calculate the standard deviation of muscles
    for ( constErodeMaskIterator.GoToBegin(),constDataIterator.GoToBegin(); !constErodeMaskIterator.IsAtEnd(); ++constErodeMaskIterator, ++constDataIterator ) {
        if ( constErodeMaskIterator.Get() == SEMIT ) {
            if (constDataIterator.Get() > 0) {
	        skewSemit += pow(constDataIterator.Get() - meanSemit, 3);
	        kortSemit += pow(constDataIterator.Get() - meanSemit, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == REC_FEM ) {
            if (constDataIterator.Get() > 0) {
	        skewRect += pow(constDataIterator.Get() - meanRect, 3);
	        kortRect += pow(constDataIterator.Get() - meanRect, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == CRAN_SART ) {
            if (constDataIterator.Get() > 0) {
	        skewCran += pow(constDataIterator.Get() - meanCran, 3);
	        kortCran += pow(constDataIterator.Get() - meanCran, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == ADDUCTOR ) {
            if (constDataIterator.Get() > 0) {
	        skewAdductor += pow(constDataIterator.Get() - meanAdductor, 3);
	        kortAdductor += pow(constDataIterator.Get() - meanAdductor, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == GRACILIS ) {
            if (constDataIterator.Get() > 0) {
	        skewGracilis += pow(constDataIterator.Get() - meanGracilis, 3);
	        kortGracilis += pow(constDataIterator.Get() - meanGracilis, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == BIC_FEM ) {
            if (constDataIterator.Get() > 0) {
	        skewBic += pow(constDataIterator.Get() - meanBic, 3);
	        kortBic += pow(constDataIterator.Get() - meanBic, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == SEMI_MEM ) {
            if (constDataIterator.Get() > 0) {
	        skewSemimem += pow(constDataIterator.Get() - meanSemimem, 3);
	        kortSemimem += pow(constDataIterator.Get() - meanSemimem, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == VAS_LAT ) {
            if (constDataIterator.Get() > 0) {
	        skewVaslat += pow(constDataIterator.Get() - meanVaslat, 3);
	        kortVaslat += pow(constDataIterator.Get() - meanVaslat, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == VAS_MED ) {
            if (constDataIterator.Get() > 0) {
	        skewVasmed += pow(constDataIterator.Get() - meanVasmed, 3);
	        kortVasmed += pow(constDataIterator.Get() - meanVasmed, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == VAS_INT ) {
            if (constDataIterator.Get() > 0) {
	        skewVasint += pow(constDataIterator.Get() - meanVasint, 3);
	        kortVasint += pow(constDataIterator.Get() - meanVasint, 4);
            }
        }
        if ( constErodeMaskIterator.Get() == CAUD_SART ) {
            if (constDataIterator.Get() > 0) {
	        skewCaud += pow(constDataIterator.Get() - meanCaud, 3);
	        kortCaud += pow(constDataIterator.Get() - meanCaud, 4);
            }
        }
    }
    if(sdSemit)
        skewSemit = skewSemit / (pow(sdSemit, 3) * (volSemitReal - 1));
    else
        skewSemit = 0;
    if(sdSemit)
        kortSemit = kortSemit / (pow(sdSemit, 4) * (volSemitReal - 1)) - 3;
    else
        kortSemit = 0;
    if(sdRect)
        skewRect = skewRect / (pow(sdRect, 3) * (volRectReal - 1));
    else
        skewRect = 0;
    if(sdRect)
        kortRect = kortRect / (pow(sdRect, 4) * (volRectReal - 1)) - 3;
    else
        kortRect = 0;

    if(sdCran)
        skewCran = skewCran / (pow(sdCran, 3) * (volCranReal - 1));
    else
        skewCran = 0;
    if(sdCran)
        kortCran = kortCran / (pow(sdCran, 4) * (volCranReal - 1)) - 3;
    else
        kortCran = 0;

    if(sdAdductor)
        skewAdductor = skewAdductor / (pow(sdAdductor, 3) * (volAdductorReal - 1));
    else
        skewAdductor = 0;
    if(sdAdductor)
        kortAdductor = kortAdductor / (pow(sdAdductor, 4) * (volAdductorReal - 1)) - 3;
    else
        kortAdductor = 0;

    if(sdGracilis)
        skewGracilis = skewGracilis / (pow(sdGracilis, 3) * (volGracilisReal - 1));
    else
        skewGracilis = 0;
    if(sdGracilis)
        kortGracilis = kortGracilis / (pow(sdGracilis, 4) * (volGracilisReal - 1)) - 3;
    else
        kortGracilis = 0;

    if(sdBic)
        skewBic = skewBic / (pow(sdBic, 3) * (volBicReal - 1));
    else
        skewBic = 0;
    if(sdBic)
        kortBic = kortBic / (pow(sdBic, 4) * (volBicReal - 1)) - 3;
    else
        kortBic = 0;

    if(sdSemimem)
        skewSemimem = skewSemimem / (pow(sdSemimem, 3) * (volSemimemReal - 1));
    else
        skewSemimem = 0;
    if(sdSemimem)
        kortSemimem = kortSemimem / (pow(sdSemimem, 4) * (volSemimemReal - 1)) - 3;
    else
        kortSemimem = 0;

    if(sdVaslat)
        skewVaslat = skewVaslat / (pow(sdVaslat, 3) * (volVaslatReal - 1));
    else
        skewVaslat = 0;
    if(sdVaslat)
        kortVaslat = kortVaslat / (pow(sdVaslat, 4) * (volVaslatReal - 1)) - 3;
    else
        kortVaslat = 0;

    if(sdVasmed)
        skewVasmed = skewVasmed / (pow(sdVasmed, 3) * (volVasmedReal - 1));
    else
        skewVasmed = 0;
    if(sdVasmed)
        kortVasmed = kortVasmed / (pow(sdVasmed, 4) * (volVasmedReal - 1)) - 3;
    else
        kortVasmed = 0;

    if(sdVasint)
        skewVasint = skewVasint / (pow(sdVasint, 3) * (volVasintReal - 1));
    else
        skewVasint = 0;
    if(sdVasint)
        kortVasint = kortVasint / (pow(sdVasint, 4) * (volVasintReal - 1)) - 3;
    else
        kortVasint = 0;

    if(sdCaud)
        skewCaud = skewCaud / (pow(sdCaud, 3) * (volCaudReal - 1));
    else
        skewCaud = 0;
    if(sdCaud)
        kortCaud = kortCaud / (pow(sdCaud, 4) * (volCaudReal - 1)) - 3;
    else
        kortCaud = 0;
/*
    if(volOutputMark){
        efile << "                   CS     RF     BF     GR     SM     AD     ST    VL    VI    VM    CA" << "\n";
        efile << "muscle_volume_in_voxel       " ;
        efile << volCran << "  " << volRect << "  " << volBic << "  " << volGracilis << "  " << volSemimem << "  " << volAdductor << "  " << volSemit << "  " << volVaslat << "  " << volVasint << "  " << volVasmed << "  " << volCaud << "\n\n";
        efile << "muscle_volume       " ;
        efile << volCranNonErod * voxelSize << "  " << volRectNonErod * voxelSize << "  " << volBicNonErod * voxelSize << "  " << volGracilisNonErod * voxelSize << "  " << volSemimemNonErod * voxelSize << "  " << volAdductorNonErod * voxelSize << "  " << volSemitNonErod * voxelSize << "  " << volVaslatNonErod * voxelSize << "  " << volVasintNonErod * voxelSize << "  " << volVasmedNonErod * voxelSize << "  " << volCaudNonErod * voxelSize << "\n";
        efile << "t2_weight          ";
    }
    efile << meanCran << "  " << meanRect << "  " << meanBic << "  " << meanGracilis << "  " << meanSemimem << "  " << meanAdductor << "  " << meanSemit << "  " << meanVaslat << "  " << meanVasint << "  " << meanVasmed << "  " << meanCaud << "\n";

#ifdef HISTOGRAM_TEXTURE
    efile << "SD        ";
    efile << sdCran << "  " << sdRect << "  " << sdBic << "  " << sdGracilis << "  " << sdAdductor << "  " << sdSemit << "  " << sdVaslat << "  " << sdVasint << "  " << sdVasmed << "  " << sdSemimem << "  " << sdCaud << "\n";
    efile << "skewness        ";
    efile << skewCran << "  " << skewRect << "  " << skewBic << "  " << skewGracilis << "  " << skewSemimem << "  " << skewAdductor << "  " << skewSemit << "  " << skewVaslat << "  " << skewVasint << "  " << skewVasmed << "  " << skewCaud << "\n";
    efile << "kortness                                                                                                                                                                                                             ";
    efile << kortCran << "  " << kortRect << "  " << kortBic << "  " << kortGracilis << "  " << kortSemimem << "  " << kortAdductor << "  " << kortSemit << "  " << kortVaslat << "  " << kortVasint << "  " << kortVasmed << "  " << kortCaud << "\n";
    efile << "entropy        ";
    efile << entropyCran << "  " << entropyRect << "  " << entropyBic << "  " << entropyGracilis << "  " << entropySemimem << "  " << entropyAdductor << "  " << entropySemit << "  " << entropyVaslat << "  " << entropyVasint << "  " << entropyVasmed << "  " << entropyCaud << "\n";
#endif
*/
    if(volOutputMark){
        efile << "                   CS     RF     BF     GR     SM     AD     ST    VL    VI    VM    CA" << "\n";
        efile << "muscle_volume_in_voxel       " ;
        efile << volCran << "  " << volRect << "  " << volBic << "  " << volGracilis << "  " << volSemimem << "  " << volAdductor << "  " << volSemit << "  " << volVaslat << "  " << volVasint << "  " << volVasmed << "  " << volCaud << "\n\n";
        efile << "muscle_volume       " ;
        efile << volCranNonErod * voxelSize << "  " << volRectNonErod * voxelSize << "  " << volBicNonErod * voxelSize << "  " << volGracilisNonErod * voxelSize << "  " << volSemimemNonErod * voxelSize << "  " << volAdductorNonErod * voxelSize << "  " << volSemitNonErod * voxelSize << "  " << volVaslatNonErod * voxelSize << "  " << volVasintNonErod * voxelSize << "  " << volVasmedNonErod * voxelSize << "  " << volCaudNonErod * voxelSize << "\n";
        efile << "t2_weight          ";
    }
    efile << meanCran << "  " << meanRect << "  " << meanBic << "  " << meanGracilis << "  " << meanSemimem << "  " << meanAdductor << "  " << meanSemit << "  " << meanVaslat << "  " << meanVasint << "  " << meanVasmed << "  " << meanCaud << "\n";

#ifdef HISTOGRAM_TEXTURE
    efile << "SD        ";
    efile << sdCran << "  " << sdRect << "  " << sdBic << "  " << sdGracilis << "  " << sdAdductor << "  " << sdSemit << "  " << sdVaslat << "  " << sdVasint << "  " << sdVasmed << "  " << sdSemimem << "  " << sdCaud << "\n";
    efile << "skewness        ";
    efile << skewCran << "  " << skewRect << "  " << skewBic << "  " << skewGracilis << "  " << skewSemimem << "  " << skewAdductor << "  " << skewSemit << "  " << skewVaslat << "  " << skewVasint << "  " << skewVasmed << "  " << skewCaud << "\n";
    efile << "kortness                                                                                                                                                                                                             ";
    efile << kortCran << "  " << kortRect << "  " << kortBic << "  " << kortGracilis << "  " << kortSemimem << "  " << kortAdductor << "  " << kortSemit << "  " << kortVaslat << "  " << kortVasint << "  " << kortVasmed << "  " << kortCaud << "\n";
    efile << "entropy        ";
    efile << entropyCran << "  " << entropyRect << "  " << entropyBic << "  " << entropyGracilis << "  " << entropySemimem << "  " << entropyAdductor << "  " << entropySemit << "  " << entropyVaslat << "  " << entropyVasint << "  " << entropyVasmed << "  " << entropyCaud << "\n";
#endif
    efile.close();

    std::cout << "the mean and sd of Semitendinosus is:  " << meanSemit << "  " << sdSemit << std::endl;
    std::cout << "the mean and sd of Rectus femoris is:  " << meanRect << "  " << sdRect << std::endl;
    std::cout << "the mean and sd of Cranial Sartorius is:  " << meanCran << "  " << sdCran << std::endl;
    std::cout << "the mean and sd of Adductor is:  " << meanAdductor << "  " << sdAdductor << std::endl;
    std::cout << "the mean and sd of Gracilis is:  " << meanGracilis << "  " << sdGracilis << std::endl;
    std::cout << "the mean and sd of Biceps femoris is:  " << meanBic << "  " << sdBic << std::endl;

    std::cout << "the mean and sd of Semimem femoris is:  " << meanSemimem << "  " << sdSemimem << std::endl;
    std::cout << "the mean and sd of Vaslat femoris is:  " << meanVaslat << "  " << sdVaslat << std::endl;
    std::cout << "the mean and sd of Vasmed femoris is:  " << meanVasmed << "  " << sdVasmed << std::endl;
    std::cout << "the mean and sd of Vasint femoris is:  " << meanVasint << "  " << sdVasint << std::endl;
    std::cout << "the mean and sd of Caud femoris is:  " << meanCaud << "  " << sdCaud << std::endl;

}

#endif
