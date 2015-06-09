/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef __vtkTestMain_h
#define __vtkTestMain_h

// This file is used to create TestDriver executables
// These executables are able to register a function pointer to a string name
// in a lookup table.   By including this file, it creates a main function
// that calls RegisterTests() then looks up the function pointer for the test
// specified on the command line.
#include "itkWin32Header.h"
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "itksys/SystemTools.hxx"
#include <limits>
#include <sstream>

// VTK testing
#include <vtkVersion.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSphereSource.h>
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataReader.h>
#define ITK_TEST_DIMENSION_MAX 6

typedef int ( *MainFuncPointer )( int , char *[] ) ;
std::map<std::string, MainFuncPointer> StringToTestFunctionMap ;

#define REGISTER_TEST( test )       \
    extern int test( int , char *[] ) ; \
    StringToTestFunctionMap[ #test ] = test

int RegressionTestPolyData( const char *testImageFilename ,
                           const char *baselineImageFilename
                           ) ;

std::map<std::string, int> RegressionTestBaselines( char * ) ;

void RegisterTests() ;

void PrintAvailableTests()
{
    std::cout << "Available tests:\n" ;
    std::map<std::string, MainFuncPointer>::iterator j = StringToTestFunctionMap.begin() ;
    int                                              i = 0 ;
    while( j != StringToTestFunctionMap.end() )
    {
        std::cout << i << ". " << j->first << "\n" ;
        ++i ;
        ++j ;
    }
}

int main( int ac , char *av[] )
{
    typedef std::pair< char * , char * > ComparePairType ;
    std::vector< ComparePairType > compareList ;
    RegisterTests() ;
    std::string testToRun ;
    if( ac < 2 )
    {
        PrintAvailableTests() ;
        std::cout << "To run a test, enter the test number: " ;
        int testNum = 0 ;
        std::cin >> testNum ;
        std::map<std::string, MainFuncPointer>::iterator j = StringToTestFunctionMap.begin() ;
        int                                              i = 0 ;
        while( j != StringToTestFunctionMap.end() && i < testNum )
        {
            ++i ;
            ++j ;
        }

        if( j == StringToTestFunctionMap.end() )
        {
            std::cerr << testNum << " is an invalid test number\n" ;
            return -1 ;
        }
        testToRun = j->first ;
    }
    else
    {
        while( ac > 0 && testToRun.empty() )
        {
            if( ac > 3 && strcmp(av[1], "--compare") == 0 )
            {
                compareList.push_back( ComparePairType( av[ 2 ] , av[ 3 ] ) ) ;
                av += 3 ;
                ac -= 3 ;
            }
            else
            {
                testToRun = av[ 1 ] ;
            }
        }
    }
    std::map< std::string , MainFuncPointer >::iterator j = StringToTestFunctionMap.find( testToRun ) ;
    if( j != StringToTestFunctionMap.end() )
    {
        MainFuncPointer f = j->second ;
        int             result ;
        try
        {
            // Invoke the test's "main" function.
            result = ( *f )( ac - 1, av + 1 ) ;
            // Make a list of possible baselines
            for( int i = 0 ; i < static_cast< int >( compareList.size() ) ; i++ )
            {
                char *                               baselineFilename = compareList[ i ].first ;
                char *                               testFilename = compareList[ i ].second ;
                std::map<std::string, int>           baselines = RegressionTestBaselines( baselineFilename ) ;
                std::map<std::string, int>::iterator baseline = baselines.begin() ;
                std::string                          bestBaseline ;
                int                                  bestBaselineStatus = std::numeric_limits< int >::max() ;
                while( baseline != baselines.end() )
                {
                    baseline->second = RegressionTestPolyData( testFilename ,
                                                              ( baseline->first ).c_str()
                                                              );
                    if( baseline->second < bestBaselineStatus )
                    {
                        bestBaseline = baseline->first ;
                        bestBaselineStatus = baseline->second ;
                    }
                    if( baseline->second == 0 )
                    {
                        break ;
                    }
                    ++baseline ;
                }

                // if the best we can do still has errors, generate the error images
                if( bestBaselineStatus )
                {
                    std::cout << "None of the baseline correspond to the output" << std::endl ;
                }

                // output the matching baseline
                std::cout << "<DartMeasurement name=\"BaselineImageName\" type=\"text/string\">" ;
                std::cout << itksys::SystemTools::GetFilenameName( bestBaseline ) ;
                std::cout << "</DartMeasurement>" << std::endl ;
                result += bestBaselineStatus ;
            }
        }
        catch( ... )
        {
            std::cerr << "VTK test driver caught an unknown exception!!!\n" ;
            result = -1 ;
        }
        return result ;
    }
    PrintAvailableTests() ;
    std::cerr << "Failed: " << testToRun << ": No test registered with name " << testToRun << "\n" ;
    return -1 ;
}

// Regression Testing Code

int ReadVTKVTPFile( const char* fileName , std::string extension , vtkSmartPointer< vtkPolyData > &polyData )
{
    if( extension == ".vtk" )
    {
        vtkSmartPointer< vtkPolyDataReader > reader = vtkSmartPointer< vtkPolyDataReader >::New() ;
        reader->SetFileName( fileName ) ;
        reader->Update() ;
        polyData = reader->GetOutput() ;
        return reader->GetErrorCode() ;
    }
    else if( extension == ".vtp" )
    {
        vtkSmartPointer< vtkXMLPolyDataReader > reader = vtkSmartPointer< vtkXMLPolyDataReader >::New() ;
        reader->SetFileName( fileName ) ;
        reader->Update() ;
        polyData = reader->GetOutput() ;
        return reader->GetErrorCode() ;
    }
    else
    {
        return 1 ;
    }
}

int RegressionTestPolyData( const char *testImageFilename ,
                           const char *baselineImageFilename
                           )
{
    vtkSmartPointer<vtkPolyData> vtkpolydata1 =
            vtkSmartPointer<vtkPolyData>::New() ;
    vtkSmartPointer<vtkPolyData> vtkpolydata2 =
            vtkSmartPointer<vtkPolyData>::New() ;
    std::ifstream infile( testImageFilename ) ;
    if( !infile )
    {
        std::cerr << "Unable to find: " << testImageFilename << std::endl ;
        return -1 ;
    }
    std::ifstream infileBaseline( baselineImageFilename ) ;
    if( !infileBaseline )
    {
        std::cerr << "Unable to find: " << baselineImageFilename << std::endl ;
        return -1 ;
    }
    if( ReadVTKVTPFile( testImageFilename , itksys::SystemTools::GetFilenameExtension( testImageFilename ) , vtkpolydata1 ) != 0 )
    {
        return -1 ;
    }
    if( ReadVTKVTPFile( baselineImageFilename , itksys::SystemTools::GetFilenameExtension( baselineImageFilename) , vtkpolydata2 ) !=0 )
    {
            return -1 ;
    }
    if( vtkpolydata2->GetNumberOfCells() != vtkpolydata1->GetNumberOfCells() )
    {
        return 1 ;
    }
    if( vtkpolydata2->GetNumberOfStrips() != vtkpolydata1->GetNumberOfStrips() )
    {
        return 2 ;
    }
    if( vtkpolydata2->GetNumberOfPolys() != vtkpolydata1->GetNumberOfPolys() )
    {
        return 3 ;
    }
    if( vtkpolydata2->GetNumberOfLines() != vtkpolydata1->GetNumberOfLines() )
    {
        return 4 ;
    }
    if( vtkpolydata2->GetNumberOfVerts() != vtkpolydata1->GetNumberOfVerts() )
    {
        return 5 ;
    }
    if( vtkpolydata2->GetNumberOfPieces() != vtkpolydata1->GetNumberOfPieces() )
    {
        return 6 ;
    }
    if( vtkpolydata2->GetMaxCellSize() != vtkpolydata1->GetMaxCellSize() )
    {
        return 7 ;
    }
    if( vtkpolydata2->GetCellData()->GetNumberOfArrays() != vtkpolydata1->GetCellData()->GetNumberOfArrays() )
    {
        std::cout << "8" << std::endl ;
        return 8;
    }
    for( int i = 0 ; i < vtkpolydata2->GetCellData()->GetNumberOfArrays() ; i++ )
    {
        if( vtkpolydata2->GetCellData()->GetArray( i )->GetNumberOfTuples() != vtkpolydata1->GetCellData()->GetArray( i )->GetNumberOfTuples() )
        {
            std::cout << "9" << std::endl ;
            return 9 ;
        }
        for( int j = 0 ; j < vtkpolydata2->GetCellData()->GetArray( i )->GetNumberOfTuples() ; j++ )
        {
            if( vtkpolydata2->GetCellData()->GetArray( i )->GetNumberOfComponents() != vtkpolydata1->GetCellData()->GetArray( i )->GetNumberOfComponents() )
            {
                std::cout << "10"<<std::endl ;
                return 10 ;
            }
            for( int k = 0 ; k < vtkpolydata2->GetCellData()->GetArray( i )->GetNumberOfComponents() ; k++ )
            {
                if( vtkpolydata2->GetCellData()->GetArray( i )->GetComponent( j , k ) != vtkpolydata1->GetCellData()->GetArray( i )->GetComponent( j , k ) )
                {
                    std::cout << vtkpolydata2->GetCellData()->GetArray( i )->GetComponent( j , k ) << std::endl<<vtkpolydata1->GetCellData()->GetArray( i )->GetComponent( j , k ) << std::endl ;
                    std::cout << "11" << std::endl ;
                    return 11 ;
                }
            }
        }
    }
        return 0;
}

//
// Generate all of the possible baselines
// The possible baselines are generated fromn the baselineFilename using the
// following algorithm:
// 1) strip the suffix
// 2) append a digit .x
// 3) append the original suffix.
// It the file exists, increment x and continue
//
std::map< std::string , int > RegressionTestBaselines( char *baselineFilename )
{
    std::map< std::string , int > baselines ;
    baselines[ std::string( baselineFilename ) ] = 0 ;

    std::string originalBaseline( baselineFilename ) ;

    int                    x = 0 ;
    std::string::size_type suffixPos = originalBaseline.rfind(".") ;
    std::string            suffix ;
    if( suffixPos != std::string::npos )
    {
        suffix = originalBaseline.substr( suffixPos, originalBaseline.length() ) ;
        originalBaseline.erase( suffixPos, originalBaseline.length() ) ;
    }
    while( ++x )
    {
        std::ostringstream filename ;
        filename << originalBaseline << "." << x << suffix ;
        std::ifstream filestream( filename.str().c_str() ) ;
        if( !filestream )
        {
            break ;
        }
        baselines[ filename.str() ] = 0 ;
        filestream.close() ;
    }

    return baselines ;
}

#endif
