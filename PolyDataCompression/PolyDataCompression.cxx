#include "PolyDataCompressionCLP.h"
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <itksys/SystemTools.hxx>

std::string ChangeEndOfFileName ( std::string fileName, std::string change , std::string extension )
{
    if( itksys::SystemTools::GetFilenamePath( fileName ).empty() == true )
    {
        return itksys::SystemTools::GetFilenameWithoutLastExtension( fileName ) + change + "." + extension ;
    }
    return  itksys::SystemTools::GetFilenamePath( fileName )
            + "/" + itksys::SystemTools::GetFilenameWithoutLastExtension( fileName ) + change +"." + extension ;
}

int ReadFile( const char* fileName , std::string extension , vtkSmartPointer< vtkPolyData > &polyData )
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

int WriteFile( std::string encoding , std::string extension , const char* outputFileName ,
                int compressionLevel , vtkSmartPointer< vtkPolyData > readerPolyData )
{
    if( extension == "vtp" )
    {
        vtkSmartPointer< vtkXMLPolyDataWriter > writer = vtkSmartPointer< vtkXMLPolyDataWriter >::New() ;
        writer->SetFileName( outputFileName ) ;
        writer->SetInputData( readerPolyData ) ;
        vtkZLibDataCompressor *compressor = dynamic_cast< vtkZLibDataCompressor* > ( writer->GetCompressor() ) ;
        if( compressor )
        {
            compressor->SetCompressionLevel( compressionLevel ) ;
        }
        if( encoding == "binary" )
        {
            writer->SetDataModeToBinary() ;
        }
        else if( encoding == "appended" )
        {
            writer->SetDataModeToAppended() ;
        }
        else if( encoding == "ascii" )
        {
            writer->SetDataModeToAscii() ;
        }
        else // should not arrive here. tested in main before already and program should exit at that time.
        {
            return -1 ;
        }
        writer->Update() ;
        return writer->GetErrorCode() ;
    }
    else
    {
        vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New() ;
        writer->SetFileName( outputFileName ) ;
        writer->SetInputData( readerPolyData ) ;
        if( encoding =="binary" )
        {
            writer->SetFileTypeToBinary();
        }
        else if( encoding == "ascii" )
        {
            writer->SetFileTypeToASCII() ;
        }
        else
        {
            return -1 ;
        }

        writer->Update() ;
        return writer->GetErrorCode() ;
    }
}

int main( int argc, char *argv[] )
{
    std::vector< std::string > fileNameList ;
    bool flagFound = false ;
    for( int i = 0 ; i < argc ; i++ )
    {
        if( flagFound )
        {
            if( argv[ i ][ 0 ] == '-' )
            {
                std::cout << "Please write all your flags before the -f flag. -f flag should be the last flag. "<< std::endl ;
                return EXIT_FAILURE ;
            }
            else
            {
              fileNameList.push_back( argv[ i ] ) ;
            }
        }
        if( !strcmp( argv[ i ] , "-f" ) || !strcmp( argv[ i ] , "--fileNameList" ) )
        {
            flagFound = true ;
        }
    }
    argc = argc - ( fileNameList.empty() ? 0 : fileNameList.size() - 1 ) ;
    PARSE_ARGS ;
    if( !flagFound )
    {
       std::cerr << "Missing fileNameList flag, please give at least one file." << std::endl ;
       return EXIT_FAILURE ;
    }
    std::string append = "-compressed" ;
    for( size_t i = 0 ; i < fileNameList.size() ; i++ )
    {
        std::string testFileName = itksys::SystemTools::GetFilenameName( ( fileNameList[ i ] ) ) ;
        if(  testFileName.find( append + ".vtp" ) != std::string::npos ) //looks if "-compressed" is already at the end of the filename -> if yes, the file will be skipped
        {
            std::cout << fileNameList[ i ] << " already compressed." << std::endl ;
            continue ;
        }
        else
        {
            vtkSmartPointer< vtkPolyData > readerPolyData = vtkSmartPointer< vtkPolyData >::New() ;
            std::ifstream infile( fileNameList[ i ].c_str() ) ;
            if( !infile )
            {
                std::cerr << "Unable to find input file: " << fileNameList[ i ] << std::endl ;
                return EXIT_FAILURE ;
            }
            if( ReadFile( fileNameList[ i ].c_str() , itksys::SystemTools::GetFilenameLastExtension( fileNameList[ i ] ) , readerPolyData ) != 0 )
            {
                std::cerr << "Unable to open input file: " << fileNameList[ i ] << std::endl ;
                return EXIT_FAILURE ;
            }
            if( !readerPolyData.GetPointer() )
            {
                return EXIT_FAILURE ;
            }
            std::cout << "Compressing: " << fileNameList[ i ] << std::endl ;
            if( overwrite == 0 )
            {
                std::string outputFileName = ChangeEndOfFileName( fileNameList[ i ] , append , extension ) ;
                if( WriteFile( encoding , extension , outputFileName.c_str() , compressionLevel , readerPolyData ) )
                {
                    std::cerr << "Unable to write output file." << std::endl ;
                    return EXIT_FAILURE ;
                }
            }
            else
            {
                std::string outputFileName = ChangeEndOfFileName( fileNameList[ i ] , append , extension ) ;
                if( WriteFile( encoding , extension , outputFileName.c_str() , compressionLevel , readerPolyData ) )
                {
                    std::cerr << "Unable to write temporary compressed file." << std::endl ;
                    return EXIT_FAILURE ;
                }
                if( WriteFile( encoding , extension , fileNameList[ i ].c_str() , compressionLevel , readerPolyData ) )
                {
                    std::cerr << "Unable to write output file." << std::endl ;
                    return EXIT_FAILURE ;
                }
            }
        }
    }
    return EXIT_SUCCESS ;
}
