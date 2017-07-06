#ifndef FIBERFILEIO_H
#define FIBERFILEIO_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <itkMacro.h> //For itkEception

#define BLUE "\033[0;34m"
#define BLUE_BOLD "\033[1;34m"
#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define GREEN_BOLD "\033[1;32m"
#define YELLOW "\033[0;33m"
#define CYAN "\033[0;36m"
#define CYAN_BOLD "\033[1;36m"
#define NC "\033[0m"
#define NC_BOLD "\033[1m"

/**
 * Read a VTK File
 * @param   - Location of the file to read
 * @return  - Data of the file
 */
vtkSmartPointer<vtkPolyData> readVTKFile (std::string filename);

/**
 * Write in a VTK Format (.vtp or .vtk) a file at the location thrown in parameter
 * @param filename - Location of the file to write
 * @param output   - Data of the file to write
 */
void writeVTKFile (std::string filename, vtkSmartPointer<vtkPolyData> output);

#ifndef ITK_MANUAL_INSTANTIATION
#include "fiberfileIO.hxx"
#endif
 
 
#endif