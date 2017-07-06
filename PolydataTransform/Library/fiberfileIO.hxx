#ifndef __FiberFileIO_hxx
#define __FiberFileIO_hxx

#include "fiberfileIO.h"
#include "vtkVersion.h"
#define NB_LINES_MAX 250
#define NB_WORDS_MAX 250
vtkSmartPointer<vtkPolyData> readVTKFile (std::string filename)
{           // VTK
	if (filename.rfind(".vtk") != std::string::npos)
	{
        std::cout<<"---Reading VTK file "<<GREEN<<filename.c_str()<<NC<<std::endl;
        vtkSmartPointer<vtkPolyDataReader> fiberReader = vtkPolyDataReader::New();
        fiberReader->SetFileName(filename.c_str());
        if(fiberReader->OpenVTKFile())
        {
            fiberReader->Update();
            return fiberReader->GetOutput();
        }
        else
        {
            throw itk::ExceptionObject("File Non Valid");
        }

    }
	        // XML
    else if (filename.rfind(".vtp") != std::string::npos)
    {
        std::cout<<"---Reading VTP file "<<GREEN<<filename.c_str()<<NC<<std::endl;
        vtkSmartPointer<vtkXMLPolyDataReader> fiberReader = vtkXMLPolyDataReader::New();
        fiberReader->SetFileName(filename.c_str());
        if(fiberReader->CanReadFile(filename.c_str()))
        {
            fiberReader->Update();
            return fiberReader->GetOutput();
        }
        else
        {
            throw itk::ExceptionObject("File Non Valid");
        }
    }
    else
    {
        throw itk::ExceptionObject("Unknown file format for fibers");
    }

}

void writeVTKFile (std::string filename, vtkSmartPointer<vtkPolyData> output)
{	
	if (filename.rfind(".vtk") != std::string::npos)
	{
        std::cout<<"---Writing VTK file "<<GREEN<<filename.c_str()<<NC<<std::endl;
        vtkSmartPointer<vtkPolyDataWriter> fiberWriter = vtkPolyDataWriter::New();
        fiberWriter->SetFileName(filename.c_str());

		#if (VTK_MAJOR_VERSION < 6)
            fiberWriter->SetInput(output);
		#else
            fiberWriter->SetInputData(output);
		#endif
            fiberWriter->Update();
    }
	        // XML
    else if (filename.rfind(".vtp") != std::string::npos)
    {
        std::cout<<"---Writing VTP file "<<GREEN<<filename.c_str()<<NC<<std::endl;
    	vtkSmartPointer<vtkXMLPolyDataWriter> fiberWriter = vtkXMLPolyDataWriter::New();
    	fiberWriter->SetFileName(filename.c_str());
		#if (VTK_MAJOR_VERSION < 6)
        	fiberWriter->SetInput(output);
		#else
        	fiberWriter->SetInputData(output);
		#endif
        	fiberWriter->Update();
    }
    else
    {
        throw itk::ExceptionObject("Unknown file format for fibers");
    }

	
}

vtkSmartPointer<vtkPolyData> readFCSVFile (std::string filename)
{
    if (filename.rfind(".fcsv") != std::string::npos)
    {
        std::cout<<"---Reading FCSV file "<<GREEN<<filename.c_str()<<NC<<std::endl;
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        std::fstream fcsvfile(filename.c_str());
        std::string line;
        std::string mot;
        std::string words[NB_LINES_MAX][NB_WORDS_MAX]; 
        int i,j;
        // vtkSmartPointer<vtkPoints> landmarksPtToAdd = vtkSmartPointer<vtkPoints>::New();
        vtkPoints    * points = vtkPoints::New();
        if(fcsvfile)
        {
            getline(fcsvfile, line);
            fcsvfile>>mot;
            while(mot=="#")
            {
                if(getline(fcsvfile, line))
                    fcsvfile>>mot;
                else
                    mot="#";
            }

            i=0;
            do
            {
                std::size_t pos_end;// = mot.find(",,");
                std::size_t pos1;
                j=0;
                do
                {
                    if(j>NB_WORDS_MAX)
                        throw itk::ExceptionObject("Too many words per line in the FCSV File Max=250");
                    if(i>NB_LINES_MAX)
                        throw itk::ExceptionObject("Too many lines in the FCSV File Max = 250");
                    std::size_t pos0 = 0;
                    pos1 = mot.find(',');
                    pos_end = mot.find(",,");
                    words[i][j] = mot.substr(pos0, pos1-pos0);
                    mot = mot.substr(pos1+1);
                    j++;
                }
                while(pos1+1<pos_end);
                i++;
            }
            while(fcsvfile>>mot);
            int NbPoints = i;
            for (int i = 0; i < NbPoints; ++i)
            {
                double pt[3];
                pt[0] = atof(words[i][1].c_str());
                pt[1] = atof(words[i][2].c_str());
                pt[2] = atof(words[i][3].c_str());
                points->InsertPoint(i,pt[0],pt[1],pt[2]);
            }
            polydata->SetPoints(points);
            return polydata;
        }
        else
        {
            std::cout<<"Error !";
        }
    }
    else
    {
        throw itk::ExceptionObject("Unknown file format for fibers");
    }

}

void writeFCSVFile (std::string filename, vtkSmartPointer<vtkPolyData> output)
{
    if (filename.rfind(".fcsv") != std::string::npos)
    {
        std::cout<<"---Writting FCSV file "<<GREEN<<filename.c_str()<<NC<<std::endl;
        std::ofstream fcsvfile;
        fcsvfile.open(filename.c_str());
        vtkPoints* points = output->GetPoints();
        int nbPoints = output->GetNumberOfPoints();
        fcsvfile << "# Markups fiducial file version = 4.5\n";
        fcsvfile << "# CoordinateSystem = 0\n";
        fcsvfile << "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n";

        for(int i=0; i<nbPoints; i++)
        {
            fcsvfile <<"Landmark_"<<i<<","<<points->GetPoint(i)[0]<<","<<points->GetPoint(i)[1]<<","<<points->GetPoint(i)[2];
            fcsvfile <<",0,0,0,1,1,1,0,F-"<<i+1<<",,\n";
        }
        fcsvfile.close();
    }
    else
    {
        throw itk::ExceptionObject("Unknown file format for fibers");
    }

}
#endif