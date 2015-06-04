#include "IntensityRescalerCommandLine.h"

#ifdef QT_GUI
#define QT_ALTERNATE_QTSMANIP
#include <qtextstream.h>
#include <qfile.h>
#endif

#include <fstream>
#include <iostream>
#include "ImageIntensityNormalizer.h"

IntensityRescalerCommandLine::IntensityRescalerCommandLine()
{
  m_target="";
  m_targetems="";
  m_source.clear();
  m_sourceems.clear();
  m_label.clear();
  m_targetwindowing = true;
  m_sourcewindowing = true;
  m_classmatching = true;
  m_sigma = 3;
  m_outputsuffix = "-irescaled";
  m_outputdir = "";
}

IntensityRescalerCommandLine::~IntensityRescalerCommandLine()
{
}

#ifdef QT_GUI
void IntensityRescalerCommandLine::Create(QString filename)
{
  if (QFile(filename).exists())
  {
    std::cerr << "Error: File already exists!" << std::endl;
  }
#else
void IntensityRescalerCommandLine::Create(std::string filename)
{
  std::ifstream f(filename.c_str());
  if (f.good())
  {
    f.close() ;
    std::cerr << "Error: File already exists!" << std::endl;
    return ;
  }
#endif
  else
  {
     std::ofstream m_file;
     #ifdef QT_GUI
     m_file.open(filename.latin1(),std::ofstream::binary);
     #else
     f.close() ;
     m_file.open(filename.c_str(),std::ofstream::binary);
     #endif
     m_file << "# Example script for Intensity Rescaler" << std::endl;
     m_file << std::endl;
     m_file << "# Target image: One image which is the reference" << std::endl;
     m_file << "Target=/compute/target.gipl" << std::endl;
     m_file << std::endl;
     m_file << "# Target segmentation image: EM Segmentation of the target image" << std::endl;
     m_file << "TargetEMS=/compute/EMS-target.gipl" << std::endl;
     m_file << std::endl;

     m_file << "# Source image(s): Image(s) to be intensity rescaled base on reference image" << std::endl;
     m_file << "# Source segmentation image(s): EM Segmentation image(s)" << std::endl;
     m_file << "Source=/compute/source1.gipl" << std::endl;
     m_file << "SourceEMS=/compute/EMS-source1.gipl" << std::endl;
     m_file << "Source=/compute/source2.gipl" << std::endl;
     m_file << "SourceEMS=/compute/EMS-source2.gipl" << std::endl;
     m_file << std::endl;

     m_file << "# Label List: Usually 1,2 and 3 are the labels for White, Gray and CSF pattern" << std::endl;
     m_file << "Label=1" << std::endl;
     m_file << "Label=2" << std::endl;
     m_file << "Label=3" << std::endl;
     m_file << std::endl;

     m_file << "# Target instensity windowing: Do you want to adjust min/max intensity of the target ? [ON/OFF]" << std::endl;
     m_file << "TargetWindowing=ON" << std::endl;
     m_file << std::endl;

     m_file << "# Source instensity windowing: Do you want to adjust min/max intensity of source image(s) ? [ON/OFF]" << std::endl;
     m_file << "SourceWindowing=ON" << std::endl;
     m_file << std::endl;

     m_file << "# Class matching: Do you want to iteratively adjust classes ? [ON/OFF]" << std::endl;
     m_file << "ClassMatching=ON" << std::endl;
     m_file << std::endl;

     m_file << "# Sigma for class matching: Standard deviation for min/max adjustment" << std::endl;
     m_file << "Sigma=3" << std::endl;
     m_file << std::endl;

     m_file << "# OutputSuffix: Suffix for output filename" << std::endl;
     m_file << "OutputSuffix=-irescaled" << std::endl;
     m_file << std::endl;

     m_file << "# OutputDir: Output Directory (Not necessary)" << std::endl;
     m_file << "# OutputDir=C:/IntensityRescaling" << std::endl;
     m_file << std::endl;

     m_file << "# Intensity rescaler script end" << std::endl;
     m_file << std::endl;
     m_file.close() ;
  }
}

#ifdef QT_GUI
void IntensityRescalerCommandLine::Load(QString filename)
{
  QFile f(filename);
  if ( f.open(IO_ReadOnly) )
  {
    QTextStream t( &f );
    QString s;
    while (!t.eof())
    {
      s = t.readLine();
      s = s.stripWhiteSpace();
      if (!s.startsWith("#") && (s.length()>2))
      {
        QString m_name = s.mid(0,s.find("="));
        m_name = m_name.stripWhiteSpace();
        QString m_value = s.mid(s.find("=")+1);
        m_value = m_value.stripWhiteSpace();
        AddOption(m_name,m_value);
      }
    }
    f.close();
  }
}
#else
// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

void IntensityRescalerCommandLine::Load(std::string filename)
{
  std::ifstream f(filename.c_str());
  if ( f.good() )
  {
    std::string s;
    while (!f.eof())
    {
      getline(f,s);
      s = trim(s);
      if ( (s.size()>2) &&  !s[0] == '#' )
      {
        std::string m_name = s.substr(0,s.find("="));
        m_name = trim(m_name);
        std::string m_value = s.substr(s.find("=")+1);
        m_value = trim(m_value);
        AddOption(m_name,m_value);
      }
    }
  }
    f.close();
}
#endif


#ifdef QT_GUI
void IntensityRescalerCommandLine::AddOption(QString name,QString value)
{
  if (name.lower() == "target")     m_target=value;
  if (name.lower() == "targetems")  m_targetems=value;
  if (name.lower() == "source")     m_source.push_back(value);
  if (name.lower() == "sourceems")  m_sourceems.push_back(value);
  if (name.lower() == "label")      m_label.push_back(value.toInt());
  if (name.lower() == "targetwindowing")  {if (value.lower()=="off") m_targetwindowing=false;};
  if (name.lower() == "sourcewindowing")  {if (value.lower()=="off") m_sourcewindowing=false;};
  if (name.lower() == "classmatching")  {if (value.lower()=="off") m_classmatching=false;};
  if (name.lower() == "sigma")  m_sigma=value.toFloat();
  if (name.lower() == "outputsuffix")  m_outputsuffix=value;
  if (name.lower() == "outputdir")  m_outputdir=value;
}
#else
void IntensityRescalerCommandLine::AddOption(std::string name,std::string value)
{
  std::string nameLower ;
  std::transform( name.begin(), name.end(), nameLower.begin(), ::tolower) ;
  std::string valueLower ;
  std::transform( value.begin(), value.end(), valueLower.begin() , ::tolower ) ;
  if (nameLower == "target")     m_target=value;
  if (nameLower == "targetems")  m_targetems=value;
  if (nameLower == "source")     m_source.push_back(value);
  if (nameLower == "sourceems")  m_sourceems.push_back(value);
  if (nameLower == "label")      m_label.push_back(std::stoi(value));
  if (nameLower == "targetwindowing")  {if (valueLower=="off") m_targetwindowing=false;};
  if (nameLower == "sourcewindowing")  {if (valueLower=="off") m_sourcewindowing=false;};
  if (nameLower == "classmatching")  {if (valueLower=="off") m_classmatching=false;};
  if (nameLower == "sigma")  m_sigma=std::stof(value);
  if (nameLower == "outputsuffix")  m_outputsuffix=value;
  if (nameLower == "outputdir")  m_outputdir=value;
}
#endif


void IntensityRescalerCommandLine::DisplayOptions()
{
  std::cout << "--------------------------------------" << std::endl;
  std::cout << "Selected config file" << std::endl;
  std::cout << "--------------------------------------" << std::endl;
  #ifdef QT_GUI
  std::cout << "Target Image: " << m_target.latin1() << std::endl;
  std::cout << "Target Segmentation: " << m_targetems.latin1() << std::endl;
  #else
  std::cout << "Target Image: " << m_target << std::endl;
  std::cout << "Target Segmentation: " << m_targetems << std::endl;
  #endif
  if (m_source.size() != m_sourceems.size())
  {
    std::cerr << "Error: Number of source image(s) and segmentation(s) are not compatible! " << std::endl;
    return;
  }

  std::cout << "Source image(s):" << std::endl;
  for (unsigned int i=0;i<m_source.size();i++)
  {
    #ifdef QT_GUI
    std::cout << "\t#" << i+1 << " - Source Image : " << m_source[i].latin1() << std::endl;
    std::cout << "\t#" << i+1 << " - Source Segmentation: " << m_sourceems[i].latin1() << std::endl;
    #else
    std::cout << "\t#" << i+1 << " - Source Image : " << m_source[i] << std::endl;
    std::cout << "\t#" << i+1 << " - Source Segmentation: " << m_sourceems[i] << std::endl;
    #endif
  }

  std::cout << "Label(s):" << std::endl;
  for (unsigned int j=0;j<m_label.size();j++)
  {
    std::cout << "\t#" << j+1 << ": " << m_label[j] << std::endl;
  }

  if (m_targetwindowing)
    std::cout << "Target windowing: ON"<< std::endl;
  else
    std::cout << "Target windowing: OFF"<< std::endl;

  if (m_sourcewindowing)
    std::cout << "Source windowing: ON"<< std::endl;
  else
    std::cout << "Source windowing: OFF"<< std::endl;

  if (m_classmatching)
    std::cout << "Class matching: ON"<< std::endl;
  else
    std::cout << "Class matching: OFF"<< std::endl;

  std::cout << "Sigma: " << m_sigma << std::endl;
  #ifdef QT_GUI
  std::cout << "Output suffix: " << m_outputsuffix.latin1() << std::endl;
  std::cout << "Output directory: " << m_outputdir.latin1() << std::endl;
  #else
  std::cout << "Output suffix: " << m_outputsuffix << std::endl;
  std::cout << "Output directory: " << m_outputdir << std::endl;
  #endif
  std::cout << "--------------------------------------" << std::endl;
}

void IntensityRescalerCommandLine::Run(bool m_verbose)
{
  if (m_source.size() != m_sourceems.size())
  {
    std::cerr << "Error: Number of source image(s) and segmentation(s) are not compatible! " << std::endl;
    return;
  }

  if (m_label.size() == 0)
  {
    std::cerr << "Error: Number of labels should be greater than 0!" << std::endl;
    return;
  }

  if (m_target.size() == 0)
  {
    std::cerr << "Error: No target image!" << std::endl;
    return;
  }

  if (m_targetems.size() == 0)
  {
    std::cerr << "Error: No target segmentation!" << std::endl;
    return;
  }

  if (m_source.size() == 0)
  {
    std::cerr << "Warning: No source images! Only target intensity windowing can be performed!" << std::endl;
  }

  ImageIntensityNormalizer m_rescaler;
  double newmax;

  ImageIntensityNormalizer::ImagePointer m_image;
  ImageIntensityNormalizer::VectorList* m_vectorlist;

  ImageIntensityNormalizer::ImagePointer m_sourceimg;
  ImageIntensityNormalizer::ImagePointer m_sourcesegimg;
  ImageIntensityNormalizer::ImagePointer m_targetimg;
  ImageIntensityNormalizer::ImagePointer m_targetsegimg;

  //Add labellist
  for (unsigned k=0;k<m_label.size();k++)
    m_rescaler.AddLabel(m_label[k]);

  #ifdef QT_GUI
  m_targetimg = m_rescaler.ReadImage(m_target);
  #else
  m_targetimg = m_rescaler.ReadImage(m_target.c_str());
  #endif
  if (m_target.size() == 0)
  {
    #ifdef QT_GUI
    std::cerr << "Error: Target image " << m_target.latin1() << " could not be read!" << std::endl;
    #else
    std::cerr << "Error: Target image " << m_target << " could not be read!" << std::endl;
    #endif
    return;
  }
  #ifdef QT_GUI
  m_targetsegimg = m_rescaler.ReadImage(m_targetems);
  #else
  m_targetsegimg = m_rescaler.ReadImage(m_targetems.c_str());
  #endif
  if (m_targetsegimg.IsNull())
  {
    #ifdef QT_GUI
    std::cerr << "Error: Target segmentation " << m_targetems.latin1() << " could not be read!" << std::endl;
    #else
    std::cerr << "Error: Target segmentation " << m_targetems << " could not be read!" << std::endl;
    #endif
    return;
  }

  //Target intensity windowing
  if (m_targetwindowing)
  {
    if (m_verbose) std::cout << "Target intensity windowing ...";
     //Compute Mean values for target
    m_vectorlist = m_rescaler.ComputeMax(m_targetimg,m_targetsegimg);
    //Intensity windowing
    m_targetimg = m_rescaler.TargetIntensityWindowing(m_targetimg,*m_vectorlist,m_sigma,&newmax);
    if (m_verbose) std::cout << " done!" << std::endl;
    //Save target
    #ifdef QT_GUI
    QString m_extension = m_target.mid(m_target.findRev("."));
    QString m_targetfilename = m_target.mid(0,m_target.findRev("."));
    if (m_extension == ".gz")
    {
      m_extension = m_target.mid(m_targetfilename.findRev("."));
      m_targetfilename = m_targetfilename.mid(0,m_targetfilename.findRev("."));
    }

    if (m_outputdir.length() != 0)
      m_targetfilename = m_outputdir + "/" + m_targetfilename.mid(m_targetfilename.findRev("/")+1);
    QString m_classfilename = m_targetfilename + m_outputsuffix + "-class.txt";
    #else
    std::string m_extension = m_target.substr(m_target.find_last_of("."));
    std::string m_targetfilename = m_target.substr(0,m_target.find_last_of("."));
    if (m_extension == ".gz")
    {
      m_extension = m_target.substr(m_targetfilename.find_last_of("."));
      m_targetfilename = m_targetfilename.substr(0,m_targetfilename.find_last_of("."));
    }

    if (m_outputdir.size() != 0)
      m_targetfilename = m_outputdir + "/" + m_targetfilename.substr(m_targetfilename.find_last_of("/")+1);
    std::string m_classfilename = m_targetfilename + m_outputsuffix + "-class.txt";
    #endif
    m_targetfilename = m_targetfilename + m_outputsuffix + m_extension;
    #ifdef QT_GUI
    if (m_verbose) std::cout << "Save target irescaled: " << m_targetfilename.latin1() << "...";
    m_rescaler.SaveImage(m_targetimg,m_targetfilename);
    #else
    if (m_verbose) std::cout << "Save target irescaled: " << m_targetfilename << "...";
    m_rescaler.SaveImage(m_targetimg,m_targetfilename.c_str() );
    #endif

    if (m_verbose) std::cout << " done!" << std::endl;
    #ifdef QT_GUI
    if (m_verbose) std::cout << "Save target class model info : " << m_classfilename.latin1() << " ...";
    #else
    if (m_verbose) std::cout << "Save target class model info : " << m_classfilename << " ...";
    #endif

    // recompute class means and write the model to disk
    m_vectorlist = m_rescaler.ComputeMax(m_targetimg,m_targetsegimg);
    #ifdef QT_GUI
    ofstream efile(m_classfilename.latin1(), ios::out);
    #else
    ofstream efile(m_classfilename.c_str(), ios::out);
    #endif
    if (!efile) {
      #ifdef QT_GUI
      cerr << "Error: open of file \"" << m_classfilename.latin1() << "\" failed." << endl;
      #else
      cerr << "Error: open of file \"" << m_classfilename << "\" failed." << endl;
      #endif
      exit(-1);
    }
    for (unsigned int i = 0; i < m_vectorlist->size(); i++) {
      efile << "class " << ((*m_vectorlist)[i])[0] << " " << ((*m_vectorlist)[i])[1]
            << " " << ((*m_vectorlist)[i])[2] << endl;
    }
    efile.close();
    if (m_verbose) std::cout << " done!" << std::endl;
  }
  else
  {
    typedef itk::MinimumMaximumImageCalculator< ImageIntensityNormalizer::ImageType > CalculatorType;
    CalculatorType::Pointer calculator  =  CalculatorType::New();
    calculator->SetImage(m_targetimg);
    calculator->Compute();
    newmax = calculator->GetMaximum();
  }


  for (unsigned int i=0;i<m_source.size();i++)
  {
    #ifdef QT_QUI
    if (m_verbose) std::cout << "Processing source: " << m_source[i].latin1() << ": ";
    m_sourceimg = m_rescaler.ReadImage(m_source[i].latin1());
    m_sourcesegimg = m_rescaler.ReadImage(m_sourceems[i].latin1());
    if (m_sourceimg.IsNull())
    {
      std::cerr << "Warning: Source segmentation " << m_source[i].latin1() << " could not be read!" << std::endl;
    }
    else
    if (m_sourcesegimg.IsNull())
    {
      std::cerr << "Warning: Source segmentation " << m_sourceems[i].latin1() << " could not be read!" << std::endl;
    }
    #else
    if (m_verbose) std::cout << "Processing source: " << m_source[i] << ": ";
    m_sourceimg = m_rescaler.ReadImage(m_source[i].c_str());
    m_sourcesegimg = m_rescaler.ReadImage(m_sourceems[i].c_str());
    if (m_sourceimg.IsNull())
    {
      std::cerr << "Warning: Source segmentation " << m_source[i] << " could not be read!" << std::endl;
    }
    else
    if (m_sourcesegimg.IsNull())
    {
      std::cerr << "Warning: Source segmentation " << m_sourceems[i] << " could not be read!" << std::endl;
    }
    #endif
    else
    {
      if (m_sourcewindowing)
      {
        if (m_verbose) std::cout << "Source intensity windowing ...";
        //Compute Mean Values for source
        m_vectorlist = m_rescaler.ComputeMax(m_sourceimg,m_sourcesegimg);
        //Copmute new image for source
        m_sourceimg = m_rescaler.IntensityWindowing(m_sourceimg,*m_vectorlist,m_sigma,newmax);
        if (m_verbose) std::cout << " done!" << std::endl;
      }

      if (m_classmatching)
      {
        if (m_verbose) std::cout << "Source class matching ...";
        m_sourceimg = m_rescaler.AdjustClasses(m_targetimg,m_targetsegimg,m_sourceimg,m_sourcesegimg);
        if (m_verbose) std::cout << " done!" << std::endl;
      }

      #ifdef QT_GUI
      QString m_extension = m_source[i].mid(m_source[i].findRev("."));
      QString m_sourcefilename = m_source[i].mid(0,m_source[i].findRev("."));

      if (m_extension == ".gz")
      {
        m_extension = m_source[i].mid(m_sourcefilename.findRev("."));
        m_sourcefilename = m_sourcefilename.mid(0,m_sourcefilename.findRev("."));
      }

      if (m_outputdir.length() != 0)
        m_sourcefilename = m_outputdir + "/" + m_sourcefilename.mid(m_sourcefilename.findRev("/")+1);
      QString m_classfilename = m_sourcefilename + m_outputsuffix + "-class.txt";
      #else
      std::string m_extension = m_source[i].substr(m_source[i].find_last_of("."));
      std::string m_sourcefilename = m_source[i].substr(0,m_source[i].find_last_of("."));

      if (m_extension == ".gz")
      {
        m_extension = m_source[i].substr(m_sourcefilename.find_last_of("."));
        m_sourcefilename = m_sourcefilename.substr(0,m_sourcefilename.find_last_of("."));
      }

      if (m_outputdir.length() != 0)
        m_sourcefilename = m_outputdir + "/" + m_sourcefilename.substr(m_sourcefilename.find_last_of("/")+1);
      std::string m_classfilename = m_sourcefilename + m_outputsuffix + "-class.txt";
      #endif
      m_sourcefilename = m_sourcefilename + m_outputsuffix + m_extension;
      #ifdef QT_GUI
      if (m_verbose) std::cout << "Save source irescaled: " << m_sourcefilename.latin1() << " ...";
      m_rescaler.SaveImage(m_sourceimg,m_sourcefilename.latin1());
      if (m_verbose) std::cout << " done!" << std::endl;

      if (m_verbose) std::cout << "Save source class model info: " << m_classfilename.latin1() << " ...";

      // recompute class means and write the model to disk
      m_vectorlist = m_rescaler.ComputeMax(m_sourceimg,m_sourcesegimg);
      ofstream efile(m_classfilename.latin1(), ios::out);
      if (!efile) {
        cerr << "Error: open of file \"" << m_classfilename.latin1() << "\" failed." << endl;
        exit(-1);
      }
      #else
      if (m_verbose) std::cout << "Save source irescaled: " << m_sourcefilename << " ...";
      m_rescaler.SaveImage(m_sourceimg,m_sourcefilename.c_str());
      if (m_verbose) std::cout << " done!" << std::endl;

      if (m_verbose) std::cout << "Save source class model info: " << m_classfilename << " ...";

      // recompute class means and write the model to disk
      m_vectorlist = m_rescaler.ComputeMax(m_sourceimg,m_sourcesegimg);
      ofstream efile(m_classfilename.c_str(), ios::out);
      if (!efile) {
        cerr << "Error: open of file \"" << m_classfilename << "\" failed." << endl;
        exit(-1);
      }
      #endif
      for (unsigned int i = 0; i < m_vectorlist->size(); i++) {
        efile << "class " << ((*m_vectorlist)[i])[0] << " " << ((*m_vectorlist)[i])[1]
              << " " << ((*m_vectorlist)[i])[2] << endl;
      }
      efile.close();
      if (m_verbose) std::cout << " done!" << std::endl;
    }
  }
}
