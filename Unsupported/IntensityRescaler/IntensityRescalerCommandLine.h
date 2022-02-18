#ifndef IntensityRescalerCommandLine_h
#define IntensityRescalerCommandLine_h

#ifdef QT_GUI
#include <qstring.h>
#endif
#include <string>
#include <vector>
#include <cstdlib>

class IntensityRescalerCommandLine
{
public:

  IntensityRescalerCommandLine();
  ~IntensityRescalerCommandLine();
  #ifdef QT_GUI
  void Create(QString filename);
  void Load(QString filename);
  void AddOption(QString name,QString value);
  #else
  void Create(std::string filename);
  void Load(std::string filename);
  void AddOption(std::string name,std::string value);
  #endif
  void DisplayOptions();
  void Run(bool m_verbose=false);

private:
  std::vector<int> m_label;
  bool m_targetwindowing;
  bool m_sourcewindowing;
  bool m_classmatching;
  float m_sigma;
  #ifdef QT_GUI
  QString m_target;
  QString m_targetems;
  QString m_outputsuffix;
  QString m_outputdir;
  std::vector<QString> m_source;
  std::vector<QString> m_sourceems;
  #else
  std::string m_target;
  std::string m_targetems;
  std::string m_outputsuffix;
  std::string m_outputdir;
  std::vector<std::string> m_source;
  std::vector<std::string> m_sourceems;
  #endif

};

#endif

