#ifndef IntensityRescalerCommandLine_h
#define IntensityRescalerCommandLine_h
  
#include <qstring.h>
#include <vector>
#include <stdlib.h>

class IntensityRescalerCommandLine
{
public:

  IntensityRescalerCommandLine();
  ~IntensityRescalerCommandLine();
  void Create(QString filename);
  void Load(QString filename);
  void AddOption(QString name,QString value);
  void DisplayOptions();
  void Run(bool m_verbose=false);

private:
  QString m_target;
  QString m_targetems;
  std::vector<QString> m_source;
  std::vector<QString> m_sourceems;
  std::vector<int> m_label;
  bool m_targetwindowing;
  bool m_sourcewindowing;
  bool m_classmatching;
  float m_sigma;
  QString m_outputsuffix;
  QString m_outputdir;

};

#endif

