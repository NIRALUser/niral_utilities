#include "XMLIniIO.h"

XMLIniIO::XMLIniIO()
{
  m_reader = new neurolib::XMLReader();
  m_writer = new neurolib::XMLWriter();
}


XMLIniIO::~XMLIniIO()
{
}

void XMLIniIO::SetFileName(const char* filename)
{
   m_filename = filename;
}


void XMLIniIO::Read()
{
   if (m_reader->Open(m_filename) == -1)
     return;
   
   balisestruct m_value;
   m_value.balise = "none";
   while(m_value.balise.length() != 0)
   {
      m_value.balise =  m_reader->GetBalise();
      m_value.value =  m_reader->GetValue();
      if (m_value.balise.length() != 0)
         m_list.push_back(m_value);
   }
   m_reader->Close();
   
}

void XMLIniIO::Write()
{
   m_writer->Open(m_filename);
   balisestruct m_value;
   for (unsigned int i=0;i<m_list.size();i++)
   {
      m_value =  m_list[i];
      m_writer->Write(m_value.balise,m_value.value);
   }
   m_writer->Close();
}

QString XMLIniIO::Find(QString balise)
{
  for(vector<balisestruct>::iterator itNum = m_list.begin(); itNum < m_list.end(); itNum++)
      if ((*itNum).balise == balise)
          return (*itNum).value;

  return "";
}

void XMLIniIO::Update(QString balise,QString value)
{
  balisestruct m_value;
  bool m_append = true;
  for (unsigned int i=0;i<m_list.size();i++)
   {
    if (m_list[i].balise == balise)
    {
      m_list[i].value = value;
      m_append = false;
    }
  }

  if (m_append)
  {
     m_value.balise = balise;
     m_value.value = value;
     m_list.push_back(m_value);
  }
}

void XMLIniIO::Update(QString balise,bool value)
{
  balisestruct m_value;
  bool m_append = true;
  for (unsigned int i=0;i<m_list.size();i++)
   {
    if (m_list[i].balise == balise)
    {
      m_list[i].value = QString("%1").arg(value);
      m_append = false;
    }
  }

  if (m_append)
  {
     m_value.balise = balise;
     m_value.value = QString("%1").arg(value);
     m_list.push_back(m_value);
  }
}

void XMLIniIO::Remove(QString balise)
{
  for(vector<balisestruct>::iterator itNum = m_list.begin(); itNum < m_list.end(); itNum++)
      if ((*itNum).balise == balise)
          m_list.erase(itNum);


}
