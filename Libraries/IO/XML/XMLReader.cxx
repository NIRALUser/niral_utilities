#include "XMLReader.h"
#include "qmessagebox.h"

namespace neurolib {

XMLReader::XMLReader()
{
}


XMLReader::~XMLReader()
{
}

int XMLReader::Open(const char* filename)
{
   /** Open file for reading */
   networkfile.open(filename,ifstream::binary);
   if (networkfile == NULL) return -1;

   
   return 0;
} 

QString XMLReader::GetValue()
{
   return m_value;
}

QString XMLReader::GetBalise()
{
   char* data = (char*)malloc(1000);
   networkfile.getline(data,1000);
   QString line = data;
   int begin_balise_start = line.find("<");
   int begin_balise_end =  line.find(">",begin_balise_start+1);
   
   int end_balise_begin = line.find("<",begin_balise_end+1);

   m_balise = line.mid(begin_balise_start+1,begin_balise_end-begin_balise_start-1);

   if (end_balise_begin != -1)
      m_value = line.mid(begin_balise_end+1,end_balise_begin-begin_balise_end-1);
   else
       {
           m_value = "";
       }

   return m_balise;
}


QString XMLReader::GetCurrentBalise()
{
   return m_balise;
}


void XMLReader::Close()
{
   networkfile.close();
}

}
