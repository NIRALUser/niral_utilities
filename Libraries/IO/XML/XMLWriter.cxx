#include "XMLWriter.h"

namespace neurolib {

XMLWriter::XMLWriter()
{
   /** default values */
   moduleid = 0;
   treeoffset = 0;
}


XMLWriter::~XMLWriter()
{
}

int XMLWriter::Open(const char* filename)
{
   /** Open file for writing */
   networkfile.open(filename,std::ofstream::binary);
   if (networkfile == NULL) return -1;
   
   return 0;
}


void XMLWriter::Start(char* name)
{
   startTab();
   networkfile << "<" << name << ">" << std::endl;
}

void XMLWriter::End(char* name)
{
   endTab();
   networkfile << "</" << name << ">" << std::endl;
   
}

void XMLWriter::startTab()
{
   for (int i=0;i<treeoffset;i++)
   networkfile << "\t";
   treeoffset++;
}

void XMLWriter::endTab()
{
   treeoffset--;
   for (int i=0;i<treeoffset;i++)
      networkfile << "\t";
}

void XMLWriter::Write(QString balise,QString value)
{
   /** Write Identification number which is necessary for connections */
   startTab();
   networkfile << "<" << balise.latin1() << ">" << value.latin1() << "</" << balise.latin1() << ">" << std::endl;
   treeoffset--;
}

void XMLWriter::Write(QString balise,int value)
{
   /** Write Identification number which is necessary for connections */
   startTab();
   networkfile << "<" << balise.latin1() << ">" << value << "</" << balise.latin1() << ">" << std::endl;
   treeoffset--;
}

void XMLWriter::Write(QString balise,float value)
{
   /** Write Identification number which is necessary for connections */
   startTab();
//   std::cout << "Value" << balise.latin1() << " : " << value <<std::endl;
   networkfile << "<" << balise.latin1() << ">" << value << "</" << balise.latin1() << ">" << std::endl;
   treeoffset--;
}


void XMLWriter::Close()
{
   networkfile.close();   
}

}
