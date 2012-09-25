#ifndef _XMLREADER_H
#define _XMLREADER_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <qstring.h>

using namespace std;


namespace neurolib {
  
class XMLReader
{
public:
   XMLReader();
   ~XMLReader();
   int Open(const char* filename);
   QString GetBalise();
   QString GetCurrentBalise();
   QString GetValue();
   void Close();

private:
   std::ifstream networkfile;
   QString m_value;
   QString m_balise;


};
}

#endif
