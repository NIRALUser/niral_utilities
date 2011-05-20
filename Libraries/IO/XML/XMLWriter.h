#ifndef _XMLWRITER_H
#define _XMLWRITER_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <qstring.h>

//using namespace std;

namespace neurolib {
  
class XMLWriter
{
public:
   XMLWriter();
   ~XMLWriter();
   void Start(char *);
   void End(char *);
   void Write(QString balise,QString value);
   void Write(QString balise,int value);
   void Write(QString balise,float value);
   int Open(const char* filename);
   void Close();
   void startTab();
   void endTab();

protected:

private:
   char filename[300];
   int moduleid;
   std::ofstream networkfile;
   int treeoffset;
};

}

#endif
