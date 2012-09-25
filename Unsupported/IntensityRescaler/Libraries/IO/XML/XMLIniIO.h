#ifndef _XMLINIIO_H
#define _XMLINIIO_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "XMLWriter.h"
#include "XMLReader.h"

class XMLIniIO
{
public:
   XMLIniIO();
   ~XMLIniIO();
   void SetFileName(const char* filename);
   void Read();
   void Write();
   QString Find(QString balise);
   void Update(QString balise,QString value);
   void Update(QString balise,bool value);
   void Remove(QString balise);

protected:

private:
   QString m_filename;
   struct balisestruct
   {
      QString balise;
      QString value;
   };
   
   std::vector<balisestruct> m_list;
   neurolib::XMLReader* m_reader;
   neurolib::XMLWriter* m_writer;   
};

#endif
