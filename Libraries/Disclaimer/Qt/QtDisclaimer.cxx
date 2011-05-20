/*=========================================================================
  
  Program:   QtDisclaimer
  Module:    $RCSfile: QtDisclaimer.cxx,v $
  Language:  C++
  Date:      $Date: 2008/06/17 20:35:08 $
  Version:   $Revision: 1.2 $
  Author:    Matthieu Jomier

  Copyright (c) 2004 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "QtDisclaimer.h"
#include <qlabel.h>
#include <stdlib.h>

QtDisclaimer::QtDisclaimer( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : QtDisclaimerGUI( parent, name, modal,  fl )
{
}

QtDisclaimer::~QtDisclaimer()
{
}


void QtDisclaimer::Show()
{
  Update();
  setCaption( QString(m_application.c_str()) + " - Disclaimer (Qt)");
  QString m_devtext = m_title.c_str();
  m_devtext += "\n\n";
  m_devtext += m_author.c_str();
  g_developer->setText(m_devtext);
  g_disclaimer->setText(m_disclaimer.c_str());
  if (!exec())
    exit(0);
}

