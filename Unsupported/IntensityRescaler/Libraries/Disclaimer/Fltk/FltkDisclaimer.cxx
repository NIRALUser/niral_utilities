/*=========================================================================
  
  Program:   FltkDisclaimer
  Module:    $RCSfile: FltkDisclaimer.cxx,v $
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

#include "FltkDisclaimer.h"
#include <stdlib.h>

FltkDisclaimer::FltkDisclaimer( )
    : FltkDisclaimerGUI()
{
  MakeWindow();
  m_buffer = new Fl_Text_Buffer();
  g_developer->buffer(m_buffer);
  m_buffer2 = new Fl_Text_Buffer();
  g_disclaimer->buffer(m_buffer2);
}

FltkDisclaimer::~FltkDisclaimer()
{
}

void FltkDisclaimer::Show()
{
  Update();
  std::string m_caption =  m_application;
  m_caption += " - Disclaimer";
  FltkDisclaimerGUIWindow->label(m_caption.c_str());
  g_developer->insert(m_title.c_str());
  g_developer->insert("\n");
  g_developer->insert(m_author.c_str());
  g_disclaimer->insert(m_disclaimer.c_str());
  FltkDisclaimerGUIWindow->show();
  Fl::run();
}

void FltkDisclaimer::Accept()
{
  FltkDisclaimerGUIWindow->hide();
  delete m_buffer;
  delete m_buffer2;
}

void FltkDisclaimer::Reject()
{
  delete m_buffer;
  delete m_buffer2;
  exit(0);
}
