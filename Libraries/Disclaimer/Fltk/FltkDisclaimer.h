/*=========================================================================

  Program:   FltkDisclaimer
  Module:    $RCSfile: FltkDisclaimer.h,v $
  Language:  C++
  Date:      $Date: 2005/09/17 08:00:59 $
  Version:   $Revision: 1.1 $
  Author:    Matthieu Jomier

  Copyright (c) 2005 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef FltkDisclaimer_h
#define FltkDisclaimer_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "FltkDisclaimerGUI.h"
#include "Disclaimer.h"
#include <vector>
#include <string.h>

class FltkDisclaimer : public FltkDisclaimerGUI, public Disclaimer
{ 
public:
    
  FltkDisclaimer();
  ~FltkDisclaimer();

  void Show();
  void Accept();
  void Reject();

private:
  Fl_Text_Buffer* m_buffer;
  Fl_Text_Buffer* m_buffer2;
};

#endif
