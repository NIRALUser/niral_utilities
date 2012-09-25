/*=========================================================================

  Program:   QtDisclaimer
  Module:    $RCSfile: QtDisclaimer.h,v $
  Language:  C++
  Date:      $Date: 2005/09/17 08:00:26 $
  Version:   $Revision: 1.1 $
  Author:    Matthieu Jomier

  Copyright (c) 2005 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef QtDisclaimer_h
#define QtDisclaimer_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "Disclaimer.h"
#include "QtDisclaimerGUI.h"
#include <vector>

class QtDisclaimer : public QtDisclaimerGUI, public Disclaimer
{ 
public:
    
  QtDisclaimer( QWidget* parent = 0, const char* name = 0, bool modal= true, WFlags fl = 0 );
  ~QtDisclaimer();

  void Show();



};

#endif
