/*=========================================================================

  Program:   QtDisclaimer
  Module:    $RCSfile: Disclaimer.h,v $
  Language:  C++
  Date:      $Date: 2005/09/20 14:31:38 $
  Version:   $Revision: 1.3 $
  Author:    Matthieu Jomier

  Copyright (c) 2005 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef Disclaimer_h
#define Disclaimer_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <stdio.h>
#include <vector>
#include <string>

class Disclaimer 
{ 
public:
    
  Disclaimer();
  virtual ~Disclaimer();

  virtual void SetApplication(std::string application);
  virtual void SetDirector(std::string director);
  virtual void SetDeveloper(std::string developer);
  virtual void Update();
  virtual void Show();

protected:
  std::string m_application;
  std::vector<std::string> m_director;
  std::vector<std::string> m_developer;
  std::string m_title;
  std::string m_author;
  std::string m_disclaimer;
};

#endif
