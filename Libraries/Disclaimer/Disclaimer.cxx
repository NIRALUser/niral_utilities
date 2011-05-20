/*=========================================================================
  
  Program:   Disclaimer
  Module:    $RCSfile: Disclaimer.cxx,v $
  Language:  C++
  Date:      $Date: 2005/09/17 07:54:13 $
  Version:   $Revision: 1.1 $
  Author:    Matthieu Jomier

  Copyright (c) 2004 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "Disclaimer.h"
#include <iostream>

Disclaimer::Disclaimer()
{
}

Disclaimer::~Disclaimer()
{
}

void Disclaimer::SetApplication(std::string application)
{
  m_application = application;
}

void Disclaimer::SetDirector(std::string director)
{
  m_director.push_back(director);
}

void Disclaimer::SetDeveloper(std::string developer)
{
  m_developer.push_back(developer);
}

void Disclaimer::Show()
{
  Update();
  std::cout << m_title.c_str() << std::endl;
  std::cout << std::endl;
  std::cout <<  m_author.c_str() << std::endl; 
  std::cout << std::endl;
  std::cout <<  m_disclaimer.c_str() << std::endl; 
}

void Disclaimer::Update()
{
  m_title= "UNC Neuroimage Analysis Group";

  m_author = "";

 if (m_director.size() == 1)
 {
   m_author += "Director: ";
   m_author += m_director[0];
   m_author += "\n";
 }
 else
 {
   for (unsigned int i=0;i<m_developer.size();i++)
   {
     if (i==0)
     {
      m_author +="Directors: ";
      m_author += m_director[i];
      m_author += "\n";
     }
      else
     {
       m_author += "          ";
       m_author += m_director[i];
       m_author += "\n";
     }
   }
 }



 if (m_developer.size() == 1)
 {
    m_author += "Developer: ";
    m_author += m_developer[0];
    m_author += "\n";
 }
 else
 {
   for (unsigned int i=0;i<m_developer.size();i++)
   {
     if (i==0)
     {
        m_author += "Developers: ";
        m_author += m_developer[i];
        m_author += "\n";
     }
     else
     {
       m_author += "            ";
       m_author += m_developer[i];
       m_author += "\n";
     }
   }
 }


 m_disclaimer = m_application;
 m_disclaimer += ", and the contents of the documentation, are intended for educational, research, and\n";
 m_disclaimer += "informational purposes only.\n";
 m_disclaimer += m_application;
 m_disclaimer += " or information derived from ";
 m_disclaimer += m_application;
 m_disclaimer += "may be used only for these purposes and may not\n";
 m_disclaimer += "under any circumstances whatsoever be used for clinical purposes.\n";
 m_disclaimer += "This software has not been approved by the FDA and is not intended for treating or diagnosing\n";
 m_disclaimer += "human subjects and the recipient and user will not use the ";
 m_disclaimer += m_application;
 m_disclaimer += " software for such purposes.\n\n";
 m_disclaimer += "The ";
 m_disclaimer +=  m_application; 
 m_disclaimer += " copyright holders and contributors and all affiliated organizations shall not be liable\n";
 m_disclaimer += "for any damages arising out of the use of ";
 m_disclaimer += m_application;
 m_disclaimer += " by any party for any purpose.\n";
 m_disclaimer += "\n";
 m_disclaimer += "Every publication using this software needs to properly acknowledge the source\n\n";
 m_disclaimer += "More information available on the NeuroLib website: http://www.ia.unc.edu/dev/license\n";

}
