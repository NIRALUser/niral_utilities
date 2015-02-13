/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $HeadURL$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// This file is intended to be compiled and linked against a shared
// library CLP to prevent the need to compile twice.

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#ifdef WIN32
#ifdef BUILD_SHARED_LIBS
#define MODULE_IMPORT __declspec(dllimport)
#else
#define MODULE_IMPORT
#endif
#else
#define MODULE_IMPORT
#endif

extern "C" MODULE_IMPORT int ModuleEntryPoint(const int, const char **);

int main(const int argc, const char** argv)
{
  return ModuleEntryPoint(argc, argv);
}
