/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: NeuroLibOpenGLSliceTexture.txx,v $
  Language:  C++
  Date:      $Date: 2005/04/02 12:07:26 $
  Version:   $Revision: 1.3 $
  Copyright (c) 2003 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
template<class TPixel>
NeuroLibOpenGLSliceTexture<TPixel>
::NeuroLibOpenGLSliceTexture()
{
  // Set to -1 to force a call to 'generate'
  m_IsTextureInitalized = false;

  // Create the filter
  m_PadFilter = FilterType::New();

  // Set the update time to -1
  m_UpdateTime = 0;

  // Init the GL settings to uchar, luminance defautls, which are harmless
  m_GlComponents = 1;
  m_GlFormat = GL_LUMINANCE;
  m_GlType = GL_UNSIGNED_BYTE;
  m_Image = 0;
}

template<class TPixel>
NeuroLibOpenGLSliceTexture<TPixel>
::~NeuroLibOpenGLSliceTexture()
{
  if(m_IsTextureInitalized)
    glDeleteTextures(1,&m_TextureIndex);
}

template<class TPixel>
void
NeuroLibOpenGLSliceTexture<TPixel>
::SetImage(ImagePointer inImage)
{
  if(m_Image != inImage)
  {
    m_Image = inImage;
    m_UpdateTime = 0;
  }
}

template<class TPixel>
void
NeuroLibOpenGLSliceTexture<TPixel>
::Update()
{
  // Better have an image
    //Previous Assert
  if(!m_Image)
    return;

  // Update the image (necessary?)
  if(m_Image->GetSource())
  {
     m_Image->GetSource()->UpdateLargestPossibleRegion();
  }

  // Check if everything is up-to-date and no computation is needed
 // if (m_IsTextureInitalized && m_UpdateTime == m_Image->GetPipelineMTime())
  //  return;

  // Promote the image dimensions to powers of 2
  itk::Size<2> szImage = m_Image->GetLargestPossibleRegion().GetSize();
  m_TextureSize[0] = 1; //Vector2ui(1);
  m_TextureSize[1] = 1;

  // Use shift to quickly double the coordinates
  for (unsigned int i=0;i<2;i++)
    while (m_TextureSize[i] < szImage[i])
      m_TextureSize[i] <<= 1;

  // Compute the pad offset
  vnl_vector_fixed<unsigned long,2> 
    offset(m_TextureSize[0] - szImage[0],m_TextureSize[1] - szImage[1]);

  // Set the parameters of the pad filter
   m_PadFilter->SetInput(m_Image);
   m_PadFilter->SetPadUpperBound(offset.data_block());

  // Apply the padding
   m_PadFilter->UpdateLargestPossibleRegion();

  // Create the texture index if necessary
  if(!m_IsTextureInitalized)
    {
      // Generate one texture
      glGenTextures(1,&m_TextureIndex);
      m_IsTextureInitalized = true;
    }

  // Select the texture for pixel pumping
  glBindTexture(GL_TEXTURE_2D,m_TextureIndex);

  // Properties for the texture
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  
  // Pump the pixels into the texture
  /* glTexImage2D(GL_TEXTURE_2D,0,m_GlComponents,
               m_TextureSize(0),m_TextureSize(1),
               0,m_GlFormat,m_GlType,
               m_PadFilter->GetOutput()->GetBufferPointer()); */


  //m_TextureSize[0] = szImage[0]; 
  //m_TextureSize[1] = szImage[1];
  gluBuild2DMipmaps(
    GL_TEXTURE_2D, m_GlComponents, 
    m_TextureSize[0], m_TextureSize[1],
    m_GlFormat,m_GlType,
//  m_Image->GetBufferPointer());

    m_PadFilter->GetOutput()->GetBufferPointer());

  // Remember the image's timestamp
   m_UpdateTime = m_Image->GetPipelineMTime();
}

template<class TPixel>
void
NeuroLibOpenGLSliceTexture<TPixel>
::Draw()
{
  // Update the texture
  Update();

  // Should have a texture number
  //Previous Assert
  if(!m_IsTextureInitalized)
   return;

  // GL settings
  glPushAttrib(GL_TEXTURE_BIT);
  glEnable(GL_TEXTURE_2D);

  // Select our texture
  glBindTexture(GL_TEXTURE_2D,m_TextureIndex);

  // Set the color to the background color
  double clrBackground[3];
  clrBackground[0] = 255;
  clrBackground[1] = 255;
  clrBackground[2] = 255;
  glColor3dv(clrBackground);

  int w = m_Image->GetBufferedRegion().GetSize()[0];
  int h = m_Image->GetBufferedRegion().GetSize()[1];
  double tx = w * 1.0 / m_TextureSize[0];
  double ty = h * 1.0 / m_TextureSize[1];
    
  // Draw quad 
  glBegin(GL_QUADS);
  glTexCoord2d(0,0);
  glVertex2d(0,0);
  glTexCoord2d(0,ty);
  glVertex2d(0,h);
  glTexCoord2d(tx,ty);
  glVertex2d(w,h);
  glTexCoord2d(tx,0);
  glVertex2d(w,0);
  glEnd();

  glPopAttrib();
}

template<class TPixel>
void
NeuroLibOpenGLSliceTexture<TPixel>
::DrawTransparent(unsigned char alpha)
{
  // Update the texture
  Update();

  // Should have a texture number
  //Previous Assert
  if(!m_IsTextureInitalized)
  return;

  // GL settings
  glPushAttrib(GL_TEXTURE_BIT | GL_COLOR_BUFFER_BIT);
  glEnable(GL_TEXTURE_2D);

  // Turn on alpha blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  // Select our texture
  glBindTexture(GL_TEXTURE_2D,m_TextureIndex);

  // Set the color to white
  glColor4ub(255,255,255,alpha);
    
  // Draw quad 
  glBegin(GL_QUADS);
  glTexCoord2d(0,0);
  glVertex2d(0,0);
  glTexCoord2d(0,1);
  glVertex2d(0,m_TextureSize[1]);
  glTexCoord2d(1,1);
  glVertex2d(m_TextureSize[0],m_TextureSize[1]);
  glTexCoord2d(1,0);
  glVertex2d(m_TextureSize[0],0);
  glEnd();

  glPopAttrib();
}

