#ifndef _PALETTE_H
#define _PALETTE_H

#include <math.h>

class Palette
{
public:

  struct colorstruct
  {
    int r;
    int g;
    int b;
  };

  enum PaletteType {PALETTE_GRAY, 
                    PALETTE_RGB, 
                    PALETTE_RB,
                    PALETTE_MATLAB};

  static colorstruct* Generate(PaletteType scale)
  {
    colorstruct* colormap = new colorstruct[65536];

    unsigned short red, green, blue;
    int i;
    float a = (float)65535/(float)(32768*(32768-65535));
    float b = -65535*a;
    switch (scale) {
    //Classical gray palette
    case PALETTE_GRAY:
      for (i=0;i<65536;i++) {
      colormap[i].r = i;
      colormap[i].g = i;	
      colormap[i].b = i;
      }

      break;

    case PALETTE_RGB:
      for (i=0;i<65536;i++) {
      if (i <32768)
        red = (int)(a*(32768+i)*(32768+i)+b*(32768+i));
      else
        red = 0;
      green = (int)(a*(i)*(i)+b*(i));


      if (i > 32768)
        blue = (int)(a*(i-32768)*(i-32768)+b*(i-32768));
      else
        blue = 0;
      colormap[i].r = red;
      colormap[i].g = green;	
      colormap[i].b = blue;
      }

      break;

    case PALETTE_RB:
      for (i=0;i<65536;i++) {
      if (i <32768)
        blue = (int)(a*(32768+i)*(32768+i)+b*(32768+i));
      else
        blue = 0;

      red = (int)(a*(i)*(i)+b*(i));


      if (i > 32768)
        {
        red = 65535; //a*(i-32768)*(i-32768)+b*(i-32768);
        green = int(a*(i-32768)*(i-32768)+b*(i-32768));
        }
      else
        green = 0;

      colormap[i].r = red;
      colormap[i].g = green;	
      colormap[i].b = blue;
      }
      break;

      //Matlab palette
    case PALETTE_MATLAB:
      for (i=0;i<65536;i++) {
				
      //blue
      if (i<8192)
        blue = 4*(i+8192);
      else
        if ((i>=8192) && (i<24576))
          blue = 65535;
        else
          if ((i>=24576) && (i<40960))
            blue = 65535-4*(i-24576);
          else
            blue = 0;
				
      //green
      if ((i>8192) && (i<24576))
        green = 4*(i-8192);
      else
        if ((i>=24576) && (i<40960))
          green = 65535;
        else
          if ((i>=40960) && (i<57344))
            green = 65535-4*(i-40960);
          else
            green = 0;


      //red 
      if ((i>24576) && (i<40960))
        red = 4*(i-24576);
      else
        if ((i>=40960) && (i<57344))
          red = 65535;
        else
          if ((i>=57344))
            red = 65535-4*(i-57344);
          else
            red = 0;



      colormap[i].r = red;
      colormap[i].g = green;	
      colormap[i].b = blue;
      }
      break;
    }

    return colormap;
  }

};
  

#endif

