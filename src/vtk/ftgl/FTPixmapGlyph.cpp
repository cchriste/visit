#include  "FTPixmapGlyph.h"
#ifdef FTGL_DEBUG
  #include "mmgr.h"
#endif


FTPixmapGlyph::FTPixmapGlyph( FT_Glyph _glyph)
:  FTGlyph(),
  destWidth(0),
  destHeight(0),
  numGreys(0),
  data(0)
{
  this->glyph = _glyph;
  bBox = FTBBox(this->glyph);
  advance = (float)(this->glyph->advance.x >> 16);
}

void FTPixmapGlyph::ConvertGlyph( const int horizontal,
                                  const FTGLRenderContext *context)
{
  // This function will always fail if the glyph's format isn't scalable????
  err = FT_Glyph_To_Bitmap( &this->glyph, ft_render_mode_normal, 0, 1);
  if( err || ft_glyph_format_bitmap != this->glyph->format)
  {
    return;
  }

  FT_BitmapGlyph  bitmap = reinterpret_cast<FT_BitmapGlyph>(this->glyph);
  FT_Bitmap*      source = &bitmap->bitmap;

  //check the pixel mode
  //ft_pixel_mode_grays
      
  int srcWidth;
  int srcHeight;
  int srcPitch;
  if (horizontal)
    {
    srcWidth = source->width;
    srcHeight = source->rows;
    srcPitch = source->pitch;
    }
  else
    {
    srcWidth = source->rows;
    srcHeight = source->width;
    srcPitch = 1;
    }
  
  // FIXME What about dest alignment?
  destWidth = srcWidth;
  destHeight = srcHeight;
    
  if( destWidth && destHeight)
    {
    data = new unsigned char[destWidth * destHeight * 4];
    
    // Get the current glColor.
    float ftglColour[4];
#ifdef FTGL_SUPPORT_MANGLE_MESA
    if (context && context->UseMangleMesa)
      {
      this->GetCurrentColorMesa(ftglColour, context);
      }
    else
#endif
      {
      this->GetCurrentColorOpenGL(ftglColour, context);
      }
      
#if 1
    unsigned char red = static_cast<unsigned char>(ftglColour[0]*255.0f);
    unsigned char green = static_cast<unsigned char>(ftglColour[1]*255.0f);
    unsigned char blue = static_cast<unsigned char>(ftglColour[2]*255.0f);

    unsigned char *src = source->buffer;
    unsigned char *src_row;

    unsigned char *dest;
    size_t dest_step;
    if (horizontal)
      {
      dest = data + ((destHeight - 1) * destWidth) * 4;
      dest_step = destWidth * 4 * 2;
      }
    else
      {
      dest = data;
      dest_step = 0;
      }

    if (ftglColour[3] == 1.0f)
      {
      for(int y = 0; y < srcHeight; ++y)
        {
        src_row = src;
        for(int x = 0; x < srcWidth; ++x)
          {
          *dest++ = red;
          *dest++ = green;
          *dest++ = blue;
          if (horizontal)
            {
            *dest++ = *src_row++;
            }
          else
            {
            *dest++ = *src_row;
            src_row += srcHeight;
            }
          }
        src += srcPitch;
        dest -= dest_step;
        }
      }
    else
      {
      for(int y = 0; y < srcHeight; ++y)
        {
        src_row = src;
        for(int x = 0; x < srcWidth; ++x)
          {
          *dest++ = red;
          *dest++ = green;
          *dest++ = blue;
          if (horizontal)
            {
            *dest++ = static_cast<unsigned char>(ftglColour[3] * *src_row++);
            }
          else
            {
            *dest++ = static_cast<unsigned char>(ftglColour[3] * *src_row);
            src_row += srcHeight;
            }
          }
        src += srcPitch;
        dest -= dest_step;
        }
      }
#else
      
    for(int y = 0; y < srcHeight; ++y)
      {
      --destHeight;
      for(int x = 0; x < srcWidth; ++x)
        {
        *( data + ( destHeight * destWidth  + x) * 4 + 0) = static_cast<unsigned char>( ftglColour[0] * 255.0f);
        *( data + ( destHeight * destWidth  + x) * 4 + 1) = static_cast<unsigned char>( ftglColour[1] * 255.0f);
        *( data + ( destHeight * destWidth  + x) * 4 + 2) = static_cast<unsigned char>( ftglColour[2] * 255.0f);
        *( data + ( destHeight * destWidth  + x) * 4 + 3) = static_cast<unsigned char>( ftglColour[3] * (*( source->buffer + ( y * srcPitch) + x)));
        }      
      }

#endif  

    destHeight = srcHeight;
    }
  
  numGreys = source->num_grays;
  if (horizontal)
    {
    pos.x = bitmap->left;
    pos.y = srcHeight - bitmap->top;
    }
  else
    {
    pos.x = - bitmap->top;
    pos.y = - bitmap->left;
    }
  
  this->glyphHasBeenConverted = 1;
}


FTPixmapGlyph::~FTPixmapGlyph()
{
  if( data)
    delete [] data;
}


float FTPixmapGlyph::Render( const FT_Vector& pen,
                             const int horizontal,
                             const FTGLRenderContext *context)
{
  if (!this->glyphHasBeenConverted)
    {
    this->ConvertGlyph(horizontal, context);
    }

  if( data)
    {
#ifdef FTGL_SUPPORT_MANGLE_MESA
    if (context && context->UseMangleMesa)
      {
      this->RenderMesa(pen, context);
      }
    else
#endif
      {
      this->RenderOpenGL(pen, context);
      }
    }

  return advance;
}
