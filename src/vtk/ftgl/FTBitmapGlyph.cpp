#include  "FTBitmapGlyph.h"
#ifdef FTGL_DEBUG
  #include "mmgr.h"
#endif


FTBitmapGlyph::FTBitmapGlyph( FT_Glyph _glyph)
:  FTGlyph(),
  destWidth(0),
  destHeight(0),
  data(0)
{
  this->glyph = _glyph;
  bBox = FTBBox(this->glyph);
  advance = (float)(this->glyph->advance.x >> 16);
}

void FTBitmapGlyph::ConvertGlyph( const int horizontal)
{
  // This function will always fail if the glyph's format isn't scalable????
  err = FT_Glyph_To_Bitmap( &glyph, ft_render_mode_mono, 0, 1);
  if( err || ft_glyph_format_bitmap != this->glyph->format)
  {return;}

  FT_BitmapGlyph  bitmap = (FT_BitmapGlyph)this->glyph;
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
    data = new unsigned char[srcPitch * destHeight];
    
#if 1
    unsigned char *src = source->buffer;
    unsigned char *src_row;
    
    unsigned char *dest;
    size_t dest_step;
    if (horizontal)
      {
      dest = data + ((destHeight - 1) * srcPitch);
      dest_step = srcPitch * 2;
      }
    else
      {
      dest = data;
      dest_step = 0;
      }
    
    for(int y = 0; y < srcHeight; ++y)
      {
      src_row = src;
      for(int x = 0; x < srcPitch; ++x)
        {
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

#else

    for(int y = 0; y < srcHeight; ++y)
      {
      --destHeight;
      for(int x = 0; x < srcPitch; ++x)
        {
        *( data + ( destHeight * srcPitch + x)) = *( source->buffer + ( y * srcPitch) + x);
        }      
      }

#endif    
    destHeight = srcHeight;
    }
  
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


FTBitmapGlyph::~FTBitmapGlyph()
{
  if( data)
    delete [] data;
}


float FTBitmapGlyph::Render( const FT_Vector& pen,
                             const int horizontal,
                             const FTGLRenderContext *context)
{
  if (!this->glyphHasBeenConverted)
    {
    this->ConvertGlyph(horizontal);
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
