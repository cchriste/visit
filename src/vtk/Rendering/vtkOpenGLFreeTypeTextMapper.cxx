/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOpenGLFreeTypeTextMapper.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOpenGLFreeTypeTextMapper.h"

#include "vtkActor2D.h"
#include "vtkObjectFactory.h"
#include "vtkProperty2D.h"
#include "vtkTextProperty.h"
#include "vtkViewport.h"
#include "vtkWindow.h"
#include "vtkToolkits.h"  // for VTK_USE_GL2PS

#include "vtkFreeTypeFontCache.h"
#include "vtkfreetypeConfig.h"
#include "vtkftglConfig.h"

#include "vtkgluPickMatrix.h"

#include "FTFont.h"

#ifdef VTK_USE_GL2PS
#include "gl2ps.h"
#endif // VTK_USE_GL2PS


//----------------------------------------------------------------------------
// Print debug info

#define VTK_FTTM_DEBUG 0
#define VTK_FTTM_DEBUG_CD 0

//----------------------------------------------------------------------------
// GL2PS related internal helper functions.

#ifdef VTK_USE_GL2PS
static void
vtkOpenGLFreeTypeTextMapper_GetGL2PSFontName(vtkTextProperty *tprop,
                                             char *ps_font)
{
 // For speed we use ARIAL == 0, COURIER == 1, TIMES == 2
  static char const *family[] = {"Helvetica", "Courier", "Times"};
  static char const *italic[] = {"Oblique", "Oblique", "Italic"};
  static char const *base[] = {"", "", "-Roman"};

  int font = tprop->GetFontFamily();

  if (font > 2)
    {
    sprintf(ps_font, "%s", tprop->GetFontFamilyAsString());
    if (tprop->GetBold())
      {
      sprintf(ps_font, "%s%s", ps_font, "Bold");
      }
    if (tprop->GetItalic())
      {
      sprintf(ps_font, "%s%s", ps_font, "Italic");
      }
      return;
    }

  if (tprop->GetBold())
    {
    sprintf(ps_font, "%s-%s", family[font], "Bold");
    if (tprop->GetItalic())
      {
      sprintf(ps_font, "%s%s", ps_font, italic[font]);
      }
    }
  else if (tprop->GetItalic())
    {
    sprintf(ps_font, "%s-%s", family[font], italic[font]);
    }
  else
    {
    sprintf(ps_font, "%s%s", family[font], base[font]);
    }
}
#endif

//----------------------------------------------------------------------------
#ifndef VTK_IMPLEMENT_MESA_CXX
vtkCxxRevisionMacro(vtkOpenGLFreeTypeTextMapper, "$Revision: 1.39 $");
vtkStandardNewMacro(vtkOpenGLFreeTypeTextMapper);
#endif

//----------------------------------------------------------------------------
vtkOpenGLFreeTypeTextMapper::vtkOpenGLFreeTypeTextMapper()
{
  this->LastSize[0] = 0;
  this->LastSize[1] = 0;
}

//----------------------------------------------------------------------------
vtkOpenGLFreeTypeTextMapper::~vtkOpenGLFreeTypeTextMapper()
{
  if (this->LastWindow)
    {
    this->ReleaseGraphicsResources(this->LastWindow);
    }  
}

//----------------------------------------------------------------------------
void vtkOpenGLFreeTypeTextMapper::ReleaseGraphicsResources(vtkWindow *vtkNotUsed(win))
{
#if VTK_FTTM_DEBUG
    printf("vtkOpenGLFreeTypeTextMapper::ReleaseGraphicsResources\n");
#endif
  
  this->LastWindow = NULL;
  
  // Very important
  // the release of graphics resources indicates that significant changes have
  // occurred. Old fonts, cached sizes etc are all no longer valid, so we send
  // ourselves a general modified message.

  // this->Modified();
}

//----------------------------------------------------------------------------
void vtkOpenGLFreeTypeTextMapper::GetSize(vtkViewport* viewport, int *size)
{
  // Check for multiline

  if (this->NumberOfLines > 1)
    {
    this->GetMultiLineSize(viewport, size);
    return;
    }

  // Check for input

  if (this->Input == NULL || this->Input[0] == '\0') 
    {
    size[0] = size[1] = 0;
    return;
    }

  vtkTextProperty *tprop = this->GetTextProperty();
  if (!tprop)
    {
    vtkErrorMacro(<< "Need a text property to get size");
    size[0] = size[1] = 0;
    return;
    }

  // Check to see whether we have to rebuild anything

  if (this->GetMTime() < this->SizeBuildTime &&
      tprop->GetMTime() < this->SizeBuildTime)
    {
#if VTK_FTTM_DEBUG
  printf("vtkOpenGLFreeTypeTextMapper::GetSize: In cache!\n");
#endif

    size[0] = this->LastSize[0];
    size[1] = this->LastSize[1];
    return;
    }

  // Check for font and try to set the size

  vtkFreeTypeFontCache::Entry *entry = vtkFreeTypeFontCache::GetInstance()->GetFont(tprop);
  FTFont *font = entry->Font;
  if (!font) 
    {
    vtkErrorMacro(<< "Render - No font");
    size[0] = size[1] = 0;
    return;
    }
  
  // The font global ascender and descender might just be too high
  // for given a face. Let's get a compromise by computing these values
  // from some usual ascii chars.
  
  if (entry->LargestAscender < 0 || entry->LargestDescender < 0)
    {
    float llx, lly, llz, urx, ury, urz;
    font->BBox("_/7Agfy", llx, lly, llz, urx, ury, urz);
    entry->LargestAscender = ury;
    entry->LargestDescender = lly;
    }
  
  this->LastSize[0] = size[0] = (int)font->Advance(this->Input);
  this->LastSize[1] = size[1] =
    (int)(entry->LargestAscender - entry->LargestDescender);
  this->LastLargestDescender = (int)entry->LargestDescender;

  this->SizeBuildTime.Modified();
}

//----------------------------------------------------------------------------
void vtkOpenGLFreeTypeTextMapper::RenderOverlay(vtkViewport* viewport, 
                                                vtkActor2D* actor)
{
  vtkDebugMacro (<< "RenderOverlay");

  // Check for input

  if (this->Input == NULL || this->Input[0] == '\0') 
    {
    return;
    }

  // Check for multi-lines

  if (this->NumberOfLines > 1)
    {
    this->RenderOverlayMultipleLines(viewport, actor);
    return;
    }

  // Get text property

  vtkTextProperty *tprop = this->GetTextProperty();
  if (!tprop)
    {
    vtkErrorMacro(<< "Need a text property to render mapper");
    return;
    }

  // Get the window information for display

  vtkWindow* window = viewport->GetVTKWindow();
  if (this->LastWindow && this->LastWindow != window)
    {
    this->ReleaseGraphicsResources(this->LastWindow);
    }
  this->LastWindow = window;

  // Get size of text

  int size[2];
  this->GetSize(viewport, size);

  // Get the position of the text actor

  int* actorPos;
  actorPos= 
    actor->GetActualPositionCoordinate()->GetComputedViewportValue(viewport);
  
  // Define bounding rectangle
  int orientVector[2];
  int upVector[2];
  if (tprop->GetOrientation() == VTK_TEXT_HORIZONTAL)
    {
    orientVector[0] = 1;
    orientVector[1] = 0;
    upVector[0] = 0;
    upVector[1] = 1;
    }
  else
    {
    orientVector[0] = 0;
    orientVector[1] = 1;
    upVector[0] = -1;
    upVector[1] =  0;
    }

  int pos[2];
  pos[0] = (int)(actorPos[0] - upVector[0] * tprop->GetLineOffset());
  pos[1] = (int)(actorPos[1] - upVector[1] * tprop->GetLineOffset());

  switch (tprop->GetJustification())
    {
    case VTK_TEXT_LEFT: 
      break;
    case VTK_TEXT_CENTERED:
      pos[0] = pos[0] - orientVector[0] * size[0] / 2;
      pos[1] = pos[1] - orientVector[1] * size[0] / 2;
      break;
    case VTK_TEXT_RIGHT: 
      pos[0] = pos[0] - orientVector[0] * size[0];
      pos[1] = pos[1] - orientVector[1] * size[0];
      break;
    }

  switch (tprop->GetVerticalJustification())
    {
    case VTK_TEXT_TOP: 
      pos[0] = pos[0] - upVector[0] * (size[1] - this->LastLargestDescender);
      pos[1] = pos[1] - upVector[1] * (size[1] - this->LastLargestDescender);
      break;
    case VTK_TEXT_CENTERED:
      pos[0] = pos[0] - upVector[0] - (size[1] / 2 + this->LastLargestDescender / 2);
      pos[1] = pos[1] - upVector[1] - (size[1] / 2 + this->LastLargestDescender / 2);
      break;
    case VTK_TEXT_BOTTOM: 
      break;
    }
  
  // Push a 2D matrix on the stack

  int *vsize = viewport->GetSize();
  double *vport = viewport->GetViewport();
  double *tileViewport = viewport->GetVTKWindow()->GetTileViewport();
  double visVP[4];

  visVP[0] = (vport[0] >= tileViewport[0]) ? vport[0] : tileViewport[0];
  visVP[1] = (vport[1] >= tileViewport[1]) ? vport[1] : tileViewport[1];
  visVP[2] = (vport[2] <= tileViewport[2]) ? vport[2] : tileViewport[2];
  visVP[3] = (vport[3] <= tileViewport[3]) ? vport[3] : tileViewport[3];

  if (visVP[0] == visVP[2] || visVP[1] == visVP[3])
    {
    return;
    }

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  if(viewport->GetIsPicking())
    {
    vtkgluPickMatrix(viewport->GetPickX(), viewport->GetPickY(),
                     1, 1, viewport->GetOrigin(), viewport->GetSize());
    }
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_LIGHTING);
  glDepthFunc(GL_ALWAYS);

  if (actor->GetProperty()->GetDisplayLocation() == VTK_FOREGROUND_LOCATION)
    {
    glOrtho(0, vsize[0] - 1, 0, vsize[1] - 1, 0, 1);
    }
  else
    {
    glOrtho(0, vsize[0] - 1, 0, vsize[1] - 1, -1, 0);
    }

  int *winSize = viewport->GetVTKWindow()->GetSize();

  int xoff = static_cast<int>(pos[0] - winSize[0] * (visVP[0] - vport[0]));
  int yoff = static_cast<int>(pos[1] - winSize[1] * (visVP[1] - vport[1]));
  
  // When picking draw the bounds of the text as a rectangle,
  // as text only picks when the pick point is exactly on the
  // origin of the text 

  if(viewport->GetIsPicking())
    {
    float x1 = 2.0 * (float)actorPos[0] / vsize[0] - 1;
    float y1 = 2.0 * ((float)actorPos[1] - tprop->GetLineOffset())/vsize[1] - 1;
    float width = 2.0 * (float)size[0] / vsize[0];
    float height = 2.0 * (float)size[1] / vsize[1];
    glRectf(x1, y1, x1 + width, y1 + height);

    // Clean up and return after drawing the rectangle

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glEnable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    
    return;
    }

  // Get the font color from the text actor

  unsigned char red, green, blue, alpha;
  
  // TOFIX: the default text prop color is set to a special (-1, -1, -1) value
  // to maintain backward compatibility for a while. Text mapper classes will
  // use the Actor2D color instead of the text prop color if this value is 
  // found (i.e. if the text prop color has not been set).
  double* tpropColor = tprop->GetColor();
  if (tpropColor[0] < 0.0 && tpropColor[1] < 0.0 && tpropColor[2] < 0.0)
    {
    tpropColor = actor->GetProperty()->GetColor();
    }

  // TOFIX: same goes for opacity

  float opacity = tprop->GetOpacity();
  if (opacity < 0.0)
    {
    opacity = actor->GetProperty()->GetOpacity();
    }

  red   = (unsigned char) (tpropColor[0] * 255.0);
  green = (unsigned char) (tpropColor[1] * 255.0);
  blue  = (unsigned char) (tpropColor[2] * 255.0);
  alpha = (unsigned char) (opacity       * 255.0);

  // Get the font
  
  FTFont *font = 
    vtkFreeTypeFontCache::GetInstance()->GetFont(tprop, 1, red,green,blue)->Font;
  if (!font) 
    {
    vtkErrorMacro(<< "Render - No font");
    return;
    }

  struct FTGLRenderContext *ftgl_context = 0;

#ifdef VTK_IMPLEMENT_MESA_CXX
  // If we support Mangle Mesa, VTK_IMPLEMENT_MESA_CXX will be defined to
  // compile this unit as a Mesa text mapper. In that case, provide a
  // context to FTGL to switch dynamically to Mangle Mesa rendering.
  struct FTGLRenderContext ftgl_context_mesa;
  ftgl_context_mesa.UseMangleMesa = 1;
  ftgl_context = &ftgl_context_mesa;
#endif

  // Setup the fonts for GL2PS output.

#ifdef VTK_USE_GL2PS
  char ps_font[64];
  vtkOpenGLFreeTypeTextMapper_GetGL2PSFontName(tprop, ps_font);
#endif // VTK_USE_GL2PS

  // Set up the shadow color

  int antialiasing_requested = 
    (tprop->GetGlobalAntiAliasing() == VTK_TEXT_GLOBAL_ANTIALIASING_ALL || 
     (tprop->GetGlobalAntiAliasing() == VTK_TEXT_GLOBAL_ANTIALIASING_SOME 
      && tprop->GetAntiAliasing())) ? 1 : 0;

  if (tprop->GetShadow())
    {
    unsigned char rgb = (red + green + blue) / 3.0 > 128.0 ? 0 : 255;
    unsigned char shadow_red = rgb, shadow_green = rgb, shadow_blue = rgb; 

    // Get the shadow font
  
#if VTK_FTFC_CACHE_BY_RGBA
    FTFont *shadow_font;
    if (antialiasing_requested)
      {
      shadow_font = vtkFreeTypeFontCache::GetInstance()->GetFont(
        tprop, 1, shadow_red, shadow_green, shadow_blue)->Font;
      if (!shadow_font) 
        {
        vtkErrorMacro(<< "Render - No shadow font");
        return;
        }
      } 
    else 
      {
      shadow_font = font;
      }
#endif
    
    // Set the color here since load/render glyphs is done
    // on demand and this color has to be consistent for a given font entry.
    
    glColor4ub(shadow_red, shadow_green, shadow_blue, alpha);

    // Required for clipping to work correctly

    glRasterPos2i(0, 0);
    glBitmap(0, 0, 0, 0, xoff + 1, yoff - 1, NULL);
    
    // Draw the shadow text
    
#if VTK_FTFC_CACHE_BY_RGBA
    shadow_font->render(this->Input, ftgl_context);

    // Get the font again, Duh, since it may have been freed from the 
    // cache by the shadow font

    if (antialiasing_requested)
      {
      font = 
        vtkFreeTypeFontCache::GetInstance()->GetFont(tprop, 1, red, green, blue)->Font;
      if (!font) 
        {
        vtkErrorMacro(<< "Render - No font");
        return;
        }
      }
#else
    font->render(this->Input, ftgl_context);
#endif

    // Shadow text for GL2PS.

#ifdef VTK_USE_GL2PS
    gl2psText(this->Input, ps_font, tprop->GetFontSize());
#endif // VTK_USE_GL2PS
    }

  // Set the color here since load/render glyphs is done
  // on demand and this color has to be consistent for a given font entry.

  glColor4ub(red, green, blue, alpha);

  // Required for clipping to work correctly

  glRasterPos2i(0, 0);
  glBitmap(0, 0, 0, 0, xoff, yoff, NULL);

  // Display a string

  font->render(this->Input, ftgl_context);

  glFlush();

  // Normal text for GL2PS.

#ifdef VTK_USE_GL2PS
  gl2psText(this->Input, ps_font, tprop->GetFontSize());
#endif // VTK_USE_GL2PS

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_LIGHTING);
  glDepthFunc(GL_LEQUAL);
}

void vtkOpenGLFreeTypeTextMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
