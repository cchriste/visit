/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFreeTypeUtilities.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFreeTypeUtilities.h"

#include "vtkTextProperty.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"
#include "vtkImageData.h"
//#include "vtkDebugLeaks.h"

// FTGL

#include "vtkftglConfig.h"
#include "FTLibrary.h"
#include "FTGLPixmapFont.h"

// The embedded fonts

#include "fonts/vtkEmbeddedFonts.h"

#include <sys/stat.h>

// Print debug info

#define VTK_FTFC_DEBUG 0
#define VTK_FTFC_DEBUG_CD 0

//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkFreeTypeUtilities, "$Revision: 1.12 $");
vtkInstantiatorNewMacro(vtkFreeTypeUtilities);

//----------------------------------------------------------------------------
// The singleton, and the singleton cleanup

vtkFreeTypeUtilities* vtkFreeTypeUtilities::Instance = NULL;
vtkFreeTypeUtilitiesCleanup vtkFreeTypeUtilities::Cleanup;

//----------------------------------------------------------------------------
// The embedded fonts
// Create a lookup table between the text mapper attributes 
// and the font buffers.

struct EmbeddedFontStruct
{
  size_t length;
  unsigned char *ptr;
};

//----------------------------------------------------------------------------
// This callback will be called by the FTGLibrary singleton cleanup destructor
// if it happens to be destroyed before our singleton (this order is not 
// deterministic). It will destroy our singleton, if needed.

void vtkFreeTypeUtilitiesCleanupCallback ()
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilitiesCleanupCallback\n");
#endif
  vtkFreeTypeUtilities::SetInstance(NULL);
}

//----------------------------------------------------------------------------
// Create the singleton cleanup
// Register our singleton cleanup callback against the FTLibrary so that
// it might be called before the FTLibrary singleton is destroyed.

vtkFreeTypeUtilitiesCleanup::vtkFreeTypeUtilitiesCleanup()
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilitiesCleanup::vtkFreeTypeUtilitiesCleanup\n");
#endif
  FTLibraryCleanup::AddDependency(&vtkFreeTypeUtilitiesCleanupCallback);
}

//----------------------------------------------------------------------------
// Delete the singleton cleanup 
// The callback called here might have been called by the FTLibrary singleton
// cleanup first (depending on the destruction order), but in case ours is
// destroyed first, let's call it too.

vtkFreeTypeUtilitiesCleanup::~vtkFreeTypeUtilitiesCleanup()
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilitiesCleanup::~vtkFreeTypeUtilitiesCleanup\n");
#endif
  vtkFreeTypeUtilitiesCleanupCallback();
}

//----------------------------------------------------------------------------
vtkFreeTypeUtilities* vtkFreeTypeUtilities::GetInstance()
{
  if (!vtkFreeTypeUtilities::Instance)
    {
    vtkFreeTypeUtilities::Instance = (vtkFreeTypeUtilities*)
      vtkObjectFactory::CreateInstance("vtkFreeTypeUtilities");
    if (!vtkFreeTypeUtilities::Instance)
      {
#ifdef VTK_DEBUG_LEAKS
      vtkDebugLeaks::DestructClass("vtkFreeTypeUtilities");
#endif
      vtkFreeTypeUtilities::Instance = new vtkFreeTypeUtilities;
      }
    }
  return vtkFreeTypeUtilities::Instance;
}

//----------------------------------------------------------------------------
void vtkFreeTypeUtilities::SetInstance(vtkFreeTypeUtilities* instance)
{
  if (vtkFreeTypeUtilities::Instance == instance)
    {
    return;
    }

  if (vtkFreeTypeUtilities::Instance)
    {
    vtkFreeTypeUtilities::Instance->Delete();
    }

  vtkFreeTypeUtilities::Instance = instance;

  // User will call ->Delete() after setting instance

  if (instance)
    {
    instance->Register(NULL);
    }
}

//----------------------------------------------------------------------------
vtkFreeTypeUtilities* vtkFreeTypeUtilities::New()
{
  vtkFreeTypeUtilities* ret = vtkFreeTypeUtilities::GetInstance();
  ret->Register(NULL);
  return ret;
}

//----------------------------------------------------------------------------
vtkFreeTypeUtilities::vtkFreeTypeUtilities() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::vtkFreeTypeUtilities\n");
#endif

  this->MaximumNumberOfFaces = 30; // combinations of family+bold+italic
  this->MaximumNumberOfSizes = this->MaximumNumberOfFaces * 20; // sizes
  this->MaximumNumberOfBytes = 300000UL * this->MaximumNumberOfSizes;

  this->CacheManager = NULL;
  this->ImageCache   = NULL;
  this->CMapCache    = NULL;

  this->NumberOfEntries = 0;
  this->InitializeCache();
}

//----------------------------------------------------------------------------
vtkFreeTypeUtilities::~vtkFreeTypeUtilities() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::~vtkFreeTypeUtilities\n");
#endif

  this->ReleaseCache();
  this->ReleaseCacheManager();
}

//----------------------------------------------------------------------------
FT_Library* vtkFreeTypeUtilities::GetLibrary() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::GetLibrary\n");
#endif

  FTLibrary *ftgl_lib = FTLibrary::GetInstance();
  if (ftgl_lib)
    {
    return ftgl_lib->GetLibrary();
    }

  return NULL;
}

//----------------------------------------------------------------------------
FTC_Manager* vtkFreeTypeUtilities::GetCacheManager() 
{
  if (!this->CacheManager)
    {
    this->InitializeCacheManager();
    }

  return this->CacheManager;
}

//----------------------------------------------------------------------------
FTC_ImageCache* vtkFreeTypeUtilities::GetImageCache() 
{
  if (!this->ImageCache)
    {
    this->InitializeCacheManager();
    }

  return this->ImageCache;
}

//----------------------------------------------------------------------------
FTC_CMapCache* vtkFreeTypeUtilities::GetCMapCache() 
{
  if (!this->CMapCache)
    {
    this->InitializeCacheManager();
    }

  return this->CMapCache;
}

//----------------------------------------------------------------------------
void vtkFreeTypeUtilities::MapTextPropertyToId(vtkTextProperty *tprop, 
                                               unsigned long *id) 
{
  if (!tprop || !id)
    {
    vtkErrorMacro(<< "Wrong parameters, one of them is NULL");
    return;
    }

  // Set the first bit to avoid id = 0 
  // (the id will be mapped to a pointer, FTC_FaceID, so let's avoid NULL)

  *id = 1;
  int bits = 1;

  // The font family is in 4 bits (= 5 bits so far)
  // (2 would be enough right now, but who knows, it might grow)

  int fam = (tprop->GetFontFamily() - tprop->GetFontFamilyMinValue()) << bits;
  bits += 4;

  // Bold is in 1 bit (= 6 bits so far)

  int bold = (tprop->GetBold() ? 1 : 0) << bits;
  bits++;

  // Italic is in 1 bit (= 7 bits so far)

  int italic = (tprop->GetItalic() ? 1 : 0) << bits;
  bits++;

  // Orientation (in degrees)
  // We need 9 bits for 0 to 360. What do we need for more precisions:
  // - 1/10th degree: 12 bits (11.8)

  int angle = (vtkMath::Round(tprop->GetOrientation() * 10.0) % 3600) << bits;
  
  // We really should not use more than 32 bits

  // Now final id

  *id |= fam | bold | italic | angle;
}

//----------------------------------------------------------------------------
void vtkFreeTypeUtilities::MapIdToTextProperty(unsigned long id,
                                               vtkTextProperty *tprop) 
{
  if (!tprop)
    {
    vtkErrorMacro(<< "Wrong parameters, one of them is NULL");
    return;
    }

  // The first was set to avoid id = 0

  int bits = 1;

  // The font family is in 4 bits
  // (2 would be enough right now, but who knows, it might grow)

  int fam = id >> bits;
  bits += 4;
  tprop->SetFontFamily((fam & ((1 << 4) - 1))+ tprop->GetFontFamilyMinValue());

  // Bold is in 1 bit

  int bold = id >> bits;
  bits++;
  tprop->SetBold(bold & 0x1);

  // Italic is in 1 bit

  int italic = id >> bits;
  bits++;
  tprop->SetItalic(italic & 0x1);

  // Orientation (in degrees)
  // We need 9 bits for 0 to 360. What do we need for more precisions:
  // - 1/10th degree: 12 bits (11.8)

  int angle = id >> bits;
  tprop->SetOrientation((float)(angle & ((1 << 12) - 1)) / 10.0);

  // We really should not use more than 32 bits
}

//----------------------------------------------------------------------------
FT_CALLBACK_DEF(FT_Error)
vtkFreeTypeUtilitiesFaceRequester(FTC_FaceID face_id,
                                  FT_Library lib,
                                  FT_Pointer request_data,
                                  FT_Face* face)
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilitiesFaceRequester()\n");
#endif

  FT_UNUSED(request_data);

  vtkFreeTypeUtilities *self = 
    reinterpret_cast<vtkFreeTypeUtilities*>(request_data);

  // Map the ID to a text property

  vtkTextProperty *tprop = vtkTextProperty::New();
  self->MapIdToTextProperty(reinterpret_cast<unsigned long>(face_id), tprop);

  // Fonts, organized by [Family][Bold][Italic]
    
  static EmbeddedFontStruct EmbeddedFonts[3][2][2] = 
    {
      {
        {
          { // VTK_ARIAL: Bold [ ] Italic [ ]
            
            face_arial_buffer_length, face_arial_buffer
          },
          { // VTK_ARIAL: Bold [ ] Italic [x]
            face_arial_italic_buffer_length, face_arial_italic_buffer
          }
        },
        {
          { // VTK_ARIAL: Bold [x] Italic [ ]
            face_arial_bold_buffer_length, face_arial_bold_buffer
          },
          { // VTK_ARIAL: Bold [x] Italic [x]
            face_arial_bold_italic_buffer_length, face_arial_bold_italic_buffer
          }
        }
      },
      {
        {
          { // VTK_COURIER: Bold [ ] Italic [ ]
            face_courier_buffer_length, face_courier_buffer
          },
          { // VTK_COURIER: Bold [ ] Italic [x]
            face_courier_italic_buffer_length, face_courier_italic_buffer
          }
        },
        {
          { // VTK_COURIER: Bold [x] Italic [ ]
            face_courier_bold_buffer_length, face_courier_bold_buffer
          },
          { // VTK_COURIER: Bold [x] Italic [x]
            face_courier_bold_italic_buffer_length, 
            face_courier_bold_italic_buffer
          }
        }
      },
      {
        {
          { // VTK_TIMES: Bold [ ] Italic [ ]
            face_times_buffer_length, face_times_buffer
          },
          { // VTK_TIMES: Bold [ ] Italic [x]
            face_times_italic_buffer_length, face_times_italic_buffer
          }
        },
        {
          { // VTK_TIMES: Bold [x] Italic [ ]
            face_times_bold_buffer_length, face_times_bold_buffer
          },
          { // VTK_TIMES: Bold [x] Italic [x]
            face_times_bold_italic_buffer_length, face_times_bold_italic_buffer
          }
        }
      }
    };
    
  FT_Long length = EmbeddedFonts
    [tprop->GetFontFamily()][tprop->GetBold()][tprop->GetItalic()].length;
  FT_Byte *ptr = EmbeddedFonts
    [tprop->GetFontFamily()][tprop->GetBold()][tprop->GetItalic()].ptr;

  // Create a new face

  FT_Error error = FT_New_Memory_Face(lib, ptr, length, 0, face);
  if (error)
    {
    vtkErrorWithObjectMacro(
      tprop,
      << "Unable to create font !" << " (family: " << tprop->GetFontFamily() 
      << ", bold: " << tprop->GetBold() << ", italic: " << tprop->GetItalic()
      << ", length: " << length << ")");
    }
  else
    {
#if VTK_FTFC_DEBUG
    cout << "Requested: " << *face
         << " (F: " << tprop->GetFontFamily() 
         << ", B: " << tprop->GetBold() 
         << ", I: " << tprop->GetItalic()
         << ", O: " << tprop->GetOrientation() << ")" << endl;
#endif    
    if (tprop->GetOrientation() != 0.0)
      {
      // FreeType documentation says that the transform should not be set
      // but we cache faces also by transform, so that there is a unique
      // (face, orientation) cache entry

      FT_Matrix matrix;
      float angle = vtkMath::DegreesToRadians() * tprop->GetOrientation(); 
      matrix.xx = (FT_Fixed)( cos(angle) * 0x10000L);
      matrix.xy = (FT_Fixed)(-sin(angle) * 0x10000L);
      matrix.yx = (FT_Fixed)( sin(angle) * 0x10000L);
      matrix.yy = (FT_Fixed)( cos(angle) * 0x10000L);
      FT_Set_Transform(*face, &matrix, NULL);
      }
    }

  tprop->Delete();

  return error;
}

//----------------------------------------------------------------------------
void vtkFreeTypeUtilities::InitializeCacheManager() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::InitializeCacheManager()\n");
#endif

  this->ReleaseCacheManager();

  FT_Error error;

  // Create the cache manager itself

  this->CacheManager = new FTC_Manager;

  error = FTC_Manager_New(*this->GetLibrary(), 
                          this->MaximumNumberOfFaces,
                          this->MaximumNumberOfSizes,
                          this->MaximumNumberOfBytes,
                          vtkFreeTypeUtilitiesFaceRequester,
                          (FT_Pointer)this, 
                          this->CacheManager);

  if (error)
    {
    vtkErrorMacro(<< "Failed allocating a new FreeType Cache Manager");
    }

  // The image cache

  this->ImageCache = new FTC_ImageCache;

  error = FTC_ImageCache_New(*this->CacheManager, this->ImageCache);

  if (error)
    {
    vtkErrorMacro(<< "Failed allocating a new FreeType Image Cache");
    }

  // The charmap cache

  this->CMapCache = new FTC_CMapCache;

  error = FTC_CMapCache_New(*this->CacheManager, this->CMapCache);

  if (error)
    {
    vtkErrorMacro(<< "Failed allocating a new FreeType CMap Cache");
    }
}

//----------------------------------------------------------------------------
void vtkFreeTypeUtilities::ReleaseCacheManager() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::ReleaseCacheManager()\n");
#endif

  if (this->CacheManager)
    {
    FTC_Manager_Done(*this->CacheManager);

    delete this->CacheManager;
    this->CacheManager = NULL;
    }

  if (this->ImageCache)
    {
    delete this->ImageCache;
    this->ImageCache = NULL;
    }

  if (this->CMapCache)
    {
    delete this->CMapCache;
    this->CMapCache = NULL;
    }
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetSize(unsigned long tprop_cache_id,
                                  int font_size,
                                  FT_Size *size)
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::GetSize()\n");
#endif

  if (!size || font_size <= 0)
    {
    vtkErrorMacro(<< "Wrong parameters, size is NULL or invalid font size");
    return 0;
    }

  FTC_Manager *manager = this->GetCacheManager();
  if (!manager)
    {
    vtkErrorMacro(<< "Failed querying the cache manager !");
    return 0;
    }

  // Map the id of a text property in the cache to a FTC_FaceID

#if (FREETYPE_MAJOR >=2 && FREETYPE_MINOR >= 1 && FREETYPE_PATCH >= 9)
  FTC_FaceID face_id = reinterpret_cast<FTC_FaceID>(tprop_cache_id);

  FTC_ScalerRec scaler_rec;
  scaler_rec.face_id = face_id;
  scaler_rec.width = font_size;
  scaler_rec.height = font_size;
  scaler_rec.pixel = 1;

  FT_Error error = FTC_Manager_LookupSize(*manager, &scaler_rec, size);
  if (error)
    {
    vtkErrorMacro(<< "Failed looking up a FreeType Size");
    }

  return error ? 0 : 1;
#else
  return 0;
#endif
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetSize(vtkTextProperty *tprop, 
                                  FT_Size *size)
{
  if (!tprop)
    {
    vtkErrorMacro(<< "Wrong parameters, text property is NULL");
    return 0;
    }

  // Map the text property to a unique id that will be used as face id

  unsigned long tprop_cache_id;
  this->MapTextPropertyToId(tprop, &tprop_cache_id);

  return this->GetSize(tprop_cache_id, tprop->GetFontSize(), size);
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetFace(unsigned long tprop_cache_id,
                                  FT_Face *face)
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::GetFace()\n");
#endif

  if (!face)
    {
    vtkErrorMacro(<< "Wrong parameters, face is NULL");
    return 0;
    }

  FTC_Manager *manager = this->GetCacheManager();
  if (!manager)
    {
    vtkErrorMacro(<< "Failed querying the cache manager !");
    return 0;
    }

  // Map the id of a text property in the cache to a FTC_FaceID

  FTC_FaceID face_id = reinterpret_cast<FTC_FaceID>(tprop_cache_id);

  FT_Error error = FTC_Manager_LookupFace(*manager, face_id, face);
  if (error)
    {
    vtkErrorMacro(<< "Failed looking up a FreeType Face");
    }

  return error ? 0 : 1;
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetFace(vtkTextProperty *tprop,
                                  FT_Face *face)
{
  if (!tprop)
    {
    vtkErrorMacro(<< "Wrong parameters, face is NULL");
    return 0;
    }

  // Map the text property to a unique id that will be used as face id

  unsigned long tprop_cache_id;
  this->MapTextPropertyToId(tprop, &tprop_cache_id);

  return this->GetFace(tprop_cache_id, face);
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetGlyphIndex(unsigned long tprop_cache_id,
                                        char c, 
                                        FT_UInt *gindex)
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::GetGlyphIndex()\n");
#endif

  if (!gindex)
    {
    vtkErrorMacro(<< "Wrong parameters, gindex is NULL");
    return 0;
    }

  FTC_CMapCache *cmap_cache = this->GetCMapCache();
  if (!cmap_cache)
    {
    vtkErrorMacro(<< "Failed querying the charmap cache manager !");
    return 0;
    }

  // Map the id of a text property in the cache to a FTC_FaceID

  FTC_FaceID face_id = reinterpret_cast<FTC_FaceID>(tprop_cache_id);

  // Lookup the glyph index

#if (FREETYPE_MAJOR >=2 && FREETYPE_MINOR >= 1 && FREETYPE_PATCH >= 9)
  *gindex = FTC_CMapCache_Lookup(*cmap_cache, face_id, 0, c);

  return *gindex ? 1 : 0;
#else
  return 0;
#endif
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetGlyphIndex(vtkTextProperty *tprop,
                                        char c, 
                                        FT_UInt *gindex)
{
  if (!tprop)
    {
    vtkErrorMacro(<< "Wrong parameters, text property is NULL");
    return 0;
    }

  // Map the text property to a unique id that will be used as face id

  unsigned long tprop_cache_id;
  this->MapTextPropertyToId(tprop, &tprop_cache_id);

  return this->GetGlyphIndex(tprop_cache_id, c, gindex);
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetGlyph(unsigned long tprop_cache_id, 
                                   int font_size,
                                   FT_UInt gindex,
                                   FT_Glyph *glyph,
                                   int request)
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::GetGlyph()\n");
#endif

  if (!glyph)
    {
    vtkErrorMacro(<< "Wrong parameters, one of them is NULL");
    return 0;
    }

  FTC_ImageCache *image_cache = this->GetImageCache();
  if (!image_cache)
    {
    vtkErrorMacro(<< "Failed querying the image cache manager !");
    return 0;
    }

  // Map the id of a text property in the cache to a FTC_FaceID

  FTC_FaceID face_id = reinterpret_cast<FTC_FaceID>(tprop_cache_id);

  // Which font are we looking for

#if (FREETYPE_MAJOR >=2 && FREETYPE_MINOR >= 1 && FREETYPE_PATCH >= 9)
  FTC_ImageTypeRec image_type_rec;
  image_type_rec.face_id = face_id;
  image_type_rec.width = font_size;
  image_type_rec.height = font_size;
  image_type_rec.flags = FT_LOAD_DEFAULT;
  if (request == GLYPH_REQUEST_BITMAP)
    {
    image_type_rec.flags |= FT_LOAD_RENDER;
    }
  else if (request == GLYPH_REQUEST_OUTLINE)
    {
    image_type_rec.flags |= FT_LOAD_NO_BITMAP;
    }

  // Lookup the glyph

  FT_Error error = FTC_ImageCache_Lookup(
    *image_cache, &image_type_rec, gindex, glyph, NULL);

  return error ? 0 : 1;
#else
  return 0;
#endif
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetGlyph(vtkTextProperty *tprop,
                                   char c, 
                                   FT_Glyph *glyph,
                                   int request)
{
  if (!tprop)
    {
    vtkErrorMacro(<< "Wrong parameters, text property is NULL");
    return 0;
    }

  // Map the text property to a unique id that will be used as face id

  unsigned long tprop_cache_id;
  this->MapTextPropertyToId(tprop, &tprop_cache_id);

  // Get the character/glyph index

  FT_UInt gindex;
  if (!this->GetGlyphIndex(tprop_cache_id, c, &gindex))
    {
    vtkErrorMacro(<< "Failed querying a glyph index");
    return 0;
    }
  
  // Get the glyph

  return this->GetGlyph(
    tprop_cache_id, tprop->GetFontSize(), gindex, glyph, request);
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::IsBoundingBoxValid(int bbox[4])
{
  return (!bbox ||
          bbox[0] == VTK_INT_MAX || bbox[1] == VTK_INT_MIN ||
          bbox[2] == VTK_INT_MAX || bbox[3] == VTK_INT_MIN) ? 0 : 1;
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::GetBoundingBox(vtkTextProperty *tprop, 
                                         const char *str,
                                         int bbox[4])
{
  // We need the tprop and bbox

  if (!tprop || !bbox)
    {
    vtkErrorMacro(<< "Wrong parameters, one of them is NULL or zero");
    return 0;
    }
  
  // Initialize bbox to some large values

  bbox[0] = bbox[2] = VTK_INT_MAX;
  bbox[1] = bbox[3] = VTK_INT_MIN;

  // No string to render, bail out now

  if (!str)
    {
    return 1;
    }

  // Map the text property to a unique id that will be used as face id

  unsigned long tprop_cache_id;
  this->MapTextPropertyToId(tprop, &tprop_cache_id);

  // Get the face

  FT_Face face;
  if (!this->GetFace(tprop_cache_id, &face))
    {
    vtkErrorMacro(<< "Failed retrieving the face");
    return 0;
    }

  int face_has_kerning = FT_HAS_KERNING(face);

  // Iterate char by char

  FT_Glyph glyph;
  FT_BitmapGlyph bitmap_glyph;
  FT_Bitmap *bitmap;
  FT_UInt gindex, previous_gindex = 0;
  FT_Vector kerning_delta;

  int x = 0, y = 0;

  for (; *str; str++)
    {
    // Get the glyph index

    if (!this->GetGlyphIndex(tprop_cache_id, *str, &gindex))
      {
      continue;
      }

    // Get the glyph as a bitmap

    if (!this->GetGlyph(tprop_cache_id, 
                        tprop->GetFontSize(), 
                        gindex, 
                        &glyph, 
                        vtkFreeTypeUtilities::GLYPH_REQUEST_BITMAP) ||
        glyph->format != FT_GLYPH_FORMAT_BITMAP)
      {
      continue;
      }

    bitmap_glyph = reinterpret_cast<FT_BitmapGlyph>(glyph);
    bitmap = &bitmap_glyph->bitmap;

    if (bitmap->width && bitmap->rows)
      {
      // Starting position given the bearings

      // Substract 1 to the bearing Y, because this is the vertical distance
      // from the glyph origin (0,0) to the topmost pixel of the glyph bitmap
      // (more precisely, to the pixel just above the bitmap). This distance is
      // expressed in integer pixels, and is positive for upwards y.

      int pen_x = x + bitmap_glyph->left;
      int pen_y = y + bitmap_glyph->top - 1;

      // Add the kerning

      if (face_has_kerning && previous_gindex && gindex)
        {
        FT_Get_Kerning(
          face, previous_gindex, gindex, FT_KERNING_DEFAULT, &kerning_delta);
        pen_x += kerning_delta.x >> 6;
        pen_y += kerning_delta.y >> 6;
        }

      previous_gindex = gindex;

      // Update bounding box

      if (pen_x < bbox[0])
        {
        bbox[0] = pen_x;
        }
      if (pen_x > bbox[1])
        {
        bbox[1] = pen_x;
        }
      if (pen_y < bbox[2])
        {
        bbox[2] = pen_y;
        }
      if (pen_y > bbox[3])
        {
        bbox[3] = pen_y;
        }

      pen_x += (bitmap->width - 1);
      pen_y -= (bitmap->rows - 1);

      if (pen_x < bbox[0])
        {
        bbox[0] = pen_x;
        }
      if (pen_x > bbox[1])
        {
        bbox[1] = pen_x;
        }
      if (pen_y < bbox[2])
        {
        bbox[2] = pen_y;
        }
      if (pen_y > bbox[3])
        {
        bbox[3] = pen_y;
        }
      }

    // Advance to next char

    x += (bitmap_glyph->root.advance.x + 0x8000) >> 16;
    y += (bitmap_glyph->root.advance.y + 0x8000) >> 16;
    }

  // Margin for shadow (x + 1, y - 1)

  if (tprop->GetShadow() && this->IsBoundingBoxValid(bbox))
    {
    bbox[1]++;
    bbox[2]--;
    }

  return 1;
}

//----------------------------------------------------------------------------
template <class T>
int vtkFreeTypeUtilitiesRenderString(
  vtkFreeTypeUtilities *self, 
  vtkTextProperty *tprop, 
  const char *str,
  int x, int y,
  vtkImageData *data,
  T *vtkNotUsed(ptr),
  int use_shadow_color)
{
  // Map the text property to a unique id that will be used as face id

  unsigned long tprop_cache_id;
  self->MapTextPropertyToId(tprop, &tprop_cache_id);

  // Get the face

  FT_Face face;
  if (!self->GetFace(tprop_cache_id, &face))
    {
    vtkErrorWithObjectMacro(self, << "Failed retrieving the face");
    return 0;
    }

  int face_has_kerning = FT_HAS_KERNING(face);

  // Text property size and opacity

  int tprop_font_size = tprop->GetFontSize();

  float tprop_opacity = tprop->GetOpacity();

  // Text color (get the shadow color if we are actually drawing the shadow)
  // Also compute the luminance, if we are drawing to a grayscale image

  double color[3];
  if (use_shadow_color)
    {
    tprop->GetShadowColor(color);
    }
  else
    {
    tprop->GetColor(color);
    }
  float tprop_r = color[0];
  float tprop_g = color[1];
  float tprop_b = color[2];

  float tprop_l = 0.3 * tprop_r + 0.59 * tprop_g + 0.11 * tprop_b;

  // Image params (comps, increments, range)

  int data_nb_comp = data->GetNumberOfScalarComponents();

  vtkIdType data_inc_x, data_inc_y, data_inc_z;
  data->GetIncrements(data_inc_x, data_inc_y, data_inc_z);

  double data_min, data_max;
  if (data->GetScalarType() == VTK_DOUBLE ||
      data->GetScalarType() == VTK_FLOAT)
    {
    data_min = 0.0;
    data_max = 1.0;
    }
  else
    {
    data_min = data->GetScalarTypeMin();
    data_max = data->GetScalarTypeMax();
    }
  double data_range = (data_max - data_min);

  // Render char by char

  FT_Glyph glyph;
  FT_BitmapGlyph bitmap_glyph;
  FT_Bitmap *bitmap;
  FT_UInt gindex, previous_gindex = 0;
  FT_Vector kerning_delta;

  for (; *str; str++)
    {
    // Get the glyph index

    if (!self->GetGlyphIndex(tprop_cache_id, *str, &gindex))
      {
      continue;
      }

    // Get the glyph as a bitmap

    if (!self->GetGlyph(tprop_cache_id, 
                        tprop_font_size, 
                        gindex, 
                        &glyph, 
                        vtkFreeTypeUtilities::GLYPH_REQUEST_BITMAP) ||
        glyph->format != FT_GLYPH_FORMAT_BITMAP)
      {
      continue;
      }

    bitmap_glyph = reinterpret_cast<FT_BitmapGlyph>(glyph);
    bitmap = &bitmap_glyph->bitmap;

    if (bitmap->pixel_mode != FT_PIXEL_MODE_GRAY)
      {
      continue;
      }

    if (bitmap->width && bitmap->rows)
      {
      // Starting position given the bearings

      // Substract 1 to the bearing Y, because this is the vertical distance
      // from the glyph origin (0,0) to the topmost pixel of the glyph bitmap
      // (more precisely, to the pixel just above the bitmap). This distance is
      // expressed in integer pixels, and is positive for upwards y.

#if VTK_FTFC_DEBUG
      cout << *str << ", orient: " << tprop->GetOrientation() << ", x: " << x << ", y: " << y << ", left: " << bitmap_glyph->left << ", top: " << bitmap_glyph->top << ", width: " << bitmap->width << ", rows: " << bitmap->rows << endl;
#endif

      int pen_x = x + bitmap_glyph->left;
      int pen_y = y + bitmap_glyph->top - 1;

      // Add the kerning

      if (face_has_kerning && previous_gindex && gindex)
        {
        FT_Get_Kerning(
          face, previous_gindex, gindex, FT_KERNING_DEFAULT, &kerning_delta);
        pen_x += kerning_delta.x >> 6;
        pen_y += kerning_delta.y >> 6;
        }

      previous_gindex = gindex;

      // Render

      T *data_ptr = (T*)data->GetScalarPointer(pen_x, pen_y, 0);

      int data_pitch = (-data->GetDimensions()[0] - bitmap->width) * data_inc_x;

      unsigned char *glyph_ptr_row = bitmap->buffer;
      unsigned char *glyph_ptr;

      float t_alpha, t_1_m_alpha, data_alpha;

      int i, j;
      for (j = 0; j < bitmap->rows; j++)
        {
        glyph_ptr = glyph_ptr_row;

        // That loop should probably be located inside the switch for efficency

        for (i = 0; i < bitmap->width; i++)
          {
          t_alpha = tprop_opacity * (*glyph_ptr / 255.0);
          t_1_m_alpha = 1.0 - t_alpha;

          switch (data_nb_comp)
            {
            // L

            case 1:
              *data_ptr = (T)(
                (data_min + data_range * tprop_l * t_alpha) + 
                *data_ptr * t_1_m_alpha);
              glyph_ptr++;
              data_ptr++;
              break;
            
              // L,A
              // TOFIX: that code is so wrong (alpha)  

            case 2:
              data_alpha = (data_ptr[1] - data_min) / data_range;
              *data_ptr = (T)(
                (data_min + data_range * tprop_l * t_alpha) + 
                (*data_ptr * data_alpha) * t_1_m_alpha);
              data_ptr++;
              *data_ptr = (T)(
                data_min + data_range * (t_alpha + data_alpha * t_1_m_alpha));
              data_ptr++;
              glyph_ptr++;
              break;

              // RGB

            case 3:
              *data_ptr = (T)(
                (data_min + data_range * tprop_r * t_alpha) + 
                *data_ptr * t_1_m_alpha);
              data_ptr++;
              *data_ptr = (T)(
                (data_min + data_range * tprop_g * t_alpha) + 
                *data_ptr * t_1_m_alpha);
              data_ptr++;
              *data_ptr = (T)(
                (data_min + data_range * tprop_b * t_alpha) + 
                *data_ptr * t_1_m_alpha);
              data_ptr++;
              glyph_ptr++;
              break;

              // RGB,A
              // TOFIX: that code is so wrong (alpha)  

            case 4:
              data_alpha = (data_ptr[1] - data_min) / data_range;
              *data_ptr = (T)(
                (data_min + data_range * tprop_r * t_alpha) + 
                (*data_ptr * data_alpha) * t_1_m_alpha);
              data_ptr++;
              *data_ptr = (T)(
                (data_min + data_range * tprop_g * t_alpha) + 
                (*data_ptr * data_alpha) * t_1_m_alpha);
              data_ptr++;
              *data_ptr = (T)(
                (data_min + data_range * tprop_b * t_alpha) + 
                (*data_ptr * data_alpha) * t_1_m_alpha);
              data_ptr++;
              *data_ptr = (T)(
                data_min + data_range * (t_alpha + data_alpha * t_1_m_alpha));
              data_ptr++;
              glyph_ptr++;
              break;
            }
          }
        glyph_ptr_row += bitmap->pitch;
        data_ptr += data_pitch;
        }
      }

    // Advance to next char

    x += (bitmap_glyph->root.advance.x + 0x8000) >> 16;
    y += (bitmap_glyph->root.advance.y + 0x8000) >> 16;
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkFreeTypeUtilities::RenderString(vtkTextProperty *tprop, 
                                       const char *str,
                                       int x, int y,
                                       vtkImageData *data)
{
  // Check

  if (!tprop || !str || !data)
    {
    vtkErrorMacro(<< "Wrong parameters, one of them is NULL or zero");
    return 0;
    }

  if (data->GetNumberOfScalarComponents() > 4)
    {
    vtkErrorMacro("The image data must have a maximum of four components");
    return 0;
    }

  // Execute shadow

  int res = 1;

  if (tprop->GetShadow())
    {
    switch (data->GetScalarType())
      {
      vtkTemplateMacro(res &= vtkFreeTypeUtilitiesRenderString( 
                         this, 
                         tprop,
                         str,
                         x + 1, y - 1,
                         data, 
                         (VTK_TT *)(NULL),
                         1));
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
        return 0;
      }
    }

  // Execute text

  switch (data->GetScalarType())
    {
    vtkTemplateMacro(res &= vtkFreeTypeUtilitiesRenderString( 
                       this, 
                       tprop,
                       str,
                       x, y,
                       data, 
                       (VTK_TT *)(NULL),
                       0));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return 0;
    }

  return res;
}

//----------------------------------------------------------------------------
void vtkFreeTypeUtilities::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "MaximumNumberOfFaces: " 
     << this->MaximumNumberOfFaces << endl;

  os << indent << "MaximumNumberOfSizes: " 
     << this->MaximumNumberOfSizes << endl;

  os << indent << "MaximumNumberOfBytes: " 
     << this->MaximumNumberOfBytes << endl;
}

//----------------------------------------------------------------------------
// Print entry

void vtkFreeTypeUtilities::PrintEntry(int i, char *msg) 
{
  if (!this->Entries[i])
    {
    return;
    }

  printf("%s: [%2d] =", msg, i);

  vtkTextProperty *tprop = this->Entries[i]->TextProperty;
  if (tprop)
    {
    printf(" [S: %2d]", tprop->GetFontSize());

    printf(" [S: %d]", tprop->GetOrientation());

    double *color = tprop->GetColor();
    if (color)
      {
      printf(" [RGBA: %.2f/%.2f/%.2f (%.2f)]", 
             color[0], color[1], color[2], tprop->GetOpacity());
      }

    /*
    if (tprop->GetFaceFileName())
      {
      printf(" [F: %s]", tprop->GetFaceFileName());
      }
    else
    */
      {
      printf(" [F: %d] [I: %d] [B: %d]", 
             tprop->GetFontFamily(), tprop->GetItalic(), tprop->GetBold());
      }
    }
    
  if (this->Entries[i]->Font)
    {
    printf(" [F: %p]", (void*)this->Entries[i]->Font);
    printf("\n                                                [f: %p]", 
           (void*)*(this->Entries[i]->Font->Face()->Face()));
    }
  
  printf("\n");
  fflush(stdout);
}

//----------------------------------------------------------------------------
// Release entry

void vtkFreeTypeUtilities::ReleaseEntry(int i) 
{
  if (!this->Entries[i])
    {
    return;
    }

#if VTK_FTFC_DEBUG
  this->PrintEntry(this->NumberOfEntries, "Rl");
#endif

  if (this->Entries[i]->TextProperty)
    {
    this->Entries[i]->TextProperty->Delete();
    this->Entries[i]->TextProperty = NULL;
    }

  if (this->Entries[i]->Font)
    {
    delete this->Entries[i]->Font;
    this->Entries[i]->Font = NULL;
    }
  
  delete this->Entries[i];
  this->Entries[i] = NULL;
}

//----------------------------------------------------------------------------
// Initialize cache

void vtkFreeTypeUtilities::InitializeCache() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::InitializeCache()\n");
#endif  
  this->ReleaseCache();

  int i;
  for (i = 0; i < VTK_FTFC_CACHE_CAPACITY; i++)
    {
    this->Entries[i] = NULL;
    }
}

//----------------------------------------------------------------------------
// Release cache

void vtkFreeTypeUtilities::ReleaseCache() 
{
#if VTK_FTFC_DEBUG_CD
  printf("vtkFreeTypeUtilities::ReleaseCache()\n");
#endif  

  int i;
  for (i = 0; i < this->NumberOfEntries; i++)
    {
#if VTK_FTFC_DEBUG
    this->PrintEntry(i, "Rl");
#endif
    this->ReleaseEntry(i);
    }

  this->NumberOfEntries = 0;
}

//----------------------------------------------------------------------------
// Get a font from the cache given the text property. If no font is
// found in the cache, one is created and stored with the given color
// parameters. If override_color is true, then red, green, blue are used as
// text color instead of the colors found in the vtkTextProperty.

vtkFreeTypeUtilities::Entry* 
vtkFreeTypeUtilities::GetFont(vtkTextProperty *tprop, 
                              double override_color[3])
{
  int i, j;

  // Get the requested color and opacity

  double tprop_color[3];
  for (i = 0; i < 3; i++)
    {
    tprop_color[i] = override_color ? override_color[i] : tprop->GetColor()[i];
    if (tprop_color[i] < 0.0)
      {
      tprop_color[i] = 0.0;
      }
    }

  float tprop_opacity = 
    (tprop->GetOpacity() < 0.0) ? 1.0 : tprop->GetOpacity();
  
  // Has the font been cached ?
  
  for (i = 0; i < this->NumberOfEntries; i++)
    {
    vtkTextProperty *entry_tprop = this->Entries[i]->TextProperty;
    double *entry_tprop_color = entry_tprop->GetColor();

    if (
      // If a face file name has been specified, it overrides the 
      // font family as well as italic and bold attributes

      ((
        /*
          !entry_tprop->GetFaceFileName() && !tprop->GetFaceFileName() && 
        */
        entry_tprop->GetFontFamily() == tprop->GetFontFamily() &&
        entry_tprop->GetItalic()     == tprop->GetItalic() &&
        entry_tprop->GetBold()       == tprop->GetBold())
       /*
        || (entry_tprop->GetFaceFileName() && tprop->GetFaceFileName() && 
        !strcmp(entry_tprop->GetFaceFileName(), tprop->GetFaceFileName()))
       */
        )

      && (entry_tprop_color[0]      == tprop_color[0] &&
          entry_tprop_color[1]      == tprop_color[1] &&
          entry_tprop_color[2]      == tprop_color[2] &&
          entry_tprop->GetOpacity() == tprop_opacity)

      && entry_tprop->GetFontSize() == tprop->GetFontSize()
      && entry_tprop->GetOrientation() == tprop->GetOrientation())
      {
      // Make this the most recently used
      if (i != 0)
        {
        vtkFreeTypeUtilities::Entry *tmp = this->Entries[i];
        for (j = i - 1; j >= 0; j--)
          {
          this->Entries[j+1] = this->Entries[j];
          }
        this->Entries[0] = tmp;
        }
      return this->Entries[0];
      }
    }

  // OK the font is not cached, try to create one

  FTFont *font = new FTGLPixmapFont;

  // A face file name has been provided, try to load it, otherwise
  // just use the embedded fonts (i.e. font family, bold and italic attribs)

  /*
  if (tprop->GetFaceFileName())
    {
    if (!font->Open(tprop->GetFaceFileName(), false))
      {
      vtkErrorWithObjectMacro(
        tprop, 
        << "Unable to load font " << tprop->GetFaceFileName());
      delete font;
      return 0;
      }

    // Try to load an AFM metrics file for the PFB/PFA Postscript fonts

    int length = strlen(tprop->GetFaceFileName());
    if (length > 4 &&
        (!strcmp(tprop->GetFaceFileName() + length - 4, ".pfb") ||
         !strcmp(tprop->GetFaceFileName() + length - 4, ".pfa")))
      {
      char *metrics = new char[length + 1];
      strncpy(metrics, tprop->GetFaceFileName(), length - 3);
      strcpy(metrics + length - 3, "afm");
      struct stat fs;
      if (stat(metrics, &fs) == 0)
        {
        font->Attach(metrics); 
        }
      delete [] metrics;
      }
    }
  else
    */
    {
    // Fonts, organized by [Family][Bold][Italic]
    
    static EmbeddedFontStruct EmbeddedFonts[3][2][2] = 
    {
      {
        {
          { // VTK_ARIAL: Bold [ ] Italic [ ]
        
            face_arial_buffer_length, face_arial_buffer
          },
          { // VTK_ARIAL: Bold [ ] Italic [x]
            face_arial_italic_buffer_length, face_arial_italic_buffer
          }
        },
        {
          { // VTK_ARIAL: Bold [x] Italic [ ]
            face_arial_bold_buffer_length, face_arial_bold_buffer
          },
          { // VTK_ARIAL: Bold [x] Italic [x]
            face_arial_bold_italic_buffer_length, face_arial_bold_italic_buffer
          }
        }
      },
      {
        {
          { // VTK_COURIER: Bold [ ] Italic [ ]
            face_courier_buffer_length, face_courier_buffer
          },
          { // VTK_COURIER: Bold [ ] Italic [x]
            face_courier_italic_buffer_length, face_courier_italic_buffer
          }
        },
        {
          { // VTK_COURIER: Bold [x] Italic [ ]
            face_courier_bold_buffer_length, face_courier_bold_buffer
          },
          { // VTK_COURIER: Bold [x] Italic [x]
            face_courier_bold_italic_buffer_length, 
            face_courier_bold_italic_buffer
          }
        }
      },
      {
        {
          { // VTK_TIMES: Bold [ ] Italic [ ]
            face_times_buffer_length, face_times_buffer
          },
          { // VTK_TIMES: Bold [ ] Italic [x]
            face_times_italic_buffer_length, face_times_italic_buffer
          }
        },
        {
          { // VTK_TIMES: Bold [x] Italic [ ]
            face_times_bold_buffer_length, face_times_bold_buffer
          },
          { // VTK_TIMES: Bold [x] Italic [x]
            face_times_bold_italic_buffer_length, face_times_bold_italic_buffer
          }
        }
      }
    };
    
    size_t length = EmbeddedFonts
      [tprop->GetFontFamily()][tprop->GetBold()][tprop->GetItalic()].length;
    unsigned char *ptr = EmbeddedFonts
      [tprop->GetFontFamily()][tprop->GetBold()][tprop->GetItalic()].ptr;

    if (!font->Open(ptr, length, false))
      {
      vtkErrorWithObjectMacro(
        tprop,
        << "Unable to create font !" << " (family: " << tprop->GetFontFamily() 
        << ", bold: " << tprop->GetBold() << ", italic: " << tprop->GetItalic()
        << ", length: " << length << ")");
      delete font;
      return 0;
      }
    }

  // Set orientation

  font->Horizontal(tprop->GetOrientation() == 0. ? 1 : 0);
  
  // Set face size

  font->FaceSize(tprop->GetFontSize());
  
  // We need to make room for a new font

  if (this->NumberOfEntries == VTK_FTFC_CACHE_CAPACITY)
    {
#if VTK_FTFC_DEBUG
    printf("Cache is full, deleting last!\n");
#endif
    this->NumberOfEntries--;
    }

  // Add the new font

  if (this->Entries[this->NumberOfEntries])
    {
    this->ReleaseEntry(this->NumberOfEntries);
    }
  this->Entries[this->NumberOfEntries] = new vtkFreeTypeUtilities::Entry;
  
  // Set the other info

  this->Entries[this->NumberOfEntries]->TextProperty = vtkTextProperty::New();

  vtkTextProperty *entry_tprop = 
    this->Entries[this->NumberOfEntries]->TextProperty;

  entry_tprop->ShallowCopy(tprop);
  entry_tprop->SetOpacity(tprop_opacity);
  entry_tprop->SetColor(tprop_color);

  this->Entries[this->NumberOfEntries]->Font = font;

  this->Entries[this->NumberOfEntries]->LargestAscender  =
  this->Entries[this->NumberOfEntries]->LargestDescender = -1;

#if VTK_FTFC_DEBUG
  this->PrintEntry(this->NumberOfEntries, "Cr");
#endif

  vtkFreeTypeUtilities::Entry *tmp = this->Entries[this->NumberOfEntries];

  this->NumberOfEntries++;
  return tmp;
}
