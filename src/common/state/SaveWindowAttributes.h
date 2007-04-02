/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef SAVEWINDOWATTRIBUTES_H
#define SAVEWINDOWATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: SaveWindowAttributes
//
// Purpose:
//    This class contains the attributes used for saving windows.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Sun Jan 22 19:51:42 PST 2006
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SaveWindowAttributes : public AttributeSubject
{
public:
    enum FileFormat
    {
        BMP,
        CURVE,
        JPEG,
        OBJ,
        PNG,
        POSTSCRIPT,
        PPM,
        RGB,
        STL,
        TIFF,
        ULTRA,
        VTK
    };
    enum CompressionType
    {
        None,
        PackBits,
        Jpeg,
        Deflate
    };

    SaveWindowAttributes();
    SaveWindowAttributes(const SaveWindowAttributes &obj);
    virtual ~SaveWindowAttributes();

    virtual SaveWindowAttributes& operator = (const SaveWindowAttributes &obj);
    virtual bool operator == (const SaveWindowAttributes &obj) const;
    virtual bool operator != (const SaveWindowAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectOutputDirectory();
    void SelectFileName();
    void SelectLastRealFilename();

    // Property setting methods
    void SetOutputToCurrentDirectory(bool outputToCurrentDirectory_);
    void SetOutputDirectory(const std::string &outputDirectory_);
    void SetFileName(const std::string &fileName_);
    void SetFamily(bool family_);
    void SetFormat(FileFormat format_);
    void SetMaintainAspect(bool maintainAspect_);
    void SetWidth(int width_);
    void SetHeight(int height_);
    void SetScreenCapture(bool screenCapture_);
    void SetSaveTiled(bool saveTiled_);
    void SetQuality(int quality_);
    void SetProgressive(bool progressive_);
    void SetBinary(bool binary_);
    void SetLastRealFilename(const std::string &lastRealFilename_);
    void SetStereo(bool stereo_);
    void SetCompression(CompressionType compression_);

    // Property getting methods
    bool              GetOutputToCurrentDirectory() const;
    const std::string &GetOutputDirectory() const;
          std::string &GetOutputDirectory();
    const std::string &GetFileName() const;
          std::string &GetFileName();
    bool              GetFamily() const;
    FileFormat        GetFormat() const;
    bool              GetMaintainAspect() const;
    int               GetWidth() const;
    int               GetHeight() const;
    bool              GetScreenCapture() const;
    bool              GetSaveTiled() const;
    int               GetQuality() const;
    bool              GetProgressive() const;
    bool              GetBinary() const;
    const std::string &GetLastRealFilename() const;
          std::string &GetLastRealFilename();
    bool              GetStereo() const;
    CompressionType   GetCompression() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string FileFormat_ToString(FileFormat);
    static bool FileFormat_FromString(const std::string &, FileFormat &);
protected:
    static std::string FileFormat_ToString(int);
public:
    static std::string CompressionType_ToString(CompressionType);
    static bool CompressionType_FromString(const std::string &, CompressionType &);
protected:
    static std::string CompressionType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool CurrentFormatIsImageFormat(void);
private:
    bool        outputToCurrentDirectory;
    std::string outputDirectory;
    std::string fileName;
    bool        family;
    int         format;
    bool        maintainAspect;
    int         width;
    int         height;
    bool        screenCapture;
    bool        saveTiled;
    int         quality;
    bool        progressive;
    bool        binary;
    std::string lastRealFilename;
    bool        stereo;
    int         compression;
};

#endif
