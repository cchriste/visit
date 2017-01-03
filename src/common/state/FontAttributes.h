/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
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
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef FONTATTRIBUTES_H
#define FONTATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>

// ****************************************************************************
// Class: FontAttributes
//
// Purpose:
//    Describes font properties that we can set through VTK
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API FontAttributes : public AttributeSubject
{
public:
    enum FontName
    {
        Arial,
        Courier,
        Times
    };

    // These constructors are for objects of this class
    FontAttributes();
    FontAttributes(const FontAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    FontAttributes(private_tmfs_t tmfs);
    FontAttributes(const FontAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~FontAttributes();

    virtual FontAttributes& operator = (const FontAttributes &obj);
    virtual bool operator == (const FontAttributes &obj) const;
    virtual bool operator != (const FontAttributes &obj) const;
private:
    void Init();
    void Copy(const FontAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectColor();

    // Property setting methods
    void SetFont(FontName font_);
    void SetScale(double scale_);
    void SetUseForegroundColor(bool useForegroundColor_);
    void SetColor(const ColorAttribute &color_);
    void SetBold(bool bold_);
    void SetItalic(bool italic_);

    // Property getting methods
    FontName             GetFont() const;
    double               GetScale() const;
    bool                 GetUseForegroundColor() const;
    const ColorAttribute &GetColor() const;
          ColorAttribute &GetColor();
    bool                 GetBold() const;
    bool                 GetItalic() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string FontName_ToString(FontName);
    static bool FontName_FromString(const std::string &, FontName &);
protected:
    static std::string FontName_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_font = 0,
        ID_scale,
        ID_useForegroundColor,
        ID_color,
        ID_bold,
        ID_italic,
        ID__LAST
    };

private:
    int            font;
    double         scale;
    bool           useForegroundColor;
    ColorAttribute color;
    bool           bold;
    bool           italic;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define FONTATTRIBUTES_TMFS "idbabb"

#endif
