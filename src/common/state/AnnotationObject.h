/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#ifndef ANNOTATIONOBJECT_H
#define ANNOTATIONOBJECT_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>

// ****************************************************************************
// Class: AnnotationObject
//
// Purpose:
//    This class defines a general set of attributes that are used to set the attributes for all annotation objects.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API AnnotationObject : public AttributeSubject
{
public:
    enum AnnotationType
    {
        Text2D,
        Text3D,
        TimeSlider,
        Line2D,
        Line3D,
        Arrow2D,
        Arrow3D,
        Box,
        Image,
        LegendAttributes,
        MaxAnnotationType
    };
    enum FontFamily
    {
        Arial,
        Courier,
        Times
    };

    // These constructors are for objects of this class
    AnnotationObject();
    AnnotationObject(const AnnotationObject &obj);
protected:
    // These constructors are for objects derived from this class
    AnnotationObject(private_tmfs_t tmfs);
    AnnotationObject(const AnnotationObject &obj, private_tmfs_t tmfs);
public:
    virtual ~AnnotationObject();

    virtual AnnotationObject& operator = (const AnnotationObject &obj);
    virtual bool operator == (const AnnotationObject &obj) const;
    virtual bool operator != (const AnnotationObject &obj) const;
private:
    void Init();
    void Copy(const AnnotationObject &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectObjectName();
    void SelectPosition();
    void SelectPosition2();
    void SelectTextColor();
    void SelectColor1();
    void SelectColor2();
    void SelectText();
    void SelectDoubleVector1();
    void SelectStringVector1();
    void SelectStringVector2();

    // Property setting methods
    void SetObjectName(const std::string &objectName_);
    void SetObjectType(AnnotationType objectType_);
    void SetVisible(bool visible_);
    void SetActive(bool active_);
    void SetPosition(const double *position_);
    void SetPosition2(const double *position2_);
    void SetTextColor(const ColorAttribute &textColor_);
    void SetUseForegroundForTextColor(bool useForegroundForTextColor_);
    void SetColor1(const ColorAttribute &color1_);
    void SetColor2(const ColorAttribute &color2_);
    void SetText(const stringVector &text_);
    void SetFontFamily(FontFamily fontFamily_);
    void SetFontBold(bool fontBold_);
    void SetFontItalic(bool fontItalic_);
    void SetFontShadow(bool fontShadow_);
    void SetDoubleAttribute1(double doubleAttribute1_);
    void SetIntAttribute1(int intAttribute1_);
    void SetIntAttribute2(int intAttribute2_);
    void SetIntAttribute3(int intAttribute3_);
    void SetDoubleVector1(const doubleVector &doubleVector1_);
    void SetStringVector1(const stringVector &stringVector1_);
    void SetStringVector2(const stringVector &stringVector2_);

    // Property getting methods
    const std::string    &GetObjectName() const;
          std::string    &GetObjectName();
    AnnotationType       GetObjectType() const;
    bool                 GetVisible() const;
    bool                 GetActive() const;
    const double         *GetPosition() const;
          double         *GetPosition();
    const double         *GetPosition2() const;
          double         *GetPosition2();
    const ColorAttribute &GetTextColor() const;
          ColorAttribute &GetTextColor();
    bool                 GetUseForegroundForTextColor() const;
    const ColorAttribute &GetColor1() const;
          ColorAttribute &GetColor1();
    const ColorAttribute &GetColor2() const;
          ColorAttribute &GetColor2();
    const stringVector   &GetText() const;
          stringVector   &GetText();
    FontFamily           GetFontFamily() const;
    bool                 GetFontBold() const;
    bool                 GetFontItalic() const;
    bool                 GetFontShadow() const;
    double               GetDoubleAttribute1() const;
    int                  GetIntAttribute1() const;
    int                  GetIntAttribute2() const;
    int                  GetIntAttribute3() const;
    const doubleVector   &GetDoubleVector1() const;
          doubleVector   &GetDoubleVector1();
    const stringVector   &GetStringVector1() const;
          stringVector   &GetStringVector1();
    const stringVector   &GetStringVector2() const;
          stringVector   &GetStringVector2();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string AnnotationType_ToString(AnnotationType);
    static bool AnnotationType_FromString(const std::string &, AnnotationType &);
protected:
    static std::string AnnotationType_ToString(int);
public:
    static std::string FontFamily_ToString(FontFamily);
    static bool FontFamily_FromString(const std::string &, FontFamily &);
protected:
    static std::string FontFamily_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual void ProcessOldVersions(DataNode *parentNode, const char *configVersion);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_objectName = 0,
        ID_objectType,
        ID_visible,
        ID_active,
        ID_position,
        ID_position2,
        ID_textColor,
        ID_useForegroundForTextColor,
        ID_color1,
        ID_color2,
        ID_text,
        ID_fontFamily,
        ID_fontBold,
        ID_fontItalic,
        ID_fontShadow,
        ID_doubleAttribute1,
        ID_intAttribute1,
        ID_intAttribute2,
        ID_intAttribute3,
        ID_doubleVector1,
        ID_stringVector1,
        ID_stringVector2,
        ID__LAST
    };

private:
    std::string    objectName;
    int            objectType;
    bool           visible;
    bool           active;
    double         position[3];
    double         position2[3];
    ColorAttribute textColor;
    bool           useForegroundForTextColor;
    ColorAttribute color1;
    ColorAttribute color2;
    stringVector   text;
    int            fontFamily;
    bool           fontBold;
    bool           fontItalic;
    bool           fontShadow;
    double         doubleAttribute1;
    int            intAttribute1;
    int            intAttribute2;
    int            intAttribute3;
    doubleVector   doubleVector1;
    stringVector   stringVector1;
    stringVector   stringVector2;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define ANNOTATIONOBJECT_TMFS "sibbDDabaas*ibbbdiiid*s*s*"

#endif
