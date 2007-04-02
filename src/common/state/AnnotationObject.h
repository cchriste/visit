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
// Creation:   Fri Mar 31 14:20:56 PST 2006
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
        Arrow2D,
        Arrow3D,
        Box,
        Image
    };
    enum FontFamily
    {
        Arial,
        Courier,
        Times
    };

    AnnotationObject();
    AnnotationObject(const AnnotationObject &obj);
    virtual ~AnnotationObject();

    virtual AnnotationObject& operator = (const AnnotationObject &obj);
    virtual bool operator == (const AnnotationObject &obj) const;
    virtual bool operator != (const AnnotationObject &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectPosition();
    void SelectPosition2();
    void SelectTextColor();
    void SelectColor1();
    void SelectColor2();
    void SelectText();

    // Property setting methods
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

    // Property getting methods
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

private:
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
};

#endif
