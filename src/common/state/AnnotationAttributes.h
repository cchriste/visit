/*****************************************************************************
*
* Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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

#ifndef ANNOTATIONATTRIBUTES_H
#define ANNOTATIONATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

#include <Axes2D.h>
#include <Axes3D.h>
#include <FontAttributes.h>
#include <ColorAttribute.h>
#include <AxesArray.h>

// ****************************************************************************
// Class: AnnotationAttributes
//
// Purpose:
//    This class contains the attributes controlling annotations.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API AnnotationAttributes : public AttributeSubject
{
public:
    enum GradientStyle
    {
        TopToBottom,
        BottomToTop,
        LeftToRight,
        RightToLeft,
        Radial
    };
    enum BackgroundMode
    {
        Solid,
        Gradient,
        Image,
        ImageSphere
    };
    enum PathExpansionMode
    {
        File,
        Directory,
        Full,
        Smart,
        SmartDirectory
    };

    // These constructors are for objects of this class
    AnnotationAttributes();
    AnnotationAttributes(const AnnotationAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    AnnotationAttributes(private_tmfs_t tmfs);
    AnnotationAttributes(const AnnotationAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~AnnotationAttributes();

    virtual AnnotationAttributes& operator = (const AnnotationAttributes &obj);
    virtual bool operator == (const AnnotationAttributes &obj) const;
    virtual bool operator != (const AnnotationAttributes &obj) const;
private:
    void Init();
    void Copy(const AnnotationAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectAxes2D();
    void SelectAxes3D();
    void SelectUserInfoFont();
    void SelectDatabaseInfoFont();
    void SelectBackgroundColor();
    void SelectForegroundColor();
    void SelectGradientColor1();
    void SelectGradientColor2();
    void SelectBackgroundImage();
    void SelectAxesArray();

    // Property setting methods
    void SetAxes2D(const Axes2D &axes2D_);
    void SetAxes3D(const Axes3D &axes3D_);
    void SetUserInfoFlag(bool userInfoFlag_);
    void SetUserInfoFont(const FontAttributes &userInfoFont_);
    void SetDatabaseInfoFlag(bool databaseInfoFlag_);
    void SetTimeInfoFlag(bool timeInfoFlag_);
    void SetDatabaseInfoFont(const FontAttributes &databaseInfoFont_);
    void SetDatabaseInfoExpansionMode(PathExpansionMode databaseInfoExpansionMode_);
    void SetDatabaseInfoTimeScale(double databaseInfoTimeScale_);
    void SetDatabaseInfoTimeOffset(double databaseInfoTimeOffset_);
    void SetLegendInfoFlag(bool legendInfoFlag_);
    void SetBackgroundColor(const ColorAttribute &backgroundColor_);
    void SetForegroundColor(const ColorAttribute &foregroundColor_);
    void SetGradientBackgroundStyle(GradientStyle gradientBackgroundStyle_);
    void SetGradientColor1(const ColorAttribute &gradientColor1_);
    void SetGradientColor2(const ColorAttribute &gradientColor2_);
    void SetBackgroundMode(BackgroundMode backgroundMode_);
    void SetBackgroundImage(const std::string &backgroundImage_);
    void SetImageRepeatX(int imageRepeatX_);
    void SetImageRepeatY(int imageRepeatY_);
    void SetAxesArray(const AxesArray &axesArray_);

    // Property getting methods
    const Axes2D         &GetAxes2D() const;
          Axes2D         &GetAxes2D();
    const Axes3D         &GetAxes3D() const;
          Axes3D         &GetAxes3D();
    bool                 GetUserInfoFlag() const;
    const FontAttributes &GetUserInfoFont() const;
          FontAttributes &GetUserInfoFont();
    bool                 GetDatabaseInfoFlag() const;
    bool                 GetTimeInfoFlag() const;
    const FontAttributes &GetDatabaseInfoFont() const;
          FontAttributes &GetDatabaseInfoFont();
    PathExpansionMode    GetDatabaseInfoExpansionMode() const;
    double               GetDatabaseInfoTimeScale() const;
    double               GetDatabaseInfoTimeOffset() const;
    bool                 GetLegendInfoFlag() const;
    const ColorAttribute &GetBackgroundColor() const;
          ColorAttribute &GetBackgroundColor();
    const ColorAttribute &GetForegroundColor() const;
          ColorAttribute &GetForegroundColor();
    GradientStyle        GetGradientBackgroundStyle() const;
    const ColorAttribute &GetGradientColor1() const;
          ColorAttribute &GetGradientColor1();
    const ColorAttribute &GetGradientColor2() const;
          ColorAttribute &GetGradientColor2();
    BackgroundMode       GetBackgroundMode() const;
    const std::string    &GetBackgroundImage() const;
          std::string    &GetBackgroundImage();
    int                  GetImageRepeatX() const;
    int                  GetImageRepeatY() const;
    const AxesArray      &GetAxesArray() const;
          AxesArray      &GetAxesArray();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string GradientStyle_ToString(GradientStyle);
    static bool GradientStyle_FromString(const std::string &, GradientStyle &);
protected:
    static std::string GradientStyle_ToString(int);
public:
    static std::string BackgroundMode_ToString(BackgroundMode);
    static bool BackgroundMode_FromString(const std::string &, BackgroundMode &);
protected:
    static std::string BackgroundMode_ToString(int);
public:
    static std::string PathExpansionMode_ToString(PathExpansionMode);
    static bool PathExpansionMode_FromString(const std::string &, PathExpansionMode &);
protected:
    static std::string PathExpansionMode_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    const ColorAttribute GetDiscernibleBackgroundColor() const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_axes2D = 0,
        ID_axes3D,
        ID_userInfoFlag,
        ID_userInfoFont,
        ID_databaseInfoFlag,
        ID_timeInfoFlag,
        ID_databaseInfoFont,
        ID_databaseInfoExpansionMode,
        ID_databaseInfoTimeScale,
        ID_databaseInfoTimeOffset,
        ID_legendInfoFlag,
        ID_backgroundColor,
        ID_foregroundColor,
        ID_gradientBackgroundStyle,
        ID_gradientColor1,
        ID_gradientColor2,
        ID_backgroundMode,
        ID_backgroundImage,
        ID_imageRepeatX,
        ID_imageRepeatY,
        ID_axesArray,
        ID__LAST
    };

private:
    Axes2D         axes2D;
    Axes3D         axes3D;
    bool           userInfoFlag;
    FontAttributes userInfoFont;
    bool           databaseInfoFlag;
    bool           timeInfoFlag;
    FontAttributes databaseInfoFont;
    int            databaseInfoExpansionMode;
    double         databaseInfoTimeScale;
    double         databaseInfoTimeOffset;
    bool           legendInfoFlag;
    ColorAttribute backgroundColor;
    ColorAttribute foregroundColor;
    int            gradientBackgroundStyle;
    ColorAttribute gradientColor1;
    ColorAttribute gradientColor2;
    int            backgroundMode;
    std::string    backgroundImage;
    int            imageRepeatX;
    int            imageRepeatY;
    AxesArray      axesArray;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define ANNOTATIONATTRIBUTES_TMFS "aababbaiddbaaiaaisiia"

#endif
