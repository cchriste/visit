/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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

#ifndef CONTOURATTRIBUTES_H
#define CONTOURATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorControlPointList.h>
#include <ColorAttribute.h>
#include <ColorAttributeList.h>

// ****************************************************************************
// Class: ContourAttributes
//
// Purpose:
//    This class contains the plot attributes for the contour plot.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class ContourAttributes : public AttributeSubject
{
public:
    enum Select_by
    {
        Level,
        Value,
        Percent
    };
    enum Scaling
    {
        Linear,
        Log
    };
    enum ColoringMethod
    {
        ColorBySingleColor,
        ColorByMultipleColors,
        ColorByColorTable
    };
    static const int MAX_CONTOURS;

    ContourAttributes();
    ContourAttributes(const ContourAttributes &obj);
    virtual ~ContourAttributes();

    virtual ContourAttributes& operator = (const ContourAttributes &obj);
    virtual bool operator == (const ContourAttributes &obj) const;
    virtual bool operator != (const ContourAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectDefaultPalette();
    void SelectChangedColors();
    void SelectColorTableName();
    void SelectSingleColor();
    void SelectMultiColor();
    void SelectContourValue();
    void SelectContourPercent();

    // Property setting methods
    void SetDefaultPalette(const ColorControlPointList &defaultPalette_);
    void SetChangedColors(const unsignedCharVector &changedColors_);
    void SetColorType(ColoringMethod colorType_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetLegendFlag(bool legendFlag_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetMultiColor(const ColorAttributeList &multiColor_);
    void SetContourNLevels(int contourNLevels_);
    void SetContourValue(const doubleVector &contourValue_);
    void SetContourPercent(const doubleVector &contourPercent_);
    void SetContourMethod(Select_by contourMethod_);
    void SetMinFlag(bool minFlag_);
    void SetMaxFlag(bool maxFlag_);
    void SetMin(double min_);
    void SetMax(double max_);
    void SetScaling(Scaling scaling_);
    void SetWireframe(bool wireframe_);

    // Property getting methods
    const ColorControlPointList &GetDefaultPalette() const;
          ColorControlPointList &GetDefaultPalette();
    const unsignedCharVector    &GetChangedColors() const;
          unsignedCharVector    &GetChangedColors();
    ColoringMethod              GetColorType() const;
    const std::string           &GetColorTableName() const;
          std::string           &GetColorTableName();
    bool                        GetLegendFlag() const;
    int                         GetLineStyle() const;
    int                         GetLineWidth() const;
    const ColorAttribute        &GetSingleColor() const;
          ColorAttribute        &GetSingleColor();
    const ColorAttributeList    &GetMultiColor() const;
          ColorAttributeList    &GetMultiColor();
    int                         GetContourNLevels() const;
    const doubleVector          &GetContourValue() const;
          doubleVector          &GetContourValue();
    const doubleVector          &GetContourPercent() const;
          doubleVector          &GetContourPercent();
    Select_by                   GetContourMethod() const;
    bool                        GetMinFlag() const;
    bool                        GetMaxFlag() const;
    double                      GetMin() const;
    double                      GetMax() const;
    Scaling                     GetScaling() const;
    bool                        GetWireframe() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Select_by_ToString(Select_by);
    static bool Select_by_FromString(const std::string &, Select_by &);
protected:
    static std::string Select_by_ToString(int);
public:
    static std::string Scaling_ToString(Scaling);
    static bool Scaling_FromString(const std::string &, Scaling &);
protected:
    static std::string Scaling_ToString(int);
public:
    static std::string ColoringMethod_ToString(ColoringMethod);
    static bool ColoringMethod_FromString(const std::string &, ColoringMethod &);
protected:
    static std::string ColoringMethod_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const ContourAttributes &obj);
    void SetContourValue(int i, double d);
    void SetContourPercent(int i, double d);
    void EnlargeMultiColor(int newSize);
    bool ColorIsChanged(int index) const;
    void MarkColorAsChanged(int index);
    virtual bool SetValue(const std::string &name, const int &value);
    virtual bool SetValue(const std::string &name, const doubleVector &value);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_defaultPalette = 0,
        ID_changedColors,
        ID_colorType,
        ID_colorTableName,
        ID_legendFlag,
        ID_lineStyle,
        ID_lineWidth,
        ID_singleColor,
        ID_multiColor,
        ID_contourNLevels,
        ID_contourValue,
        ID_contourPercent,
        ID_contourMethod,
        ID_minFlag,
        ID_maxFlag,
        ID_min,
        ID_max,
        ID_scaling,
        ID_wireframe
    };

private:
    ColorControlPointList defaultPalette;
    unsignedCharVector    changedColors;
    int                   colorType;
    std::string           colorTableName;
    bool                  legendFlag;
    int                   lineStyle;
    int                   lineWidth;
    ColorAttribute        singleColor;
    ColorAttributeList    multiColor;
    int                   contourNLevels;
    doubleVector          contourValue;
    doubleVector          contourPercent;
    int                   contourMethod;
    bool                  minFlag;
    bool                  maxFlag;
    double                min;
    double                max;
    int                   scaling;
    bool                  wireframe;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
