/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

#ifndef MULTICURVEATTRIBUTES_H
#define MULTICURVEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorControlPointList.h>
#include <ColorAttribute.h>
#include <ColorAttributeList.h>

// ****************************************************************************
// Class: MultiCurveAttributes
//
// Purpose:
//    This class contains the plot attributes for the MultiCurve plot.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class MultiCurveAttributes : public AttributeSubject
{
public:
    enum ColoringMethod
    {
        ColorBySingleColor,
        ColorByMultipleColors
    };

    // These constructors are for objects of this class
    MultiCurveAttributes();
    MultiCurveAttributes(const MultiCurveAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    MultiCurveAttributes(private_tmfs_t tmfs);
    MultiCurveAttributes(const MultiCurveAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~MultiCurveAttributes();

    virtual MultiCurveAttributes& operator = (const MultiCurveAttributes &obj);
    virtual bool operator == (const MultiCurveAttributes &obj) const;
    virtual bool operator != (const MultiCurveAttributes &obj) const;
private:
    void Init();
    void Copy(const MultiCurveAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectDefaultPalette();
    void SelectChangedColors();
    void SelectSingleColor();
    void SelectMultiColor();
    void SelectYAxisTitleFormat();
    void SelectMarkerVariable();
    void SelectIdVariable();

    // Property setting methods
    void SetDefaultPalette(const ColorControlPointList &defaultPalette_);
    void SetChangedColors(const unsignedCharVector &changedColors_);
    void SetColorType(ColoringMethod colorType_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetMultiColor(const ColorAttributeList &multiColor_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetYAxisTitleFormat(const std::string &yAxisTitleFormat_);
    void SetUseYAxisTickSpacing(bool useYAxisTickSpacing_);
    void SetYAxisTickSpacing(double yAxisTickSpacing_);
    void SetDisplayMarkers(bool displayMarkers_);
    void SetMarkerScale(double markerScale_);
    void SetMarkerLineWidth(int markerLineWidth_);
    void SetMarkerVariable(const std::string &markerVariable_);
    void SetDisplayIds(bool displayIds_);
    void SetIdVariable(const std::string &idVariable_);
    void SetLegendFlag(bool legendFlag_);

    // Property getting methods
    const ColorControlPointList &GetDefaultPalette() const;
          ColorControlPointList &GetDefaultPalette();
    const unsignedCharVector    &GetChangedColors() const;
          unsignedCharVector    &GetChangedColors();
    ColoringMethod              GetColorType() const;
    const ColorAttribute        &GetSingleColor() const;
          ColorAttribute        &GetSingleColor();
    const ColorAttributeList    &GetMultiColor() const;
          ColorAttributeList    &GetMultiColor();
    int                         GetLineStyle() const;
    int                         GetLineWidth() const;
    const std::string           &GetYAxisTitleFormat() const;
          std::string           &GetYAxisTitleFormat();
    bool                        GetUseYAxisTickSpacing() const;
    double                      GetYAxisTickSpacing() const;
    bool                        GetDisplayMarkers() const;
    double                      GetMarkerScale() const;
    int                         GetMarkerLineWidth() const;
    const std::string           &GetMarkerVariable() const;
          std::string           &GetMarkerVariable();
    bool                        GetDisplayIds() const;
    const std::string           &GetIdVariable() const;
          std::string           &GetIdVariable();
    bool                        GetLegendFlag() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
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
    void EnlargeMultiColor(int newSize);
    bool ColorIsChanged(int index) const;
    void MarkColorAsChanged(int index);
    bool ChangesRequireRecalculation(const MultiCurveAttributes &) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_defaultPalette = 0,
        ID_changedColors,
        ID_colorType,
        ID_singleColor,
        ID_multiColor,
        ID_lineStyle,
        ID_lineWidth,
        ID_yAxisTitleFormat,
        ID_useYAxisTickSpacing,
        ID_yAxisTickSpacing,
        ID_displayMarkers,
        ID_markerScale,
        ID_markerLineWidth,
        ID_markerVariable,
        ID_displayIds,
        ID_idVariable,
        ID_legendFlag,
        ID__LAST
    };

private:
    ColorControlPointList defaultPalette;
    unsignedCharVector    changedColors;
    int                   colorType;
    ColorAttribute        singleColor;
    ColorAttributeList    multiColor;
    int                   lineStyle;
    int                   lineWidth;
    std::string           yAxisTitleFormat;
    bool                  useYAxisTickSpacing;
    double                yAxisTickSpacing;
    bool                  displayMarkers;
    double                markerScale;
    int                   markerLineWidth;
    std::string           markerVariable;
    bool                  displayIds;
    std::string           idVariable;
    bool                  legendFlag;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define MULTICURVEATTRIBUTES_TMFS "au*iaaiisbdbdisbsb"

#endif
