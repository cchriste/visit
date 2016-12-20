/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

#ifndef FILLEDBOUNDARYATTRIBUTES_H
#define FILLEDBOUNDARYATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>
#include <ColorAttributeList.h>

// ****************************************************************************
// Class: FilledBoundaryAttributes
//
// Purpose:
//    This class contains the plot attributes for the filled boundary plot.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class FilledBoundaryAttributes : public AttributeSubject
{
public:
    enum ColoringMethod
    {
        ColorBySingleColor,
        ColorByMultipleColors,
        ColorByColorTable
    };
    enum PointType
    {
        Box,
        Axis,
        Icosahedron,
        Octahedron,
        Tetrahedron,
        SphereGeometry,
        Point,
        Sphere
    };

    // These constructors are for objects of this class
    FilledBoundaryAttributes();
    FilledBoundaryAttributes(const FilledBoundaryAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    FilledBoundaryAttributes(private_tmfs_t tmfs);
    FilledBoundaryAttributes(const FilledBoundaryAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~FilledBoundaryAttributes();

    virtual FilledBoundaryAttributes& operator = (const FilledBoundaryAttributes &obj);
    virtual bool operator == (const FilledBoundaryAttributes &obj) const;
    virtual bool operator != (const FilledBoundaryAttributes &obj) const;
private:
    void Init();
    void Copy(const FilledBoundaryAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectColorTableName();
    void SelectSingleColor();
    void SelectMultiColor();
    void SelectBoundaryNames();
    void SelectMixedColor();
    void SelectPointSizeVar();

    // Property setting methods
    void SetColorType(ColoringMethod colorType_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetInvertColorTable(bool invertColorTable_);
    void SetLegendFlag(bool legendFlag_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetMultiColor(const ColorAttributeList &multiColor_);
    void SetBoundaryNames(const stringVector &boundaryNames_);
    void SetOpacity(double opacity_);
    void SetWireframe(bool wireframe_);
    void SetDrawInternal(bool drawInternal_);
    void SetSmoothingLevel(int smoothingLevel_);
    void SetCleanZonesOnly(bool cleanZonesOnly_);
    void SetMixedColor(const ColorAttribute &mixedColor_);
    void SetPointSize(double pointSize_);
    void SetPointType(PointType pointType_);
    void SetPointSizeVarEnabled(bool pointSizeVarEnabled_);
    void SetPointSizeVar(const std::string &pointSizeVar_);
    void SetPointSizePixels(int pointSizePixels_);

    // Property getting methods
    ColoringMethod           GetColorType() const;
    const std::string        &GetColorTableName() const;
          std::string        &GetColorTableName();
    bool                     GetInvertColorTable() const;
    bool                     GetLegendFlag() const;
    int                      GetLineStyle() const;
    int                      GetLineWidth() const;
    const ColorAttribute     &GetSingleColor() const;
          ColorAttribute     &GetSingleColor();
    const ColorAttributeList &GetMultiColor() const;
          ColorAttributeList &GetMultiColor();
    const stringVector       &GetBoundaryNames() const;
          stringVector       &GetBoundaryNames();
    double                   GetOpacity() const;
    bool                     GetWireframe() const;
    bool                     GetDrawInternal() const;
    int                      GetSmoothingLevel() const;
    bool                     GetCleanZonesOnly() const;
    const ColorAttribute     &GetMixedColor() const;
          ColorAttribute     &GetMixedColor();
    double                   GetPointSize() const;
    PointType                GetPointType() const;
    bool                     GetPointSizeVarEnabled() const;
    const std::string        &GetPointSizeVar() const;
          std::string        &GetPointSizeVar();
    int                      GetPointSizePixels() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string ColoringMethod_ToString(ColoringMethod);
    static bool ColoringMethod_FromString(const std::string &, ColoringMethod &);
protected:
    static std::string ColoringMethod_ToString(int);
public:
    static std::string PointType_ToString(PointType);
    static bool PointType_FromString(const std::string &, PointType &);
protected:
    static std::string PointType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const FilledBoundaryAttributes &obj);
    virtual bool VarChangeRequiresReset(void);
    virtual void ProcessOldVersions(DataNode *parentNode, const char *configVersion);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_colorType = 0,
        ID_colorTableName,
        ID_invertColorTable,
        ID_legendFlag,
        ID_lineStyle,
        ID_lineWidth,
        ID_singleColor,
        ID_multiColor,
        ID_boundaryNames,
        ID_opacity,
        ID_wireframe,
        ID_drawInternal,
        ID_smoothingLevel,
        ID_cleanZonesOnly,
        ID_mixedColor,
        ID_pointSize,
        ID_pointType,
        ID_pointSizeVarEnabled,
        ID_pointSizeVar,
        ID_pointSizePixels,
        ID__LAST
    };

private:
    int                colorType;
    std::string        colorTableName;
    bool               invertColorTable;
    bool               legendFlag;
    int                lineStyle;
    int                lineWidth;
    ColorAttribute     singleColor;
    ColorAttributeList multiColor;
    stringVector       boundaryNames;
    double             opacity;
    bool               wireframe;
    bool               drawInternal;
    int                smoothingLevel;
    bool               cleanZonesOnly;
    ColorAttribute     mixedColor;
    double             pointSize;
    int                pointType;
    bool               pointSizeVarEnabled;
    std::string        pointSizeVar;
    int                pointSizePixels;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define FILLEDBOUNDARYATTRIBUTES_TMFS "isbbiiaas*dbbibadibsi"

#endif
