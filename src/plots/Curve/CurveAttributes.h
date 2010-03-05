/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
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

#ifndef CURVEATTRIBUTES_H
#define CURVEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>
#include <visitstream.h>

// ****************************************************************************
// Class: CurveAttributes
//
// Purpose:
//    Attributes for the xy plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class CurveAttributes : public AttributeSubject
{
public:
    enum RenderMode
    {
        RenderAsLines,
        RenderAsSymbols
    };
    enum CurveColor
    {
        Cycle,
        Custom
    };
    enum SymbolTypes
    {
        TriangleUp,
        TriangleDown,
        Square,
        Circle,
        Plus,
        X
    };

    // These constructors are for objects of this class
    CurveAttributes();
    CurveAttributes(const CurveAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    CurveAttributes(private_tmfs_t tmfs);
    CurveAttributes(const CurveAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~CurveAttributes();

    virtual CurveAttributes& operator = (const CurveAttributes &obj);
    virtual bool operator == (const CurveAttributes &obj) const;
    virtual bool operator != (const CurveAttributes &obj) const;
private:
    void Init();
    void Copy(const CurveAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectCurveColor();
    void SelectDesignator();

    // Property setting methods
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetCurveColor(const ColorAttribute &curveColor_);
    void SetShowLabels(bool showLabels_);
    void SetDesignator(const std::string &designator_);
    void SetShowPoints(bool showPoints_);
    void SetPointSize(double pointSize_);
    void SetShowLegend(bool showLegend_);
    void SetCurveColorSource(CurveColor curveColorSource_);
    void SetRenderMode(RenderMode renderMode_);
    void SetSymbol(SymbolTypes symbol_);
    void SetSymbolDensity(int symbolDensity_);

    // Property getting methods
    int                  GetLineStyle() const;
    int                  GetLineWidth() const;
    const ColorAttribute &GetCurveColor() const;
          ColorAttribute &GetCurveColor();
    bool                 GetShowLabels() const;
    const std::string    &GetDesignator() const;
          std::string    &GetDesignator();
    bool                 GetShowPoints() const;
    double               GetPointSize() const;
    bool                 GetShowLegend() const;
    CurveColor           GetCurveColorSource() const;
    RenderMode           GetRenderMode() const;
    SymbolTypes          GetSymbol() const;
    int                  GetSymbolDensity() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string RenderMode_ToString(RenderMode);
    static bool RenderMode_FromString(const std::string &, RenderMode &);
protected:
    static std::string RenderMode_ToString(int);
public:
    static std::string CurveColor_ToString(CurveColor);
    static bool CurveColor_FromString(const std::string &, CurveColor &);
protected:
    static std::string CurveColor_ToString(int);
public:
    static std::string SymbolTypes_ToString(SymbolTypes);
    static bool SymbolTypes_FromString(const std::string &, SymbolTypes &);
protected:
    static std::string SymbolTypes_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const CurveAttributes &) const;
    void Print(ostream &, bool) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_lineStyle = 0,
        ID_lineWidth,
        ID_curveColor,
        ID_showLabels,
        ID_designator,
        ID_showPoints,
        ID_pointSize,
        ID_showLegend,
        ID_curveColorSource,
        ID_renderMode,
        ID_symbol,
        ID_symbolDensity,
        ID__LAST
    };

private:
    int            lineStyle;
    int            lineWidth;
    ColorAttribute curveColor;
    bool           showLabels;
    std::string    designator;
    bool           showPoints;
    double         pointSize;
    bool           showLegend;
    int            curveColorSource;
    int            renderMode;
    int            symbol;
    int            symbolDensity;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define CURVEATTRIBUTES_TMFS "iiabsbdbiiii"

#endif
