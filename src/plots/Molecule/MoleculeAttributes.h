/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

#ifndef MOLECULEATTRIBUTES_H
#define MOLECULEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: MoleculeAttributes
//
// Purpose:
//    This class contains the plot attributes for the molecule plot.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class MoleculeAttributes : public AttributeSubject
{
public:
    enum AtomRenderingMode
    {
        NoAtoms,
        SphereAtoms,
        ImposterAtoms
    };
    enum RadiusType
    {
        Fixed,
        Covalent,
        Atomic,
        Variable
    };
    enum BondColoringMode
    {
        ColorByAtom,
        SingleColor
    };
    enum BondRenderingMode
    {
        NoBonds,
        LineBonds,
        CylinderBonds
    };
    enum DetailLevel
    {
        Low,
        Medium,
        High,
        Super
    };

    MoleculeAttributes();
    MoleculeAttributes(const MoleculeAttributes &obj);
    virtual ~MoleculeAttributes();

    virtual MoleculeAttributes& operator = (const MoleculeAttributes &obj);
    virtual bool operator == (const MoleculeAttributes &obj) const;
    virtual bool operator != (const MoleculeAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectBondSingleColor();
    void SelectRadiusVariable();
    void SelectElementColorTable();
    void SelectResidueTypeColorTable();
    void SelectResidueSequenceColorTable();
    void SelectContinuousColorTable();

    // Property setting methods
    void SetDrawAtomsAs(AtomRenderingMode drawAtomsAs_);
    void SetScaleRadiusBy(RadiusType scaleRadiusBy_);
    void SetDrawBondsAs(BondRenderingMode drawBondsAs_);
    void SetColorBonds(BondColoringMode colorBonds_);
    void SetBondSingleColor(const ColorAttribute &bondSingleColor_);
    void SetRadiusVariable(const std::string &radiusVariable_);
    void SetRadiusScaleFactor(float radiusScaleFactor_);
    void SetRadiusFixed(float radiusFixed_);
    void SetAtomSphereQuality(DetailLevel atomSphereQuality_);
    void SetBondCylinderQuality(DetailLevel bondCylinderQuality_);
    void SetBondRadius(float bondRadius_);
    void SetBondLineWidth(int bondLineWidth_);
    void SetBondLineStyle(int bondLineStyle_);
    void SetElementColorTable(const std::string &elementColorTable_);
    void SetResidueTypeColorTable(const std::string &residueTypeColorTable_);
    void SetResidueSequenceColorTable(const std::string &residueSequenceColorTable_);
    void SetContinuousColorTable(const std::string &continuousColorTable_);
    void SetLegendFlag(bool legendFlag_);
    void SetMinFlag(bool minFlag_);
    void SetScalarMin(float scalarMin_);
    void SetMaxFlag(bool maxFlag_);
    void SetScalarMax(float scalarMax_);

    // Property getting methods
    AtomRenderingMode    GetDrawAtomsAs() const;
    RadiusType           GetScaleRadiusBy() const;
    BondRenderingMode    GetDrawBondsAs() const;
    BondColoringMode     GetColorBonds() const;
    const ColorAttribute &GetBondSingleColor() const;
          ColorAttribute &GetBondSingleColor();
    const std::string    &GetRadiusVariable() const;
          std::string    &GetRadiusVariable();
    float                GetRadiusScaleFactor() const;
    float                GetRadiusFixed() const;
    DetailLevel          GetAtomSphereQuality() const;
    DetailLevel          GetBondCylinderQuality() const;
    float                GetBondRadius() const;
    int                  GetBondLineWidth() const;
    int                  GetBondLineStyle() const;
    const std::string    &GetElementColorTable() const;
          std::string    &GetElementColorTable();
    const std::string    &GetResidueTypeColorTable() const;
          std::string    &GetResidueTypeColorTable();
    const std::string    &GetResidueSequenceColorTable() const;
          std::string    &GetResidueSequenceColorTable();
    const std::string    &GetContinuousColorTable() const;
          std::string    &GetContinuousColorTable();
    bool                 GetLegendFlag() const;
    bool                 GetMinFlag() const;
    float                GetScalarMin() const;
    bool                 GetMaxFlag() const;
    float                GetScalarMax() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string AtomRenderingMode_ToString(AtomRenderingMode);
    static bool AtomRenderingMode_FromString(const std::string &, AtomRenderingMode &);
protected:
    static std::string AtomRenderingMode_ToString(int);
public:
    static std::string RadiusType_ToString(RadiusType);
    static bool RadiusType_FromString(const std::string &, RadiusType &);
protected:
    static std::string RadiusType_ToString(int);
public:
    static std::string BondColoringMode_ToString(BondColoringMode);
    static bool BondColoringMode_FromString(const std::string &, BondColoringMode &);
protected:
    static std::string BondColoringMode_ToString(int);
public:
    static std::string BondRenderingMode_ToString(BondRenderingMode);
    static bool BondRenderingMode_FromString(const std::string &, BondRenderingMode &);
protected:
    static std::string BondRenderingMode_ToString(int);
public:
    static std::string DetailLevel_ToString(DetailLevel);
    static bool DetailLevel_FromString(const std::string &, DetailLevel &);
protected:
    static std::string DetailLevel_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const MoleculeAttributes &);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_drawAtomsAs = 0,
        ID_scaleRadiusBy,
        ID_drawBondsAs,
        ID_colorBonds,
        ID_bondSingleColor,
        ID_radiusVariable,
        ID_radiusScaleFactor,
        ID_radiusFixed,
        ID_atomSphereQuality,
        ID_bondCylinderQuality,
        ID_bondRadius,
        ID_bondLineWidth,
        ID_bondLineStyle,
        ID_elementColorTable,
        ID_residueTypeColorTable,
        ID_residueSequenceColorTable,
        ID_continuousColorTable,
        ID_legendFlag,
        ID_minFlag,
        ID_scalarMin,
        ID_maxFlag,
        ID_scalarMax
    };

private:
    int            drawAtomsAs;
    int            scaleRadiusBy;
    int            drawBondsAs;
    int            colorBonds;
    ColorAttribute bondSingleColor;
    std::string    radiusVariable;
    float          radiusScaleFactor;
    float          radiusFixed;
    int            atomSphereQuality;
    int            bondCylinderQuality;
    float          bondRadius;
    int            bondLineWidth;
    int            bondLineStyle;
    std::string    elementColorTable;
    std::string    residueTypeColorTable;
    std::string    residueSequenceColorTable;
    std::string    continuousColorTable;
    bool           legendFlag;
    bool           minFlag;
    float          scalarMin;
    bool           maxFlag;
    float          scalarMax;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
