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

#ifndef SCATTERATTRIBUTES_H
#define SCATTERATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: ScatterAttributes
//
// Purpose:
//    Attributes for the scatter plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Mon Dec 17 12:11:50 PDT 2007
//
// Modifications:
//   
// ****************************************************************************

class ScatterAttributes : public AttributeSubject
{
public:
    enum Scaling
    {
        Linear,
        Log,
        Skew
    };
    enum PointType
    {
        Box,
        Axis,
        Icosahedron,
        Point,
        Sphere
    };
    enum VariableRole
    {
        Coordinate0,
        Coordinate1,
        Coordinate2,
        Color,
        None
    };

    ScatterAttributes();
    ScatterAttributes(const ScatterAttributes &obj);
    virtual ~ScatterAttributes();

    virtual ScatterAttributes& operator = (const ScatterAttributes &obj);
    virtual bool operator == (const ScatterAttributes &obj) const;
    virtual bool operator != (const ScatterAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVar2();
    void SelectVar3();
    void SelectVar4();
    void SelectColorTableName();
    void SelectSingleColor();

    // Property setting methods
    void SetVar1Role(VariableRole var1Role_);
    void SetVar1MinFlag(bool var1MinFlag_);
    void SetVar1MaxFlag(bool var1MaxFlag_);
    void SetVar1Min(double var1Min_);
    void SetVar1Max(double var1Max_);
    void SetVar1Scaling(Scaling var1Scaling_);
    void SetVar1SkewFactor(double var1SkewFactor_);
    void SetVar2Role(VariableRole var2Role_);
    void SetVar2(const std::string &var2_);
    void SetVar2MinFlag(bool var2MinFlag_);
    void SetVar2MaxFlag(bool var2MaxFlag_);
    void SetVar2Min(double var2Min_);
    void SetVar2Max(double var2Max_);
    void SetVar2Scaling(Scaling var2Scaling_);
    void SetVar2SkewFactor(double var2SkewFactor_);
    void SetVar3Role(VariableRole var3Role_);
    void SetVar3(const std::string &var3_);
    void SetVar3MinFlag(bool var3MinFlag_);
    void SetVar3MaxFlag(bool var3MaxFlag_);
    void SetVar3Min(double var3Min_);
    void SetVar3Max(double var3Max_);
    void SetVar3Scaling(Scaling var3Scaling_);
    void SetVar3SkewFactor(double var3SkewFactor_);
    void SetVar4Role(VariableRole var4Role_);
    void SetVar4(const std::string &var4_);
    void SetVar4MinFlag(bool var4MinFlag_);
    void SetVar4MaxFlag(bool var4MaxFlag_);
    void SetVar4Min(double var4Min_);
    void SetVar4Max(double var4Max_);
    void SetVar4Scaling(Scaling var4Scaling_);
    void SetVar4SkewFactor(double var4SkewFactor_);
    void SetPointSize(double pointSize_);
    void SetPointSizePixels(int pointSizePixels_);
    void SetPointType(PointType pointType_);
    void SetScaleCube(bool scaleCube_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetForegroundFlag(bool foregroundFlag_);
    void SetLegendFlag(bool legendFlag_);

    // Property getting methods
    VariableRole         GetVar1Role() const;
    bool                 GetVar1MinFlag() const;
    bool                 GetVar1MaxFlag() const;
    double               GetVar1Min() const;
    double               GetVar1Max() const;
    Scaling              GetVar1Scaling() const;
    double               GetVar1SkewFactor() const;
    VariableRole         GetVar2Role() const;
    const std::string    &GetVar2() const;
          std::string    &GetVar2();
    bool                 GetVar2MinFlag() const;
    bool                 GetVar2MaxFlag() const;
    double               GetVar2Min() const;
    double               GetVar2Max() const;
    Scaling              GetVar2Scaling() const;
    double               GetVar2SkewFactor() const;
    VariableRole         GetVar3Role() const;
    const std::string    &GetVar3() const;
          std::string    &GetVar3();
    bool                 GetVar3MinFlag() const;
    bool                 GetVar3MaxFlag() const;
    double               GetVar3Min() const;
    double               GetVar3Max() const;
    Scaling              GetVar3Scaling() const;
    double               GetVar3SkewFactor() const;
    VariableRole         GetVar4Role() const;
    const std::string    &GetVar4() const;
          std::string    &GetVar4();
    bool                 GetVar4MinFlag() const;
    bool                 GetVar4MaxFlag() const;
    double               GetVar4Min() const;
    double               GetVar4Max() const;
    Scaling              GetVar4Scaling() const;
    double               GetVar4SkewFactor() const;
    double               GetPointSize() const;
    int                  GetPointSizePixels() const;
    PointType            GetPointType() const;
    bool                 GetScaleCube() const;
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    const ColorAttribute &GetSingleColor() const;
          ColorAttribute &GetSingleColor();
    bool                 GetForegroundFlag() const;
    bool                 GetLegendFlag() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Scaling_ToString(Scaling);
    static bool Scaling_FromString(const std::string &, Scaling &);
protected:
    static std::string Scaling_ToString(int);
public:
    static std::string PointType_ToString(PointType);
    static bool PointType_FromString(const std::string &, PointType &);
protected:
    static std::string PointType_ToString(int);
public:
    static std::string VariableRole_ToString(VariableRole);
    static bool VariableRole_FromString(const std::string &, VariableRole &);
protected:
    static std::string VariableRole_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const ScatterAttributes &) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_var1Role,
        ID_var1MinFlag,
        ID_var1MaxFlag,
        ID_var1Min,
        ID_var1Max,
        ID_var1Scaling,
        ID_var1SkewFactor,
        ID_var2Role,
        ID_var2,
        ID_var2MinFlag,
        ID_var2MaxFlag,
        ID_var2Min,
        ID_var2Max,
        ID_var2Scaling,
        ID_var2SkewFactor,
        ID_var3Role,
        ID_var3,
        ID_var3MinFlag,
        ID_var3MaxFlag,
        ID_var3Min,
        ID_var3Max,
        ID_var3Scaling,
        ID_var3SkewFactor,
        ID_var4Role,
        ID_var4,
        ID_var4MinFlag,
        ID_var4MaxFlag,
        ID_var4Min,
        ID_var4Max,
        ID_var4Scaling,
        ID_var4SkewFactor,
        ID_pointSize,
        ID_pointSizePixels,
        ID_pointType,
        ID_scaleCube,
        ID_colorTableName,
        ID_singleColor,
        ID_foregroundFlag,
        ID_legendFlag
    };

private:
    int            var1Role;
    bool           var1MinFlag;
    bool           var1MaxFlag;
    double         var1Min;
    double         var1Max;
    int            var1Scaling;
    double         var1SkewFactor;
    int            var2Role;
    std::string    var2;
    bool           var2MinFlag;
    bool           var2MaxFlag;
    double         var2Min;
    double         var2Max;
    int            var2Scaling;
    double         var2SkewFactor;
    int            var3Role;
    std::string    var3;
    bool           var3MinFlag;
    bool           var3MaxFlag;
    double         var3Min;
    double         var3Max;
    int            var3Scaling;
    double         var3SkewFactor;
    int            var4Role;
    std::string    var4;
    bool           var4MinFlag;
    bool           var4MaxFlag;
    double         var4Min;
    double         var4Max;
    int            var4Scaling;
    double         var4SkewFactor;
    double         pointSize;
    int            pointSizePixels;
    int            pointType;
    bool           scaleCube;
    std::string    colorTableName;
    ColorAttribute singleColor;
    bool           foregroundFlag;
    bool           legendFlag;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
