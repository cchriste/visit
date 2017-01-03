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

#ifndef PDFATTRIBUTES_H
#define PDFATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: PDFAttributes
//
// Purpose:
//    Attributes for the PDF operator
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class PDFAttributes : public AttributeSubject
{
public:
    enum Scaling
    {
        Linear,
        Log,
        Skew
    };
    enum NumAxes
    {
        Two,
        Three
    };
    enum DensityType
    {
        Probability,
        ZoneCount
    };

    // These constructors are for objects of this class
    PDFAttributes();
    PDFAttributes(const PDFAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    PDFAttributes(private_tmfs_t tmfs);
    PDFAttributes(const PDFAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~PDFAttributes();

    virtual PDFAttributes& operator = (const PDFAttributes &obj);
    virtual bool operator == (const PDFAttributes &obj) const;
    virtual bool operator != (const PDFAttributes &obj) const;
private:
    void Init();
    void Copy(const PDFAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVar1();
    void SelectVar2();
    void SelectVar3();

    // Property setting methods
    void SetVar1(const std::string &var1_);
    void SetVar1MinFlag(bool var1MinFlag_);
    void SetVar1MaxFlag(bool var1MaxFlag_);
    void SetVar1Min(double var1Min_);
    void SetVar1Max(double var1Max_);
    void SetVar1Scaling(Scaling var1Scaling_);
    void SetVar1SkewFactor(double var1SkewFactor_);
    void SetVar1NumSamples(int var1NumSamples_);
    void SetVar2(const std::string &var2_);
    void SetVar2MinFlag(bool var2MinFlag_);
    void SetVar2MaxFlag(bool var2MaxFlag_);
    void SetVar2Min(double var2Min_);
    void SetVar2Max(double var2Max_);
    void SetVar2Scaling(Scaling var2Scaling_);
    void SetVar2SkewFactor(double var2SkewFactor_);
    void SetVar2NumSamples(int var2NumSamples_);
    void SetNumAxes(NumAxes numAxes_);
    void SetVar3(const std::string &var3_);
    void SetVar3MinFlag(bool var3MinFlag_);
    void SetVar3MaxFlag(bool var3MaxFlag_);
    void SetVar3Min(double var3Min_);
    void SetVar3Max(double var3Max_);
    void SetVar3Scaling(Scaling var3Scaling_);
    void SetVar3SkewFactor(double var3SkewFactor_);
    void SetVar3NumSamples(int var3NumSamples_);
    void SetScaleCube(bool scaleCube_);
    void SetDensityType(DensityType densityType_);

    // Property getting methods
    const std::string &GetVar1() const;
          std::string &GetVar1();
    bool              GetVar1MinFlag() const;
    bool              GetVar1MaxFlag() const;
    double            GetVar1Min() const;
    double            GetVar1Max() const;
    Scaling           GetVar1Scaling() const;
    double            GetVar1SkewFactor() const;
    int               GetVar1NumSamples() const;
    const std::string &GetVar2() const;
          std::string &GetVar2();
    bool              GetVar2MinFlag() const;
    bool              GetVar2MaxFlag() const;
    double            GetVar2Min() const;
    double            GetVar2Max() const;
    Scaling           GetVar2Scaling() const;
    double            GetVar2SkewFactor() const;
    int               GetVar2NumSamples() const;
    NumAxes           GetNumAxes() const;
    const std::string &GetVar3() const;
          std::string &GetVar3();
    bool              GetVar3MinFlag() const;
    bool              GetVar3MaxFlag() const;
    double            GetVar3Min() const;
    double            GetVar3Max() const;
    Scaling           GetVar3Scaling() const;
    double            GetVar3SkewFactor() const;
    int               GetVar3NumSamples() const;
    bool              GetScaleCube() const;
    DensityType       GetDensityType() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Scaling_ToString(Scaling);
    static bool Scaling_FromString(const std::string &, Scaling &);
protected:
    static std::string Scaling_ToString(int);
public:
    static std::string NumAxes_ToString(NumAxes);
    static bool NumAxes_FromString(const std::string &, NumAxes &);
protected:
    static std::string NumAxes_ToString(int);
public:
    static std::string DensityType_ToString(DensityType);
    static bool DensityType_FromString(const std::string &, DensityType &);
protected:
    static std::string DensityType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_var1 = 0,
        ID_var1MinFlag,
        ID_var1MaxFlag,
        ID_var1Min,
        ID_var1Max,
        ID_var1Scaling,
        ID_var1SkewFactor,
        ID_var1NumSamples,
        ID_var2,
        ID_var2MinFlag,
        ID_var2MaxFlag,
        ID_var2Min,
        ID_var2Max,
        ID_var2Scaling,
        ID_var2SkewFactor,
        ID_var2NumSamples,
        ID_numAxes,
        ID_var3,
        ID_var3MinFlag,
        ID_var3MaxFlag,
        ID_var3Min,
        ID_var3Max,
        ID_var3Scaling,
        ID_var3SkewFactor,
        ID_var3NumSamples,
        ID_scaleCube,
        ID_densityType,
        ID__LAST
    };

private:
    std::string var1;
    bool        var1MinFlag;
    bool        var1MaxFlag;
    double      var1Min;
    double      var1Max;
    int         var1Scaling;
    double      var1SkewFactor;
    int         var1NumSamples;
    std::string var2;
    bool        var2MinFlag;
    bool        var2MaxFlag;
    double      var2Min;
    double      var2Max;
    int         var2Scaling;
    double      var2SkewFactor;
    int         var2NumSamples;
    int         numAxes;
    std::string var3;
    bool        var3MinFlag;
    bool        var3MaxFlag;
    double      var3Min;
    double      var3Max;
    int         var3Scaling;
    double      var3SkewFactor;
    int         var3NumSamples;
    bool        scaleCube;
    int         densityType;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define PDFATTRIBUTES_TMFS "sbbddidisbbddidiisbbddidibi"

#endif
