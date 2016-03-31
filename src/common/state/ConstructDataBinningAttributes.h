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

#ifndef CONSTRUCTDATABINNINGATTRIBUTES_H
#define CONSTRUCTDATABINNINGATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ConstructDataBinningAttributes
//
// Purpose:
//    Attributes for constructing a data binning
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ConstructDataBinningAttributes : public AttributeSubject
{
public:
    enum BinningScheme
    {
        Uniform,
        Unknown
    };
    enum ReductionOperator
    {
        Average,
        Minimum,
        Maximum,
        StandardDeviation,
        Variance,
        Sum,
        Count,
        RMS,
        PDF
    };
    enum OutOfBoundsBehavior
    {
        Clamp,
        Discard
    };
    enum BinType
    {
        Variable,
        X,
        Y,
        Z
    };

    // These constructors are for objects of this class
    ConstructDataBinningAttributes();
    ConstructDataBinningAttributes(const ConstructDataBinningAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    ConstructDataBinningAttributes(private_tmfs_t tmfs);
    ConstructDataBinningAttributes(const ConstructDataBinningAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~ConstructDataBinningAttributes();

    virtual ConstructDataBinningAttributes& operator = (const ConstructDataBinningAttributes &obj);
    virtual bool operator == (const ConstructDataBinningAttributes &obj) const;
    virtual bool operator != (const ConstructDataBinningAttributes &obj) const;
private:
    void Init();
    void Copy(const ConstructDataBinningAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectName();
    void SelectVarnames();
    void SelectBinType();
    void SelectBinBoundaries();
    void SelectVarForReductionOperator();
    void SelectNumBins();

    // Property setting methods
    void SetName(const std::string &name_);
    void SetVarnames(const stringVector &varnames_);
    void SetBinType(const unsignedCharVector &binType_);
    void SetBinBoundaries(const doubleVector &binBoundaries_);
    void SetReductionOperator(ReductionOperator reductionOperator_);
    void SetVarForReductionOperator(const std::string &varForReductionOperator_);
    void SetUndefinedValue(double undefinedValue_);
    void SetBinningScheme(BinningScheme binningScheme_);
    void SetNumBins(const intVector &numBins_);
    void SetOverTime(bool overTime_);
    void SetTimeStart(int timeStart_);
    void SetTimeEnd(int timeEnd_);
    void SetTimeStride(int timeStride_);
    void SetOutOfBoundsBehavior(OutOfBoundsBehavior outOfBoundsBehavior_);

    // Property getting methods
    const std::string        &GetName() const;
          std::string        &GetName();
    const stringVector       &GetVarnames() const;
          stringVector       &GetVarnames();
    const unsignedCharVector &GetBinType() const;
          unsignedCharVector &GetBinType();
    const doubleVector       &GetBinBoundaries() const;
          doubleVector       &GetBinBoundaries();
    ReductionOperator        GetReductionOperator() const;
    const std::string        &GetVarForReductionOperator() const;
          std::string        &GetVarForReductionOperator();
    double                   GetUndefinedValue() const;
    BinningScheme            GetBinningScheme() const;
    const intVector          &GetNumBins() const;
          intVector          &GetNumBins();
    bool                     GetOverTime() const;
    int                      GetTimeStart() const;
    int                      GetTimeEnd() const;
    int                      GetTimeStride() const;
    OutOfBoundsBehavior      GetOutOfBoundsBehavior() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string BinningScheme_ToString(BinningScheme);
    static bool BinningScheme_FromString(const std::string &, BinningScheme &);
protected:
    static std::string BinningScheme_ToString(int);
public:
    static std::string ReductionOperator_ToString(ReductionOperator);
    static bool ReductionOperator_FromString(const std::string &, ReductionOperator &);
protected:
    static std::string ReductionOperator_ToString(int);
public:
    static std::string OutOfBoundsBehavior_ToString(OutOfBoundsBehavior);
    static bool OutOfBoundsBehavior_FromString(const std::string &, OutOfBoundsBehavior &);
protected:
    static std::string OutOfBoundsBehavior_ToString(int);
public:
    static std::string BinType_ToString(BinType);
    static bool BinType_FromString(const std::string &, BinType &);
protected:
    static std::string BinType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ReductionRequiresVariable(void);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_name = 0,
        ID_varnames,
        ID_binType,
        ID_binBoundaries,
        ID_reductionOperator,
        ID_varForReductionOperator,
        ID_undefinedValue,
        ID_binningScheme,
        ID_numBins,
        ID_overTime,
        ID_timeStart,
        ID_timeEnd,
        ID_timeStride,
        ID_outOfBoundsBehavior,
        ID__LAST
    };

private:
    std::string        name;
    stringVector       varnames;
    unsignedCharVector binType;
    doubleVector       binBoundaries;
    int                reductionOperator;
    std::string        varForReductionOperator;
    double             undefinedValue;
    int                binningScheme;
    intVector          numBins;
    bool               overTime;
    int                timeStart;
    int                timeEnd;
    int                timeStride;
    int                outOfBoundsBehavior;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define CONSTRUCTDATABINNINGATTRIBUTES_TMFS "ss*u*d*isdii*biiii"

#endif
