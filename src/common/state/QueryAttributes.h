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

#ifndef QUERYATTRIBUTES_H
#define QUERYATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

#include <visitstream.h>

// ****************************************************************************
// Class: QueryAttributes
//
// Purpose:
//    This class contains attributes used for query.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API QueryAttributes : public AttributeSubject
{
public:
    enum VarType
    {
        Mesh,
        Scalar,
        Vector,
        Tensor,
        Symmetric_Tensor,
        Array,
        Label,
        Material,
        Species,
        Curve,
        Unknown
    };

    // These constructors are for objects of this class
    QueryAttributes();
    QueryAttributes(const QueryAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    QueryAttributes(private_tmfs_t tmfs);
    QueryAttributes(const QueryAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~QueryAttributes();

    virtual QueryAttributes& operator = (const QueryAttributes &obj);
    virtual bool operator == (const QueryAttributes &obj) const;
    virtual bool operator != (const QueryAttributes &obj) const;
private:
    void Init();
    void Copy(const QueryAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectResultsMessage();
    void SelectResultsValue();
    void SelectVarTypes();
    void SelectXUnits();
    void SelectYUnits();
    void SelectFloatFormat();
    void SelectXmlResult();
    void SelectQueryInputParams();
    void SelectDefaultName();
    void SelectDefaultVars();

    // Property setting methods
    void SetResultsMessage(const std::string &resultsMessage_);
    void SetResultsValue(const doubleVector &resultsValue_);
    void SetTimeStep(int timeStep_);
    void SetVarTypes(const intVector &varTypes_);
    void SetPipeIndex(int pipeIndex_);
    void SetXUnits(const std::string &xUnits_);
    void SetYUnits(const std::string &yUnits_);
    void SetFloatFormat(const std::string &floatFormat_);
    void SetXmlResult(const std::string &xmlResult_);
    void SetSuppressOutput(bool suppressOutput_);
    void SetQueryInputParams(const MapNode &queryInputParams_);
    void SetDefaultName(const std::string &defaultName_);
    void SetDefaultVars(const stringVector &defaultVars_);

    // Property getting methods
    const std::string  &GetResultsMessage() const;
          std::string  &GetResultsMessage();
    const doubleVector &GetResultsValue() const;
          doubleVector &GetResultsValue();
    int                GetTimeStep() const;
    const intVector    &GetVarTypes() const;
          intVector    &GetVarTypes();
    int                GetPipeIndex() const;
    const std::string  &GetXUnits() const;
          std::string  &GetXUnits();
    const std::string  &GetYUnits() const;
          std::string  &GetYUnits();
    const std::string  &GetFloatFormat() const;
          std::string  &GetFloatFormat();
    const std::string  &GetXmlResult() const;
          std::string  &GetXmlResult();
    bool               GetSuppressOutput() const;
    const MapNode      &GetQueryInputParams() const;
          MapNode      &GetQueryInputParams();
    const std::string  &GetDefaultName() const;
          std::string  &GetDefaultName();
    const stringVector &GetDefaultVars() const;
          stringVector &GetDefaultVars();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string VarType_ToString(VarType);
    static bool VarType_FromString(const std::string &, VarType &);
protected:
    static std::string VarType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void SetVariables(const stringVector &variables_);
    const std::string &GetName() const;
    std::string &GetName();
    const stringVector &GetVariables() const;
    stringVector &GetVariables();
    void Reset();
    void PrintSelf(ostream &os);
    void SetResultsValue(const double);
    void SetResultsValues(const double*, const int);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_resultsMessage = 0,
        ID_resultsValue,
        ID_timeStep,
        ID_varTypes,
        ID_pipeIndex,
        ID_xUnits,
        ID_yUnits,
        ID_floatFormat,
        ID_xmlResult,
        ID_suppressOutput,
        ID_queryInputParams,
        ID_defaultName,
        ID_defaultVars,
        ID__LAST
    };

private:
    std::string  resultsMessage;
    doubleVector resultsValue;
    int          timeStep;
    intVector    varTypes;
    int          pipeIndex;
    std::string  xUnits;
    std::string  yUnits;
    std::string  floatFormat;
    std::string  xmlResult;
    bool         suppressOutput;
    MapNode      queryInputParams;
    std::string  defaultName;
    stringVector defaultVars;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define QUERYATTRIBUTES_TMFS "sd*ii*issssbmss*"

#endif
