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

#ifndef STATISTICALTRENDSATTRIBUTES_H
#define STATISTICALTRENDSATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: StatisticalTrendsAttributes
//
// Purpose:
//    This class contains attributes for the StatisticalTrends operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class StatisticalTrendsAttributes : public AttributeSubject
{
public:
    enum TrendTypeEnum
    {
        Absolute,
        Relative
    };
    enum StatisticTypeEnum
    {
        Sum,
        Mean,
        Variance,
        Slope,
        Residuals
    };
    enum TrendAxisEnum
    {
        Step,
        Time,
        Cycle
    };
    enum VariableSourceEnum
    {
        Default,
        OperatorExpression
    };

    // These constructors are for objects of this class
    StatisticalTrendsAttributes();
    StatisticalTrendsAttributes(const StatisticalTrendsAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    StatisticalTrendsAttributes(private_tmfs_t tmfs);
    StatisticalTrendsAttributes(const StatisticalTrendsAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~StatisticalTrendsAttributes();

    virtual StatisticalTrendsAttributes& operator = (const StatisticalTrendsAttributes &obj);
    virtual bool operator == (const StatisticalTrendsAttributes &obj) const;
    virtual bool operator != (const StatisticalTrendsAttributes &obj) const;
private:
    void Init();
    void Copy(const StatisticalTrendsAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetStartIndex(int startIndex_);
    void SetStopIndex(int stopIndex_);
    void SetStride(int stride_);
    void SetStartTrendType(TrendTypeEnum startTrendType_);
    void SetStopTrendType(TrendTypeEnum stopTrendType_);
    void SetStatisticType(StatisticTypeEnum statisticType_);
    void SetTrendAxis(TrendAxisEnum trendAxis_);
    void SetVariableSource(VariableSourceEnum variableSource_);

    // Property getting methods
    int GetStartIndex() const;
    int GetStopIndex() const;
    int GetStride() const;
    TrendTypeEnum GetStartTrendType() const;
    TrendTypeEnum GetStopTrendType() const;
    StatisticTypeEnum GetStatisticType() const;
    TrendAxisEnum GetTrendAxis() const;
    VariableSourceEnum GetVariableSource() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string TrendTypeEnum_ToString(TrendTypeEnum);
    static bool TrendTypeEnum_FromString(const std::string &, TrendTypeEnum &);
protected:
    static std::string TrendTypeEnum_ToString(int);
public:
    static std::string StatisticTypeEnum_ToString(StatisticTypeEnum);
    static bool StatisticTypeEnum_FromString(const std::string &, StatisticTypeEnum &);
protected:
    static std::string StatisticTypeEnum_ToString(int);
public:
    static std::string TrendAxisEnum_ToString(TrendAxisEnum);
    static bool TrendAxisEnum_FromString(const std::string &, TrendAxisEnum &);
protected:
    static std::string TrendAxisEnum_ToString(int);
public:
    static std::string VariableSourceEnum_ToString(VariableSourceEnum);
    static bool VariableSourceEnum_FromString(const std::string &, VariableSourceEnum &);
protected:
    static std::string VariableSourceEnum_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_startIndex = 0,
        ID_stopIndex,
        ID_stride,
        ID_startTrendType,
        ID_stopTrendType,
        ID_statisticType,
        ID_trendAxis,
        ID_variableSource,
        ID__LAST
    };

private:
    int startIndex;
    int stopIndex;
    int stride;
    int startTrendType;
    int stopTrendType;
    int statisticType;
    int trendAxis;
    int variableSource;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define STATISTICALTRENDSATTRIBUTES_TMFS "iiiiiiii"

#endif
