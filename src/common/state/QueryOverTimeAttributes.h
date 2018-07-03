/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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

#ifndef QUERYOVERTIMEATTRIBUTES_H
#define QUERYOVERTIMEATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

#include <QueryAttributes.h>
#include <PickAttributes.h>

// ****************************************************************************
// Class: QueryOverTimeAttributes
//
// Purpose:
//    Attributes for queries over time.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API QueryOverTimeAttributes : public AttributeSubject
{
public:
    enum TimeType
    {
        Cycle,
        DTime,
        Timestep
    };

    // These constructors are for objects of this class
    QueryOverTimeAttributes();
    QueryOverTimeAttributes(const QueryOverTimeAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    QueryOverTimeAttributes(private_tmfs_t tmfs);
    QueryOverTimeAttributes(const QueryOverTimeAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~QueryOverTimeAttributes();

    virtual QueryOverTimeAttributes& operator = (const QueryOverTimeAttributes &obj);
    virtual bool operator == (const QueryOverTimeAttributes &obj) const;
    virtual bool operator != (const QueryOverTimeAttributes &obj) const;
private:
    void Init();
    void Copy(const QueryOverTimeAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectQueryAtts();
    void SelectPickAtts();
    void SelectCachedCurvePts();

    // Property setting methods
    void SetTimeType(TimeType timeType_);
    void SetStartTimeFlag(bool startTimeFlag_);
    void SetStartTime(int startTime_);
    void SetEndTimeFlag(bool endTimeFlag_);
    void SetEndTime(int endTime_);
    void SetStrideFlag(bool strideFlag_);
    void SetStride(int stride_);
    void SetCreateWindow(bool createWindow_);
    void SetWindowId(int windowId_);
    void SetQueryAtts(const QueryAttributes &queryAtts_);
    void SetPickAtts(const PickAttributes &pickAtts_);
    void SetCachedCurvePts(const doubleVector &cachedCurvePts_);
    void SetUseCachedPts(bool useCachedPts_);

    // Property getting methods
    TimeType              GetTimeType() const;
    bool                  GetStartTimeFlag() const;
    int                   GetStartTime() const;
    bool                  GetEndTimeFlag() const;
    int                   GetEndTime() const;
    bool                  GetStrideFlag() const;
    int                   GetStride() const;
    bool                  GetCreateWindow() const;
    int                   GetWindowId() const;
    const QueryAttributes &GetQueryAtts() const;
          QueryAttributes &GetQueryAtts();
    const PickAttributes  &GetPickAtts() const;
          PickAttributes  &GetPickAtts();
    const doubleVector    &GetCachedCurvePts() const;
          doubleVector    &GetCachedCurvePts();
    bool                  GetUseCachedPts() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string TimeType_ToString(TimeType);
    static bool TimeType_FromString(const std::string &, TimeType &);
protected:
    static std::string TimeType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_timeType = 0,
        ID_startTimeFlag,
        ID_startTime,
        ID_endTimeFlag,
        ID_endTime,
        ID_strideFlag,
        ID_stride,
        ID_createWindow,
        ID_windowId,
        ID_queryAtts,
        ID_pickAtts,
        ID_cachedCurvePts,
        ID_useCachedPts,
        ID__LAST
    };

protected:
    doubleVector    cachedCurvePts;
    bool            useCachedPts;
private:
    int             timeType;
    bool            startTimeFlag;
    int             startTime;
    bool            endTimeFlag;
    int             endTime;
    bool            strideFlag;
    int             stride;
    bool            createWindow;
    int             windowId;
    QueryAttributes queryAtts;
    PickAttributes  pickAtts;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define QUERYOVERTIMEATTRIBUTES_TMFS "ibibibibiaad*b"

#endif
