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

#ifndef FTLEATTRIBUTES_H
#define FTLEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: FTLEAttributes
//
// Purpose:
//    Attributes for FTLE
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class FTLEAttributes : public AttributeSubject
{
public:
    enum Region
    {
        NativeResolutionOfMesh,
        RegularGrid
    };
    enum Direction
    {
        Forward,
        Backward
    };
    enum FlowType
    {
        Unsteady,
        Steady
    };

    // These constructors are for objects of this class
    FTLEAttributes();
    FTLEAttributes(const FTLEAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    FTLEAttributes(private_tmfs_t tmfs);
    FTLEAttributes(const FTLEAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~FTLEAttributes();

    virtual FTLEAttributes& operator = (const FTLEAttributes &obj);
    virtual bool operator == (const FTLEAttributes &obj) const;
    virtual bool operator != (const FTLEAttributes &obj) const;
private:
    void Init();
    void Copy(const FTLEAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectResolution();
    void SelectStartPosition();
    void SelectEndPosition();

    // Property setting methods
    void SetIntegrationTime(double integrationTime_);
    void SetRegionType(Region regionType_);
    void SetResolution(const int *Resolution_);
    void SetUseDataSetStart(bool UseDataSetStart_);
    void SetStartPosition(const double *StartPosition_);
    void SetUseDataSetEnd(bool UseDataSetEnd_);
    void SetEndPosition(const double *EndPosition_);
    void SetDirection(Direction direction_);
    void SetFlowType(FlowType flowType_);

    // Property getting methods
    double       GetIntegrationTime() const;
    Region       GetRegionType() const;
    const int    *GetResolution() const;
          int    *GetResolution();
    bool         GetUseDataSetStart() const;
    const double *GetStartPosition() const;
          double *GetStartPosition();
    bool         GetUseDataSetEnd() const;
    const double *GetEndPosition() const;
          double *GetEndPosition();
    Direction    GetDirection() const;
    FlowType     GetFlowType() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Region_ToString(Region);
    static bool Region_FromString(const std::string &, Region &);
protected:
    static std::string Region_ToString(int);
public:
    static std::string Direction_ToString(Direction);
    static bool Direction_FromString(const std::string &, Direction &);
protected:
    static std::string Direction_ToString(int);
public:
    static std::string FlowType_ToString(FlowType);
    static bool FlowType_FromString(const std::string &, FlowType &);
protected:
    static std::string FlowType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_integrationTime = 0,
        ID_regionType,
        ID_Resolution,
        ID_UseDataSetStart,
        ID_StartPosition,
        ID_UseDataSetEnd,
        ID_EndPosition,
        ID_direction,
        ID_flowType,
        ID__LAST
    };

private:
    double integrationTime;
    int    regionType;
    int    Resolution[3];
    bool   UseDataSetStart;
    double StartPosition[3];
    bool   UseDataSetEnd;
    double EndPosition[3];
    int    direction;
    int    flowType;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define FTLEATTRIBUTES_TMFS "diIbDbDii"

#endif
