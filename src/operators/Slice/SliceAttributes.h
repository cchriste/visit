/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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

#ifndef SLICEATTRIBUTES_H
#define SLICEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <cmath>

// ****************************************************************************
// Class: SliceAttributes
//
// Purpose:
//    This class contains attributes for the arbitrary slice.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class SliceAttributes : public AttributeSubject
{
public:
    enum AxisType
    {
        XAxis,
        YAxis,
        ZAxis,
        Arbitrary,
        ThetaPhi
    };
    enum OriginType
    {
        Point,
        Intercept,
        Percent,
        Zone,
        Node
    };

    SliceAttributes();
    SliceAttributes(const SliceAttributes &obj);
    virtual ~SliceAttributes();

    virtual SliceAttributes& operator = (const SliceAttributes &obj);
    virtual bool operator == (const SliceAttributes &obj) const;
    virtual bool operator != (const SliceAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectOriginPoint();
    void SelectNormal();
    void SelectUpAxis();
    void SelectMeshName();

    // Property setting methods
    void SetOriginType(OriginType originType_);
    void SetOriginPoint(const double *originPoint_);
    void SetOriginIntercept(double originIntercept_);
    void SetOriginPercent(double originPercent_);
    void SetOriginZone(int originZone_);
    void SetOriginNode(int originNode_);
    void SetNormal(const double *normal_);
    void SetAxisType(AxisType axisType_);
    void SetUpAxis(const double *upAxis_);
    void SetProject2d(bool project2d_);
    void SetInteractive(bool interactive_);
    void SetFlip(bool flip_);
    void SetOriginZoneDomain(int originZoneDomain_);
    void SetOriginNodeDomain(int originNodeDomain_);
    void SetMeshName(const std::string &meshName_);
    void SetTheta(double theta_);
    void SetPhi(double phi_);

    // Property getting methods
    OriginType        GetOriginType() const;
    const double      *GetOriginPoint() const;
          double      *GetOriginPoint();
    double            GetOriginIntercept() const;
    double            GetOriginPercent() const;
    int               GetOriginZone() const;
    int               GetOriginNode() const;
    const double      *GetNormal() const;
          double      *GetNormal();
    AxisType          GetAxisType() const;
    const double      *GetUpAxis() const;
          double      *GetUpAxis();
    bool              GetProject2d() const;
    bool              GetInteractive() const;
    bool              GetFlip() const;
    int               GetOriginZoneDomain() const;
    int               GetOriginNodeDomain() const;
    const std::string &GetMeshName() const;
          std::string &GetMeshName();
    double            GetTheta() const;
    double            GetPhi() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string AxisType_ToString(AxisType);
    static bool AxisType_FromString(const std::string &, AxisType &);
protected:
    static std::string AxisType_ToString(int);
public:
    static std::string OriginType_ToString(OriginType);
    static bool OriginType_FromString(const std::string &, OriginType &);
protected:
    static std::string OriginType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void UpdateOrthogonalAxes();
    virtual bool EqualTo(const AttributeGroup *atts) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_originType = 0,
        ID_originPoint,
        ID_originIntercept,
        ID_originPercent,
        ID_originZone,
        ID_originNode,
        ID_normal,
        ID_axisType,
        ID_upAxis,
        ID_project2d,
        ID_interactive,
        ID_flip,
        ID_originZoneDomain,
        ID_originNodeDomain,
        ID_meshName,
        ID_theta,
        ID_phi
    };

private:
    int         originType;
    double      originPoint[3];
    double      originIntercept;
    double      originPercent;
    int         originZone;
    int         originNode;
    double      normal[3];
    int         axisType;
    double      upAxis[3];
    bool        project2d;
    bool        interactive;
    bool        flip;
    int         originZoneDomain;
    int         originNodeDomain;
    std::string meshName;
    double      theta;
    double      phi;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
