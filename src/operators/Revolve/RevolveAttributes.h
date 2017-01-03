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

#ifndef REVOLVEATTRIBUTES_H
#define REVOLVEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: RevolveAttributes
//
// Purpose:
//    This class contains attributes for the revolve operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class RevolveAttributes : public AttributeSubject
{
public:
    enum MeshType
    {
        Auto,
        XY,
        RZ,
        ZR
    };

    // These constructors are for objects of this class
    RevolveAttributes();
    RevolveAttributes(const RevolveAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    RevolveAttributes(private_tmfs_t tmfs);
    RevolveAttributes(const RevolveAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~RevolveAttributes();

    virtual RevolveAttributes& operator = (const RevolveAttributes &obj);
    virtual bool operator == (const RevolveAttributes &obj) const;
    virtual bool operator != (const RevolveAttributes &obj) const;
private:
    void Init();
    void Copy(const RevolveAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectAxis();

    // Property setting methods
    void SetMeshType(MeshType meshType_);
    void SetAutoAxis(bool autoAxis_);
    void SetAxis(const double *axis_);
    void SetStartAngle(double startAngle_);
    void SetStopAngle(double stopAngle_);
    void SetSteps(int steps_);

    // Property getting methods
    MeshType     GetMeshType() const;
    bool         GetAutoAxis() const;
    const double *GetAxis() const;
          double *GetAxis();
    double       GetStartAngle() const;
    double       GetStopAngle() const;
    int          GetSteps() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string MeshType_ToString(MeshType);
    static bool MeshType_FromString(const std::string &, MeshType &);
protected:
    static std::string MeshType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_meshType = 0,
        ID_autoAxis,
        ID_axis,
        ID_startAngle,
        ID_stopAngle,
        ID_steps,
        ID__LAST
    };

private:
    int    meshType;
    bool   autoAxis;
    double axis[3];
    double startAngle;
    double stopAngle;
    int    steps;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define REVOLVEATTRIBUTES_TMFS "ibDddi"

#endif
