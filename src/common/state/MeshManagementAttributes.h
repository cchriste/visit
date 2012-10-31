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

#ifndef MESHMANAGEMENTATTRIBUTES_H
#define MESHMANAGEMENTATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: MeshManagementAttributes
//
// Purpose:
//    Global variables controlling reading and conversion of non-standard meshes
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API MeshManagementAttributes : public AttributeSubject
{
public:
    enum DiscretizationModes
    {
        Uniform,
        Adaptive,
        MultiPass
    };

    // These constructors are for objects of this class
    MeshManagementAttributes();
    MeshManagementAttributes(const MeshManagementAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    MeshManagementAttributes(private_tmfs_t tmfs);
    MeshManagementAttributes(const MeshManagementAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~MeshManagementAttributes();

    virtual MeshManagementAttributes& operator = (const MeshManagementAttributes &obj);
    virtual bool operator == (const MeshManagementAttributes &obj) const;
    virtual bool operator != (const MeshManagementAttributes &obj) const;
private:
    void Init();
    void Copy(const MeshManagementAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectDiscretizationTolerance();
    void SelectDiscretizationToleranceX();
    void SelectDiscretizationToleranceY();
    void SelectDiscretizationToleranceZ();

    // Property setting methods
    void SetDiscretizationTolerance(const doubleVector &discretizationTolerance_);
    void SetDiscretizationToleranceX(const doubleVector &discretizationToleranceX_);
    void SetDiscretizationToleranceY(const doubleVector &discretizationToleranceY_);
    void SetDiscretizationToleranceZ(const doubleVector &discretizationToleranceZ_);
    void SetDiscretizationMode(DiscretizationModes discretizationMode_);
    void SetDiscretizeBoundaryOnly(bool discretizeBoundaryOnly_);
    void SetPassNativeCSG(bool passNativeCSG_);

    // Property getting methods
    const doubleVector &GetDiscretizationTolerance() const;
          doubleVector &GetDiscretizationTolerance();
    const doubleVector &GetDiscretizationToleranceX() const;
          doubleVector &GetDiscretizationToleranceX();
    const doubleVector &GetDiscretizationToleranceY() const;
          doubleVector &GetDiscretizationToleranceY();
    const doubleVector &GetDiscretizationToleranceZ() const;
          doubleVector &GetDiscretizationToleranceZ();
    DiscretizationModes GetDiscretizationMode() const;
    bool               GetDiscretizeBoundaryOnly() const;
    bool               GetPassNativeCSG() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string DiscretizationModes_ToString(DiscretizationModes);
    static bool DiscretizationModes_FromString(const std::string &, DiscretizationModes &);
protected:
    static std::string DiscretizationModes_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_discretizationTolerance = 0,
        ID_discretizationToleranceX,
        ID_discretizationToleranceY,
        ID_discretizationToleranceZ,
        ID_discretizationMode,
        ID_discretizeBoundaryOnly,
        ID_passNativeCSG,
        ID__LAST
    };

private:
    doubleVector discretizationTolerance;
    doubleVector discretizationToleranceX;
    doubleVector discretizationToleranceY;
    doubleVector discretizationToleranceZ;
    int          discretizationMode;
    bool         discretizeBoundaryOnly;
    bool         passNativeCSG;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define MESHMANAGEMENTATTRIBUTES_TMFS "d*d*d*d*ibb"

#endif
