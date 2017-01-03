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

#ifndef REFLECTATTRIBUTES_H
#define REFLECTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ReflectAttributes
//
// Purpose:
//    This class contains attributes for the reflect operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class ReflectAttributes : public AttributeSubject
{
public:
    enum Octant
    {
        PXPYPZ,
        NXPYPZ,
        PXNYPZ,
        NXNYPZ,
        PXPYNZ,
        NXPYNZ,
        PXNYNZ,
        NXNYNZ
    };

    // These constructors are for objects of this class
    ReflectAttributes();
    ReflectAttributes(const ReflectAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    ReflectAttributes(private_tmfs_t tmfs);
    ReflectAttributes(const ReflectAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~ReflectAttributes();

    virtual ReflectAttributes& operator = (const ReflectAttributes &obj);
    virtual bool operator == (const ReflectAttributes &obj) const;
    virtual bool operator != (const ReflectAttributes &obj) const;
private:
    void Init();
    void Copy(const ReflectAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectReflections();

    // Property setting methods
    void SetOctant(Octant octant_);
    void SetUseXBoundary(bool useXBoundary_);
    void SetSpecifiedX(double specifiedX_);
    void SetUseYBoundary(bool useYBoundary_);
    void SetSpecifiedY(double specifiedY_);
    void SetUseZBoundary(bool useZBoundary_);
    void SetSpecifiedZ(double specifiedZ_);
    void SetReflections(const int *reflections_);

    // Property getting methods
    Octant    GetOctant() const;
    bool      GetUseXBoundary() const;
    double    GetSpecifiedX() const;
    bool      GetUseYBoundary() const;
    double    GetSpecifiedY() const;
    bool      GetUseZBoundary() const;
    double    GetSpecifiedZ() const;
    const int *GetReflections() const;
          int *GetReflections();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Octant_ToString(Octant);
    static bool Octant_FromString(const std::string &, Octant &);
protected:
    static std::string Octant_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual void ProcessOldVersions(DataNode *node, const char *configVersion);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_octant = 0,
        ID_useXBoundary,
        ID_specifiedX,
        ID_useYBoundary,
        ID_specifiedY,
        ID_useZBoundary,
        ID_specifiedZ,
        ID_reflections,
        ID__LAST
    };

private:
    int    octant;
    bool   useXBoundary;
    double specifiedX;
    bool   useYBoundary;
    double specifiedY;
    bool   useZBoundary;
    double specifiedZ;
    int    reflections[8];

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define REFLECTATTRIBUTES_TMFS "ibdbdbdI"

#endif
