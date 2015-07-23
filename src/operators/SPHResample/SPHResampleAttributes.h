/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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

#ifndef SPHRESAMPLEATTRIBUTES_H
#define SPHRESAMPLEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: SPHResampleAttributes
//
// Purpose:
//    
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class SPHResampleAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    SPHResampleAttributes();
    SPHResampleAttributes(const SPHResampleAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    SPHResampleAttributes(private_tmfs_t tmfs);
    SPHResampleAttributes(const SPHResampleAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~SPHResampleAttributes();

    virtual SPHResampleAttributes& operator = (const SPHResampleAttributes &obj);
    virtual bool operator == (const SPHResampleAttributes &obj) const;
    virtual bool operator != (const SPHResampleAttributes &obj) const;
private:
    void Init();
    void Copy(const SPHResampleAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectTensorSupportVariable();
    void SelectWeightVariable();

    // Property setting methods
    void SetMinX(float minX_);
    void SetMaxX(float maxX_);
    void SetXnum(int xnum_);
    void SetMinY(float minY_);
    void SetMaxY(float maxY_);
    void SetYnum(int ynum_);
    void SetMinZ(float minZ_);
    void SetMaxZ(float maxZ_);
    void SetZnum(int znum_);
    void SetTensorSupportVariable(const std::string &tensorSupportVariable_);
    void SetWeightVariable(const std::string &weightVariable_);
    void SetRK(bool RK_);

    // Property getting methods
    float             GetMinX() const;
    float             GetMaxX() const;
    int               GetXnum() const;
    float             GetMinY() const;
    float             GetMaxY() const;
    int               GetYnum() const;
    float             GetMinZ() const;
    float             GetMaxZ() const;
    int               GetZnum() const;
    const std::string &GetTensorSupportVariable() const;
          std::string &GetTensorSupportVariable();
    const std::string &GetWeightVariable() const;
          std::string &GetWeightVariable();
    bool              GetRK() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_minX = 0,
        ID_maxX,
        ID_xnum,
        ID_minY,
        ID_maxY,
        ID_ynum,
        ID_minZ,
        ID_maxZ,
        ID_znum,
        ID_tensorSupportVariable,
        ID_weightVariable,
        ID_RK,
        ID__LAST
    };

private:
    float       minX;
    float       maxX;
    int         xnum;
    float       minY;
    float       maxY;
    int         ynum;
    float       minZ;
    float       maxZ;
    int         znum;
    std::string tensorSupportVariable;
    std::string weightVariable;
    bool        RK;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define SPHRESAMPLEATTRIBUTES_TMFS "ffiffiffissb"

#endif
