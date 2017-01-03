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

#ifndef LIGHTLIST_H
#define LIGHTLIST_H
#include <state_exports.h>
#include <AttributeSubject.h>

#include <LightAttributes.h>

// ****************************************************************************
// Class: LightList
//
// Purpose:
//    This class contains a list of LightAttributes.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API LightList : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    LightList();
    LightList(const LightList &obj);
protected:
    // These constructors are for objects derived from this class
    LightList(private_tmfs_t tmfs);
    LightList(const LightList &obj, private_tmfs_t tmfs);
public:
    virtual ~LightList();

    virtual LightList& operator = (const LightList &obj);
    virtual bool operator == (const LightList &obj) const;
    virtual bool operator != (const LightList &obj) const;
private:
    void Init();
    void Copy(const LightList &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectLight0();
    void SelectLight1();
    void SelectLight2();
    void SelectLight3();
    void SelectLight4();
    void SelectLight5();
    void SelectLight6();
    void SelectLight7();

    // Property setting methods
    void SetLight0(const LightAttributes &light0_);
    void SetLight1(const LightAttributes &light1_);
    void SetLight2(const LightAttributes &light2_);
    void SetLight3(const LightAttributes &light3_);
    void SetLight4(const LightAttributes &light4_);
    void SetLight5(const LightAttributes &light5_);
    void SetLight6(const LightAttributes &light6_);
    void SetLight7(const LightAttributes &light7_);

    // Property getting methods
    const LightAttributes &GetLight0() const;
          LightAttributes &GetLight0();
    const LightAttributes &GetLight1() const;
          LightAttributes &GetLight1();
    const LightAttributes &GetLight2() const;
          LightAttributes &GetLight2();
    const LightAttributes &GetLight3() const;
          LightAttributes &GetLight3();
    const LightAttributes &GetLight4() const;
          LightAttributes &GetLight4();
    const LightAttributes &GetLight5() const;
          LightAttributes &GetLight5();
    const LightAttributes &GetLight6() const;
          LightAttributes &GetLight6();
    const LightAttributes &GetLight7() const;
          LightAttributes &GetLight7();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    const LightAttributes &GetLight(int i) const;
    LightAttributes &GetLight(int i);
    int NumLights() const;
    void SelectLight(int i);
    void SetLight(int i, const LightAttributes &l);
    void SetDefaultEnabledStates();

    // IDs that can be used to identify fields in case statements
    enum {
        ID_light0 = 0,
        ID_light1,
        ID_light2,
        ID_light3,
        ID_light4,
        ID_light5,
        ID_light6,
        ID_light7,
        ID__LAST
    };

private:
    LightAttributes light0;
    LightAttributes light1;
    LightAttributes light2;
    LightAttributes light3;
    LightAttributes light4;
    LightAttributes light5;
    LightAttributes light6;
    LightAttributes light7;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define LIGHTLIST_TMFS "aaaaaaaa"

#endif
