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

#ifndef CLIENTMETHOD_H
#define CLIENTMETHOD_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ClientMethod
//
// Purpose:
//    This class contains the attributes needed to make one client execute a method on another client.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ClientMethod : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    ClientMethod();
    ClientMethod(const ClientMethod &obj);
protected:
    // These constructors are for objects derived from this class
    ClientMethod(private_tmfs_t tmfs);
    ClientMethod(const ClientMethod &obj, private_tmfs_t tmfs);
public:
    virtual ~ClientMethod();

    virtual ClientMethod& operator = (const ClientMethod &obj);
    virtual bool operator == (const ClientMethod &obj) const;
    virtual bool operator != (const ClientMethod &obj) const;
private:
    void Init();
    void Copy(const ClientMethod &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectMethodName();
    void SelectIntArgs();
    void SelectDoubleArgs();
    void SelectStringArgs();

    // Property setting methods
    void SetMethodName(const std::string &methodName_);
    void SetIntArgs(const intVector &intArgs_);
    void SetDoubleArgs(const doubleVector &doubleArgs_);
    void SetStringArgs(const stringVector &stringArgs_);

    // Property getting methods
    const std::string  &GetMethodName() const;
          std::string  &GetMethodName();
    const intVector    &GetIntArgs() const;
          intVector    &GetIntArgs();
    const doubleVector &GetDoubleArgs() const;
          doubleVector &GetDoubleArgs();
    const stringVector &GetStringArgs() const;
          stringVector &GetStringArgs();


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void ClearArgs();
    void AddArgument(int);
    void AddArgument(double);
    void AddArgument(const std::string &);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_methodName = 0,
        ID_intArgs,
        ID_doubleArgs,
        ID_stringArgs,
        ID__LAST
    };

private:
    std::string  methodName;
    intVector    intArgs;
    doubleVector doubleArgs;
    stringVector stringArgs;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define CLIENTMETHOD_TMFS "si*d*s*"

#endif
