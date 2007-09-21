/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef PICKVARINFO_H
#define PICKVARINFO_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
#include "snprintf.h"

// ****************************************************************************
// Class: PickVarInfo
//
// Purpose:
//    This class contains PickVarInfo.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Mon Sep 17 16:05:41 PST 2007
//
// Modifications:
//   
// ****************************************************************************

class STATE_API PickVarInfo : public AttributeSubject
{
public:
    enum Centering
    {
        Nodal,
        Zonal,
        None
    };

    PickVarInfo();
    PickVarInfo(const PickVarInfo &obj);
    virtual ~PickVarInfo();

    virtual PickVarInfo& operator = (const PickVarInfo &obj);
    virtual bool operator == (const PickVarInfo &obj) const;
    virtual bool operator != (const PickVarInfo &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVariableName();
    void SelectVariableType();
    void SelectNames();
    void SelectValues();
    void SelectMixNames();
    void SelectMixValues();
    void SelectMiscMessage();
    void SelectNumMatsPerZone();
    void SelectMatNames();
    void SelectNumSpecsPerMat();
    void SelectFloatFormat();

    // Property setting methods
    void SetVariableName(const std::string &variableName_);
    void SetVariableType(const std::string &variableType_);
    void SetNames(const stringVector &names_);
    void SetValues(const doubleVector &values_);
    void SetMixNames(const stringVector &mixNames_);
    void SetMixValues(const doubleVector &mixValues_);
    void SetMixVar(bool mixVar_);
    void SetCentering(Centering centering_);
    void SetMiscMessage(const std::string &miscMessage_);
    void SetNumMatsPerZone(const intVector &numMatsPerZone_);
    void SetMatNames(const stringVector &matNames_);
    void SetNumSpecsPerMat(const intVector &numSpecsPerMat_);
    void SetTreatAsASCII(bool treatAsASCII_);
    void SetFloatFormat(const std::string &floatFormat_);

    // Property getting methods
    const std::string  &GetVariableName() const;
          std::string  &GetVariableName();
    const std::string  &GetVariableType() const;
          std::string  &GetVariableType();
    const stringVector &GetNames() const;
          stringVector &GetNames();
    const doubleVector &GetValues() const;
          doubleVector &GetValues();
    const stringVector &GetMixNames() const;
          stringVector &GetMixNames();
    const doubleVector &GetMixValues() const;
          doubleVector &GetMixValues();
    bool               GetMixVar() const;
    Centering          GetCentering() const;
    const std::string  &GetMiscMessage() const;
          std::string  &GetMiscMessage();
    const intVector    &GetNumMatsPerZone() const;
          intVector    &GetNumMatsPerZone();
    const stringVector &GetMatNames() const;
          stringVector &GetMatNames();
    const intVector    &GetNumSpecsPerMat() const;
          intVector    &GetNumSpecsPerMat();
    bool               GetTreatAsASCII() const;
    const std::string  &GetFloatFormat() const;
          std::string  &GetFloatFormat();

    // Enum conversion functions
    static std::string Centering_ToString(Centering);
    static bool Centering_FromString(const std::string &, Centering &);
protected:
    static std::string Centering_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void PrintSelf(ostream &os);
    void CreateOutputString(std::string &, const std::string &);
    void PrintTensor(std::string &, const std::vector<double> &, int, int, int);
    void PrintSymmetricTensor(std::string &, const std::vector<double> &, int, int, int);
    bool HasInfo(void);
    void PrintArray(std::string &, const std::vector<double> &, int, int, int);
private:
    std::string  variableName;
    std::string  variableType;
    stringVector names;
    doubleVector values;
    stringVector mixNames;
    doubleVector mixValues;
    bool         mixVar;
    int          centering;
    std::string  miscMessage;
    intVector    numMatsPerZone;
    stringVector matNames;
    intVector    numSpecsPerMat;
    bool         treatAsASCII;
    std::string  floatFormat;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
