/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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

#ifndef COLORTABLEATTRIBUTES_H
#define COLORTABLEATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

class ColorControlPointList;

// ****************************************************************************
// Class: ColorTableAttributes
//
// Purpose:
//    This class contains the list of colortables.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ColorTableAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    ColorTableAttributes();
    ColorTableAttributes(const ColorTableAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    ColorTableAttributes(private_tmfs_t tmfs);
    ColorTableAttributes(const ColorTableAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~ColorTableAttributes();

    virtual ColorTableAttributes& operator = (const ColorTableAttributes &obj);
    virtual bool operator == (const ColorTableAttributes &obj) const;
    virtual bool operator != (const ColorTableAttributes &obj) const;
private:
    void Init();
    void Copy(const ColorTableAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectNames();
    void SelectColorTables();
    void SelectActiveContinuous();
    void SelectActiveDiscrete();

    // Property setting methods
    void SetNames(const stringVector &names_);
    void SetActiveContinuous(const std::string &activeContinuous_);
    void SetActiveDiscrete(const std::string &activeDiscrete_);

    // Property getting methods
    const stringVector &GetNames() const;
          stringVector &GetNames();
    const AttributeGroupVector &GetColorTables() const;
          AttributeGroupVector &GetColorTables();
    const std::string  &GetActiveContinuous() const;
          std::string  &GetActiveContinuous();
    const std::string  &GetActiveDiscrete() const;
          std::string  &GetActiveDiscrete();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddColorTables(const ColorControlPointList &);
    void ClearColorTables();
    void RemoveColorTables(int i);
    int  GetNumColorTables() const;
    ColorControlPointList &GetColorTables(int i);
    const ColorControlPointList &GetColorTables(int i) const;

    ColorControlPointList &operator [] (int i);
    const ColorControlPointList &operator [] (int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    int GetColorTableIndex(const std::string &name) const;
    const ColorControlPointList *GetColorControlPoints(int index) const;
    const ColorControlPointList *GetColorControlPoints(const std::string &name) const;
    void AddColorTable(const std::string &name, const ColorControlPointList &cpts);
    void RemoveColorTable(const std::string &name);
    void RemoveColorTable(int index);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_names = 0,
        ID_colorTables,
        ID_activeContinuous,
        ID_activeDiscrete,
        ID__LAST
    };

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    stringVector         names;
    AttributeGroupVector colorTables;
    std::string          activeContinuous;
    std::string          activeDiscrete;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define COLORTABLEATTRIBUTES_TMFS "s*a*ss"

#endif
