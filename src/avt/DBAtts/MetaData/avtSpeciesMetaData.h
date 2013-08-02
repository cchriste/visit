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

#ifndef AVTSPECIESMETADATA_H
#define AVTSPECIESMETADATA_H
#include <dbatts_exports.h>
#include <string>
#include <AttributeSubject.h>

class avtMatSpeciesMetaData;

// ****************************************************************************
// Class: avtSpeciesMetaData
//
// Purpose:
//    Contains species metadata attributes
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class DBATTS_API avtSpeciesMetaData : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    avtSpeciesMetaData();
    avtSpeciesMetaData(const avtSpeciesMetaData &obj);
protected:
    // These constructors are for objects derived from this class
    avtSpeciesMetaData(private_tmfs_t tmfs);
    avtSpeciesMetaData(const avtSpeciesMetaData &obj, private_tmfs_t tmfs);
public:
    virtual ~avtSpeciesMetaData();

    virtual avtSpeciesMetaData& operator = (const avtSpeciesMetaData &obj);
    virtual bool operator == (const avtSpeciesMetaData &obj) const;
    virtual bool operator != (const avtSpeciesMetaData &obj) const;
private:
    void Init();
    void Copy(const avtSpeciesMetaData &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSpecies();

    // Property setting methods

    // Property getting methods
    const AttributeGroupVector &GetSpecies() const;
          AttributeGroupVector &GetSpecies();


    // Attributegroup convenience methods
    void AddSpecies(const avtMatSpeciesMetaData &);
    void ClearSpecies();
    void RemoveSpecies(int i);
    int  GetNumSpecies() const;
    avtMatSpeciesMetaData &GetSpecies(int i);
    const avtMatSpeciesMetaData &GetSpecies(int i) const;

    avtMatSpeciesMetaData &operator [] (int i);
    const avtMatSpeciesMetaData &operator [] (int i) const;

    // User-defined methods
    avtSpeciesMetaData(const std::string &n, const std::string &meshn, const std::string &matn, int nummat, const intVector &ns, const std::vector<stringVector> &sn);
    void Print(ostream &, int = 0) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_name = 0,
        ID_originalName,
        ID_validVariable,
        ID_meshName,
        ID_materialName,
        ID_numMaterials,
        ID_species,
        ID__LAST
    };

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
public:
    std::string          name;
    std::string          originalName;
    bool                 validVariable;
    std::string          meshName;
    std::string          materialName;
    int                  numMaterials;
private:
    AttributeGroupVector species;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define AVTSPECIESMETADATA_TMFS "ssbssia*"

#endif
