/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

#ifndef REPLICATEATTRIBUTES_H
#define REPLICATEATTRIBUTES_H
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ReplicateAttributes
//
// Purpose:
//    This class contains attributes for the replicate operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class ReplicateAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    ReplicateAttributes();
    ReplicateAttributes(const ReplicateAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    ReplicateAttributes(private_tmfs_t tmfs);
    ReplicateAttributes(const ReplicateAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~ReplicateAttributes();

    virtual ReplicateAttributes& operator = (const ReplicateAttributes &obj);
    virtual bool operator == (const ReplicateAttributes &obj) const;
    virtual bool operator != (const ReplicateAttributes &obj) const;
private:
    void Init();
    void Copy(const ReplicateAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectXVector();
    void SelectYVector();
    void SelectZVector();
    void SelectNewPeriodicOrigin();

    // Property setting methods
    void SetUseUnitCellVectors(bool useUnitCellVectors_);
    void SetXVector(const double *xVector_);
    void SetYVector(const double *yVector_);
    void SetZVector(const double *zVector_);
    void SetXReplications(int xReplications_);
    void SetYReplications(int yReplications_);
    void SetZReplications(int zReplications_);
    void SetMergeResults(bool mergeResults_);
    void SetReplicateUnitCellAtoms(bool replicateUnitCellAtoms_);
    void SetShiftPeriodicAtomOrigin(bool shiftPeriodicAtomOrigin_);
    void SetNewPeriodicOrigin(const double *newPeriodicOrigin_);

    // Property getting methods
    bool         GetUseUnitCellVectors() const;
    const double *GetXVector() const;
          double *GetXVector();
    const double *GetYVector() const;
          double *GetYVector();
    const double *GetZVector() const;
          double *GetZVector();
    int          GetXReplications() const;
    int          GetYReplications() const;
    int          GetZReplications() const;
    bool         GetMergeResults() const;
    bool         GetReplicateUnitCellAtoms() const;
    bool         GetShiftPeriodicAtomOrigin() const;
    const double *GetNewPeriodicOrigin() const;
          double *GetNewPeriodicOrigin();

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
        ID_useUnitCellVectors = 0,
        ID_xVector,
        ID_yVector,
        ID_zVector,
        ID_xReplications,
        ID_yReplications,
        ID_zReplications,
        ID_mergeResults,
        ID_replicateUnitCellAtoms,
        ID_shiftPeriodicAtomOrigin,
        ID_newPeriodicOrigin,
        ID__LAST
    };

private:
    bool   useUnitCellVectors;
    double xVector[3];
    double yVector[3];
    double zVector[3];
    int    xReplications;
    int    yReplications;
    int    zReplications;
    bool   mergeResults;
    bool   replicateUnitCellAtoms;
    bool   shiftPeriodicAtomOrigin;
    double newPeriodicOrigin[3];

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define REPLICATEATTRIBUTES_TMFS "bDDDiiibbbD"

#endif
