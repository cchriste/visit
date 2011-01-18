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

#ifndef AVTMESHMETADATA_H
#define AVTMESHMETADATA_H
#include <dbatts_exports.h>
#include <string>
#include <avtTypes.h>
#include <AttributeSubject.h>

#include <NameschemeAttributes.h>

// ****************************************************************************
// Class: avtMeshMetaData
//
// Purpose:
//    Contains mesh metadata attributes
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class DBATTS_API avtMeshMetaData : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    avtMeshMetaData();
    avtMeshMetaData(const avtMeshMetaData &obj);
protected:
    // These constructors are for objects derived from this class
    avtMeshMetaData(private_tmfs_t tmfs);
    avtMeshMetaData(const avtMeshMetaData &obj, private_tmfs_t tmfs);
public:
    virtual ~avtMeshMetaData();

    virtual avtMeshMetaData& operator = (const avtMeshMetaData &obj);
    virtual bool operator == (const avtMeshMetaData &obj) const;
    virtual bool operator != (const avtMeshMetaData &obj) const;
private:
    void Init();
    void Copy(const avtMeshMetaData &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // User-defined methods
    avtMeshMetaData(const double *, std::string, int, int, int, int, int, int, avtMeshType);
    avtMeshMetaData(std::string, int, int, int, int, int, int, avtMeshType);
    void SetExtents(const double *);
    void UnsetExtents();
    void Print(ostream &, int = 0) const;
    void SetAMRInfo(const std::string &levelName, const std::string &patchName, int origin,const std::vector<int> &patchesPerLevel);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_name = 0,
        ID_originalName,
        ID_validVariable,
        ID_meshType,
        ID_meshCoordType,
        ID_cellOrigin,
        ID_spatialDimension,
        ID_topologicalDimension,
        ID_xUnits,
        ID_yUnits,
        ID_zUnits,
        ID_xLabel,
        ID_yLabel,
        ID_zLabel,
        ID_hasSpatialExtents,
        ID_minSpatialExtents,
        ID_maxSpatialExtents,
        ID_numBlocks,
        ID_blockOrigin,
        ID_blockPieceName,
        ID_blockTitle,
        ID_blockNames,
        ID_blockNameScheme,
        ID_numGroups,
        ID_groupOrigin,
        ID_groupPieceName,
        ID_groupTitle,
        ID_groupIds,
        ID_groupIdsBasedOnRange,
        ID_disjointElements,
        ID_containsGhostZones,
        ID_containsOriginalCells,
        ID_containsOriginalNodes,
        ID_containsGlobalNodeIds,
        ID_containsGlobalZoneIds,
        ID_loadBalanceScheme,
        ID_nodesAreCritical,
        ID_unitCellVectors,
        ID_unitCellOrigin,
        ID_rectilinearGridHasTransform,
        ID_rectilinearGridTransform,
        ID_nodeOrigin,
        ID_containsExteriorBoundaryGhosts,
        ID_hideFromGUI,
        ID_LODs,
        ID__LAST
    };

public:
    std::string          name;
    std::string          originalName;
    bool                 validVariable;
    avtMeshType          meshType;
    avtMeshCoordType     meshCoordType;
    int                  cellOrigin;
    int                  spatialDimension;
    int                  topologicalDimension;
    std::string          xUnits;
    std::string          yUnits;
    std::string          zUnits;
    std::string          xLabel;
    std::string          yLabel;
    std::string          zLabel;
    bool                 hasSpatialExtents;
    double               minSpatialExtents[3];
    double               maxSpatialExtents[3];
    int                  numBlocks;
    int                  blockOrigin;
    std::string          blockPieceName;
    std::string          blockTitle;
    stringVector         blockNames;
    NameschemeAttributes blockNameScheme;
    int                  numGroups;
    int                  groupOrigin;
    std::string          groupPieceName;
    std::string          groupTitle;
    intVector            groupIds;
    intVector            groupIdsBasedOnRange;
    bool                 disjointElements;
    avtGhostType         containsGhostZones;
    bool                 containsOriginalCells;
    bool                 containsOriginalNodes;
    bool                 containsGlobalNodeIds;
    bool                 containsGlobalZoneIds;
    LoadBalanceScheme    loadBalanceScheme;
    bool                 nodesAreCritical;
    float                unitCellVectors[9];
    float                unitCellOrigin[3];
    bool                 rectilinearGridHasTransform;
    double               rectilinearGridTransform[16];
    int                  nodeOrigin;
    bool                 containsExteriorBoundaryGhosts;
    bool                 hideFromGUI;
    int                  LODs;

private:
    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define AVTMESHMETADATA_TMFS "ssbiiiiissssssbDDiisss*aiissi*i*bibbbbibFFbDibbi"

#endif
