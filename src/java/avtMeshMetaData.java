// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: avtMeshMetaData
//
// Purpose:
//    Contains mesh metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class avtMeshMetaData extends AttributeSubject
{
    public avtMeshMetaData()
    {
        super(41);

        name = new String("mesh");
        originalName = new String("");
        validVariable = true;
        meshType = 0;
        meshCoordType = 0;
        cellOrigin = 0;
        spatialDimension = 3;
        topologicalDimension = 3;
        xUnits = new String("");
        yUnits = new String("");
        zUnits = new String("");
        xLabel = new String("X-Axis");
        yLabel = new String("Y-Axis");
        zLabel = new String("Z-Axis");
        hasSpatialExtents = false;
        minSpatialExtents = new double[3];
        minSpatialExtents[0] = 0;
        minSpatialExtents[1] = 0;
        minSpatialExtents[2] = 0;
        maxSpatialExtents = new double[3];
        maxSpatialExtents[0] = 0;
        maxSpatialExtents[1] = 0;
        maxSpatialExtents[2] = 0;
        numBlocks = 1;
        blockOrigin = 0;
        blockPieceName = new String("domain");
        blockTitle = new String("domains");
        blockNames = new Vector();
        numGroups = 0;
        groupOrigin = 0;
        groupPieceName = new String("group");
        groupTitle = new String("groups");
        groupIds = new Vector();
        disjointElements = false;
        containsGhostZones = 0;
        containsOriginalCells = false;
        containsOriginalNodes = false;
        containsGlobalNodeIds = false;
        containsGlobalZoneIds = false;
        loadBalanceScheme = 0;
        nodesAreCritical = false;
        unitCellVectors = new float[9];
        unitCellVectors[0] = 1f;
        unitCellVectors[1] = 0f;
        unitCellVectors[2] = 0f;
        unitCellVectors[3] = 0f;
        unitCellVectors[4] = 1f;
        unitCellVectors[5] = 0f;
        unitCellVectors[6] = 0f;
        unitCellVectors[7] = 0f;
        unitCellVectors[8] = 1f;
        rectilinearGridHasTransform = false;
        rectilinearGridTransform = new double[16];
        rectilinearGridTransform[0] = 1;
        rectilinearGridTransform[1] = 0;
        rectilinearGridTransform[2] = 0;
        rectilinearGridTransform[3] = 0;
        rectilinearGridTransform[4] = 0;
        rectilinearGridTransform[5] = 1;
        rectilinearGridTransform[6] = 0;
        rectilinearGridTransform[7] = 0;
        rectilinearGridTransform[8] = 0;
        rectilinearGridTransform[9] = 0;
        rectilinearGridTransform[10] = 1;
        rectilinearGridTransform[11] = 0;
        rectilinearGridTransform[12] = 0;
        rectilinearGridTransform[13] = 0;
        rectilinearGridTransform[14] = 0;
        rectilinearGridTransform[15] = 1;
        nodeOrigin = 0;
        containsExteriorBoundaryGhosts = false;
        hideFromGUI = false;
    }

    public avtMeshMetaData(avtMeshMetaData obj)
    {
        super(41);

        int i;

        name = new String(obj.name);
        originalName = new String(obj.originalName);
        validVariable = obj.validVariable;
        meshType = obj.meshType;
        meshCoordType = obj.meshCoordType;
        cellOrigin = obj.cellOrigin;
        spatialDimension = obj.spatialDimension;
        topologicalDimension = obj.topologicalDimension;
        xUnits = new String(obj.xUnits);
        yUnits = new String(obj.yUnits);
        zUnits = new String(obj.zUnits);
        xLabel = new String(obj.xLabel);
        yLabel = new String(obj.yLabel);
        zLabel = new String(obj.zLabel);
        hasSpatialExtents = obj.hasSpatialExtents;
        minSpatialExtents = new double[3];
        minSpatialExtents[0] = obj.minSpatialExtents[0];
        minSpatialExtents[1] = obj.minSpatialExtents[1];
        minSpatialExtents[2] = obj.minSpatialExtents[2];

        maxSpatialExtents = new double[3];
        maxSpatialExtents[0] = obj.maxSpatialExtents[0];
        maxSpatialExtents[1] = obj.maxSpatialExtents[1];
        maxSpatialExtents[2] = obj.maxSpatialExtents[2];

        numBlocks = obj.numBlocks;
        blockOrigin = obj.blockOrigin;
        blockPieceName = new String(obj.blockPieceName);
        blockTitle = new String(obj.blockTitle);
        blockNames = new Vector(obj.blockNames.size());
        for(i = 0; i < obj.blockNames.size(); ++i)
            blockNames.addElement(new String((String)obj.blockNames.elementAt(i)));

        numGroups = obj.numGroups;
        groupOrigin = obj.groupOrigin;
        groupPieceName = new String(obj.groupPieceName);
        groupTitle = new String(obj.groupTitle);
        groupIds = new Vector();
        for(i = 0; i < obj.groupIds.size(); ++i)
        {
            Integer iv = (Integer)obj.groupIds.elementAt(i);
            groupIds.addElement(new Integer(iv.intValue()));
        }
        disjointElements = obj.disjointElements;
        containsGhostZones = obj.containsGhostZones;
        containsOriginalCells = obj.containsOriginalCells;
        containsOriginalNodes = obj.containsOriginalNodes;
        containsGlobalNodeIds = obj.containsGlobalNodeIds;
        containsGlobalZoneIds = obj.containsGlobalZoneIds;
        loadBalanceScheme = obj.loadBalanceScheme;
        nodesAreCritical = obj.nodesAreCritical;
        unitCellVectors = new float[9];
        for(i = 0; i < obj.unitCellVectors.length; ++i)
            unitCellVectors[i] = obj.unitCellVectors[i];

        rectilinearGridHasTransform = obj.rectilinearGridHasTransform;
        rectilinearGridTransform = new double[16];
        for(i = 0; i < obj.rectilinearGridTransform.length; ++i)
            rectilinearGridTransform[i] = obj.rectilinearGridTransform[i];

        nodeOrigin = obj.nodeOrigin;
        containsExteriorBoundaryGhosts = obj.containsExteriorBoundaryGhosts;
        hideFromGUI = obj.hideFromGUI;

        SelectAll();
    }

    public boolean equals(avtMeshMetaData obj)
    {
        int i;

        // Compare the minSpatialExtents arrays.
        boolean minSpatialExtents_equal = true;
        for(i = 0; i < 3 && minSpatialExtents_equal; ++i)
            minSpatialExtents_equal = (minSpatialExtents[i] == obj.minSpatialExtents[i]);

        // Compare the maxSpatialExtents arrays.
        boolean maxSpatialExtents_equal = true;
        for(i = 0; i < 3 && maxSpatialExtents_equal; ++i)
            maxSpatialExtents_equal = (maxSpatialExtents[i] == obj.maxSpatialExtents[i]);

        // Compare the elements in the blockNames vector.
        boolean blockNames_equal = (obj.blockNames.size() == blockNames.size());
        for(i = 0; (i < blockNames.size()) && blockNames_equal; ++i)
        {
            // Make references to String from Object.
            String blockNames1 = (String)blockNames.elementAt(i);
            String blockNames2 = (String)obj.blockNames.elementAt(i);
            blockNames_equal = blockNames1.equals(blockNames2);
        }
        // Compare the elements in the groupIds vector.
        boolean groupIds_equal = (obj.groupIds.size() == groupIds.size());
        for(i = 0; (i < groupIds.size()) && groupIds_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer groupIds1 = (Integer)groupIds.elementAt(i);
            Integer groupIds2 = (Integer)obj.groupIds.elementAt(i);
            groupIds_equal = groupIds1.equals(groupIds2);
        }
        // Compare the unitCellVectors arrays.
        boolean unitCellVectors_equal = true;
        for(i = 0; i < 9 && unitCellVectors_equal; ++i)
            unitCellVectors_equal = (unitCellVectors[i] == obj.unitCellVectors[i]);

        // Compare the rectilinearGridTransform arrays.
        boolean rectilinearGridTransform_equal = true;
        for(i = 0; i < 16 && rectilinearGridTransform_equal; ++i)
            rectilinearGridTransform_equal = (rectilinearGridTransform[i] == obj.rectilinearGridTransform[i]);

        // Create the return value
        return ((name.equals(obj.name)) &&
                (originalName.equals(obj.originalName)) &&
                (validVariable == obj.validVariable) &&
                (meshType == obj.meshType) &&
                (meshCoordType == obj.meshCoordType) &&
                (cellOrigin == obj.cellOrigin) &&
                (spatialDimension == obj.spatialDimension) &&
                (topologicalDimension == obj.topologicalDimension) &&
                (xUnits.equals(obj.xUnits)) &&
                (yUnits.equals(obj.yUnits)) &&
                (zUnits.equals(obj.zUnits)) &&
                (xLabel.equals(obj.xLabel)) &&
                (yLabel.equals(obj.yLabel)) &&
                (zLabel.equals(obj.zLabel)) &&
                (hasSpatialExtents == obj.hasSpatialExtents) &&
                minSpatialExtents_equal &&
                maxSpatialExtents_equal &&
                (numBlocks == obj.numBlocks) &&
                (blockOrigin == obj.blockOrigin) &&
                (blockPieceName.equals(obj.blockPieceName)) &&
                (blockTitle.equals(obj.blockTitle)) &&
                blockNames_equal &&
                (numGroups == obj.numGroups) &&
                (groupOrigin == obj.groupOrigin) &&
                (groupPieceName.equals(obj.groupPieceName)) &&
                (groupTitle.equals(obj.groupTitle)) &&
                groupIds_equal &&
                (disjointElements == obj.disjointElements) &&
                (containsGhostZones == obj.containsGhostZones) &&
                (containsOriginalCells == obj.containsOriginalCells) &&
                (containsOriginalNodes == obj.containsOriginalNodes) &&
                (containsGlobalNodeIds == obj.containsGlobalNodeIds) &&
                (containsGlobalZoneIds == obj.containsGlobalZoneIds) &&
                (loadBalanceScheme == obj.loadBalanceScheme) &&
                (nodesAreCritical == obj.nodesAreCritical) &&
                unitCellVectors_equal &&
                (rectilinearGridHasTransform == obj.rectilinearGridHasTransform) &&
                rectilinearGridTransform_equal &&
                (nodeOrigin == obj.nodeOrigin) &&
                (containsExteriorBoundaryGhosts == obj.containsExteriorBoundaryGhosts) &&
                (hideFromGUI == obj.hideFromGUI));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetOriginalName(String originalName_)
    {
        originalName = originalName_;
        Select(1);
    }

    public void SetValidVariable(boolean validVariable_)
    {
        validVariable = validVariable_;
        Select(2);
    }

    public void SetMeshType(int meshType_)
    {
        meshType = meshType_;
        Select(3);
    }

    public void SetMeshCoordType(int meshCoordType_)
    {
        meshCoordType = meshCoordType_;
        Select(4);
    }

    public void SetCellOrigin(int cellOrigin_)
    {
        cellOrigin = cellOrigin_;
        Select(5);
    }

    public void SetSpatialDimension(int spatialDimension_)
    {
        spatialDimension = spatialDimension_;
        Select(6);
    }

    public void SetTopologicalDimension(int topologicalDimension_)
    {
        topologicalDimension = topologicalDimension_;
        Select(7);
    }

    public void SetXUnits(String xUnits_)
    {
        xUnits = xUnits_;
        Select(8);
    }

    public void SetYUnits(String yUnits_)
    {
        yUnits = yUnits_;
        Select(9);
    }

    public void SetZUnits(String zUnits_)
    {
        zUnits = zUnits_;
        Select(10);
    }

    public void SetXLabel(String xLabel_)
    {
        xLabel = xLabel_;
        Select(11);
    }

    public void SetYLabel(String yLabel_)
    {
        yLabel = yLabel_;
        Select(12);
    }

    public void SetZLabel(String zLabel_)
    {
        zLabel = zLabel_;
        Select(13);
    }

    public void SetHasSpatialExtents(boolean hasSpatialExtents_)
    {
        hasSpatialExtents = hasSpatialExtents_;
        Select(14);
    }

    public void SetMinSpatialExtents(double[] minSpatialExtents_)
    {
        minSpatialExtents[0] = minSpatialExtents_[0];
        minSpatialExtents[1] = minSpatialExtents_[1];
        minSpatialExtents[2] = minSpatialExtents_[2];
        Select(15);
    }

    public void SetMinSpatialExtents(double e0, double e1, double e2)
    {
        minSpatialExtents[0] = e0;
        minSpatialExtents[1] = e1;
        minSpatialExtents[2] = e2;
        Select(15);
    }

    public void SetMaxSpatialExtents(double[] maxSpatialExtents_)
    {
        maxSpatialExtents[0] = maxSpatialExtents_[0];
        maxSpatialExtents[1] = maxSpatialExtents_[1];
        maxSpatialExtents[2] = maxSpatialExtents_[2];
        Select(16);
    }

    public void SetMaxSpatialExtents(double e0, double e1, double e2)
    {
        maxSpatialExtents[0] = e0;
        maxSpatialExtents[1] = e1;
        maxSpatialExtents[2] = e2;
        Select(16);
    }

    public void SetNumBlocks(int numBlocks_)
    {
        numBlocks = numBlocks_;
        Select(17);
    }

    public void SetBlockOrigin(int blockOrigin_)
    {
        blockOrigin = blockOrigin_;
        Select(18);
    }

    public void SetBlockPieceName(String blockPieceName_)
    {
        blockPieceName = blockPieceName_;
        Select(19);
    }

    public void SetBlockTitle(String blockTitle_)
    {
        blockTitle = blockTitle_;
        Select(20);
    }

    public void SetBlockNames(Vector blockNames_)
    {
        blockNames = blockNames_;
        Select(21);
    }

    public void SetNumGroups(int numGroups_)
    {
        numGroups = numGroups_;
        Select(22);
    }

    public void SetGroupOrigin(int groupOrigin_)
    {
        groupOrigin = groupOrigin_;
        Select(23);
    }

    public void SetGroupPieceName(String groupPieceName_)
    {
        groupPieceName = groupPieceName_;
        Select(24);
    }

    public void SetGroupTitle(String groupTitle_)
    {
        groupTitle = groupTitle_;
        Select(25);
    }

    public void SetGroupIds(Vector groupIds_)
    {
        groupIds = groupIds_;
        Select(26);
    }

    public void SetDisjointElements(boolean disjointElements_)
    {
        disjointElements = disjointElements_;
        Select(27);
    }

    public void SetContainsGhostZones(int containsGhostZones_)
    {
        containsGhostZones = containsGhostZones_;
        Select(28);
    }

    public void SetContainsOriginalCells(boolean containsOriginalCells_)
    {
        containsOriginalCells = containsOriginalCells_;
        Select(29);
    }

    public void SetContainsOriginalNodes(boolean containsOriginalNodes_)
    {
        containsOriginalNodes = containsOriginalNodes_;
        Select(30);
    }

    public void SetContainsGlobalNodeIds(boolean containsGlobalNodeIds_)
    {
        containsGlobalNodeIds = containsGlobalNodeIds_;
        Select(31);
    }

    public void SetContainsGlobalZoneIds(boolean containsGlobalZoneIds_)
    {
        containsGlobalZoneIds = containsGlobalZoneIds_;
        Select(32);
    }

    public void SetLoadBalanceScheme(int loadBalanceScheme_)
    {
        loadBalanceScheme = loadBalanceScheme_;
        Select(33);
    }

    public void SetNodesAreCritical(boolean nodesAreCritical_)
    {
        nodesAreCritical = nodesAreCritical_;
        Select(34);
    }

    public void SetUnitCellVectors(float[] unitCellVectors_)
    {
        for(int i = 0; i < 9; ++i)
             unitCellVectors[i] = unitCellVectors_[i];
        Select(35);
    }

    public void SetRectilinearGridHasTransform(boolean rectilinearGridHasTransform_)
    {
        rectilinearGridHasTransform = rectilinearGridHasTransform_;
        Select(36);
    }

    public void SetRectilinearGridTransform(double[] rectilinearGridTransform_)
    {
        for(int i = 0; i < 16; ++i)
             rectilinearGridTransform[i] = rectilinearGridTransform_[i];
        Select(37);
    }

    public void SetNodeOrigin(int nodeOrigin_)
    {
        nodeOrigin = nodeOrigin_;
        Select(38);
    }

    public void SetContainsExteriorBoundaryGhosts(boolean containsExteriorBoundaryGhosts_)
    {
        containsExteriorBoundaryGhosts = containsExteriorBoundaryGhosts_;
        Select(39);
    }

    public void SetHideFromGUI(boolean hideFromGUI_)
    {
        hideFromGUI = hideFromGUI_;
        Select(40);
    }

    // Property getting methods
    public String   GetName() { return name; }
    public String   GetOriginalName() { return originalName; }
    public boolean  GetValidVariable() { return validVariable; }
    public int      GetMeshType() { return meshType; }
    public int      GetMeshCoordType() { return meshCoordType; }
    public int      GetCellOrigin() { return cellOrigin; }
    public int      GetSpatialDimension() { return spatialDimension; }
    public int      GetTopologicalDimension() { return topologicalDimension; }
    public String   GetXUnits() { return xUnits; }
    public String   GetYUnits() { return yUnits; }
    public String   GetZUnits() { return zUnits; }
    public String   GetXLabel() { return xLabel; }
    public String   GetYLabel() { return yLabel; }
    public String   GetZLabel() { return zLabel; }
    public boolean  GetHasSpatialExtents() { return hasSpatialExtents; }
    public double[] GetMinSpatialExtents() { return minSpatialExtents; }
    public double[] GetMaxSpatialExtents() { return maxSpatialExtents; }
    public int      GetNumBlocks() { return numBlocks; }
    public int      GetBlockOrigin() { return blockOrigin; }
    public String   GetBlockPieceName() { return blockPieceName; }
    public String   GetBlockTitle() { return blockTitle; }
    public Vector   GetBlockNames() { return blockNames; }
    public int      GetNumGroups() { return numGroups; }
    public int      GetGroupOrigin() { return groupOrigin; }
    public String   GetGroupPieceName() { return groupPieceName; }
    public String   GetGroupTitle() { return groupTitle; }
    public Vector   GetGroupIds() { return groupIds; }
    public boolean  GetDisjointElements() { return disjointElements; }
    public int      GetContainsGhostZones() { return containsGhostZones; }
    public boolean  GetContainsOriginalCells() { return containsOriginalCells; }
    public boolean  GetContainsOriginalNodes() { return containsOriginalNodes; }
    public boolean  GetContainsGlobalNodeIds() { return containsGlobalNodeIds; }
    public boolean  GetContainsGlobalZoneIds() { return containsGlobalZoneIds; }
    public int      GetLoadBalanceScheme() { return loadBalanceScheme; }
    public boolean  GetNodesAreCritical() { return nodesAreCritical; }
    public float[]  GetUnitCellVectors() { return unitCellVectors; }
    public boolean  GetRectilinearGridHasTransform() { return rectilinearGridHasTransform; }
    public double[] GetRectilinearGridTransform() { return rectilinearGridTransform; }
    public int      GetNodeOrigin() { return nodeOrigin; }
    public boolean  GetContainsExteriorBoundaryGhosts() { return containsExteriorBoundaryGhosts; }
    public boolean  GetHideFromGUI() { return hideFromGUI; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteString(originalName);
        if(WriteSelect(2, buf))
            buf.WriteBool(validVariable);
        if(WriteSelect(3, buf))
            buf.WriteInt(meshType);
        if(WriteSelect(4, buf))
            buf.WriteInt(meshCoordType);
        if(WriteSelect(5, buf))
            buf.WriteInt(cellOrigin);
        if(WriteSelect(6, buf))
            buf.WriteInt(spatialDimension);
        if(WriteSelect(7, buf))
            buf.WriteInt(topologicalDimension);
        if(WriteSelect(8, buf))
            buf.WriteString(xUnits);
        if(WriteSelect(9, buf))
            buf.WriteString(yUnits);
        if(WriteSelect(10, buf))
            buf.WriteString(zUnits);
        if(WriteSelect(11, buf))
            buf.WriteString(xLabel);
        if(WriteSelect(12, buf))
            buf.WriteString(yLabel);
        if(WriteSelect(13, buf))
            buf.WriteString(zLabel);
        if(WriteSelect(14, buf))
            buf.WriteBool(hasSpatialExtents);
        if(WriteSelect(15, buf))
            buf.WriteDoubleArray(minSpatialExtents);
        if(WriteSelect(16, buf))
            buf.WriteDoubleArray(maxSpatialExtents);
        if(WriteSelect(17, buf))
            buf.WriteInt(numBlocks);
        if(WriteSelect(18, buf))
            buf.WriteInt(blockOrigin);
        if(WriteSelect(19, buf))
            buf.WriteString(blockPieceName);
        if(WriteSelect(20, buf))
            buf.WriteString(blockTitle);
        if(WriteSelect(21, buf))
            buf.WriteStringVector(blockNames);
        if(WriteSelect(22, buf))
            buf.WriteInt(numGroups);
        if(WriteSelect(23, buf))
            buf.WriteInt(groupOrigin);
        if(WriteSelect(24, buf))
            buf.WriteString(groupPieceName);
        if(WriteSelect(25, buf))
            buf.WriteString(groupTitle);
        if(WriteSelect(26, buf))
            buf.WriteIntVector(groupIds);
        if(WriteSelect(27, buf))
            buf.WriteBool(disjointElements);
        if(WriteSelect(28, buf))
            buf.WriteInt(containsGhostZones);
        if(WriteSelect(29, buf))
            buf.WriteBool(containsOriginalCells);
        if(WriteSelect(30, buf))
            buf.WriteBool(containsOriginalNodes);
        if(WriteSelect(31, buf))
            buf.WriteBool(containsGlobalNodeIds);
        if(WriteSelect(32, buf))
            buf.WriteBool(containsGlobalZoneIds);
        if(WriteSelect(33, buf))
            buf.WriteInt(loadBalanceScheme);
        if(WriteSelect(34, buf))
            buf.WriteBool(nodesAreCritical);
        if(WriteSelect(35, buf))
            buf.WriteFloatArray(unitCellVectors);
        if(WriteSelect(36, buf))
            buf.WriteBool(rectilinearGridHasTransform);
        if(WriteSelect(37, buf))
            buf.WriteDoubleArray(rectilinearGridTransform);
        if(WriteSelect(38, buf))
            buf.WriteInt(nodeOrigin);
        if(WriteSelect(39, buf))
            buf.WriteBool(containsExteriorBoundaryGhosts);
        if(WriteSelect(40, buf))
            buf.WriteBool(hideFromGUI);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetName(buf.ReadString());
                break;
            case 1:
                SetOriginalName(buf.ReadString());
                break;
            case 2:
                SetValidVariable(buf.ReadBool());
                break;
            case 3:
                SetMeshType(buf.ReadInt());
                break;
            case 4:
                SetMeshCoordType(buf.ReadInt());
                break;
            case 5:
                SetCellOrigin(buf.ReadInt());
                break;
            case 6:
                SetSpatialDimension(buf.ReadInt());
                break;
            case 7:
                SetTopologicalDimension(buf.ReadInt());
                break;
            case 8:
                SetXUnits(buf.ReadString());
                break;
            case 9:
                SetYUnits(buf.ReadString());
                break;
            case 10:
                SetZUnits(buf.ReadString());
                break;
            case 11:
                SetXLabel(buf.ReadString());
                break;
            case 12:
                SetYLabel(buf.ReadString());
                break;
            case 13:
                SetZLabel(buf.ReadString());
                break;
            case 14:
                SetHasSpatialExtents(buf.ReadBool());
                break;
            case 15:
                SetMinSpatialExtents(buf.ReadDoubleArray());
                break;
            case 16:
                SetMaxSpatialExtents(buf.ReadDoubleArray());
                break;
            case 17:
                SetNumBlocks(buf.ReadInt());
                break;
            case 18:
                SetBlockOrigin(buf.ReadInt());
                break;
            case 19:
                SetBlockPieceName(buf.ReadString());
                break;
            case 20:
                SetBlockTitle(buf.ReadString());
                break;
            case 21:
                SetBlockNames(buf.ReadStringVector());
                break;
            case 22:
                SetNumGroups(buf.ReadInt());
                break;
            case 23:
                SetGroupOrigin(buf.ReadInt());
                break;
            case 24:
                SetGroupPieceName(buf.ReadString());
                break;
            case 25:
                SetGroupTitle(buf.ReadString());
                break;
            case 26:
                SetGroupIds(buf.ReadIntVector());
                break;
            case 27:
                SetDisjointElements(buf.ReadBool());
                break;
            case 28:
                SetContainsGhostZones(buf.ReadInt());
                break;
            case 29:
                SetContainsOriginalCells(buf.ReadBool());
                break;
            case 30:
                SetContainsOriginalNodes(buf.ReadBool());
                break;
            case 31:
                SetContainsGlobalNodeIds(buf.ReadBool());
                break;
            case 32:
                SetContainsGlobalZoneIds(buf.ReadBool());
                break;
            case 33:
                SetLoadBalanceScheme(buf.ReadInt());
                break;
            case 34:
                SetNodesAreCritical(buf.ReadBool());
                break;
            case 35:
                SetUnitCellVectors(buf.ReadFloatArray());
                break;
            case 36:
                SetRectilinearGridHasTransform(buf.ReadBool());
                break;
            case 37:
                SetRectilinearGridTransform(buf.ReadDoubleArray());
                break;
            case 38:
                SetNodeOrigin(buf.ReadInt());
                break;
            case 39:
                SetContainsExteriorBoundaryGhosts(buf.ReadBool());
                break;
            case 40:
                SetHideFromGUI(buf.ReadBool());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + stringToString("originalName", originalName, indent) + "\n";
        str = str + boolToString("validVariable", validVariable, indent) + "\n";
        str = str + intToString("meshType", meshType, indent) + "\n";
        str = str + intToString("meshCoordType", meshCoordType, indent) + "\n";
        str = str + intToString("cellOrigin", cellOrigin, indent) + "\n";
        str = str + intToString("spatialDimension", spatialDimension, indent) + "\n";
        str = str + intToString("topologicalDimension", topologicalDimension, indent) + "\n";
        str = str + stringToString("xUnits", xUnits, indent) + "\n";
        str = str + stringToString("yUnits", yUnits, indent) + "\n";
        str = str + stringToString("zUnits", zUnits, indent) + "\n";
        str = str + stringToString("xLabel", xLabel, indent) + "\n";
        str = str + stringToString("yLabel", yLabel, indent) + "\n";
        str = str + stringToString("zLabel", zLabel, indent) + "\n";
        str = str + boolToString("hasSpatialExtents", hasSpatialExtents, indent) + "\n";
        str = str + doubleArrayToString("minSpatialExtents", minSpatialExtents, indent) + "\n";
        str = str + doubleArrayToString("maxSpatialExtents", maxSpatialExtents, indent) + "\n";
        str = str + intToString("numBlocks", numBlocks, indent) + "\n";
        str = str + intToString("blockOrigin", blockOrigin, indent) + "\n";
        str = str + stringToString("blockPieceName", blockPieceName, indent) + "\n";
        str = str + stringToString("blockTitle", blockTitle, indent) + "\n";
        str = str + stringVectorToString("blockNames", blockNames, indent) + "\n";
        str = str + intToString("numGroups", numGroups, indent) + "\n";
        str = str + intToString("groupOrigin", groupOrigin, indent) + "\n";
        str = str + stringToString("groupPieceName", groupPieceName, indent) + "\n";
        str = str + stringToString("groupTitle", groupTitle, indent) + "\n";
        str = str + intVectorToString("groupIds", groupIds, indent) + "\n";
        str = str + boolToString("disjointElements", disjointElements, indent) + "\n";
        str = str + intToString("containsGhostZones", containsGhostZones, indent) + "\n";
        str = str + boolToString("containsOriginalCells", containsOriginalCells, indent) + "\n";
        str = str + boolToString("containsOriginalNodes", containsOriginalNodes, indent) + "\n";
        str = str + boolToString("containsGlobalNodeIds", containsGlobalNodeIds, indent) + "\n";
        str = str + boolToString("containsGlobalZoneIds", containsGlobalZoneIds, indent) + "\n";
        str = str + intToString("loadBalanceScheme", loadBalanceScheme, indent) + "\n";
        str = str + boolToString("nodesAreCritical", nodesAreCritical, indent) + "\n";
        str = str + floatArrayToString("unitCellVectors", unitCellVectors, indent) + "\n";
        str = str + boolToString("rectilinearGridHasTransform", rectilinearGridHasTransform, indent) + "\n";
        str = str + doubleArrayToString("rectilinearGridTransform", rectilinearGridTransform, indent) + "\n";
        str = str + intToString("nodeOrigin", nodeOrigin, indent) + "\n";
        str = str + boolToString("containsExteriorBoundaryGhosts", containsExteriorBoundaryGhosts, indent) + "\n";
        str = str + boolToString("hideFromGUI", hideFromGUI, indent) + "\n";
        return str;
    }


    // Attributes
    private String   name;
    private String   originalName;
    private boolean  validVariable;
    private int      meshType;
    private int      meshCoordType;
    private int      cellOrigin;
    private int      spatialDimension;
    private int      topologicalDimension;
    private String   xUnits;
    private String   yUnits;
    private String   zUnits;
    private String   xLabel;
    private String   yLabel;
    private String   zLabel;
    private boolean  hasSpatialExtents;
    private double[] minSpatialExtents;
    private double[] maxSpatialExtents;
    private int      numBlocks;
    private int      blockOrigin;
    private String   blockPieceName;
    private String   blockTitle;
    private Vector   blockNames; // vector of String objects
    private int      numGroups;
    private int      groupOrigin;
    private String   groupPieceName;
    private String   groupTitle;
    private Vector   groupIds; // vector of Integer objects
    private boolean  disjointElements;
    private int      containsGhostZones;
    private boolean  containsOriginalCells;
    private boolean  containsOriginalNodes;
    private boolean  containsGlobalNodeIds;
    private boolean  containsGlobalZoneIds;
    private int      loadBalanceScheme;
    private boolean  nodesAreCritical;
    private float[]  unitCellVectors;
    private boolean  rectilinearGridHasTransform;
    private double[] rectilinearGridTransform;
    private int      nodeOrigin;
    private boolean  containsExteriorBoundaryGhosts;
    private boolean  hideFromGUI;
}

