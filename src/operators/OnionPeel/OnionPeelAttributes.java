// ***************************************************************************
//
// Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-442911
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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: OnionPeelAttributes
//
// Purpose:
//    Attributes for the onion peel operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class OnionPeelAttributes extends AttributeSubject implements Plugin
{
    private static int OnionPeelAttributes_numAdditionalAtts = 9;

    // Enum values
    public final static int NODEFACE_NODE = 0;
    public final static int NODEFACE_FACE = 1;

    public final static int SEEDIDTYPE_SEEDCELL = 0;
    public final static int SEEDIDTYPE_SEEDNODE = 1;


    public OnionPeelAttributes()
    {
        super(OnionPeelAttributes_numAdditionalAtts);

        adjacencyType = NODEFACE_NODE;
        useGlobalId = false;
        categoryName = new String("Whole");
        subsetName = new String("Whole");
        index = new Vector();
        index.addElement(new Integer(1));
        logical = false;
        requestedLayer = 0;
        seedType = SEEDIDTYPE_SEEDCELL;
        honorOriginalMesh = true;
    }

    public OnionPeelAttributes(int nMoreFields)
    {
        super(OnionPeelAttributes_numAdditionalAtts + nMoreFields);

        adjacencyType = NODEFACE_NODE;
        useGlobalId = false;
        categoryName = new String("Whole");
        subsetName = new String("Whole");
        index = new Vector();
        index.addElement(new Integer(1));
        logical = false;
        requestedLayer = 0;
        seedType = SEEDIDTYPE_SEEDCELL;
        honorOriginalMesh = true;
    }

    public OnionPeelAttributes(OnionPeelAttributes obj)
    {
        super(obj);

        int i;

        adjacencyType = obj.adjacencyType;
        useGlobalId = obj.useGlobalId;
        categoryName = new String(obj.categoryName);
        subsetName = new String(obj.subsetName);
        index = new Vector();
        for(i = 0; i < obj.index.size(); ++i)
        {
            Integer iv = (Integer)obj.index.elementAt(i);
            index.addElement(new Integer(iv.intValue()));
        }
        logical = obj.logical;
        requestedLayer = obj.requestedLayer;
        seedType = obj.seedType;
        honorOriginalMesh = obj.honorOriginalMesh;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return OnionPeelAttributes_numAdditionalAtts;
    }

    public boolean equals(OnionPeelAttributes obj)
    {
        int i;

        // Compare the elements in the index vector.
        boolean index_equal = (obj.index.size() == index.size());
        for(i = 0; (i < index.size()) && index_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer index1 = (Integer)index.elementAt(i);
            Integer index2 = (Integer)obj.index.elementAt(i);
            index_equal = index1.equals(index2);
        }
        // Create the return value
        return ((adjacencyType == obj.adjacencyType) &&
                (useGlobalId == obj.useGlobalId) &&
                (categoryName.equals(obj.categoryName)) &&
                (subsetName.equals(obj.subsetName)) &&
                index_equal &&
                (logical == obj.logical) &&
                (requestedLayer == obj.requestedLayer) &&
                (seedType == obj.seedType) &&
                (honorOriginalMesh == obj.honorOriginalMesh));
    }

    public String GetName() { return "OnionPeel"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetAdjacencyType(int adjacencyType_)
    {
        adjacencyType = adjacencyType_;
        Select(0);
    }

    public void SetUseGlobalId(boolean useGlobalId_)
    {
        useGlobalId = useGlobalId_;
        Select(1);
    }

    public void SetCategoryName(String categoryName_)
    {
        categoryName = categoryName_;
        Select(2);
    }

    public void SetSubsetName(String subsetName_)
    {
        subsetName = subsetName_;
        Select(3);
    }

    public void SetIndex(Vector index_)
    {
        index = index_;
        Select(4);
    }

    public void SetLogical(boolean logical_)
    {
        logical = logical_;
        Select(5);
    }

    public void SetRequestedLayer(int requestedLayer_)
    {
        requestedLayer = requestedLayer_;
        Select(6);
    }

    public void SetSeedType(int seedType_)
    {
        seedType = seedType_;
        Select(7);
    }

    public void SetHonorOriginalMesh(boolean honorOriginalMesh_)
    {
        honorOriginalMesh = honorOriginalMesh_;
        Select(8);
    }

    // Property getting methods
    public int     GetAdjacencyType() { return adjacencyType; }
    public boolean GetUseGlobalId() { return useGlobalId; }
    public String  GetCategoryName() { return categoryName; }
    public String  GetSubsetName() { return subsetName; }
    public Vector  GetIndex() { return index; }
    public boolean GetLogical() { return logical; }
    public int     GetRequestedLayer() { return requestedLayer; }
    public int     GetSeedType() { return seedType; }
    public boolean GetHonorOriginalMesh() { return honorOriginalMesh; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(adjacencyType);
        if(WriteSelect(1, buf))
            buf.WriteBool(useGlobalId);
        if(WriteSelect(2, buf))
            buf.WriteString(categoryName);
        if(WriteSelect(3, buf))
            buf.WriteString(subsetName);
        if(WriteSelect(4, buf))
            buf.WriteIntVector(index);
        if(WriteSelect(5, buf))
            buf.WriteBool(logical);
        if(WriteSelect(6, buf))
            buf.WriteInt(requestedLayer);
        if(WriteSelect(7, buf))
            buf.WriteInt(seedType);
        if(WriteSelect(8, buf))
            buf.WriteBool(honorOriginalMesh);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetAdjacencyType(buf.ReadInt());
            break;
        case 1:
            SetUseGlobalId(buf.ReadBool());
            break;
        case 2:
            SetCategoryName(buf.ReadString());
            break;
        case 3:
            SetSubsetName(buf.ReadString());
            break;
        case 4:
            SetIndex(buf.ReadIntVector());
            break;
        case 5:
            SetLogical(buf.ReadBool());
            break;
        case 6:
            SetRequestedLayer(buf.ReadInt());
            break;
        case 7:
            SetSeedType(buf.ReadInt());
            break;
        case 8:
            SetHonorOriginalMesh(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "adjacencyType = ";
        if(adjacencyType == NODEFACE_NODE)
            str = str + "NODEFACE_NODE";
        if(adjacencyType == NODEFACE_FACE)
            str = str + "NODEFACE_FACE";
        str = str + "\n";
        str = str + boolToString("useGlobalId", useGlobalId, indent) + "\n";
        str = str + stringToString("categoryName", categoryName, indent) + "\n";
        str = str + stringToString("subsetName", subsetName, indent) + "\n";
        str = str + intVectorToString("index", index, indent) + "\n";
        str = str + boolToString("logical", logical, indent) + "\n";
        str = str + intToString("requestedLayer", requestedLayer, indent) + "\n";
        str = str + indent + "seedType = ";
        if(seedType == SEEDIDTYPE_SEEDCELL)
            str = str + "SEEDIDTYPE_SEEDCELL";
        if(seedType == SEEDIDTYPE_SEEDNODE)
            str = str + "SEEDIDTYPE_SEEDNODE";
        str = str + "\n";
        str = str + boolToString("honorOriginalMesh", honorOriginalMesh, indent) + "\n";
        return str;
    }


    // Attributes
    private int     adjacencyType;
    private boolean useGlobalId;
    private String  categoryName;
    private String  subsetName;
    private Vector  index; // vector of Integer objects
    private boolean logical;
    private int     requestedLayer;
    private int     seedType;
    private boolean honorOriginalMesh;
}

