// ***************************************************************************
//
// Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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

package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: avtSubsetsMetaData
//
// Purpose:
//    Information about a particular category of subsets of a mesh (even for material subsets)
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class avtSubsetsMetaData extends avtVarMetaData
{
    private static int avtSubsetsMetaData_numAdditionalAtts = 12;

    // Enum values
    public final static int PARTIALCELLMODES_INCLUDE = 0;
    public final static int PARTIALCELLMODES_EXCLUDE = 1;
    public final static int PARTIALCELLMODES_DISSECT = 2;

    public final static int DECOMPMODE_NONE = 0;
    public final static int DECOMPMODE_COVER = 1;
    public final static int DECOMPMODE_PARTITION = 2;


    public avtSubsetsMetaData()
    {
        super(avtSubsetsMetaData_numAdditionalAtts);

        catName = new String("");
        catCount = 0;
        nameScheme = new NameschemeAttributes();
        colorScheme = new Vector();
        setsToChunksMaps = new Vector();
        graphEdges = new Vector();
        isChunkCat = false;
        isMaterialCat = false;
        isUnionOfChunks = false;
        hasPartialCells = false;
        decompMode = DECOMPMODE_NONE;
        maxTopoDim = 0;
    }

    public avtSubsetsMetaData(int nMoreFields)
    {
        super(avtSubsetsMetaData_numAdditionalAtts + nMoreFields);

        catName = new String("");
        catCount = 0;
        nameScheme = new NameschemeAttributes();
        colorScheme = new Vector();
        setsToChunksMaps = new Vector();
        graphEdges = new Vector();
        isChunkCat = false;
        isMaterialCat = false;
        isUnionOfChunks = false;
        hasPartialCells = false;
        decompMode = DECOMPMODE_NONE;
        maxTopoDim = 0;
    }

    public avtSubsetsMetaData(avtSubsetsMetaData obj)
    {
        super(avtSubsetsMetaData_numAdditionalAtts);

        int i;

        catName = new String(obj.catName);
        catCount = obj.catCount;
        nameScheme = new NameschemeAttributes(obj.nameScheme);
        colorScheme = new Vector(obj.colorScheme.size());
        for(i = 0; i < obj.colorScheme.size(); ++i)
            colorScheme.addElement(new String((String)obj.colorScheme.elementAt(i)));

        setsToChunksMaps = new Vector();
        for(i = 0; i < obj.setsToChunksMaps.size(); ++i)
        {
            Integer iv = (Integer)obj.setsToChunksMaps.elementAt(i);
            setsToChunksMaps.addElement(new Integer(iv.intValue()));
        }
        graphEdges = new Vector();
        for(i = 0; i < obj.graphEdges.size(); ++i)
        {
            Integer iv = (Integer)obj.graphEdges.elementAt(i);
            graphEdges.addElement(new Integer(iv.intValue()));
        }
        isChunkCat = obj.isChunkCat;
        isMaterialCat = obj.isMaterialCat;
        isUnionOfChunks = obj.isUnionOfChunks;
        hasPartialCells = obj.hasPartialCells;
        decompMode = obj.decompMode;
        maxTopoDim = obj.maxTopoDim;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return avtSubsetsMetaData_numAdditionalAtts;
    }

    public boolean equals(avtSubsetsMetaData obj)
    {
        int i;

        // Compare the elements in the colorScheme vector.
        boolean colorScheme_equal = (obj.colorScheme.size() == colorScheme.size());
        for(i = 0; (i < colorScheme.size()) && colorScheme_equal; ++i)
        {
            // Make references to String from Object.
            String colorScheme1 = (String)colorScheme.elementAt(i);
            String colorScheme2 = (String)obj.colorScheme.elementAt(i);
            colorScheme_equal = colorScheme1.equals(colorScheme2);
        }
        // Compare the elements in the setsToChunksMaps vector.
        boolean setsToChunksMaps_equal = (obj.setsToChunksMaps.size() == setsToChunksMaps.size());
        for(i = 0; (i < setsToChunksMaps.size()) && setsToChunksMaps_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer setsToChunksMaps1 = (Integer)setsToChunksMaps.elementAt(i);
            Integer setsToChunksMaps2 = (Integer)obj.setsToChunksMaps.elementAt(i);
            setsToChunksMaps_equal = setsToChunksMaps1.equals(setsToChunksMaps2);
        }
        // Compare the elements in the graphEdges vector.
        boolean graphEdges_equal = (obj.graphEdges.size() == graphEdges.size());
        for(i = 0; (i < graphEdges.size()) && graphEdges_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer graphEdges1 = (Integer)graphEdges.elementAt(i);
            Integer graphEdges2 = (Integer)obj.graphEdges.elementAt(i);
            graphEdges_equal = graphEdges1.equals(graphEdges2);
        }
        // Create the return value
        return (super.equals(obj) && (catName.equals(obj.catName)) &&
                (catCount == obj.catCount) &&
                (nameScheme.equals(obj.nameScheme)) &&
                colorScheme_equal &&
                setsToChunksMaps_equal &&
                graphEdges_equal &&
                (isChunkCat == obj.isChunkCat) &&
                (isMaterialCat == obj.isMaterialCat) &&
                (isUnionOfChunks == obj.isUnionOfChunks) &&
                (hasPartialCells == obj.hasPartialCells) &&
                (decompMode == obj.decompMode) &&
                (maxTopoDim == obj.maxTopoDim));
    }

    // Property setting methods
    public void SetCatName(String catName_)
    {
        catName = catName_;
        Select((new avtSubsetsMetaData()).Offset() + 0);
    }

    public void SetCatCount(int catCount_)
    {
        catCount = catCount_;
        Select((new avtSubsetsMetaData()).Offset() + 1);
    }

    public void SetNameScheme(NameschemeAttributes nameScheme_)
    {
        nameScheme = nameScheme_;
        Select((new avtSubsetsMetaData()).Offset() + 2);
    }

    public void SetColorScheme(Vector colorScheme_)
    {
        colorScheme = colorScheme_;
        Select((new avtSubsetsMetaData()).Offset() + 3);
    }

    public void SetSetsToChunksMaps(Vector setsToChunksMaps_)
    {
        setsToChunksMaps = setsToChunksMaps_;
        Select((new avtSubsetsMetaData()).Offset() + 4);
    }

    public void SetGraphEdges(Vector graphEdges_)
    {
        graphEdges = graphEdges_;
        Select((new avtSubsetsMetaData()).Offset() + 5);
    }

    public void SetIsChunkCat(boolean isChunkCat_)
    {
        isChunkCat = isChunkCat_;
        Select((new avtSubsetsMetaData()).Offset() + 6);
    }

    public void SetIsMaterialCat(boolean isMaterialCat_)
    {
        isMaterialCat = isMaterialCat_;
        Select((new avtSubsetsMetaData()).Offset() + 7);
    }

    public void SetIsUnionOfChunks(boolean isUnionOfChunks_)
    {
        isUnionOfChunks = isUnionOfChunks_;
        Select((new avtSubsetsMetaData()).Offset() + 8);
    }

    public void SetHasPartialCells(boolean hasPartialCells_)
    {
        hasPartialCells = hasPartialCells_;
        Select((new avtSubsetsMetaData()).Offset() + 9);
    }

    public void SetDecompMode(int decompMode_)
    {
        decompMode = decompMode_;
        Select((new avtSubsetsMetaData()).Offset() + 10);
    }

    public void SetMaxTopoDim(int maxTopoDim_)
    {
        maxTopoDim = maxTopoDim_;
        Select((new avtSubsetsMetaData()).Offset() + 11);
    }

    // Property getting methods
    public String               GetCatName() { return catName; }
    public int                  GetCatCount() { return catCount; }
    public NameschemeAttributes GetNameScheme() { return nameScheme; }
    public Vector               GetColorScheme() { return colorScheme; }
    public Vector               GetSetsToChunksMaps() { return setsToChunksMaps; }
    public Vector               GetGraphEdges() { return graphEdges; }
    public boolean              GetIsChunkCat() { return isChunkCat; }
    public boolean              GetIsMaterialCat() { return isMaterialCat; }
    public boolean              GetIsUnionOfChunks() { return isUnionOfChunks; }
    public boolean              GetHasPartialCells() { return hasPartialCells; }
    public int                  GetDecompMode() { return decompMode; }
    public int                  GetMaxTopoDim() { return maxTopoDim; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        super.WriteAtts(buf);

        int offset = (new avtSubsetsMetaData()).Offset();
        if(WriteSelect(offset + 0, buf))
            buf.WriteString(catName);
        if(WriteSelect(offset + 1, buf))
            buf.WriteInt(catCount);
        if(WriteSelect(offset + 2, buf))
            nameScheme.Write(buf);
        if(WriteSelect(offset + 3, buf))
            buf.WriteStringVector(colorScheme);
        if(WriteSelect(offset + 4, buf))
            buf.WriteIntVector(setsToChunksMaps);
        if(WriteSelect(offset + 5, buf))
            buf.WriteIntVector(graphEdges);
        if(WriteSelect(offset + 6, buf))
            buf.WriteBool(isChunkCat);
        if(WriteSelect(offset + 7, buf))
            buf.WriteBool(isMaterialCat);
        if(WriteSelect(offset + 8, buf))
            buf.WriteBool(isUnionOfChunks);
        if(WriteSelect(offset + 9, buf))
            buf.WriteBool(hasPartialCells);
        if(WriteSelect(offset + 10, buf))
            buf.WriteInt(decompMode);
        if(WriteSelect(offset + 11, buf))
            buf.WriteInt(maxTopoDim);
    }

    public void ReadAtts(int id, CommunicationBuffer buf)
    {
        int offset = (new avtSubsetsMetaData()).Offset();
        int index = id - offset;
        switch(index)
        {
        case 0:
            SetCatName(buf.ReadString());
            break;
        case 1:
            SetCatCount(buf.ReadInt());
            break;
        case 2:
            nameScheme.Read(buf);
            Select(offset + 2);
            break;
        case 3:
            SetColorScheme(buf.ReadStringVector());
            break;
        case 4:
            SetSetsToChunksMaps(buf.ReadIntVector());
            break;
        case 5:
            SetGraphEdges(buf.ReadIntVector());
            break;
        case 6:
            SetIsChunkCat(buf.ReadBool());
            break;
        case 7:
            SetIsMaterialCat(buf.ReadBool());
            break;
        case 8:
            SetIsUnionOfChunks(buf.ReadBool());
            break;
        case 9:
            SetHasPartialCells(buf.ReadBool());
            break;
        case 10:
            SetDecompMode(buf.ReadInt());
            break;
        case 11:
            SetMaxTopoDim(buf.ReadInt());
            break;
        default:
            super.ReadAtts(id, buf);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("catName", catName, indent) + "\n";
        str = str + intToString("catCount", catCount, indent) + "\n";
        str = str + indent + "nameScheme = {\n" + nameScheme.toString(indent + "    ") + indent + "}\n";
        str = str + stringVectorToString("colorScheme", colorScheme, indent) + "\n";
        str = str + intVectorToString("setsToChunksMaps", setsToChunksMaps, indent) + "\n";
        str = str + intVectorToString("graphEdges", graphEdges, indent) + "\n";
        str = str + boolToString("isChunkCat", isChunkCat, indent) + "\n";
        str = str + boolToString("isMaterialCat", isMaterialCat, indent) + "\n";
        str = str + boolToString("isUnionOfChunks", isUnionOfChunks, indent) + "\n";
        str = str + boolToString("hasPartialCells", hasPartialCells, indent) + "\n";
        str = str + indent + "decompMode = ";
        if(decompMode == DECOMPMODE_NONE)
            str = str + "DECOMPMODE_NONE";
        if(decompMode == DECOMPMODE_COVER)
            str = str + "DECOMPMODE_COVER";
        if(decompMode == DECOMPMODE_PARTITION)
            str = str + "DECOMPMODE_PARTITION";
        str = str + "\n";
        str = str + intToString("maxTopoDim", maxTopoDim, indent) + "\n";
        return super.toString(indent) + str;
    }


    // Attributes
    private String               catName;
    private int                  catCount;
    private NameschemeAttributes nameScheme;
    private Vector               colorScheme; // vector of String objects
    private Vector               setsToChunksMaps; // vector of Integer objects
    private Vector               graphEdges; // vector of Integer objects
    private boolean              isChunkCat;
    private boolean              isMaterialCat;
    private boolean              isUnionOfChunks;
    private boolean              hasPartialCells;
    private int                  decompMode;
    private int                  maxTopoDim;
}

