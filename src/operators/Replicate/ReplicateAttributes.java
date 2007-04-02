// ***************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
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

// ****************************************************************************
// Class: ReplicateAttributes
//
// Purpose:
//    This class contains attributes for the replicate operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Mar 22 14:02:25 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class ReplicateAttributes extends AttributeSubject implements Plugin
{
    public ReplicateAttributes()
    {
        super(9);

        useUnitCellVectors = false;
        xVector = new double[3];
        xVector[0] = 1;
        xVector[1] = 0;
        xVector[2] = 0;
        yVector = new double[3];
        yVector[0] = 0;
        yVector[1] = 1;
        yVector[2] = 0;
        zVector = new double[3];
        zVector[0] = 0;
        zVector[1] = 0;
        zVector[2] = 1;
        xReplications = 1;
        yReplications = 1;
        zReplications = 1;
        mergeResults = true;
        replicateUnitCellAtoms = false;
    }

    public ReplicateAttributes(ReplicateAttributes obj)
    {
        super(9);

        int i;

        useUnitCellVectors = obj.useUnitCellVectors;
        xVector = new double[3];
        xVector[0] = obj.xVector[0];
        xVector[1] = obj.xVector[1];
        xVector[2] = obj.xVector[2];

        yVector = new double[3];
        yVector[0] = obj.yVector[0];
        yVector[1] = obj.yVector[1];
        yVector[2] = obj.yVector[2];

        zVector = new double[3];
        zVector[0] = obj.zVector[0];
        zVector[1] = obj.zVector[1];
        zVector[2] = obj.zVector[2];

        xReplications = obj.xReplications;
        yReplications = obj.yReplications;
        zReplications = obj.zReplications;
        mergeResults = obj.mergeResults;
        replicateUnitCellAtoms = obj.replicateUnitCellAtoms;

        SelectAll();
    }

    public boolean equals(ReplicateAttributes obj)
    {
        int i;

        // Compare the xVector arrays.
        boolean xVector_equal = true;
        for(i = 0; i < 3 && xVector_equal; ++i)
            xVector_equal = (xVector[i] == obj.xVector[i]);

        // Compare the yVector arrays.
        boolean yVector_equal = true;
        for(i = 0; i < 3 && yVector_equal; ++i)
            yVector_equal = (yVector[i] == obj.yVector[i]);

        // Compare the zVector arrays.
        boolean zVector_equal = true;
        for(i = 0; i < 3 && zVector_equal; ++i)
            zVector_equal = (zVector[i] == obj.zVector[i]);

        // Create the return value
        return ((useUnitCellVectors == obj.useUnitCellVectors) &&
                xVector_equal &&
                yVector_equal &&
                zVector_equal &&
                (xReplications == obj.xReplications) &&
                (yReplications == obj.yReplications) &&
                (zReplications == obj.zReplications) &&
                (mergeResults == obj.mergeResults) &&
                (replicateUnitCellAtoms == obj.replicateUnitCellAtoms));
    }

    public String GetName() { return "Replicate"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetUseUnitCellVectors(boolean useUnitCellVectors_)
    {
        useUnitCellVectors = useUnitCellVectors_;
        Select(0);
    }

    public void SetXVector(double[] xVector_)
    {
        xVector[0] = xVector_[0];
        xVector[1] = xVector_[1];
        xVector[2] = xVector_[2];
        Select(1);
    }

    public void SetXVector(double e0, double e1, double e2)
    {
        xVector[0] = e0;
        xVector[1] = e1;
        xVector[2] = e2;
        Select(1);
    }

    public void SetYVector(double[] yVector_)
    {
        yVector[0] = yVector_[0];
        yVector[1] = yVector_[1];
        yVector[2] = yVector_[2];
        Select(2);
    }

    public void SetYVector(double e0, double e1, double e2)
    {
        yVector[0] = e0;
        yVector[1] = e1;
        yVector[2] = e2;
        Select(2);
    }

    public void SetZVector(double[] zVector_)
    {
        zVector[0] = zVector_[0];
        zVector[1] = zVector_[1];
        zVector[2] = zVector_[2];
        Select(3);
    }

    public void SetZVector(double e0, double e1, double e2)
    {
        zVector[0] = e0;
        zVector[1] = e1;
        zVector[2] = e2;
        Select(3);
    }

    public void SetXReplications(int xReplications_)
    {
        xReplications = xReplications_;
        Select(4);
    }

    public void SetYReplications(int yReplications_)
    {
        yReplications = yReplications_;
        Select(5);
    }

    public void SetZReplications(int zReplications_)
    {
        zReplications = zReplications_;
        Select(6);
    }

    public void SetMergeResults(boolean mergeResults_)
    {
        mergeResults = mergeResults_;
        Select(7);
    }

    public void SetReplicateUnitCellAtoms(boolean replicateUnitCellAtoms_)
    {
        replicateUnitCellAtoms = replicateUnitCellAtoms_;
        Select(8);
    }

    // Property getting methods
    public boolean  GetUseUnitCellVectors() { return useUnitCellVectors; }
    public double[] GetXVector() { return xVector; }
    public double[] GetYVector() { return yVector; }
    public double[] GetZVector() { return zVector; }
    public int      GetXReplications() { return xReplications; }
    public int      GetYReplications() { return yReplications; }
    public int      GetZReplications() { return zReplications; }
    public boolean  GetMergeResults() { return mergeResults; }
    public boolean  GetReplicateUnitCellAtoms() { return replicateUnitCellAtoms; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(useUnitCellVectors);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(xVector);
        if(WriteSelect(2, buf))
            buf.WriteDoubleArray(yVector);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(zVector);
        if(WriteSelect(4, buf))
            buf.WriteInt(xReplications);
        if(WriteSelect(5, buf))
            buf.WriteInt(yReplications);
        if(WriteSelect(6, buf))
            buf.WriteInt(zReplications);
        if(WriteSelect(7, buf))
            buf.WriteBool(mergeResults);
        if(WriteSelect(8, buf))
            buf.WriteBool(replicateUnitCellAtoms);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetUseUnitCellVectors(buf.ReadBool());
                break;
            case 1:
                SetXVector(buf.ReadDoubleArray());
                break;
            case 2:
                SetYVector(buf.ReadDoubleArray());
                break;
            case 3:
                SetZVector(buf.ReadDoubleArray());
                break;
            case 4:
                SetXReplications(buf.ReadInt());
                break;
            case 5:
                SetYReplications(buf.ReadInt());
                break;
            case 6:
                SetZReplications(buf.ReadInt());
                break;
            case 7:
                SetMergeResults(buf.ReadBool());
                break;
            case 8:
                SetReplicateUnitCellAtoms(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private boolean  useUnitCellVectors;
    private double[] xVector;
    private double[] yVector;
    private double[] zVector;
    private int      xReplications;
    private int      yReplications;
    private int      zReplications;
    private boolean  mergeResults;
    private boolean  replicateUnitCellAtoms;
}

