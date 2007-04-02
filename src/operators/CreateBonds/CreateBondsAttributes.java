// ***************************************************************************
//
// Copyright (c) 2000 - 2006, The Regents of the University of California
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
import java.lang.Integer;
import java.util.Vector;
import java.lang.Double;

// ****************************************************************************
// Class: CreateBondsAttributes
//
// Purpose:
//    Attributes for the CreateBondsOperator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Aug 29 11:09:32 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

public class CreateBondsAttributes extends AttributeSubject implements Plugin
{
    public CreateBondsAttributes()
    {
        super(11);

        elementVariable = new String("element");
        atomicNumber1 = new Vector();
        atomicNumber2 = new Vector();
        minDist = new Vector();
        maxDist = new Vector();
        simpleAlgorithmFlag = false;
        simpleMinDistWithH = 0.4;
        simpleMaxDistWithH = 1.2;
        simpleMinDistNonH = 0.4;
        simpleMaxDistNonH = 1.9;
        maxBondsClamp = 10;
    }

    public CreateBondsAttributes(CreateBondsAttributes obj)
    {
        super(11);

        int i;

        elementVariable = new String(obj.elementVariable);
        atomicNumber1 = new Vector();
        for(i = 0; i < obj.atomicNumber1.size(); ++i)
        {
            Integer iv = (Integer)obj.atomicNumber1.elementAt(i);
            atomicNumber1.addElement(new Integer(iv.intValue()));
        }
        atomicNumber2 = new Vector();
        for(i = 0; i < obj.atomicNumber2.size(); ++i)
        {
            Integer iv = (Integer)obj.atomicNumber2.elementAt(i);
            atomicNumber2.addElement(new Integer(iv.intValue()));
        }
        minDist = new Vector(obj.minDist.size());
        for(i = 0; i < obj.minDist.size(); ++i)
        {
            Double dv = (Double)obj.minDist.elementAt(i);
            minDist.addElement(new Double(dv.doubleValue()));
        }

        maxDist = new Vector(obj.maxDist.size());
        for(i = 0; i < obj.maxDist.size(); ++i)
        {
            Double dv = (Double)obj.maxDist.elementAt(i);
            maxDist.addElement(new Double(dv.doubleValue()));
        }

        simpleAlgorithmFlag = obj.simpleAlgorithmFlag;
        simpleMinDistWithH = obj.simpleMinDistWithH;
        simpleMaxDistWithH = obj.simpleMaxDistWithH;
        simpleMinDistNonH = obj.simpleMinDistNonH;
        simpleMaxDistNonH = obj.simpleMaxDistNonH;
        maxBondsClamp = obj.maxBondsClamp;

        SelectAll();
    }

    public boolean equals(CreateBondsAttributes obj)
    {
        int i;

        // Create the return value
        return ((elementVariable == obj.elementVariable) &&
                (atomicNumber1 == obj.atomicNumber1) &&
                (atomicNumber2 == obj.atomicNumber2) &&
                (minDist == obj.minDist) &&
                (maxDist == obj.maxDist) &&
                (simpleAlgorithmFlag == obj.simpleAlgorithmFlag) &&
                (simpleMinDistWithH == obj.simpleMinDistWithH) &&
                (simpleMaxDistWithH == obj.simpleMaxDistWithH) &&
                (simpleMinDistNonH == obj.simpleMinDistNonH) &&
                (simpleMaxDistNonH == obj.simpleMaxDistNonH) &&
                (maxBondsClamp == obj.maxBondsClamp));
    }

    public String GetName() { return "CreateBonds"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetElementVariable(String elementVariable_)
    {
        elementVariable = elementVariable_;
        Select(0);
    }

    public void SetAtomicNumber1(Vector atomicNumber1_)
    {
        atomicNumber1 = atomicNumber1_;
        Select(1);
    }

    public void SetAtomicNumber2(Vector atomicNumber2_)
    {
        atomicNumber2 = atomicNumber2_;
        Select(2);
    }

    public void SetMinDist(Vector minDist_)
    {
        minDist = minDist_;
        Select(3);
    }

    public void SetMaxDist(Vector maxDist_)
    {
        maxDist = maxDist_;
        Select(4);
    }

    public void SetSimpleAlgorithmFlag(boolean simpleAlgorithmFlag_)
    {
        simpleAlgorithmFlag = simpleAlgorithmFlag_;
        Select(5);
    }

    public void SetSimpleMinDistWithH(double simpleMinDistWithH_)
    {
        simpleMinDistWithH = simpleMinDistWithH_;
        Select(6);
    }

    public void SetSimpleMaxDistWithH(double simpleMaxDistWithH_)
    {
        simpleMaxDistWithH = simpleMaxDistWithH_;
        Select(7);
    }

    public void SetSimpleMinDistNonH(double simpleMinDistNonH_)
    {
        simpleMinDistNonH = simpleMinDistNonH_;
        Select(8);
    }

    public void SetSimpleMaxDistNonH(double simpleMaxDistNonH_)
    {
        simpleMaxDistNonH = simpleMaxDistNonH_;
        Select(9);
    }

    public void SetMaxBondsClamp(int maxBondsClamp_)
    {
        maxBondsClamp = maxBondsClamp_;
        Select(10);
    }

    // Property getting methods
    public String  GetElementVariable() { return elementVariable; }
    public Vector  GetAtomicNumber1() { return atomicNumber1; }
    public Vector  GetAtomicNumber2() { return atomicNumber2; }
    public Vector  GetMinDist() { return minDist; }
    public Vector  GetMaxDist() { return maxDist; }
    public boolean GetSimpleAlgorithmFlag() { return simpleAlgorithmFlag; }
    public double  GetSimpleMinDistWithH() { return simpleMinDistWithH; }
    public double  GetSimpleMaxDistWithH() { return simpleMaxDistWithH; }
    public double  GetSimpleMinDistNonH() { return simpleMinDistNonH; }
    public double  GetSimpleMaxDistNonH() { return simpleMaxDistNonH; }
    public int     GetMaxBondsClamp() { return maxBondsClamp; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(elementVariable);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(atomicNumber1);
        if(WriteSelect(2, buf))
            buf.WriteIntVector(atomicNumber2);
        if(WriteSelect(3, buf))
            buf.WriteDoubleVector(minDist);
        if(WriteSelect(4, buf))
            buf.WriteDoubleVector(maxDist);
        if(WriteSelect(5, buf))
            buf.WriteBool(simpleAlgorithmFlag);
        if(WriteSelect(6, buf))
            buf.WriteDouble(simpleMinDistWithH);
        if(WriteSelect(7, buf))
            buf.WriteDouble(simpleMaxDistWithH);
        if(WriteSelect(8, buf))
            buf.WriteDouble(simpleMinDistNonH);
        if(WriteSelect(9, buf))
            buf.WriteDouble(simpleMaxDistNonH);
        if(WriteSelect(10, buf))
            buf.WriteInt(maxBondsClamp);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetElementVariable(buf.ReadString());
                break;
            case 1:
                SetAtomicNumber1(buf.ReadIntVector());
                break;
            case 2:
                SetAtomicNumber2(buf.ReadIntVector());
                break;
            case 3:
                SetMinDist(buf.ReadDoubleVector());
                break;
            case 4:
                SetMaxDist(buf.ReadDoubleVector());
                break;
            case 5:
                SetSimpleAlgorithmFlag(buf.ReadBool());
                break;
            case 6:
                SetSimpleMinDistWithH(buf.ReadDouble());
                break;
            case 7:
                SetSimpleMaxDistWithH(buf.ReadDouble());
                break;
            case 8:
                SetSimpleMinDistNonH(buf.ReadDouble());
                break;
            case 9:
                SetSimpleMaxDistNonH(buf.ReadDouble());
                break;
            case 10:
                SetMaxBondsClamp(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private String  elementVariable;
    private Vector  atomicNumber1; // vector of Integer objects
    private Vector  atomicNumber2; // vector of Integer objects
    private Vector  minDist; // vector of Double objects
    private Vector  maxDist; // vector of Double objects
    private boolean simpleAlgorithmFlag;
    private double  simpleMinDistWithH;
    private double  simpleMaxDistWithH;
    private double  simpleMinDistNonH;
    private double  simpleMaxDistNonH;
    private int     maxBondsClamp;
}

