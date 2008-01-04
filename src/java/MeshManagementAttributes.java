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

import java.lang.Double;
import java.util.Vector;

// ****************************************************************************
// Class: MeshManagementAttributes
//
// Purpose:
//    Global variables controlling reading and conversion of non-standard meshes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:55:00 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class MeshManagementAttributes extends AttributeSubject
{
    // Enum values
    public final static int DISCRETIZATIONMODES_UNIFORM = 0;
    public final static int DISCRETIZATIONMODES_ADAPTIVE = 1;


    public MeshManagementAttributes()
    {
        super(7);

        discretizationTolerance = new Vector();
        discretizationTolerance.addElement(new Double(0.01));
        discretizationTolerance.addElement(new Double(0.025));
        discretizationTolerance.addElement(new Double(0.05));
        discretizationToleranceX = new Vector();
        discretizationToleranceY = new Vector();
        discretizationToleranceZ = new Vector();
        discretizationMode = DISCRETIZATIONMODES_ADAPTIVE;
        discretizeBoundaryOnly = false;
        passNativeCSG = false;
    }

    public MeshManagementAttributes(MeshManagementAttributes obj)
    {
        super(7);

        int i;

        discretizationTolerance = new Vector(obj.discretizationTolerance.size());
        for(i = 0; i < obj.discretizationTolerance.size(); ++i)
        {
            Double dv = (Double)obj.discretizationTolerance.elementAt(i);
            discretizationTolerance.addElement(new Double(dv.doubleValue()));
        }

        discretizationToleranceX = new Vector(obj.discretizationToleranceX.size());
        for(i = 0; i < obj.discretizationToleranceX.size(); ++i)
        {
            Double dv = (Double)obj.discretizationToleranceX.elementAt(i);
            discretizationToleranceX.addElement(new Double(dv.doubleValue()));
        }

        discretizationToleranceY = new Vector(obj.discretizationToleranceY.size());
        for(i = 0; i < obj.discretizationToleranceY.size(); ++i)
        {
            Double dv = (Double)obj.discretizationToleranceY.elementAt(i);
            discretizationToleranceY.addElement(new Double(dv.doubleValue()));
        }

        discretizationToleranceZ = new Vector(obj.discretizationToleranceZ.size());
        for(i = 0; i < obj.discretizationToleranceZ.size(); ++i)
        {
            Double dv = (Double)obj.discretizationToleranceZ.elementAt(i);
            discretizationToleranceZ.addElement(new Double(dv.doubleValue()));
        }

        discretizationMode = obj.discretizationMode;
        discretizeBoundaryOnly = obj.discretizeBoundaryOnly;
        passNativeCSG = obj.passNativeCSG;

        SelectAll();
    }

    public boolean equals(MeshManagementAttributes obj)
    {
        int i;

        // Create the return value
        return ((discretizationTolerance == obj.discretizationTolerance) &&
                (discretizationToleranceX == obj.discretizationToleranceX) &&
                (discretizationToleranceY == obj.discretizationToleranceY) &&
                (discretizationToleranceZ == obj.discretizationToleranceZ) &&
                (discretizationMode == obj.discretizationMode) &&
                (discretizeBoundaryOnly == obj.discretizeBoundaryOnly) &&
                (passNativeCSG == obj.passNativeCSG));
    }

    // Property setting methods
    public void SetDiscretizationTolerance(Vector discretizationTolerance_)
    {
        discretizationTolerance = discretizationTolerance_;
        Select(0);
    }

    public void SetDiscretizationToleranceX(Vector discretizationToleranceX_)
    {
        discretizationToleranceX = discretizationToleranceX_;
        Select(1);
    }

    public void SetDiscretizationToleranceY(Vector discretizationToleranceY_)
    {
        discretizationToleranceY = discretizationToleranceY_;
        Select(2);
    }

    public void SetDiscretizationToleranceZ(Vector discretizationToleranceZ_)
    {
        discretizationToleranceZ = discretizationToleranceZ_;
        Select(3);
    }

    public void SetDiscretizationMode(int discretizationMode_)
    {
        discretizationMode = discretizationMode_;
        Select(4);
    }

    public void SetDiscretizeBoundaryOnly(boolean discretizeBoundaryOnly_)
    {
        discretizeBoundaryOnly = discretizeBoundaryOnly_;
        Select(5);
    }

    public void SetPassNativeCSG(boolean passNativeCSG_)
    {
        passNativeCSG = passNativeCSG_;
        Select(6);
    }

    // Property getting methods
    public Vector  GetDiscretizationTolerance() { return discretizationTolerance; }
    public Vector  GetDiscretizationToleranceX() { return discretizationToleranceX; }
    public Vector  GetDiscretizationToleranceY() { return discretizationToleranceY; }
    public Vector  GetDiscretizationToleranceZ() { return discretizationToleranceZ; }
    public int     GetDiscretizationMode() { return discretizationMode; }
    public boolean GetDiscretizeBoundaryOnly() { return discretizeBoundaryOnly; }
    public boolean GetPassNativeCSG() { return passNativeCSG; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleVector(discretizationTolerance);
        if(WriteSelect(1, buf))
            buf.WriteDoubleVector(discretizationToleranceX);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(discretizationToleranceY);
        if(WriteSelect(3, buf))
            buf.WriteDoubleVector(discretizationToleranceZ);
        if(WriteSelect(4, buf))
            buf.WriteInt(discretizationMode);
        if(WriteSelect(5, buf))
            buf.WriteBool(discretizeBoundaryOnly);
        if(WriteSelect(6, buf))
            buf.WriteBool(passNativeCSG);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetDiscretizationTolerance(buf.ReadDoubleVector());
                break;
            case 1:
                SetDiscretizationToleranceX(buf.ReadDoubleVector());
                break;
            case 2:
                SetDiscretizationToleranceY(buf.ReadDoubleVector());
                break;
            case 3:
                SetDiscretizationToleranceZ(buf.ReadDoubleVector());
                break;
            case 4:
                SetDiscretizationMode(buf.ReadInt());
                break;
            case 5:
                SetDiscretizeBoundaryOnly(buf.ReadBool());
                break;
            case 6:
                SetPassNativeCSG(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private Vector  discretizationTolerance; // vector of Double objects
    private Vector  discretizationToleranceX; // vector of Double objects
    private Vector  discretizationToleranceY; // vector of Double objects
    private Vector  discretizationToleranceZ; // vector of Double objects
    private int     discretizationMode;
    private boolean discretizeBoundaryOnly;
    private boolean passNativeCSG;
}

