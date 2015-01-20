// ***************************************************************************
//
// Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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


// ****************************************************************************
// Class: SimulationUIValues
//
// Purpose:
//    Contains UI values from a simulation.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SimulationUIValues extends AttributeSubject
{
    private static int SimulationUIValues_numAdditionalAtts = 6;

    public SimulationUIValues()
    {
        super(SimulationUIValues_numAdditionalAtts);

        host = new String("");
        sim = new String("");
        name = new String("");
        ivalue = 0;
        svalue = new String("");
        enabled = true;
    }

    public SimulationUIValues(int nMoreFields)
    {
        super(SimulationUIValues_numAdditionalAtts + nMoreFields);

        host = new String("");
        sim = new String("");
        name = new String("");
        ivalue = 0;
        svalue = new String("");
        enabled = true;
    }

    public SimulationUIValues(SimulationUIValues obj)
    {
        super(SimulationUIValues_numAdditionalAtts);

        host = new String(obj.host);
        sim = new String(obj.sim);
        name = new String(obj.name);
        ivalue = obj.ivalue;
        svalue = new String(obj.svalue);
        enabled = obj.enabled;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SimulationUIValues_numAdditionalAtts;
    }

    public boolean equals(SimulationUIValues obj)
    {
        // Create the return value
        return ((host.equals(obj.host)) &&
                (sim.equals(obj.sim)) &&
                (name.equals(obj.name)) &&
                (ivalue == obj.ivalue) &&
                (svalue.equals(obj.svalue)) &&
                (enabled == obj.enabled));
    }

    // Property setting methods
    public void SetHost(String host_)
    {
        host = host_;
        Select(0);
    }

    public void SetSim(String sim_)
    {
        sim = sim_;
        Select(1);
    }

    public void SetName(String name_)
    {
        name = name_;
        Select(2);
    }

    public void SetIvalue(int ivalue_)
    {
        ivalue = ivalue_;
        Select(3);
    }

    public void SetSvalue(String svalue_)
    {
        svalue = svalue_;
        Select(4);
    }

    public void SetEnabled(boolean enabled_)
    {
        enabled = enabled_;
        Select(5);
    }

    // Property getting methods
    public String  GetHost() { return host; }
    public String  GetSim() { return sim; }
    public String  GetName() { return name; }
    public int     GetIvalue() { return ivalue; }
    public String  GetSvalue() { return svalue; }
    public boolean GetEnabled() { return enabled; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(host);
        if(WriteSelect(1, buf))
            buf.WriteString(sim);
        if(WriteSelect(2, buf))
            buf.WriteString(name);
        if(WriteSelect(3, buf))
            buf.WriteInt(ivalue);
        if(WriteSelect(4, buf))
            buf.WriteString(svalue);
        if(WriteSelect(5, buf))
            buf.WriteBool(enabled);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetHost(buf.ReadString());
            break;
        case 1:
            SetSim(buf.ReadString());
            break;
        case 2:
            SetName(buf.ReadString());
            break;
        case 3:
            SetIvalue(buf.ReadInt());
            break;
        case 4:
            SetSvalue(buf.ReadString());
            break;
        case 5:
            SetEnabled(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("host", host, indent) + "\n";
        str = str + stringToString("sim", sim, indent) + "\n";
        str = str + stringToString("name", name, indent) + "\n";
        str = str + intToString("ivalue", ivalue, indent) + "\n";
        str = str + stringToString("svalue", svalue, indent) + "\n";
        str = str + boolToString("enabled", enabled, indent) + "\n";
        return str;
    }


    // Attributes
    private String  host;
    private String  sim;
    private String  name;
    private int     ivalue;
    private String  svalue;
    private boolean enabled;
}

