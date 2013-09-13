// ***************************************************************************
//
// Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: ProcessAttributes
//
// Purpose:
//    attributes to describe a running process
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ProcessAttributes extends AttributeSubject
{
    private static int ProcessAttributes_numAdditionalAtts = 5;

    public ProcessAttributes()
    {
        super(ProcessAttributes_numAdditionalAtts);

        pids = new Vector();
        ppids = new Vector();
        hosts = new Vector();
        isParallel = false;
        memory = new Vector();
    }

    public ProcessAttributes(int nMoreFields)
    {
        super(ProcessAttributes_numAdditionalAtts + nMoreFields);

        pids = new Vector();
        ppids = new Vector();
        hosts = new Vector();
        isParallel = false;
        memory = new Vector();
    }

    public ProcessAttributes(ProcessAttributes obj)
    {
        super(ProcessAttributes_numAdditionalAtts);

        int i;

        pids = new Vector();
        for(i = 0; i < obj.pids.size(); ++i)
        {
            Integer iv = (Integer)obj.pids.elementAt(i);
            pids.addElement(new Integer(iv.intValue()));
        }
        ppids = new Vector();
        for(i = 0; i < obj.ppids.size(); ++i)
        {
            Integer iv = (Integer)obj.ppids.elementAt(i);
            ppids.addElement(new Integer(iv.intValue()));
        }
        hosts = new Vector(obj.hosts.size());
        for(i = 0; i < obj.hosts.size(); ++i)
            hosts.addElement(new String((String)obj.hosts.elementAt(i)));

        isParallel = obj.isParallel;
        memory = new Vector();
        for(i = 0; i < obj.memory.size(); ++i)
        {
            Integer iv = (Integer)obj.memory.elementAt(i);
            memory.addElement(new Integer(iv.intValue()));
        }

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ProcessAttributes_numAdditionalAtts;
    }

    public boolean equals(ProcessAttributes obj)
    {
        int i;

        // Compare the elements in the pids vector.
        boolean pids_equal = (obj.pids.size() == pids.size());
        for(i = 0; (i < pids.size()) && pids_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer pids1 = (Integer)pids.elementAt(i);
            Integer pids2 = (Integer)obj.pids.elementAt(i);
            pids_equal = pids1.equals(pids2);
        }
        // Compare the elements in the ppids vector.
        boolean ppids_equal = (obj.ppids.size() == ppids.size());
        for(i = 0; (i < ppids.size()) && ppids_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer ppids1 = (Integer)ppids.elementAt(i);
            Integer ppids2 = (Integer)obj.ppids.elementAt(i);
            ppids_equal = ppids1.equals(ppids2);
        }
        // Compare the elements in the hosts vector.
        boolean hosts_equal = (obj.hosts.size() == hosts.size());
        for(i = 0; (i < hosts.size()) && hosts_equal; ++i)
        {
            // Make references to String from Object.
            String hosts1 = (String)hosts.elementAt(i);
            String hosts2 = (String)obj.hosts.elementAt(i);
            hosts_equal = hosts1.equals(hosts2);
        }
        // Compare the elements in the memory vector.
        boolean memory_equal = (obj.memory.size() == memory.size());
        for(i = 0; (i < memory.size()) && memory_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer memory1 = (Integer)memory.elementAt(i);
            Integer memory2 = (Integer)obj.memory.elementAt(i);
            memory_equal = memory1.equals(memory2);
        }
        // Create the return value
        return (pids_equal &&
                ppids_equal &&
                hosts_equal &&
                (isParallel == obj.isParallel) &&
                memory_equal);
    }

    // Property setting methods
    public void SetPids(Vector pids_)
    {
        pids = pids_;
        Select(0);
    }

    public void SetPpids(Vector ppids_)
    {
        ppids = ppids_;
        Select(1);
    }

    public void SetHosts(Vector hosts_)
    {
        hosts = hosts_;
        Select(2);
    }

    public void SetIsParallel(boolean isParallel_)
    {
        isParallel = isParallel_;
        Select(3);
    }

    public void SetMemory(Vector memory_)
    {
        memory = memory_;
        Select(4);
    }

    // Property getting methods
    public Vector  GetPids() { return pids; }
    public Vector  GetPpids() { return ppids; }
    public Vector  GetHosts() { return hosts; }
    public boolean GetIsParallel() { return isParallel; }
    public Vector  GetMemory() { return memory; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteIntVector(pids);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(ppids);
        if(WriteSelect(2, buf))
            buf.WriteStringVector(hosts);
        if(WriteSelect(3, buf))
            buf.WriteBool(isParallel);
        if(WriteSelect(4, buf))
            buf.WriteIntVector(memory);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetPids(buf.ReadIntVector());
            break;
        case 1:
            SetPpids(buf.ReadIntVector());
            break;
        case 2:
            SetHosts(buf.ReadStringVector());
            break;
        case 3:
            SetIsParallel(buf.ReadBool());
            break;
        case 4:
            SetMemory(buf.ReadIntVector());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intVectorToString("pids", pids, indent) + "\n";
        str = str + intVectorToString("ppids", ppids, indent) + "\n";
        str = str + stringVectorToString("hosts", hosts, indent) + "\n";
        str = str + boolToString("isParallel", isParallel, indent) + "\n";
        str = str + intVectorToString("memory", memory, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector  pids; // vector of Integer objects
    private Vector  ppids; // vector of Integer objects
    private Vector  hosts; // vector of String objects
    private boolean isParallel;
    private Vector  memory; // vector of Integer objects
}

