// ***************************************************************************
//
// Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

// ****************************************************************************
// Class: VisualCueList
//
// Purpose:
//    container object for shipping vectors of visual cues
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class VisualCueList extends AttributeSubject
{
    private static int VisualCueList_numAdditionalAtts = 1;

    public VisualCueList()
    {
        super(VisualCueList_numAdditionalAtts);

        cues = new Vector();
    }

    public VisualCueList(int nMoreFields)
    {
        super(VisualCueList_numAdditionalAtts + nMoreFields);

        cues = new Vector();
    }

    public VisualCueList(VisualCueList obj)
    {
        super(VisualCueList_numAdditionalAtts);

        int i;

        // *** Copy the cues field ***
        cues = new Vector(obj.cues.size());
        for(i = 0; i < obj.cues.size(); ++i)
        {
            VisualCueInfo oldObj = (VisualCueInfo)obj.cues.elementAt(i);
            cues.addElement(new VisualCueInfo(oldObj));
        }


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return VisualCueList_numAdditionalAtts;
    }

    public boolean equals(VisualCueList obj)
    {
        int i;

        // Compare the elements in the cues vector.
        boolean cues_equal = (obj.cues.size() == cues.size());
        for(i = 0; (i < cues.size()) && cues_equal; ++i)
        {
            // Make references to VisualCueInfo from Object.
            VisualCueInfo cues1 = (VisualCueInfo)cues.elementAt(i);
            VisualCueInfo cues2 = (VisualCueInfo)obj.cues.elementAt(i);
            cues_equal = cues1.equals(cues2);
        }
        // Create the return value
        return (cues_equal);
    }

    // Property setting methods
    // Property getting methods
    public Vector GetCues() { return cues; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
        {
            buf.WriteInt(cues.size());
            for(int i = 0; i < cues.size(); ++i)
            {
                VisualCueInfo tmp = (VisualCueInfo)cues.elementAt(i);
                tmp.Write(buf);
            }
        }
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        {
            int len = buf.ReadInt();
            cues.clear();
            for(int j = 0; j < len; ++j)
            {
                VisualCueInfo tmp = new VisualCueInfo();
                tmp.Read(buf);
                cues.addElement(tmp);
            }
        }
        Select(0);
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "cues = {\n";
        for(int i = 0; i < cues.size(); ++i)
        {
            AttributeSubject s = (AttributeSubject)cues.elementAt(i);
            str = str + s.toString(indent + "    ");
            if(i < cues.size()-1)
                str = str + ", ";
            str = str + "\n";
        }
        str = str + "}\n";
        return str;
    }

    // Attributegroup convenience methods
    public void AddCues(VisualCueInfo obj)
    {
        cues.addElement(new VisualCueInfo(obj));
        Select(0);
    }

    public void ClearCues()
    {
        cues.clear();
        Select(0);
    }

    public void RemoveCues(int index)
    {
        if(index >= 0 && index < cues.size())
        {
            cues.remove(index);
            Select(0);
        }
    }

    public int GetNumCues()
    {
        return cues.size();
    }

    public VisualCueInfo GetCues(int i)
    {
        VisualCueInfo tmp = (VisualCueInfo)cues.elementAt(i);
        return tmp;
    }


    // Attributes
    private Vector cues; // vector of VisualCueInfo objects
}

