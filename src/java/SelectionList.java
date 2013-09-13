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

import java.util.Vector;

// ****************************************************************************
// Class: SelectionList
//
// Purpose:
//    Contains a list of SelectionProperties objects.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SelectionList extends AttributeSubject
{
    private static int SelectionList_numAdditionalAtts = 3;

    public SelectionList()
    {
        super(SelectionList_numAdditionalAtts);

        autoApplyUpdates = false;
        selections = new Vector();
        selectionSummary = new Vector();
    }

    public SelectionList(int nMoreFields)
    {
        super(SelectionList_numAdditionalAtts + nMoreFields);

        autoApplyUpdates = false;
        selections = new Vector();
        selectionSummary = new Vector();
    }

    public SelectionList(SelectionList obj)
    {
        super(SelectionList_numAdditionalAtts);

        int i;

        autoApplyUpdates = obj.autoApplyUpdates;
        // *** Copy the selections field ***
        selections = new Vector(obj.selections.size());
        for(i = 0; i < obj.selections.size(); ++i)
        {
            SelectionProperties oldObj = (SelectionProperties)obj.selections.elementAt(i);
            selections.addElement(new SelectionProperties(oldObj));
        }

        // *** Copy the selectionSummary field ***
        selectionSummary = new Vector(obj.selectionSummary.size());
        for(i = 0; i < obj.selectionSummary.size(); ++i)
        {
            SelectionSummary oldObj = (SelectionSummary)obj.selectionSummary.elementAt(i);
            selectionSummary.addElement(new SelectionSummary(oldObj));
        }


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SelectionList_numAdditionalAtts;
    }

    public boolean equals(SelectionList obj)
    {
        int i;

        // Compare the elements in the selections vector.
        boolean selections_equal = (obj.selections.size() == selections.size());
        for(i = 0; (i < selections.size()) && selections_equal; ++i)
        {
            // Make references to SelectionProperties from Object.
            SelectionProperties selections1 = (SelectionProperties)selections.elementAt(i);
            SelectionProperties selections2 = (SelectionProperties)obj.selections.elementAt(i);
            selections_equal = selections1.equals(selections2);
        }
        // Compare the elements in the selectionSummary vector.
        boolean selectionSummary_equal = (obj.selectionSummary.size() == selectionSummary.size());
        for(i = 0; (i < selectionSummary.size()) && selectionSummary_equal; ++i)
        {
            // Make references to SelectionSummary from Object.
            SelectionSummary selectionSummary1 = (SelectionSummary)selectionSummary.elementAt(i);
            SelectionSummary selectionSummary2 = (SelectionSummary)obj.selectionSummary.elementAt(i);
            selectionSummary_equal = selectionSummary1.equals(selectionSummary2);
        }
        // Create the return value
        return ((autoApplyUpdates == obj.autoApplyUpdates) &&
                selections_equal &&
                selectionSummary_equal);
    }

    // Property setting methods
    public void SetAutoApplyUpdates(boolean autoApplyUpdates_)
    {
        autoApplyUpdates = autoApplyUpdates_;
        Select(0);
    }

    // Property getting methods
    public boolean GetAutoApplyUpdates() { return autoApplyUpdates; }
    public Vector  GetSelections() { return selections; }
    public Vector  GetSelectionSummary() { return selectionSummary; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(autoApplyUpdates);
        if(WriteSelect(1, buf))
        {
            buf.WriteInt(selections.size());
            for(int i = 0; i < selections.size(); ++i)
            {
                SelectionProperties tmp = (SelectionProperties)selections.elementAt(i);
                tmp.Write(buf);
            }
        }
        if(WriteSelect(2, buf))
        {
            buf.WriteInt(selectionSummary.size());
            for(int i = 0; i < selectionSummary.size(); ++i)
            {
                SelectionSummary tmp = (SelectionSummary)selectionSummary.elementAt(i);
                tmp.Write(buf);
            }
        }
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetAutoApplyUpdates(buf.ReadBool());
            break;
        case 1:
            {
                int len = buf.ReadInt();
                selections.clear();
                for(int j = 0; j < len; ++j)
                {
                    SelectionProperties tmp = new SelectionProperties();
                    tmp.Read(buf);
                    selections.addElement(tmp);
                }
            }
            Select(1);
            break;
        case 2:
            {
                int len = buf.ReadInt();
                selectionSummary.clear();
                for(int j = 0; j < len; ++j)
                {
                    SelectionSummary tmp = new SelectionSummary();
                    tmp.Read(buf);
                    selectionSummary.addElement(tmp);
                }
            }
            Select(2);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("autoApplyUpdates", autoApplyUpdates, indent) + "\n";
        str = str + indent + "selections = {\n";
        for(int i = 0; i < selections.size(); ++i)
        {
            AttributeSubject s = (AttributeSubject)selections.elementAt(i);
            str = str + s.toString(indent + "    ");
            if(i < selections.size()-1)
                str = str + ", ";
            str = str + "\n";
        }
        str = str + "}\n";
        str = str + indent + "selectionSummary = {\n";
        for(int i = 0; i < selectionSummary.size(); ++i)
        {
            AttributeSubject s = (AttributeSubject)selectionSummary.elementAt(i);
            str = str + s.toString(indent + "    ");
            if(i < selectionSummary.size()-1)
                str = str + ", ";
            str = str + "\n";
        }
        str = str + "}\n";
        return str;
    }

    // Attributegroup convenience methods
    public void AddSelections(SelectionProperties obj)
    {
        selections.addElement(new SelectionProperties(obj));
        Select(1);
    }

    public void ClearSelections()
    {
        selections.clear();
        Select(1);
    }

    public void RemoveSelections(int index)
    {
        if(index >= 0 && index < selections.size())
        {
            selections.remove(index);
            Select(1);
        }
    }

    public int GetNumSelections()
    {
        return selections.size();
    }

    public SelectionProperties GetSelections(int i)
    {
        SelectionProperties tmp = (SelectionProperties)selections.elementAt(i);
        return tmp;
    }

    public void AddSelectionSummary(SelectionSummary obj)
    {
        selectionSummary.addElement(new SelectionSummary(obj));
        Select(2);
    }

    public void ClearSelectionSummarys()
    {
        selectionSummary.clear();
        Select(2);
    }

    public void RemoveSelectionSummary(int index)
    {
        if(index >= 0 && index < selectionSummary.size())
        {
            selectionSummary.remove(index);
            Select(2);
        }
    }

    public int GetNumSelectionSummarys()
    {
        return selectionSummary.size();
    }

    public SelectionSummary GetSelectionSummary(int i)
    {
        SelectionSummary tmp = (SelectionSummary)selectionSummary.elementAt(i);
        return tmp;
    }


    // Attributes
    private boolean autoApplyUpdates;
    private Vector  selections; // vector of SelectionProperties objects
    private Vector  selectionSummary; // vector of SelectionSummary objects
}

