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

import java.util.Vector;
import java.lang.Double;

// ****************************************************************************
// Class: SelectionSummary
//
// Purpose:
//    Contains attributes for a selection
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SelectionSummary extends AttributeSubject
{
    private static int SelectionSummary_numAdditionalAtts = 7;

    public SelectionSummary()
    {
        super(SelectionSummary_numAdditionalAtts);

        name = new String("");
        variables = new Vector();
        cellCount = 0;
        totalCellCount = 0;
        histogramValues = new Vector();
        histogramMinBin = 0;
        histogramMaxBin = 1;
    }

    public SelectionSummary(int nMoreFields)
    {
        super(SelectionSummary_numAdditionalAtts + nMoreFields);

        name = new String("");
        variables = new Vector();
        cellCount = 0;
        totalCellCount = 0;
        histogramValues = new Vector();
        histogramMinBin = 0;
        histogramMaxBin = 1;
    }

    public SelectionSummary(SelectionSummary obj)
    {
        super(SelectionSummary_numAdditionalAtts);

        int i;

        name = new String(obj.name);
        // *** Copy the variables field ***
        variables = new Vector(obj.variables.size());
        for(i = 0; i < obj.variables.size(); ++i)
        {
            SelectionVariableSummary oldObj = (SelectionVariableSummary)obj.variables.elementAt(i);
            variables.addElement(new SelectionVariableSummary(oldObj));
        }

        cellCount = obj.cellCount;
        totalCellCount = obj.totalCellCount;
        histogramValues = new Vector(obj.histogramValues.size());
        for(i = 0; i < obj.histogramValues.size(); ++i)
        {
            Double dv = (Double)obj.histogramValues.elementAt(i);
            histogramValues.addElement(new Double(dv.doubleValue()));
        }

        histogramMinBin = obj.histogramMinBin;
        histogramMaxBin = obj.histogramMaxBin;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SelectionSummary_numAdditionalAtts;
    }

    public boolean equals(SelectionSummary obj)
    {
        int i;

        // Compare the elements in the variables vector.
        boolean variables_equal = (obj.variables.size() == variables.size());
        for(i = 0; (i < variables.size()) && variables_equal; ++i)
        {
            // Make references to SelectionVariableSummary from Object.
            SelectionVariableSummary variables1 = (SelectionVariableSummary)variables.elementAt(i);
            SelectionVariableSummary variables2 = (SelectionVariableSummary)obj.variables.elementAt(i);
            variables_equal = variables1.equals(variables2);
        }
        // Compare the elements in the histogramValues vector.
        boolean histogramValues_equal = (obj.histogramValues.size() == histogramValues.size());
        for(i = 0; (i < histogramValues.size()) && histogramValues_equal; ++i)
        {
            // Make references to Double from Object.
            Double histogramValues1 = (Double)histogramValues.elementAt(i);
            Double histogramValues2 = (Double)obj.histogramValues.elementAt(i);
            histogramValues_equal = histogramValues1.equals(histogramValues2);
        }
        // Create the return value
        return ((name.equals(obj.name)) &&
                variables_equal &&
                (cellCount == obj.cellCount) &&
                (totalCellCount == obj.totalCellCount) &&
                histogramValues_equal &&
                (histogramMinBin == obj.histogramMinBin) &&
                (histogramMaxBin == obj.histogramMaxBin));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetCellCount(int cellCount_)
    {
        cellCount = cellCount_;
        Select(2);
    }

    public void SetTotalCellCount(int totalCellCount_)
    {
        totalCellCount = totalCellCount_;
        Select(3);
    }

    public void SetHistogramValues(Vector histogramValues_)
    {
        histogramValues = histogramValues_;
        Select(4);
    }

    public void SetHistogramMinBin(double histogramMinBin_)
    {
        histogramMinBin = histogramMinBin_;
        Select(5);
    }

    public void SetHistogramMaxBin(double histogramMaxBin_)
    {
        histogramMaxBin = histogramMaxBin_;
        Select(6);
    }

    // Property getting methods
    public String GetName() { return name; }
    public Vector GetVariables() { return variables; }
    public int    GetCellCount() { return cellCount; }
    public int    GetTotalCellCount() { return totalCellCount; }
    public Vector GetHistogramValues() { return histogramValues; }
    public double GetHistogramMinBin() { return histogramMinBin; }
    public double GetHistogramMaxBin() { return histogramMaxBin; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
        {
            buf.WriteInt(variables.size());
            for(int i = 0; i < variables.size(); ++i)
            {
                SelectionVariableSummary tmp = (SelectionVariableSummary)variables.elementAt(i);
                tmp.Write(buf);
            }
        }
        if(WriteSelect(2, buf))
            buf.WriteInt(cellCount);
        if(WriteSelect(3, buf))
            buf.WriteInt(totalCellCount);
        if(WriteSelect(4, buf))
            buf.WriteDoubleVector(histogramValues);
        if(WriteSelect(5, buf))
            buf.WriteDouble(histogramMinBin);
        if(WriteSelect(6, buf))
            buf.WriteDouble(histogramMaxBin);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetName(buf.ReadString());
            break;
        case 1:
            {
                int len = buf.ReadInt();
                variables.clear();
                for(int j = 0; j < len; ++j)
                {
                    SelectionVariableSummary tmp = new SelectionVariableSummary();
                    tmp.Read(buf);
                    variables.addElement(tmp);
                }
            }
            Select(1);
            break;
        case 2:
            SetCellCount(buf.ReadInt());
            break;
        case 3:
            SetTotalCellCount(buf.ReadInt());
            break;
        case 4:
            SetHistogramValues(buf.ReadDoubleVector());
            break;
        case 5:
            SetHistogramMinBin(buf.ReadDouble());
            break;
        case 6:
            SetHistogramMaxBin(buf.ReadDouble());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + indent + "variables = {\n";
        for(int i = 0; i < variables.size(); ++i)
        {
            AttributeSubject s = (AttributeSubject)variables.elementAt(i);
            str = str + s.toString(indent + "    ");
            if(i < variables.size()-1)
                str = str + ", ";
            str = str + "\n";
        }
        str = str + "}\n";
        str = str + intToString("cellCount", cellCount, indent) + "\n";
        str = str + intToString("totalCellCount", totalCellCount, indent) + "\n";
        str = str + doubleVectorToString("histogramValues", histogramValues, indent) + "\n";
        str = str + doubleToString("histogramMinBin", histogramMinBin, indent) + "\n";
        str = str + doubleToString("histogramMaxBin", histogramMaxBin, indent) + "\n";
        return str;
    }

    // Attributegroup convenience methods
    public void AddVariables(SelectionVariableSummary obj)
    {
        variables.addElement(new SelectionVariableSummary(obj));
        Select(1);
    }

    public void ClearVariables()
    {
        variables.clear();
        Select(1);
    }

    public void RemoveVariables(int index)
    {
        if(index >= 0 && index < variables.size())
        {
            variables.remove(index);
            Select(1);
        }
    }

    public int GetNumVariables()
    {
        return variables.size();
    }

    public SelectionVariableSummary GetVariables(int i)
    {
        SelectionVariableSummary tmp = (SelectionVariableSummary)variables.elementAt(i);
        return tmp;
    }


    // Attributes
    private String name;
    private Vector variables; // vector of SelectionVariableSummary objects
    private int    cellCount;
    private int    totalCellCount;
    private Vector histogramValues; // vector of Double objects
    private double histogramMinBin;
    private double histogramMaxBin;
}

