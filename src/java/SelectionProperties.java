// ***************************************************************************
//
// Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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
// Class: SelectionProperties
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

public class SelectionProperties extends AttributeSubject
{
    private static int SelectionProperties_numAdditionalAtts = 5;

    public SelectionProperties()
    {
        super(SelectionProperties_numAdditionalAtts);

        name = new String("");
        originatingPlot = new String("");
        rangeProperty = 0;
        histogramProperty = 0;
        statisticsProperty = 0;
    }

    public SelectionProperties(int nMoreFields)
    {
        super(SelectionProperties_numAdditionalAtts + nMoreFields);

        name = new String("");
        originatingPlot = new String("");
        rangeProperty = 0;
        histogramProperty = 0;
        statisticsProperty = 0;
    }

    public SelectionProperties(SelectionProperties obj)
    {
        super(SelectionProperties_numAdditionalAtts);

        name = new String(obj.name);
        originatingPlot = new String(obj.originatingPlot);
        rangeProperty = obj.rangeProperty;
        histogramProperty = obj.histogramProperty;
        statisticsProperty = obj.statisticsProperty;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SelectionProperties_numAdditionalAtts;
    }

    public boolean equals(SelectionProperties obj)
    {
        // Create the return value
        return ((name.equals(obj.name)) &&
                (originatingPlot.equals(obj.originatingPlot)) &&
                (rangeProperty == obj.rangeProperty) &&
                (histogramProperty == obj.histogramProperty) &&
                (statisticsProperty == obj.statisticsProperty));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetOriginatingPlot(String originatingPlot_)
    {
        originatingPlot = originatingPlot_;
        Select(1);
    }

    public void SetRangeProperty(int rangeProperty_)
    {
        rangeProperty = rangeProperty_;
        Select(2);
    }

    public void SetHistogramProperty(int histogramProperty_)
    {
        histogramProperty = histogramProperty_;
        Select(3);
    }

    public void SetStatisticsProperty(int statisticsProperty_)
    {
        statisticsProperty = statisticsProperty_;
        Select(4);
    }

    // Property getting methods
    public String GetName() { return name; }
    public String GetOriginatingPlot() { return originatingPlot; }
    public int    GetRangeProperty() { return rangeProperty; }
    public int    GetHistogramProperty() { return histogramProperty; }
    public int    GetStatisticsProperty() { return statisticsProperty; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteString(originatingPlot);
        if(WriteSelect(2, buf))
            buf.WriteInt(rangeProperty);
        if(WriteSelect(3, buf))
            buf.WriteInt(histogramProperty);
        if(WriteSelect(4, buf))
            buf.WriteInt(statisticsProperty);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetName(buf.ReadString());
            break;
        case 1:
            SetOriginatingPlot(buf.ReadString());
            break;
        case 2:
            SetRangeProperty(buf.ReadInt());
            break;
        case 3:
            SetHistogramProperty(buf.ReadInt());
            break;
        case 4:
            SetStatisticsProperty(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + stringToString("originatingPlot", originatingPlot, indent) + "\n";
        str = str + intToString("rangeProperty", rangeProperty, indent) + "\n";
        str = str + intToString("histogramProperty", histogramProperty, indent) + "\n";
        str = str + intToString("statisticsProperty", statisticsProperty, indent) + "\n";
        return str;
    }


    // Attributes
    private String name;
    private String originatingPlot;
    private int    rangeProperty;
    private int    histogramProperty;
    private int    statisticsProperty;
}

