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


// ****************************************************************************
// Class: SelectionVariableSummary
//
// Purpose:
//    Contains a summary of a variable used in a selection
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SelectionVariableSummary extends AttributeSubject
{
    private static int SelectionVariableSummary_numAdditionalAtts = 4;

    public SelectionVariableSummary()
    {
        super(SelectionVariableSummary_numAdditionalAtts);

        name = new String("");
        minimum = 0;
        maximum = 0;
        histogram = new double[256];
        for (int i = 0; i < histogram.length; ++i)
            histogram[i] = 0.;
    }

    public SelectionVariableSummary(int nMoreFields)
    {
        super(SelectionVariableSummary_numAdditionalAtts + nMoreFields);

        name = new String("");
        minimum = 0;
        maximum = 0;
        histogram = new double[256];
        for (int i = 0; i < histogram.length; ++i)
            histogram[i] = 0.;
    }

    public SelectionVariableSummary(SelectionVariableSummary obj)
    {
        super(SelectionVariableSummary_numAdditionalAtts);

        int i;

        name = new String(obj.name);
        minimum = obj.minimum;
        maximum = obj.maximum;
        histogram = new double[256];
        for(i = 0; i < obj.histogram.length; ++i)
            histogram[i] = obj.histogram[i];


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SelectionVariableSummary_numAdditionalAtts;
    }

    public boolean equals(SelectionVariableSummary obj)
    {
        int i;

        // Compare the histogram arrays.
        boolean histogram_equal = true;
        for(i = 0; i < 256 && histogram_equal; ++i)
            histogram_equal = (histogram[i] == obj.histogram[i]);

        // Create the return value
        return ((name.equals(obj.name)) &&
                (minimum == obj.minimum) &&
                (maximum == obj.maximum) &&
                histogram_equal);
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetMinimum(double minimum_)
    {
        minimum = minimum_;
        Select(1);
    }

    public void SetMaximum(double maximum_)
    {
        maximum = maximum_;
        Select(2);
    }

    public void SetHistogram(double[] histogram_)
    {
        for(int i = 0; i < 256; ++i)
             histogram[i] = histogram_[i];
        Select(3);
    }

    // Property getting methods
    public String   GetName() { return name; }
    public double   GetMinimum() { return minimum; }
    public double   GetMaximum() { return maximum; }
    public double[] GetHistogram() { return histogram; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteDouble(minimum);
        if(WriteSelect(2, buf))
            buf.WriteDouble(maximum);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(histogram);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetName(buf.ReadString());
            break;
        case 1:
            SetMinimum(buf.ReadDouble());
            break;
        case 2:
            SetMaximum(buf.ReadDouble());
            break;
        case 3:
            SetHistogram(buf.ReadDoubleArray());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + doubleToString("minimum", minimum, indent) + "\n";
        str = str + doubleToString("maximum", maximum, indent) + "\n";
        str = str + doubleArrayToString("histogram", histogram, indent) + "\n";
        return str;
    }


    // Attributes
    private String   name;
    private double   minimum;
    private double   maximum;
    private double[] histogram;
}

