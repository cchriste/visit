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
import java.lang.Double;

// ****************************************************************************
// Class: AxisRestrictionAttributes
//
// Purpose:
//    Attributes for axis restrictions
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AxisRestrictionAttributes extends AttributeSubject
{
    private static int AxisRestrictionAttributes_numAdditionalAtts = 3;

    public AxisRestrictionAttributes()
    {
        super(AxisRestrictionAttributes_numAdditionalAtts);

        names = new Vector();
        minima = new Vector();
        maxima = new Vector();
    }

    public AxisRestrictionAttributes(int nMoreFields)
    {
        super(AxisRestrictionAttributes_numAdditionalAtts + nMoreFields);

        names = new Vector();
        minima = new Vector();
        maxima = new Vector();
    }

    public AxisRestrictionAttributes(AxisRestrictionAttributes obj)
    {
        super(AxisRestrictionAttributes_numAdditionalAtts);

        int i;

        names = new Vector(obj.names.size());
        for(i = 0; i < obj.names.size(); ++i)
            names.addElement(new String((String)obj.names.elementAt(i)));

        minima = new Vector(obj.minima.size());
        for(i = 0; i < obj.minima.size(); ++i)
        {
            Double dv = (Double)obj.minima.elementAt(i);
            minima.addElement(new Double(dv.doubleValue()));
        }

        maxima = new Vector(obj.maxima.size());
        for(i = 0; i < obj.maxima.size(); ++i)
        {
            Double dv = (Double)obj.maxima.elementAt(i);
            maxima.addElement(new Double(dv.doubleValue()));
        }


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AxisRestrictionAttributes_numAdditionalAtts;
    }

    public boolean equals(AxisRestrictionAttributes obj)
    {
        int i;

        // Compare the elements in the names vector.
        boolean names_equal = (obj.names.size() == names.size());
        for(i = 0; (i < names.size()) && names_equal; ++i)
        {
            // Make references to String from Object.
            String names1 = (String)names.elementAt(i);
            String names2 = (String)obj.names.elementAt(i);
            names_equal = names1.equals(names2);
        }
        // Compare the elements in the minima vector.
        boolean minima_equal = (obj.minima.size() == minima.size());
        for(i = 0; (i < minima.size()) && minima_equal; ++i)
        {
            // Make references to Double from Object.
            Double minima1 = (Double)minima.elementAt(i);
            Double minima2 = (Double)obj.minima.elementAt(i);
            minima_equal = minima1.equals(minima2);
        }
        // Compare the elements in the maxima vector.
        boolean maxima_equal = (obj.maxima.size() == maxima.size());
        for(i = 0; (i < maxima.size()) && maxima_equal; ++i)
        {
            // Make references to Double from Object.
            Double maxima1 = (Double)maxima.elementAt(i);
            Double maxima2 = (Double)obj.maxima.elementAt(i);
            maxima_equal = maxima1.equals(maxima2);
        }
        // Create the return value
        return (names_equal &&
                minima_equal &&
                maxima_equal);
    }

    // Property setting methods
    public void SetNames(Vector names_)
    {
        names = names_;
        Select(0);
    }

    public void SetMinima(Vector minima_)
    {
        minima = minima_;
        Select(1);
    }

    public void SetMaxima(Vector maxima_)
    {
        maxima = maxima_;
        Select(2);
    }

    // Property getting methods
    public Vector GetNames() { return names; }
    public Vector GetMinima() { return minima; }
    public Vector GetMaxima() { return maxima; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(names);
        if(WriteSelect(1, buf))
            buf.WriteDoubleVector(minima);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(maxima);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetNames(buf.ReadStringVector());
            break;
        case 1:
            SetMinima(buf.ReadDoubleVector());
            break;
        case 2:
            SetMaxima(buf.ReadDoubleVector());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringVectorToString("names", names, indent) + "\n";
        str = str + doubleVectorToString("minima", minima, indent) + "\n";
        str = str + doubleVectorToString("maxima", maxima, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector names; // vector of String objects
    private Vector minima; // vector of Double objects
    private Vector maxima; // vector of Double objects
}

