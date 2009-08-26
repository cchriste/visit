// ***************************************************************************
//
// Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400124
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
// Class: BoxExtents
//
// Purpose:
//    Attributes for an axis-aligned box
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class BoxExtents extends AttributeSubject
{
    private static int numAdditionalAttributes = 1;

    public BoxExtents()
    {
        super(numAdditionalAttributes);

        extents = new double[6];
        extents[0] = 0;
        extents[1] = 0;
        extents[2] = 0;
        extents[3] = 0;
        extents[4] = 0;
        extents[5] = 0;
    }

    public BoxExtents(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

        extents = new double[6];
        extents[0] = 0;
        extents[1] = 0;
        extents[2] = 0;
        extents[3] = 0;
        extents[4] = 0;
        extents[5] = 0;
    }

    public BoxExtents(BoxExtents obj)
    {
        super(numAdditionalAttributes);

        int i;

        extents = new double[6];
        for(i = 0; i < obj.extents.length; ++i)
            extents[i] = obj.extents[i];


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return numAdditionalAttributes;
    }

    public boolean equals(BoxExtents obj)
    {
        int i;

        // Compare the extents arrays.
        boolean extents_equal = true;
        for(i = 0; i < 6 && extents_equal; ++i)
            extents_equal = (extents[i] == obj.extents[i]);

        // Create the return value
        return (extents_equal);
    }

    // Property setting methods
    public void SetExtents(double[] extents_)
    {
        for(int i = 0; i < 6; ++i)
             extents[i] = extents_[i];
        Select(0);
    }

    // Property getting methods
    public double[] GetExtents() { return extents; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(extents);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        SetExtents(buf.ReadDoubleArray());
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("extents", extents, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] extents;
}

