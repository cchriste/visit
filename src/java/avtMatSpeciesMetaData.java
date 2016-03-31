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
// Class: avtMatSpeciesMetaData
//
// Purpose:
//    Contains material species metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class avtMatSpeciesMetaData extends AttributeSubject
{
    private static int avtMatSpeciesMetaData_numAdditionalAtts = 3;

    public avtMatSpeciesMetaData()
    {
        super(avtMatSpeciesMetaData_numAdditionalAtts);

        numSpecies = 0;
        speciesNames = new Vector();
        validVariable = true;
    }

    public avtMatSpeciesMetaData(int nMoreFields)
    {
        super(avtMatSpeciesMetaData_numAdditionalAtts + nMoreFields);

        numSpecies = 0;
        speciesNames = new Vector();
        validVariable = true;
    }

    public avtMatSpeciesMetaData(avtMatSpeciesMetaData obj)
    {
        super(avtMatSpeciesMetaData_numAdditionalAtts);

        int i;

        numSpecies = obj.numSpecies;
        speciesNames = new Vector(obj.speciesNames.size());
        for(i = 0; i < obj.speciesNames.size(); ++i)
            speciesNames.addElement(new String((String)obj.speciesNames.elementAt(i)));

        validVariable = obj.validVariable;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return avtMatSpeciesMetaData_numAdditionalAtts;
    }

    public boolean equals(avtMatSpeciesMetaData obj)
    {
        int i;

        // Compare the elements in the speciesNames vector.
        boolean speciesNames_equal = (obj.speciesNames.size() == speciesNames.size());
        for(i = 0; (i < speciesNames.size()) && speciesNames_equal; ++i)
        {
            // Make references to String from Object.
            String speciesNames1 = (String)speciesNames.elementAt(i);
            String speciesNames2 = (String)obj.speciesNames.elementAt(i);
            speciesNames_equal = speciesNames1.equals(speciesNames2);
        }
        // Create the return value
        return ((numSpecies == obj.numSpecies) &&
                speciesNames_equal &&
                (validVariable == obj.validVariable));
    }

    // Property setting methods
    public void SetNumSpecies(int numSpecies_)
    {
        numSpecies = numSpecies_;
        Select(0);
    }

    public void SetSpeciesNames(Vector speciesNames_)
    {
        speciesNames = speciesNames_;
        Select(1);
    }

    public void SetValidVariable(boolean validVariable_)
    {
        validVariable = validVariable_;
        Select(2);
    }

    // Property getting methods
    public int     GetNumSpecies() { return numSpecies; }
    public Vector  GetSpeciesNames() { return speciesNames; }
    public boolean GetValidVariable() { return validVariable; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(numSpecies);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(speciesNames);
        if(WriteSelect(2, buf))
            buf.WriteBool(validVariable);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetNumSpecies(buf.ReadInt());
            break;
        case 1:
            SetSpeciesNames(buf.ReadStringVector());
            break;
        case 2:
            SetValidVariable(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("numSpecies", numSpecies, indent) + "\n";
        str = str + stringVectorToString("speciesNames", speciesNames, indent) + "\n";
        str = str + boolToString("validVariable", validVariable, indent) + "\n";
        return str;
    }


    // Attributes
    private int     numSpecies;
    private Vector  speciesNames; // vector of String objects
    private boolean validVariable;
}

