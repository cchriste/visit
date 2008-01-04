// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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
// Class: avtSpeciesMetaData
//
// Purpose:
//    Contains species metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Apr 12 14:18:18 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class avtSpeciesMetaData extends AttributeSubject
{
    public avtSpeciesMetaData()
    {
        super(7);

        name = new String("Species");
        originalName = new String("Species");
        validVariable = true;
        meshName = new String("mesh");
        materialName = new String("material");
        numMaterials = 0;
        species = new Vector();
    }

    public avtSpeciesMetaData(avtSpeciesMetaData obj)
    {
        super(7);

        int i;

        name = new String(obj.name);
        originalName = new String(obj.originalName);
        validVariable = obj.validVariable;
        meshName = new String(obj.meshName);
        materialName = new String(obj.materialName);
        numMaterials = obj.numMaterials;
        // *** Copy the species field ***
        species = new Vector(obj.species.size());
        for(i = 0; i < obj.species.size(); ++i)
        {
            avtMatSpeciesMetaData newObj = (avtMatSpeciesMetaData)species.elementAt(i);
            species.addElement(new avtMatSpeciesMetaData(newObj));
        }


        SelectAll();
    }

    public boolean equals(avtSpeciesMetaData obj)
    {
        int i;

        boolean species_equal = (obj.species.size() == species.size());
        for(i = 0; (i < species.size()) && species_equal; ++i)
        {
            // Make references to avtMatSpeciesMetaData from Object.
            avtMatSpeciesMetaData species1 = (avtMatSpeciesMetaData)species.elementAt(i);
            avtMatSpeciesMetaData species2 = (avtMatSpeciesMetaData)obj.species.elementAt(i);
            species_equal = species1.equals(species2);
        }

        // Create the return value
        return ((name == obj.name) &&
                (originalName == obj.originalName) &&
                (validVariable == obj.validVariable) &&
                (meshName == obj.meshName) &&
                (materialName == obj.materialName) &&
                (numMaterials == obj.numMaterials) &&
                species_equal);
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetOriginalName(String originalName_)
    {
        originalName = originalName_;
        Select(1);
    }

    public void SetValidVariable(boolean validVariable_)
    {
        validVariable = validVariable_;
        Select(2);
    }

    public void SetMeshName(String meshName_)
    {
        meshName = meshName_;
        Select(3);
    }

    public void SetMaterialName(String materialName_)
    {
        materialName = materialName_;
        Select(4);
    }

    public void SetNumMaterials(int numMaterials_)
    {
        numMaterials = numMaterials_;
        Select(5);
    }

    // Property getting methods
    public String  GetName() { return name; }
    public String  GetOriginalName() { return originalName; }
    public boolean GetValidVariable() { return validVariable; }
    public String  GetMeshName() { return meshName; }
    public String  GetMaterialName() { return materialName; }
    public int     GetNumMaterials() { return numMaterials; }
    public Vector  GetSpecies() { return species; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteString(originalName);
        if(WriteSelect(2, buf))
            buf.WriteBool(validVariable);
        if(WriteSelect(3, buf))
            buf.WriteString(meshName);
        if(WriteSelect(4, buf))
            buf.WriteString(materialName);
        if(WriteSelect(5, buf))
            buf.WriteInt(numMaterials);
        if(WriteSelect(6, buf))
        {
            buf.WriteInt(species.size());
            for(int i = 0; i < species.size(); ++i)
            {
                avtMatSpeciesMetaData tmp = (avtMatSpeciesMetaData)species.elementAt(i);
                tmp.Write(buf);
            }
        }
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetName(buf.ReadString());
                break;
            case 1:
                SetOriginalName(buf.ReadString());
                break;
            case 2:
                SetValidVariable(buf.ReadBool());
                break;
            case 3:
                SetMeshName(buf.ReadString());
                break;
            case 4:
                SetMaterialName(buf.ReadString());
                break;
            case 5:
                SetNumMaterials(buf.ReadInt());
                break;
            case 6:
                {
                    int len = buf.ReadInt();
                    species.clear();
                    for(int j = 0; j < len; ++j)
                    {
                        avtMatSpeciesMetaData tmp = new avtMatSpeciesMetaData();
                        tmp.Read(buf);
                        species.addElement(tmp);
                    }
                }
                Select(6);
                break;
            }
        }
    }

    // Attributegroup convenience methods
    public void AddSpecies(avtMatSpeciesMetaData obj)
    {
        species.addElement(new avtMatSpeciesMetaData(obj));
        Select(6);
    }

    public void ClearSpeciess()
    {
        species.clear();
        Select(6);
    }

    public void RemoveSpecies(int index)
    {
        if(index >= 0 && index < species.size())
        {
            species.remove(index);
            Select(6);
        }
    }

    public int GetNumSpeciess()
    {
        return species.size();
    }

    public avtMatSpeciesMetaData GetSpecies(int i)
    {
        avtMatSpeciesMetaData tmp = (avtMatSpeciesMetaData)species.elementAt(i);
        return tmp;
    }


    // Attributes
    private String  name;
    private String  originalName;
    private boolean validVariable;
    private String  meshName;
    private String  materialName;
    private int     numMaterials;
    private Vector  species; // vector of avtMatSpeciesMetaData objects
}

