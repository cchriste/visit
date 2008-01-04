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
import java.lang.Integer;

// ****************************************************************************
// Class: SILAttributes
//
// Purpose:
//    This class contains the information needed to represent a SIL.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:55:11 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class SILAttributes extends AttributeSubject
{
    public SILAttributes()
    {
        super(10);

        nSets = 0;
        setNames = new Vector();
        setIds = new Vector();
        isWhole = new Vector();
        nCollections = 0;
        category = new Vector();
        role = new Vector();
        superset = new Vector();
        nspace = new Vector();
        matrices = new Vector();
    }

    public SILAttributes(SILAttributes obj)
    {
        super(10);

        int i;

        nSets = obj.nSets;
        setNames = new Vector(obj.setNames.size());
        for(i = 0; i < obj.setNames.size(); ++i)
            setNames.addElement(new String((String)obj.setNames.elementAt(i)));

        setIds = new Vector();
        for(i = 0; i < obj.setIds.size(); ++i)
        {
            Integer iv = (Integer)obj.setIds.elementAt(i);
            setIds.addElement(new Integer(iv.intValue()));
        }
        isWhole = new Vector();
        for(i = 0; i < obj.isWhole.size(); ++i)
        {
            Integer iv = (Integer)obj.isWhole.elementAt(i);
            isWhole.addElement(new Integer(iv.intValue()));
        }
        nCollections = obj.nCollections;
        category = new Vector(obj.category.size());
        for(i = 0; i < obj.category.size(); ++i)
            category.addElement(new String((String)obj.category.elementAt(i)));

        role = new Vector();
        for(i = 0; i < obj.role.size(); ++i)
        {
            Integer iv = (Integer)obj.role.elementAt(i);
            role.addElement(new Integer(iv.intValue()));
        }
        superset = new Vector();
        for(i = 0; i < obj.superset.size(); ++i)
        {
            Integer iv = (Integer)obj.superset.elementAt(i);
            superset.addElement(new Integer(iv.intValue()));
        }
        // *** Copy the nspace field ***
        nspace = new Vector(obj.nspace.size());
        for(i = 0; i < obj.nspace.size(); ++i)
        {
            NamespaceAttributes newObj = (NamespaceAttributes)nspace.elementAt(i);
            nspace.addElement(new NamespaceAttributes(newObj));
        }

        // *** Copy the matrices field ***
        matrices = new Vector(obj.matrices.size());
        for(i = 0; i < obj.matrices.size(); ++i)
        {
            SILMatrixAttributes newObj = (SILMatrixAttributes)matrices.elementAt(i);
            matrices.addElement(new SILMatrixAttributes(newObj));
        }


        SelectAll();
    }

    public boolean equals(SILAttributes obj)
    {
        int i;

        boolean nspace_equal = (obj.nspace.size() == nspace.size());
        for(i = 0; (i < nspace.size()) && nspace_equal; ++i)
        {
            // Make references to NamespaceAttributes from Object.
            NamespaceAttributes nspace1 = (NamespaceAttributes)nspace.elementAt(i);
            NamespaceAttributes nspace2 = (NamespaceAttributes)obj.nspace.elementAt(i);
            nspace_equal = nspace1.equals(nspace2);
        }

        boolean matrices_equal = (obj.matrices.size() == matrices.size());
        for(i = 0; (i < matrices.size()) && matrices_equal; ++i)
        {
            // Make references to SILMatrixAttributes from Object.
            SILMatrixAttributes matrices1 = (SILMatrixAttributes)matrices.elementAt(i);
            SILMatrixAttributes matrices2 = (SILMatrixAttributes)obj.matrices.elementAt(i);
            matrices_equal = matrices1.equals(matrices2);
        }

        // Create the return value
        return ((nSets == obj.nSets) &&
                (setNames == obj.setNames) &&
                (setIds == obj.setIds) &&
                (isWhole == obj.isWhole) &&
                (nCollections == obj.nCollections) &&
                (category == obj.category) &&
                (role == obj.role) &&
                (superset == obj.superset) &&
                nspace_equal &&
                matrices_equal);
    }

    // Property setting methods
    public void SetNSets(int nSets_)
    {
        nSets = nSets_;
        Select(0);
    }

    public void SetSetNames(Vector setNames_)
    {
        setNames = setNames_;
        Select(1);
    }

    public void SetSetIds(Vector setIds_)
    {
        setIds = setIds_;
        Select(2);
    }

    public void SetIsWhole(Vector isWhole_)
    {
        isWhole = isWhole_;
        Select(3);
    }

    public void SetNCollections(int nCollections_)
    {
        nCollections = nCollections_;
        Select(4);
    }

    public void SetCategory(Vector category_)
    {
        category = category_;
        Select(5);
    }

    public void SetRole(Vector role_)
    {
        role = role_;
        Select(6);
    }

    public void SetSuperset(Vector superset_)
    {
        superset = superset_;
        Select(7);
    }

    // Property getting methods
    public int    GetNSets() { return nSets; }
    public Vector GetSetNames() { return setNames; }
    public Vector GetSetIds() { return setIds; }
    public Vector GetIsWhole() { return isWhole; }
    public int    GetNCollections() { return nCollections; }
    public Vector GetCategory() { return category; }
    public Vector GetRole() { return role; }
    public Vector GetSuperset() { return superset; }
    public Vector GetNspace() { return nspace; }
    public Vector GetMatrices() { return matrices; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(nSets);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(setNames);
        if(WriteSelect(2, buf))
            buf.WriteIntVector(setIds);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(isWhole);
        if(WriteSelect(4, buf))
            buf.WriteInt(nCollections);
        if(WriteSelect(5, buf))
            buf.WriteStringVector(category);
        if(WriteSelect(6, buf))
            buf.WriteIntVector(role);
        if(WriteSelect(7, buf))
            buf.WriteIntVector(superset);
        if(WriteSelect(8, buf))
        {
            buf.WriteInt(nspace.size());
            for(int i = 0; i < nspace.size(); ++i)
            {
                NamespaceAttributes tmp = (NamespaceAttributes)nspace.elementAt(i);
                tmp.Write(buf);
            }
        }
        if(WriteSelect(9, buf))
        {
            buf.WriteInt(matrices.size());
            for(int i = 0; i < matrices.size(); ++i)
            {
                SILMatrixAttributes tmp = (SILMatrixAttributes)matrices.elementAt(i);
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
                SetNSets(buf.ReadInt());
                break;
            case 1:
                SetSetNames(buf.ReadStringVector());
                break;
            case 2:
                SetSetIds(buf.ReadIntVector());
                break;
            case 3:
                SetIsWhole(buf.ReadIntVector());
                break;
            case 4:
                SetNCollections(buf.ReadInt());
                break;
            case 5:
                SetCategory(buf.ReadStringVector());
                break;
            case 6:
                SetRole(buf.ReadIntVector());
                break;
            case 7:
                SetSuperset(buf.ReadIntVector());
                break;
            case 8:
                {
                    int len = buf.ReadInt();
                    nspace.clear();
                    for(int j = 0; j < len; ++j)
                    {
                        NamespaceAttributes tmp = new NamespaceAttributes();
                        tmp.Read(buf);
                        nspace.addElement(tmp);
                    }
                }
                Select(8);
                break;
            case 9:
                {
                    int len = buf.ReadInt();
                    matrices.clear();
                    for(int j = 0; j < len; ++j)
                    {
                        SILMatrixAttributes tmp = new SILMatrixAttributes();
                        tmp.Read(buf);
                        matrices.addElement(tmp);
                    }
                }
                Select(9);
                break;
            }
        }
    }

    // Attributegroup convenience methods
    public void AddNspace(NamespaceAttributes obj)
    {
        nspace.addElement(new NamespaceAttributes(obj));
        Select(8);
    }

    public void ClearNspaces()
    {
        nspace.clear();
        Select(8);
    }

    public void RemoveNspace(int index)
    {
        if(index >= 0 && index < nspace.size())
        {
            nspace.remove(index);
            Select(8);
        }
    }

    public int GetNumNspaces()
    {
        return nspace.size();
    }

    public NamespaceAttributes GetNspace(int i)
    {
        NamespaceAttributes tmp = (NamespaceAttributes)nspace.elementAt(i);
        return tmp;
    }

    public void AddMatrices(SILMatrixAttributes obj)
    {
        matrices.addElement(new SILMatrixAttributes(obj));
        Select(9);
    }

    public void ClearMatrices()
    {
        matrices.clear();
        Select(9);
    }

    public void RemoveMatrices(int index)
    {
        if(index >= 0 && index < matrices.size())
        {
            matrices.remove(index);
            Select(9);
        }
    }

    public int GetNumMatrices()
    {
        return matrices.size();
    }

    public SILMatrixAttributes GetMatrices(int i)
    {
        SILMatrixAttributes tmp = (SILMatrixAttributes)matrices.elementAt(i);
        return tmp;
    }


    // Attributes
    private int    nSets;
    private Vector setNames; // vector of String objects
    private Vector setIds; // vector of Integer objects
    private Vector isWhole; // vector of Integer objects
    private int    nCollections;
    private Vector category; // vector of String objects
    private Vector role; // vector of Integer objects
    private Vector superset; // vector of Integer objects
    private Vector nspace; // vector of NamespaceAttributes objects
    private Vector matrices; // vector of SILMatrixAttributes objects
}

