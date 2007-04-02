// ***************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
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
// Class: DBPluginInfoAttributes
//
// Purpose:
//    This class contains the attributes for all the database plugins.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:54:48 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class DBPluginInfoAttributes extends AttributeSubject
{
    public DBPluginInfoAttributes()
    {
        super(4);

        types = new Vector();
        hasWriter = new Vector();
        dbOptions = new Vector();
        typesFullNames = new Vector();
    }

    public DBPluginInfoAttributes(DBPluginInfoAttributes obj)
    {
        super(4);

        int i;

        types = new Vector(obj.types.size());
        for(i = 0; i < obj.types.size(); ++i)
            types.addElement(new String((String)obj.types.elementAt(i)));

        hasWriter = new Vector();
        for(i = 0; i < obj.hasWriter.size(); ++i)
        {
            Integer iv = (Integer)obj.hasWriter.elementAt(i);
            hasWriter.addElement(new Integer(iv.intValue()));
        }
        // *** Copy the dbOptions field ***
        dbOptions = new Vector(obj.dbOptions.size());
        for(i = 0; i < obj.dbOptions.size(); ++i)
        {
            DBOptionsAttributes newObj = (DBOptionsAttributes)dbOptions.elementAt(i);
            dbOptions.addElement(new DBOptionsAttributes(newObj));
        }

        typesFullNames = new Vector(obj.typesFullNames.size());
        for(i = 0; i < obj.typesFullNames.size(); ++i)
            typesFullNames.addElement(new String((String)obj.typesFullNames.elementAt(i)));


        SelectAll();
    }

    public boolean equals(DBPluginInfoAttributes obj)
    {
        int i;

        boolean dbOptions_equal = (obj.dbOptions.size() == dbOptions.size());
        for(i = 0; (i < dbOptions.size()) && dbOptions_equal; ++i)
        {
            // Make references to DBOptionsAttributes from Object.
            DBOptionsAttributes dbOptions1 = (DBOptionsAttributes)dbOptions.elementAt(i);
            DBOptionsAttributes dbOptions2 = (DBOptionsAttributes)obj.dbOptions.elementAt(i);
            dbOptions_equal = dbOptions1.equals(dbOptions2);
        }

        // Create the return value
        return ((types == obj.types) &&
                (hasWriter == obj.hasWriter) &&
                dbOptions_equal &&
                (typesFullNames == obj.typesFullNames));
    }

    // Property setting methods
    public void SetTypes(Vector types_)
    {
        types = types_;
        Select(0);
    }

    public void SetHasWriter(Vector hasWriter_)
    {
        hasWriter = hasWriter_;
        Select(1);
    }

    public void SetTypesFullNames(Vector typesFullNames_)
    {
        typesFullNames = typesFullNames_;
        Select(3);
    }

    // Property getting methods
    public Vector GetTypes() { return types; }
    public Vector GetHasWriter() { return hasWriter; }
    public Vector GetDbOptions() { return dbOptions; }
    public Vector GetTypesFullNames() { return typesFullNames; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(types);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(hasWriter);
        if(WriteSelect(2, buf))
        {
            buf.WriteInt(dbOptions.size());
            for(int i = 0; i < dbOptions.size(); ++i)
            {
                DBOptionsAttributes tmp = (DBOptionsAttributes)dbOptions.elementAt(i);
                tmp.Write(buf);
            }
        }
        if(WriteSelect(3, buf))
            buf.WriteStringVector(typesFullNames);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetTypes(buf.ReadStringVector());
                break;
            case 1:
                SetHasWriter(buf.ReadIntVector());
                break;
            case 2:
                {
                    int len = buf.ReadInt();
                    dbOptions.clear();
                    for(int j = 0; j < len; ++j)
                    {
                        DBOptionsAttributes tmp = new DBOptionsAttributes();
                        tmp.Read(buf);
                        dbOptions.addElement(tmp);
                    }
                }
                Select(2);
                break;
            case 3:
                SetTypesFullNames(buf.ReadStringVector());
                break;
            }
        }
    }

    // Attributegroup convenience methods
    public void AddDbOptions(DBOptionsAttributes obj)
    {
        dbOptions.addElement(new DBOptionsAttributes(obj));
        Select(2);
    }

    public void ClearDbOptions()
    {
        dbOptions.clear();
        Select(2);
    }

    public void RemoveDbOptions(int index)
    {
        if(index >= 0 && index < dbOptions.size())
        {
            dbOptions.remove(index);
            Select(2);
        }
    }

    public int GetNumDbOptions()
    {
        return dbOptions.size();
    }

    public DBOptionsAttributes GetDbOptions(int i)
    {
        DBOptionsAttributes tmp = (DBOptionsAttributes)dbOptions.elementAt(i);
        return tmp;
    }


    // Attributes
    private Vector types; // vector of String objects
    private Vector hasWriter; // vector of Integer objects
    private Vector dbOptions; // vector of DBOptionsAttributes objects
    private Vector typesFullNames; // vector of String objects
}

