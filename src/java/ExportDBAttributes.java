// ***************************************************************************
//
// Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

import java.util.Vector;

// ****************************************************************************
// Class: ExportDBAttributes
//
// Purpose:
//    The attributes for export a database
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ExportDBAttributes extends AttributeSubject
{
    private static int numAdditionalAttributes = 6;

    public ExportDBAttributes()
    {
        super(numAdditionalAttributes);

        db_type = new String("");
        db_type_fullname = new String("");
        filename = new String("visit_ex_db");
        dirname = new String(".");
        variables = new Vector();
        opts = new DBOptionsAttributes();
    }

    public ExportDBAttributes(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

        db_type = new String("");
        db_type_fullname = new String("");
        filename = new String("visit_ex_db");
        dirname = new String(".");
        variables = new Vector();
        opts = new DBOptionsAttributes();
    }

    public ExportDBAttributes(ExportDBAttributes obj)
    {
        super(numAdditionalAttributes);

        int i;

        db_type = new String(obj.db_type);
        db_type_fullname = new String(obj.db_type_fullname);
        filename = new String(obj.filename);
        dirname = new String(obj.dirname);
        variables = new Vector(obj.variables.size());
        for(i = 0; i < obj.variables.size(); ++i)
            variables.addElement(new String((String)obj.variables.elementAt(i)));

        opts = new DBOptionsAttributes(obj.opts);

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

    public boolean equals(ExportDBAttributes obj)
    {
        int i;

        // Compare the elements in the variables vector.
        boolean variables_equal = (obj.variables.size() == variables.size());
        for(i = 0; (i < variables.size()) && variables_equal; ++i)
        {
            // Make references to String from Object.
            String variables1 = (String)variables.elementAt(i);
            String variables2 = (String)obj.variables.elementAt(i);
            variables_equal = variables1.equals(variables2);
        }
        // Create the return value
        return ((db_type.equals(obj.db_type)) &&
                (db_type_fullname.equals(obj.db_type_fullname)) &&
                (filename.equals(obj.filename)) &&
                (dirname.equals(obj.dirname)) &&
                variables_equal &&
                (opts.equals(obj.opts)));
    }

    // Property setting methods
    public void SetDb_type(String db_type_)
    {
        db_type = db_type_;
        Select(0);
    }

    public void SetDb_type_fullname(String db_type_fullname_)
    {
        db_type_fullname = db_type_fullname_;
        Select(1);
    }

    public void SetFilename(String filename_)
    {
        filename = filename_;
        Select(2);
    }

    public void SetDirname(String dirname_)
    {
        dirname = dirname_;
        Select(3);
    }

    public void SetVariables(Vector variables_)
    {
        variables = variables_;
        Select(4);
    }

    public void SetOpts(DBOptionsAttributes opts_)
    {
        opts = opts_;
        Select(5);
    }

    // Property getting methods
    public String              GetDb_type() { return db_type; }
    public String              GetDb_type_fullname() { return db_type_fullname; }
    public String              GetFilename() { return filename; }
    public String              GetDirname() { return dirname; }
    public Vector              GetVariables() { return variables; }
    public DBOptionsAttributes GetOpts() { return opts; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(db_type);
        if(WriteSelect(1, buf))
            buf.WriteString(db_type_fullname);
        if(WriteSelect(2, buf))
            buf.WriteString(filename);
        if(WriteSelect(3, buf))
            buf.WriteString(dirname);
        if(WriteSelect(4, buf))
            buf.WriteStringVector(variables);
        if(WriteSelect(5, buf))
            opts.Write(buf);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetDb_type(buf.ReadString());
            break;
        case 1:
            SetDb_type_fullname(buf.ReadString());
            break;
        case 2:
            SetFilename(buf.ReadString());
            break;
        case 3:
            SetDirname(buf.ReadString());
            break;
        case 4:
            SetVariables(buf.ReadStringVector());
            break;
        case 5:
            opts.Read(buf);
            Select(5);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("db_type", db_type, indent) + "\n";
        str = str + stringToString("db_type_fullname", db_type_fullname, indent) + "\n";
        str = str + stringToString("filename", filename, indent) + "\n";
        str = str + stringToString("dirname", dirname, indent) + "\n";
        str = str + stringVectorToString("variables", variables, indent) + "\n";
        str = str + indent + "opts = {\n" + opts.toString(indent + "    ") + indent + "}\n";
        return str;
    }


    // Attributes
    private String              db_type;
    private String              db_type_fullname;
    private String              filename;
    private String              dirname;
    private Vector              variables; // vector of String objects
    private DBOptionsAttributes opts;
}

