/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef DATANODE_H
#define DATANODE_H
#include <state_exports.h>
#include <AttributeGroup.h>
#include <visitstream.h>

typedef enum
{
    INTERNAL_NODE = 0,
    CHAR_NODE, UNSIGNED_CHAR_NODE, INT_NODE, LONG_NODE, FLOAT_NODE,
    DOUBLE_NODE, STRING_NODE, BOOL_NODE,
    CHAR_ARRAY_NODE, UNSIGNED_CHAR_ARRAY_NODE, INT_ARRAY_NODE,
    LONG_ARRAY_NODE, FLOAT_ARRAY_NODE,
    DOUBLE_ARRAY_NODE, STRING_ARRAY_NODE, BOOL_ARRAY_NODE,
    CHAR_VECTOR_NODE, UNSIGNED_CHAR_VECTOR_NODE, INT_VECTOR_NODE,
    LONG_VECTOR_NODE, FLOAT_VECTOR_NODE,
    DOUBLE_VECTOR_NODE, STRING_VECTOR_NODE, BOOL_VECTOR_NODE,
    MAP_NODE_NODE
} NodeTypeEnum;

// ****************************************************************************
// Class: DataNode
//
// Purpose:
//   This class contains a tree structure that is the basis of the
//   configuration file when it is in memory.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Sep 28 10:14:49 PDT 2000
//
// Modifications:
//   Jeremy Meredith, Mon Feb 26 16:01:09 PST 2001
//   Added unsigned chars.
//
//   Brad Whitlock, Fri Mar 21 10:43:22 PDT 2003
//   I added methods to set the DataNode's data so it can be
//   transformed to a new type.
//
//   Brad Whitlock, Mon Feb 2 15:36:49 PST 2004
//   Added an optional argument to RemoveNode.
//
//   Eric Brugger, Tue Mar 27 15:57:03 PDT 2007
//   Added an additional RemoveNode method that takes a DataNode as an
//   argument.
//
//   Brad Whitlock, Thu Feb 14 14:49:19 PST 2008
//   Added SearchForNode and make GetNode not recurse.
//
//   Brad Whitlock, Tue Apr 22 12:14:38 PDT 2008
//   Added a Print method for debugging.
//
//    Kathleen Biagas, Fri Jun 17 16:41:27 PDT 2011
//    Add MapNode node.
//
// ****************************************************************************

class STATE_API DataNode
{
public:
    // Internal node constructor
    DataNode(const std::string &name);

    // Single element constructors
    DataNode(const std::string &name, char val);
    DataNode(const std::string &name, unsigned char val);
    DataNode(const std::string &name, int val);
    DataNode(const std::string &name, long val);
    DataNode(const std::string &name, float val);
    DataNode(const std::string &name, double val);
    DataNode(const std::string &name, const std::string &val);
    DataNode(const std::string &name, bool val);

    // Array constructors
    DataNode(const std::string &name, const char *vals, int len);
    DataNode(const std::string &name, const unsigned char *vals, int len);
    DataNode(const std::string &name, const int *vals, int len);
    DataNode(const std::string &name, const long *vals, int len);
    DataNode(const std::string &name, const float *vals, int len);
    DataNode(const std::string &name, const double *val, int len);
    DataNode(const std::string &name, const std::string *vals, int len);
    DataNode(const std::string &name, const bool *vals, int len);

    // Vector constructors
    DataNode(const std::string &name, const charVector &vec);
    DataNode(const std::string &name, const unsignedCharVector &vec);
    DataNode(const std::string &name, const intVector &vec);
    DataNode(const std::string &name, const longVector &vec);
    DataNode(const std::string &name, const floatVector &vec);
    DataNode(const std::string &name, const doubleVector &vec);
    DataNode(const std::string &name, const stringVector &vec);

    DataNode(const std::string &name, const MapNode &val);

    ~DataNode();

    // Functions to return the data as a particular type.
    char               AsChar() const;
    unsigned char      AsUnsignedChar() const;
    int                AsInt()  const;
    long               AsLong() const;
    float              AsFloat() const;
    double             AsDouble() const;
    const std::string &AsString() const;
    bool               AsBool() const;

    const char          *AsCharArray() const;
    const unsigned char *AsUnsignedCharArray() const;
    const int           *AsIntArray()  const;
    const long          *AsLongArray() const;
    const float         *AsFloatArray() const;
    const double        *AsDoubleArray() const;
    const std::string   *AsStringArray() const;
    const bool          *AsBoolArray() const;

    const charVector          &AsCharVector() const;
    const unsignedCharVector  &AsUnsignedCharVector() const;
    const intVector           &AsIntVector()  const;
    const longVector          &AsLongVector() const;
    const floatVector         &AsFloatVector() const;
    const doubleVector        &AsDoubleVector() const;
    const stringVector        &AsStringVector() const;

    const MapNode             &AsMapNode() const;

    void SetChar(char val);
    void SetUnsignedChar(unsigned char val);
    void SetInt(int val);
    void SetLong(long val);
    void SetFloat(float val);
    void SetDouble(double val);
    void SetString(const std::string &val);
    void SetBool(bool val);
    void SetCharArray(const char *vals, int len);
    void SetUnsignedCharArray(const unsigned char *vals, int len);
    void SetIntArray(const int *vals, int len);
    void SetLongArray(const long *vals, int len);
    void SetFloatArray(const float *vals, int len);
    void SetDoubleArray(const double *vals, int len);
    void SetStringArray(const std::string *vals, int len);
    void SetBoolArray(const bool *vals, int len);
    void SetCharVector(const charVector &vec);
    void SetUnsignedCharVector(const unsignedCharVector &vec);
    void SetIntVector(const intVector &vec);
    void SetLongVector(const longVector &vec);
    void SetFloatVector(const floatVector &vec);
    void SetDoubleVector(const doubleVector &vec);
    void SetStringVector(const stringVector &vec);

    void SetMapNode(const MapNode &val);

    // Node operations
    DataNode *GetNode(const std::string &key);
    DataNode *SearchForNode(const std::string &key, DataNode *parentNode = 0);
    void AddNode(DataNode *node);
    void RemoveNode(DataNode *node, bool deleteNode = true);
    void RemoveNode(const std::string &key, bool deleteNode = true);
    void Print(ostream &os, int indent = 0);

    // Functions to return private members.
    const std::string &GetKey() const;
    void SetKey(const std::string &);
    NodeTypeEnum GetNodeType() const;
    int GetLength() const;
    int GetNumChildren() const;
    int GetNumChildObjects() const;
    DataNode **GetChildren();
private:
    void FreeData();

    std::string   Key;
    NodeTypeEnum  NodeType;
    int           Length;
    void          *Data;

    // Static members
    static std::string         bogusString;
    static charVector          bogusCharVector;
    static unsignedCharVector  bogusUnsignedCharVector;
    static intVector           bogusIntVector;
    static longVector          bogusLongVector;
    static floatVector         bogusFloatVector;
    static doubleVector        bogusDoubleVector;
    static stringVector        bogusStringVector;

    static MapNode             bogusMapNode;
};

// Utility functions.
const char *NodeTypeName(NodeTypeEnum e);
NodeTypeEnum GetNodeType(const char *str);

#endif
