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

#ifndef MAP_NODE_H
#define MAP_NODE_H

#include <state_exports.h>
#include <XMLNode.h>
#include <JSONNode.h>
#include <Variant.h>
#include <map>

// ****************************************************************************
//  Class:  MapNode
//
//  Purpose:
//    Provides a nested map type. 
//
//  Programmer:  Cyrus Harrison
//  Creation:    December 14, 2007
//
//  Modifications:
//    Brad Whitlock, Tue Jan  6 15:32:50 PST 2009
//    I added methods so the MapNode can read/write itself using Connection.
//
//    Brad Whitlock, Fri Jan 16 11:36:10 PST 2009
//    I added a Merge function.
//
// ****************************************************************************

class STATE_API MapNode : public Variant
{
  public:
    MapNode();
    MapNode(const MapNode&);
    MapNode(const XMLNode&,bool decodeString = true);
    MapNode(const XMLNode*,bool decodeString = true);
    explicit MapNode(const JSONNode&, bool decodeString = true);
    explicit MapNode(const JSONNode*, bool decodeString = true);
    explicit MapNode(const JSONNode&, const JSONNode& metadata, bool decodeString = true);
    explicit MapNode(const JSONNode*,const JSONNode *metadata, bool decodeString = true);
    MapNode  &operator=(const MapNode&);
    MapNode  &operator=(bool);
    MapNode  &operator=(char);
    MapNode  &operator=(unsigned char);
    MapNode  &operator=(const char *); // interp as string
    MapNode  &operator=(int);
    MapNode  &operator=(long);
    MapNode  &operator=(float);
    MapNode  &operator=(double);
    MapNode  &operator=(const std::string &);
    MapNode  &operator=(const boolVector &);
    MapNode  &operator=(const charVector &);
    MapNode  &operator=(const unsignedCharVector &);
    MapNode  &operator=(const intVector &);
    MapNode  &operator=(const longVector &);
    MapNode  &operator=(const floatVector &);
    MapNode  &operator=(const doubleVector &);
    MapNode  &operator=(const stringVector &);
    MapNode  &operator=(const Variant &);
    virtual  ~MapNode();

    bool                 operator ==(const MapNode &obj) const;

    MapNode             &operator[](const std::string &);
    MapNode             *GetEntry(const std::string &);
    const MapNode       *GetEntry(const std::string &) const;

    void                 Merge(const MapNode &);
    
    void                 RemoveEntry(const std::string &);
    bool                 HasEntry(const std::string &) const;
    bool                 HasNumericEntry(const std::string &) const;
    bool                 HasNumericVectorEntry(const std::string &) const;
    void                 GetEntryNames(stringVector &) const;
    int                  GetNumEntries() const {return (int)entries.size();}
    void                 Reset();

    virtual std::string  ToXML(bool encodeString = true) const;
    virtual XMLNode      ToXMLNode(bool encodeString = true) const;

    virtual std::string  ToJSON(bool encodeString = true) const;
    virtual JSONNode     ToJSONNode(bool encodeString = true, bool id = true) const;

    int                  CalculateMessageSize(Connection &conn) const;
    int                  CalculateMessageSize(Connection *conn) const;
    void                 Write(Connection &conn) const;
    void                 Write(Connection *conn) const;
    void                 Read(Connection &conn);

    static int MapNodeType;
 private:
    virtual JSONNode ToJSONNodeData(bool encodeString) const;
    virtual JSONNode ToJSONNodeMetaData(bool id) const;
    void  SetValue(const XMLNode &, bool decodeString = true);
    void  SetValue(const JSONNode &, bool decodeString = true);
    void  SetValue(const JSONNode& data, const JSONNode& metadata,bool decodeString);
    std::map<std::string,MapNode> entries;
};

#endif

