/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

#ifndef SAVESUBWINDOWATTRIBUTES_H
#define SAVESUBWINDOWATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: SaveSubWindowAttributes
//
// Purpose:
//    The attributes for saving a sub window
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SaveSubWindowAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    SaveSubWindowAttributes();
    SaveSubWindowAttributes(const SaveSubWindowAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    SaveSubWindowAttributes(private_tmfs_t tmfs);
    SaveSubWindowAttributes(const SaveSubWindowAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~SaveSubWindowAttributes();

    virtual SaveSubWindowAttributes& operator = (const SaveSubWindowAttributes &obj);
    virtual bool operator == (const SaveSubWindowAttributes &obj) const;
    virtual bool operator != (const SaveSubWindowAttributes &obj) const;
private:
    void Init();
    void Copy(const SaveSubWindowAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectPosition();
    void SelectSize();

    // Property setting methods
    void SetPosition(const int *position_);
    void SetSize(const int *size_);
    void SetLayer(int layer_);
    void SetTransparency(double transparency_);
    void SetOmitWindow(bool omitWindow_);

    // Property getting methods
    const int *GetPosition() const;
          int *GetPosition();
    const int *GetSize() const;
          int *GetSize();
    int       GetLayer() const;
    double    GetTransparency() const;
    bool      GetOmitWindow() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_position = 0,
        ID_size,
        ID_layer,
        ID_transparency,
        ID_omitWindow,
        ID__LAST
    };

private:
    int    position[2];
    int    size[2];
    int    layer;
    double transparency;
    bool   omitWindow;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define SAVESUBWINDOWATTRIBUTES_TMFS "IIidb"

#endif
