/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#include <MessageAttributes.h>
#include <DataNode.h>

//
// Enum conversion methods for MessageAttributes::Severity
//

static const char *Severity_strings[] = {
"Error", "Warning", "Message", 
"ErrorClear", "Information"};

std::string
MessageAttributes::Severity_ToString(MessageAttributes::Severity t)
{
    int index = int(t);
    if(index < 0 || index >= 5) index = 0;
    return Severity_strings[index];
}

std::string
MessageAttributes::Severity_ToString(int t)
{
    int index = (t < 0 || t >= 5) ? 0 : t;
    return Severity_strings[index];
}

bool
MessageAttributes::Severity_FromString(const std::string &s, MessageAttributes::Severity &val)
{
    val = MessageAttributes::Error;
    for(int i = 0; i < 5; ++i)
    {
        if(s == Severity_strings[i])
        {
            val = (Severity)i;
            return true;
        }
    }
    return false;
}

// ****************************************************************************
// Method: MessageAttributes::MessageAttributes
//
// Purpose: 
//   Init utility for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void MessageAttributes::Init()
{
    hasUnicode = false;
    severity = Message;

    MessageAttributes::SelectAll();
}

// ****************************************************************************
// Method: MessageAttributes::MessageAttributes
//
// Purpose: 
//   Copy utility for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void MessageAttributes::Copy(const MessageAttributes &obj)
{
    text = obj.text;
    unicode = obj.unicode;
    hasUnicode = obj.hasUnicode;
    severity = obj.severity;

    MessageAttributes::SelectAll();
}

// Type map format string
const char *MessageAttributes::TypeMapFormatString = MESSAGEATTRIBUTES_TMFS;
const AttributeGroup::private_tmfs_t MessageAttributes::TmfsStruct = {MESSAGEATTRIBUTES_TMFS};


// ****************************************************************************
// Method: MessageAttributes::MessageAttributes
//
// Purpose: 
//   Default constructor for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

MessageAttributes::MessageAttributes() : 
    AttributeSubject(MessageAttributes::TypeMapFormatString)
{
    MessageAttributes::Init();
}

// ****************************************************************************
// Method: MessageAttributes::MessageAttributes
//
// Purpose: 
//   Constructor for the derived classes of MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

MessageAttributes::MessageAttributes(private_tmfs_t tmfs) : 
    AttributeSubject(tmfs.tmfs)
{
    MessageAttributes::Init();
}

// ****************************************************************************
// Method: MessageAttributes::MessageAttributes
//
// Purpose: 
//   Copy constructor for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

MessageAttributes::MessageAttributes(const MessageAttributes &obj) : 
    AttributeSubject(MessageAttributes::TypeMapFormatString)
{
    MessageAttributes::Copy(obj);
}

// ****************************************************************************
// Method: MessageAttributes::MessageAttributes
//
// Purpose: 
//   Copy constructor for derived classes of the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

MessageAttributes::MessageAttributes(const MessageAttributes &obj, private_tmfs_t tmfs) : 
    AttributeSubject(tmfs.tmfs)
{
    MessageAttributes::Copy(obj);
}

// ****************************************************************************
// Method: MessageAttributes::~MessageAttributes
//
// Purpose: 
//   Destructor for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

MessageAttributes::~MessageAttributes()
{
    // nothing here
}

// ****************************************************************************
// Method: MessageAttributes::operator = 
//
// Purpose: 
//   Assignment operator for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

MessageAttributes& 
MessageAttributes::operator = (const MessageAttributes &obj)
{
    if (this == &obj) return *this;

    MessageAttributes::Copy(obj);

    return *this;
}

// ****************************************************************************
// Method: MessageAttributes::operator == 
//
// Purpose: 
//   Comparison operator == for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
MessageAttributes::operator == (const MessageAttributes &obj) const
{
    // Create the return value
    return ((text == obj.text) &&
            (unicode == obj.unicode) &&
            (hasUnicode == obj.hasUnicode) &&
            (severity == obj.severity));
}

// ****************************************************************************
// Method: MessageAttributes::operator != 
//
// Purpose: 
//   Comparison operator != for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
MessageAttributes::operator != (const MessageAttributes &obj) const
{
    return !(this->operator == (obj));
}

// ****************************************************************************
// Method: MessageAttributes::TypeName
//
// Purpose: 
//   Type name method for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

const std::string
MessageAttributes::TypeName() const
{
    return "MessageAttributes";
}

// ****************************************************************************
// Method: MessageAttributes::CopyAttributes
//
// Purpose: 
//   CopyAttributes method for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
MessageAttributes::CopyAttributes(const AttributeGroup *atts)
{
    if(TypeName() != atts->TypeName())
        return false;

    // Call assignment operator.
    const MessageAttributes *tmp = (const MessageAttributes *)atts;
    *this = *tmp;

    return true;
}

// ****************************************************************************
// Method: MessageAttributes::CreateCompatible
//
// Purpose: 
//   CreateCompatible method for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

AttributeSubject *
MessageAttributes::CreateCompatible(const std::string &tname) const
{
    AttributeSubject *retval = 0;
    if(TypeName() == tname)
        retval = new MessageAttributes(*this);
    // Other cases could go here too. 

    return retval;
}

// ****************************************************************************
// Method: MessageAttributes::NewInstance
//
// Purpose: 
//   NewInstance method for the MessageAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

AttributeSubject *
MessageAttributes::NewInstance(bool copy) const
{
    AttributeSubject *retval = 0;
    if(copy)
        retval = new MessageAttributes(*this);
    else
        retval = new MessageAttributes;

    return retval;
}

// ****************************************************************************
// Method: MessageAttributes::SelectAll
//
// Purpose: 
//   Selects all attributes.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void
MessageAttributes::SelectAll()
{
    Select(ID_text,       (void *)&text);
    Select(ID_unicode,    (void *)&unicode);
    Select(ID_hasUnicode, (void *)&hasUnicode);
    Select(ID_severity,   (void *)&severity);
}

///////////////////////////////////////////////////////////////////////////////
// Persistence methods
///////////////////////////////////////////////////////////////////////////////

// ****************************************************************************
// Method: MessageAttributes::CreateNode
//
// Purpose: 
//   This method creates a DataNode representation of the object so it can be saved to a config file.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
MessageAttributes::CreateNode(DataNode *parentNode, bool completeSave, bool forceAdd)
{
    if(parentNode == 0)
        return false;

    MessageAttributes defaultObject;
    bool addToParent = false;
    // Create a node for MessageAttributes.
    DataNode *node = new DataNode("MessageAttributes");

    if(completeSave || !FieldsEqual(ID_text, &defaultObject))
    {
        addToParent = true;
        node->AddNode(new DataNode("text", text));
    }

    if(completeSave || !FieldsEqual(ID_unicode, &defaultObject))
    {
        addToParent = true;
        node->AddNode(new DataNode("unicode", unicode));
    }

    if(completeSave || !FieldsEqual(ID_hasUnicode, &defaultObject))
    {
        addToParent = true;
        node->AddNode(new DataNode("hasUnicode", hasUnicode));
    }

    if(completeSave || !FieldsEqual(ID_severity, &defaultObject))
    {
        addToParent = true;
        node->AddNode(new DataNode("severity", Severity_ToString(severity)));
    }


    // Add the node to the parent node.
    if(addToParent || forceAdd)
        parentNode->AddNode(node);
    else
        delete node;

    return (addToParent || forceAdd);
}

// ****************************************************************************
// Method: MessageAttributes::SetFromNode
//
// Purpose: 
//   This method sets attributes in this object from values in a DataNode representation of the object.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void
MessageAttributes::SetFromNode(DataNode *parentNode)
{
    if(parentNode == 0)
        return;

    DataNode *searchNode = parentNode->GetNode("MessageAttributes");
    if(searchNode == 0)
        return;

    DataNode *node;
    if((node = searchNode->GetNode("text")) != 0)
        SetText(node->AsString());
    if((node = searchNode->GetNode("unicode")) != 0)
        SetUnicode(node->AsUnsignedCharVector());
    if((node = searchNode->GetNode("hasUnicode")) != 0)
        SetHasUnicode(node->AsBool());
    if((node = searchNode->GetNode("severity")) != 0)
    {
        // Allow enums to be int or string in the config file
        if(node->GetNodeType() == INT_NODE)
        {
            int ival = node->AsInt();
            if(ival >= 0 && ival < 5)
                SetSeverity(Severity(ival));
        }
        else if(node->GetNodeType() == STRING_NODE)
        {
            Severity value;
            if(Severity_FromString(node->AsString(), value))
                SetSeverity(value);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Set property methods
///////////////////////////////////////////////////////////////////////////////

void
MessageAttributes::SetText(const std::string &text_)
{
    text = text_;
    Select(ID_text, (void *)&text);
}

void
MessageAttributes::SetUnicode(const unsignedCharVector &unicode_)
{
    unicode = unicode_;
    Select(ID_unicode, (void *)&unicode);
}

void
MessageAttributes::SetHasUnicode(bool hasUnicode_)
{
    hasUnicode = hasUnicode_;
    Select(ID_hasUnicode, (void *)&hasUnicode);
}

void
MessageAttributes::SetSeverity(MessageAttributes::Severity severity_)
{
    severity = severity_;
    Select(ID_severity, (void *)&severity);
}

///////////////////////////////////////////////////////////////////////////////
// Get property methods
///////////////////////////////////////////////////////////////////////////////

const std::string &
MessageAttributes::GetText() const
{
    return text;
}

std::string &
MessageAttributes::GetText()
{
    return text;
}

const unsignedCharVector &
MessageAttributes::GetUnicode() const
{
    return unicode;
}

unsignedCharVector &
MessageAttributes::GetUnicode()
{
    return unicode;
}

bool
MessageAttributes::GetHasUnicode() const
{
    return hasUnicode;
}

MessageAttributes::Severity
MessageAttributes::GetSeverity() const
{
    return Severity(severity);
}

///////////////////////////////////////////////////////////////////////////////
// Select property methods
///////////////////////////////////////////////////////////////////////////////

void
MessageAttributes::SelectText()
{
    Select(ID_text, (void *)&text);
}

void
MessageAttributes::SelectUnicode()
{
    Select(ID_unicode, (void *)&unicode);
}

///////////////////////////////////////////////////////////////////////////////
// Keyframing methods
///////////////////////////////////////////////////////////////////////////////

// ****************************************************************************
// Method: MessageAttributes::GetFieldName
//
// Purpose: 
//   This method returns the name of a field given its index.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

std::string
MessageAttributes::GetFieldName(int index) const
{
    switch (index)
    {
    case ID_text:       return "text";
    case ID_unicode:    return "unicode";
    case ID_hasUnicode: return "hasUnicode";
    case ID_severity:   return "severity";
    default:  return "invalid index";
    }
}

// ****************************************************************************
// Method: MessageAttributes::GetFieldType
//
// Purpose: 
//   This method returns the type of a field given its index.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

AttributeGroup::FieldType
MessageAttributes::GetFieldType(int index) const
{
    switch (index)
    {
    case ID_text:       return FieldType_string;
    case ID_unicode:    return FieldType_ucharVector;
    case ID_hasUnicode: return FieldType_bool;
    case ID_severity:   return FieldType_enum;
    default:  return FieldType_unknown;
    }
}

// ****************************************************************************
// Method: MessageAttributes::GetFieldTypeName
//
// Purpose: 
//   This method returns the name of a field type given its index.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

std::string
MessageAttributes::GetFieldTypeName(int index) const
{
    switch (index)
    {
    case ID_text:       return "string";
    case ID_unicode:    return "ucharVector";
    case ID_hasUnicode: return "bool";
    case ID_severity:   return "enum";
    default:  return "invalid index";
    }
}

// ****************************************************************************
// Method: MessageAttributes::FieldsEqual
//
// Purpose: 
//   This method compares two fields and return true if they are equal.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
MessageAttributes::FieldsEqual(int index_, const AttributeGroup *rhs) const
{
    const MessageAttributes &obj = *((const MessageAttributes*)rhs);
    bool retval = false;
    switch (index_)
    {
    case ID_text:
        {  // new scope
        retval = (text == obj.text);
        }
        break;
    case ID_unicode:
        {  // new scope
        retval = (unicode == obj.unicode);
        }
        break;
    case ID_hasUnicode:
        {  // new scope
        retval = (hasUnicode == obj.hasUnicode);
        }
        break;
    case ID_severity:
        {  // new scope
        retval = (severity == obj.severity);
        }
        break;
    default: retval = false;
    }

    return retval;
}

///////////////////////////////////////////////////////////////////////////////
// User-defined methods.
///////////////////////////////////////////////////////////////////////////////

