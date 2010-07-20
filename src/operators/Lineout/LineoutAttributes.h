/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

#ifndef LINEOUTATTRIBUTES_H
#define LINEOUTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <visitstream.h>

// ****************************************************************************
// Class: LineoutAttributes
//
// Purpose:
//    Attributes for the Lineout operator 
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class LineoutAttributes : public AttributeSubject
{
public:
    LineoutAttributes();
    LineoutAttributes(const LineoutAttributes &obj);
    virtual ~LineoutAttributes();

    virtual LineoutAttributes& operator = (const LineoutAttributes &obj);
    virtual bool operator == (const LineoutAttributes &obj) const;
    virtual bool operator != (const LineoutAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectPoint1();
    void SelectPoint2();
    void SelectDesignator();

    // Property setting methods
    void SetPoint1(const double *point1_);
    void SetPoint2(const double *point2_);
    void SetInteractive(bool interactive_);
    void SetIgnoreGlobal(bool ignoreGlobal_);
    void SetSamplingOn(bool samplingOn_);
    void SetNumberOfSamplePoints(int numberOfSamplePoints_);
    void SetReflineLabels(bool reflineLabels_);
    void SetDesignator(const std::string &designator_);

    // Property getting methods
    const double      *GetPoint1() const;
          double      *GetPoint1();
    const double      *GetPoint2() const;
          double      *GetPoint2();
    bool              GetInteractive() const;
    bool              GetIgnoreGlobal() const;
    bool              GetSamplingOn() const;
    int               GetNumberOfSamplePoints() const;
    bool              GetReflineLabels() const;
    const std::string &GetDesignator() const;
          std::string &GetDesignator();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void Print(ostream &, bool) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_point1 = 0,
        ID_point2,
        ID_interactive,
        ID_ignoreGlobal,
        ID_samplingOn,
        ID_numberOfSamplePoints,
        ID_reflineLabels,
        ID_designator
    };

private:
    double      point1[3];
    double      point2[3];
    bool        interactive;
    bool        ignoreGlobal;
    bool        samplingOn;
    int         numberOfSamplePoints;
    bool        reflineLabels;
    std::string designator;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
