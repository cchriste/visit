/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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

#ifndef INDEXSELECTATTRIBUTES_H
#define INDEXSELECTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: IndexSelectAttributes
//
// Purpose:
//    This class contains attributes for the index select operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class IndexSelectAttributes : public AttributeSubject
{
public:
    enum Dimension
    {
        OneD,
        TwoD,
        ThreeD
    };

    IndexSelectAttributes();
    IndexSelectAttributes(const IndexSelectAttributes &obj);
    virtual ~IndexSelectAttributes();

    virtual IndexSelectAttributes& operator = (const IndexSelectAttributes &obj);
    virtual bool operator == (const IndexSelectAttributes &obj) const;
    virtual bool operator != (const IndexSelectAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectCategoryName();
    void SelectSubsetName();

    // Property setting methods
    void SetDim(Dimension dim_);
    void SetXMin(int xMin_);
    void SetXMax(int xMax_);
    void SetXIncr(int xIncr_);
    void SetYMin(int yMin_);
    void SetYMax(int yMax_);
    void SetYIncr(int yIncr_);
    void SetZMin(int zMin_);
    void SetZMax(int zMax_);
    void SetZIncr(int zIncr_);
    void SetUseWholeCollection(bool useWholeCollection_);
    void SetCategoryName(const std::string &categoryName_);
    void SetSubsetName(const std::string &subsetName_);

    // Property getting methods
    Dimension         GetDim() const;
    int               GetXMin() const;
    int               GetXMax() const;
    int               GetXIncr() const;
    int               GetYMin() const;
    int               GetYMax() const;
    int               GetYIncr() const;
    int               GetZMin() const;
    int               GetZMax() const;
    int               GetZIncr() const;
    bool              GetUseWholeCollection() const;
    const std::string &GetCategoryName() const;
          std::string &GetCategoryName();
    const std::string &GetSubsetName() const;
          std::string &GetSubsetName();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Dimension_ToString(Dimension);
    static bool Dimension_FromString(const std::string &, Dimension &);
protected:
    static std::string Dimension_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual void ProcessOldVersions(DataNode *node, const char *configVersion);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_dim = 0,
        ID_xMin,
        ID_xMax,
        ID_xIncr,
        ID_yMin,
        ID_yMax,
        ID_yIncr,
        ID_zMin,
        ID_zMax,
        ID_zIncr,
        ID_useWholeCollection,
        ID_categoryName,
        ID_subsetName
    };

private:
    int         dim;
    int         xMin;
    int         xMax;
    int         xIncr;
    int         yMin;
    int         yMax;
    int         yIncr;
    int         zMin;
    int         zMax;
    int         zIncr;
    bool        useWholeCollection;
    std::string categoryName;
    std::string subsetName;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
