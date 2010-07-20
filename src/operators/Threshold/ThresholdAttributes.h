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

#ifndef THRESHOLDATTRIBUTES_H
#define THRESHOLDATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <AxisRestrictionAttributes.h>
#include <DebugStream.h>

// ****************************************************************************
// Class: ThresholdAttributes
//
// Purpose:
//    This class contains attributes for the threshold operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class ThresholdAttributes : public AttributeSubject
{
public:
    enum OutputMeshType
    {
        InputZones,
        PointMesh
    };
    enum ZonePortion
    {
        EntireZone,
        PartOfZone
    };

    ThresholdAttributes();
    ThresholdAttributes(const ThresholdAttributes &obj);
    virtual ~ThresholdAttributes();

    virtual ThresholdAttributes& operator = (const ThresholdAttributes &obj);
    virtual bool operator == (const ThresholdAttributes &obj) const;
    virtual bool operator != (const ThresholdAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectListedVarNames();
    void SelectZonePortions();
    void SelectLowerBounds();
    void SelectUpperBounds();
    void SelectDefaultVarName();

    // Property setting methods
    void SetOutputMeshType(int outputMeshType_);
    void SetListedVarNames(const stringVector &listedVarNames_);
    void SetZonePortions(const intVector &zonePortions_);
    void SetLowerBounds(const doubleVector &lowerBounds_);
    void SetUpperBounds(const doubleVector &upperBounds_);
    void SetDefaultVarName(const std::string &defaultVarName_);
    void SetDefaultVarIsScalar(bool defaultVarIsScalar_);

    // Property getting methods
    int                GetOutputMeshType() const;
    const stringVector &GetListedVarNames() const;
          stringVector &GetListedVarNames();
    const intVector    &GetZonePortions() const;
          intVector    &GetZonePortions();
    const doubleVector &GetLowerBounds() const;
          doubleVector &GetLowerBounds();
    const doubleVector &GetUpperBounds() const;
          doubleVector &GetUpperBounds();
    const std::string  &GetDefaultVarName() const;
          std::string  &GetDefaultVarName();
    bool               GetDefaultVarIsScalar() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string OutputMeshType_ToString(OutputMeshType);
    static bool OutputMeshType_FromString(const std::string &, OutputMeshType &);
protected:
    static std::string OutputMeshType_ToString(int);
public:
    static std::string ZonePortion_ToString(ZonePortion);
    static bool ZonePortion_FromString(const std::string &, ZonePortion &);
protected:
    static std::string ZonePortion_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void SupplyMissingDefaultsIfAppropriate();
    bool AttributesAreConsistent() const;
    void ForceAttributeConsistency();
    void SwitchDefaultVariableNameToTrueName();

    // IDs that can be used to identify fields in case statements
    enum {
        ID_outputMeshType = 0,
        ID_listedVarNames,
        ID_zonePortions,
        ID_lowerBounds,
        ID_upperBounds,
        ID_defaultVarName,
        ID_defaultVarIsScalar
    };

private:
    int          outputMeshType;
    stringVector listedVarNames;
    intVector    zonePortions;
    doubleVector lowerBounds;
    doubleVector upperBounds;
    std::string  defaultVarName;
    bool         defaultVarIsScalar;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
