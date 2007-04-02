/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
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

// ****************************************************************************
// Class: ThresholdAttributes
//
// Purpose:
//    This class contains attributes for the threshold operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Sep 13 08:54:28 PDT 2005
//
// Modifications:
//   
//   Mark Blair, Tue Mar  7 13:25:00 PST 2006
//   Upgraded to support multiple threshold variables.
//
//   Mark Blair, Tue Aug  8 17:47:00 PDT 2006
//   Now accommodates an empty list of threshold variables.
//
//   Mark Blair, Thu Sep 28 12:07:05 PDT 2006
//   Added SupplyMissingDefaultsIfAppropriate.
//
//   Mark Blair, Tue Oct  3 13:19:11 PDT 2006
//   Added SwitchDefaultToTrueVariableNameIfScalar and flag that indicates
//   whether default plot variable is scalar.
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

    virtual void operator = (const ThresholdAttributes &obj);
    virtual bool operator == (const ThresholdAttributes &obj) const;
    virtual bool operator != (const ThresholdAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property setting methods
    void SetOutputMeshType(OutputMeshType outputMeshType_);
    void SetOutputMeshType(int outputMeshType_);
    void SetListedVariables(const stringVector &listedVarNames_);
    void SetShownVariablePosition(int shownVarPosition_);
    void SetZonePortions(const std::vector<ZonePortion> &zonePortions_);
    void SetZonePortions(const intVector &zonePortions_);
    void SetLowerBounds(const doubleVector &lowerBounds_);
    void SetUpperBounds(const doubleVector &upperBounds_);
    void SetDefaultVarName(const std::string &defaultVarName_);
    void SetDefaultVarIsScalar(bool defaultVarIsScalar_);

    void SupplyMissingDefaultsIfAppropriate();
    bool AttributesAreConsistent() const;
    void SwitchDefaultToTrueVariableNameIfScalar();

    // Property getting methods
    OutputMeshType         GetOutputMeshType() const;
    const std::string     &GetShownVariable() const;
          std::string     &GetShownVariable();
    ZonePortion            GetZonePortion() const;
    double                 GetLowerBound() const;
    double                 GetUpperBound() const;
    const std::string     &GetDefaultVarName() const;
          std::string     &GetDefaultVarName();
    bool                   GetDefaultVarIsScalar() const;

    const stringVector    &GetListedVariables() const;
    const intVector       &GetZonePortions() const;
    const doubleVector    &GetLowerBounds() const;
    const doubleVector    &GetUpperBounds() const;

    // Property changing methods
    void                   ChangeZonePortion(ZonePortion newZonePortion_);
    void                   ChangeZonePortion(int newZonePortion_);
    void                   ChangeLowerBound(double newLowerBound_);
    void                   ChangeUpperBound(double newUpperBound_);

    void                   InsertVariable(const std::string &variable_);
    void                   DeleteVariable(const std::string &variable_);
    void                   SwapVariable(const std::string &variable_);
    void                   ShowPreviousVariable();
    void                   ShowNextVariable();

    // Property selection methods
    virtual void SelectAll();
    void SelectVariable();

    // Python compatibility methods
    void                   SetListedVarNames(const stringVector &listedVarNames_);
    void                   SetShownVarPosition(int shownVarPosition_);

    const stringVector    &GetListedVarNames() const;
          stringVector    &GetListedVarNames();
    int                    GetShownVarPosition() const;
    intVector             &GetZonePortions();
    doubleVector          &GetLowerBounds();
    doubleVector          &GetUpperBounds();

    void                   SelectOutputMeshType();
    void                   SelectListedVarNames();
    void                   SelectShownVarPosition();
    void                   SelectZonePortions();
    void                   SelectLowerBounds();
    void                   SelectUpperBounds();
    void                   SelectDefaultVarName();
    void                   SelectDefaultVarIsScalar();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string OutputMeshType_ToString(OutputMeshType);
    static bool OutputMeshType_FromString(const std::string &, OutputMeshType &);

    static std::string ZonePortion_ToString(ZonePortion);
    static bool ZonePortion_FromString(const std::string &, ZonePortion &);

protected:
    static std::string OutputMeshType_ToString(int);
    static std::string ZonePortion_ToString(int);

public:
    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    int                 outputMeshType;
    stringVector        listedVarNames;
    int                 shownVarPosition;
    intVector           zonePortions;
    doubleVector        lowerBounds;
    doubleVector        upperBounds;
    std::string         defaultVarName;
    bool                defaultVarIsScalar;
};

#endif
