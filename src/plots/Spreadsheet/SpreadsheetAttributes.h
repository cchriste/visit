/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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

#ifndef SPREADSHEETATTRIBUTES_H
#define SPREADSHEETATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>
#include <PlaneAttributes.h>
#include <PickAttributes.h>
#include <math.h>

// ****************************************************************************
// Class: SpreadsheetAttributes
//
// Purpose:
//    Contains the attributes for the visual spreadsheet.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Nov 28 15:15:50 PST 2007
//
// Modifications:
//   
// ****************************************************************************

class SpreadsheetAttributes : public AttributeSubject
{
public:
    enum NormalAxis
    {
        X,
        Y,
        Z
    };

    SpreadsheetAttributes();
    SpreadsheetAttributes(const SpreadsheetAttributes &obj);
    virtual ~SpreadsheetAttributes();

    virtual SpreadsheetAttributes& operator = (const SpreadsheetAttributes &obj);
    virtual bool operator == (const SpreadsheetAttributes &obj) const;
    virtual bool operator != (const SpreadsheetAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSubsetName();
    void SelectFormatString();
    void SelectColorTableName();
    void SelectTracerColor();
    void SelectCurrentPick();
    void SelectPastPicks();
    void SelectCurrentPickLetter();
    void SelectPastPickLetters();
    void SelectSpreadsheetFont();

    // Property setting methods
    void SetSubsetName(const std::string &subsetName_);
    void SetFormatString(const std::string &formatString_);
    void SetUseColorTable(bool useColorTable_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetShowTracerPlane(bool showTracerPlane_);
    void SetTracerColor(const ColorAttribute &tracerColor_);
    void SetNormal(NormalAxis normal_);
    void SetSliceIndex(int sliceIndex_);
    void SetCurrentPick(const double *currentPick_);
    void SetCurrentPickValid(bool currentPickValid_);
    void SetPastPicks(const doubleVector &pastPicks_);
    void SetCurrentPickLetter(const std::string &currentPickLetter_);
    void SetPastPickLetters(const stringVector &pastPickLetters_);
    void SetSpreadsheetFont(const std::string &spreadsheetFont_);
    void SetShowPatchOutline(bool showPatchOutline_);
    void SetShowCurrentCellOutline(bool showCurrentCellOutline_);

    // Property getting methods
    const std::string    &GetSubsetName() const;
          std::string    &GetSubsetName();
    const std::string    &GetFormatString() const;
          std::string    &GetFormatString();
    bool                 GetUseColorTable() const;
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    bool                 GetShowTracerPlane() const;
    const ColorAttribute &GetTracerColor() const;
          ColorAttribute &GetTracerColor();
    NormalAxis           GetNormal() const;
    int                  GetSliceIndex() const;
    const double         *GetCurrentPick() const;
          double         *GetCurrentPick();
    bool                 GetCurrentPickValid() const;
    const doubleVector   &GetPastPicks() const;
          doubleVector   &GetPastPicks();
    const std::string    &GetCurrentPickLetter() const;
          std::string    &GetCurrentPickLetter();
    const stringVector   &GetPastPickLetters() const;
          stringVector   &GetPastPickLetters();
    const std::string    &GetSpreadsheetFont() const;
          std::string    &GetSpreadsheetFont();
    bool                 GetShowPatchOutline() const;
    bool                 GetShowCurrentCellOutline() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string NormalAxis_ToString(NormalAxis);
    static bool NormalAxis_FromString(const std::string &, NormalAxis &);
protected:
    static std::string NormalAxis_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const SpreadsheetAttributes &) const;
private:
    std::string    subsetName;
    std::string    formatString;
    bool           useColorTable;
    std::string    colorTableName;
    bool           showTracerPlane;
    ColorAttribute tracerColor;
    int            normal;
    int            sliceIndex;
    double         currentPick[3];
    bool           currentPickValid;
    doubleVector   pastPicks;
    std::string    currentPickLetter;
    stringVector   pastPickLetters;
    std::string    spreadsheetFont;
    bool           showPatchOutline;
    bool           showCurrentCellOutline;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
