#ifndef LABELATTRIBUTES_H
#define LABELATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: LabelAttributes
//
// Purpose:
//    This class contains the fields that we need to set the attributes for the Label plot.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Oct 21 18:18:05 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class LabelAttributes : public AttributeSubject
{
public:
    enum LabelIndexDisplay
    {
        Natural,
        LogicalIndex,
        Index
    };
    enum LabelHorizontalAlignment
    {
        HCenter,
        Left,
        Right
    };
    enum LabelVerticalAlignment
    {
        VCenter,
        Top,
        Bottom
    };
    enum LabelDrawFacing
    {
        Front,
        Back,
        FrontAndBack
    };

    LabelAttributes();
    LabelAttributes(const LabelAttributes &obj);
    virtual ~LabelAttributes();

    virtual void operator = (const LabelAttributes &obj);
    virtual bool operator == (const LabelAttributes &obj) const;
    virtual bool operator != (const LabelAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectTextColor();
    void SelectTextLabel();

    // Property setting methods
    void SetLegendFlag(bool legendFlag_);
    void SetShowNodes(bool showNodes_);
    void SetShowCells(bool showCells_);
    void SetRestrictNumberOfLabels(bool restrictNumberOfLabels_);
    void SetDrawLabelsFacing(LabelDrawFacing drawLabelsFacing_);
    void SetShowSingleNode(bool showSingleNode_);
    void SetShowSingleCell(bool showSingleCell_);
    void SetUseForegroundTextColor(bool useForegroundTextColor_);
    void SetLabelDisplayFormat(LabelIndexDisplay labelDisplayFormat_);
    void SetNumberOfLabels(int numberOfLabels_);
    void SetTextColor(const ColorAttribute &textColor_);
    void SetTextHeight(float textHeight_);
    void SetTextLabel(const std::string &textLabel_);
    void SetHorizontalJustification(LabelHorizontalAlignment horizontalJustification_);
    void SetVerticalJustification(LabelVerticalAlignment verticalJustification_);
    void SetSingleNodeIndex(int singleNodeIndex_);
    void SetSingleCellIndex(int singleCellIndex_);

    // Property getting methods
    bool                 GetLegendFlag() const;
    bool                 GetShowNodes() const;
    bool                 GetShowCells() const;
    bool                 GetRestrictNumberOfLabels() const;
    LabelDrawFacing      GetDrawLabelsFacing() const;
    bool                 GetShowSingleNode() const;
    bool                 GetShowSingleCell() const;
    bool                 GetUseForegroundTextColor() const;
    LabelIndexDisplay    GetLabelDisplayFormat() const;
    int                  GetNumberOfLabels() const;
    const ColorAttribute &GetTextColor() const;
          ColorAttribute &GetTextColor();
    float                GetTextHeight() const;
    const std::string    &GetTextLabel() const;
          std::string    &GetTextLabel();
    LabelHorizontalAlignment GetHorizontalJustification() const;
    LabelVerticalAlignment GetVerticalJustification() const;
    int                  GetSingleNodeIndex() const;
    int                  GetSingleCellIndex() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string LabelIndexDisplay_ToString(LabelIndexDisplay);
    static bool LabelIndexDisplay_FromString(const std::string &, LabelIndexDisplay &);
protected:
    static std::string LabelIndexDisplay_ToString(int);
public:
    static std::string LabelHorizontalAlignment_ToString(LabelHorizontalAlignment);
    static bool LabelHorizontalAlignment_FromString(const std::string &, LabelHorizontalAlignment &);
protected:
    static std::string LabelHorizontalAlignment_ToString(int);
public:
    static std::string LabelVerticalAlignment_ToString(LabelVerticalAlignment);
    static bool LabelVerticalAlignment_FromString(const std::string &, LabelVerticalAlignment &);
protected:
    static std::string LabelVerticalAlignment_ToString(int);
public:
    static std::string LabelDrawFacing_ToString(LabelDrawFacing);
    static bool LabelDrawFacing_FromString(const std::string &, LabelDrawFacing &);
protected:
    static std::string LabelDrawFacing_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual bool ChangesRequireRecalculation(const LabelAttributes &) const;
private:
    bool           legendFlag;
    bool           showNodes;
    bool           showCells;
    bool           restrictNumberOfLabels;
    int            drawLabelsFacing;
    bool           showSingleNode;
    bool           showSingleCell;
    bool           useForegroundTextColor;
    int            labelDisplayFormat;
    int            numberOfLabels;
    ColorAttribute textColor;
    float          textHeight;
    std::string    textLabel;
    int            horizontalJustification;
    int            verticalJustification;
    int            singleNodeIndex;
    int            singleCellIndex;
};

#endif
