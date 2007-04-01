#ifndef QVIS_PREFERENCES_WINDOW_H
#define QVIS_PREFERENCES_WINDOW_H
#include <gui_exports.h>
#include <QvisPostableWindowObserver.h>
#include <TimeFormat.h>

class GlobalAttributes;
class QButtonGroup;
class QCheckBox;
class QSpinBox;

// ****************************************************************************
// Class: QvisPreferencesWindow
//
// Purpose: 
//   Defines QvisPreferencesWindow class.
//
// Programmer: Eric Brugger
// Creation:   Thu Mar 13 11:13:18 PST 2003
//
// Modifications:
//   Brad Whitlock, Fri Sep 5 15:41:31 PST 2003
//   Added a toggle for post windows when shown.
//
//   Brad Whitlock, Mon Oct 13 16:52:39 PST 2003
//   Added a toggle for timeState display mode.
//
//   Brad Whitlock, Fri Jan 30 14:23:26 PST 2004
//   I added a toggle for showing the selected files list.
//
//   Brad Whitlock, Fri Apr 9 14:12:16 PST 2004
//   I added a toggle for highlighting the selected files.
//
//   Brad Whitlock, Fri Aug 6 09:20:21 PDT 2004
//   I added toggles that let you set the prompting behavior for "make default"
//   and "automatically apply operator".
//
// ****************************************************************************

class GUI_API QvisPreferencesWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisPreferencesWindow(GlobalAttributes *subj,
                          const char *caption = 0,
                          const char *shortName = 0,
                          QvisNotepadArea *notepad = 0);
    virtual ~QvisPreferencesWindow();
    virtual void CreateWindowContents();

    void SetTimeStateFormat(const TimeFormat &fmt);
    void SetShowSelectedFiles(bool val);
    void SetAllowFileSelectionChange(bool val);
signals:
    void changeTimeFormat(const TimeFormat &);
    void showSelectedFiles(bool);
    void allowFileSelectionChange(bool);
public slots:
    virtual void apply();
protected:
    void UpdateWindow(bool doAll);
    void Apply(bool ignore = false);
private slots:
    void cloneWindowOnFirstRefToggled(bool val);
    void postWindowsWhenShownToggled(bool val);
    void makeDefaultConfirmToggled(bool val);
    void automaticallyApplyOperatorToggled(bool val);
    void handleTimeStateDisplayModeChange(int val);
    void timeStateNDigitsChanged(int val);
    void selectedFilesToggled(bool);
    void allowFileSelectionChangeToggled(bool);
private:
    QCheckBox        *cloneWindowOnFirstRefToggle;
    QCheckBox        *postWindowsWhenShownToggle;
    QCheckBox        *makeDefaultConfirmToggle;
    QCheckBox        *automaticallyApplyOperatorToggle;
    QCheckBox        *selectedFilesToggle;
    QCheckBox        *allowFileSelectionChangeToggle;
    QButtonGroup     *timeStateDisplayMode;
    QSpinBox         *timeStateNDigits;
    GlobalAttributes *atts;

    TimeFormat        tsFormat;
    bool              showSelFiles;
    bool              allowFileSelChange;
};

#endif
