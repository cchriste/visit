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

#ifndef QVISSEEDMEWINDOW_H
#define QVISSEEDMEWINDOW_H

#include <SeedMeAttributes.h>
#include <QvisPostableWindowSimpleObserver.h>

class QButtonGroup;
class QCheckBox;
class QLabel;
class QTextBrowser;
class QLineEdit;
class QSpinBox;
class QTabWidget;
class QVBox;
class QScrollBar;
class QPushButton;
class QFileSystemWatcher;
class QComboBox;

class QvisColorTableButton;
class QvisOpacitySlider;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisVariableButton;

// ****************************************************************************
// Class: QvisSeedMeWindow
//
// Purpose:
//    Defines QvisSeedMeWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class QvisSeedMeWindow : public QvisPostableWindowSimpleObserver
{
    Q_OBJECT
public:
    QvisSeedMeWindow(SeedMeAttributes *subj,
                     const QString &caption = QString::null,
                     const QString &shortName = QString::null,
                     QvisNotepadArea *notepad = 0);
    virtual ~QvisSeedMeWindow();
    virtual void CreateWindowContents();

    virtual void SubjectRemoved(Subject *subj);
signals:
    void runCommand(const QString &);
public slots:
    virtual void apply();
protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);
private slots:
    void collectionModeChanged(int val);
    void collectionIDProcessText();
    void sharingChanged(int val);
    void collectionTitleProcessText();
    void collectionDescriptionProcessText();
    void overwriteFilesChanged(bool val);
    void keyValueProcessText();
    void collectionEmailsProcessText();
    void uploadSequenceFileChanged(bool val);
    void sequenceTitleProcessText();
    void sequenceDescriptionProcessText();
    void createVideoChanged(bool val);
    void frameRateProcessText();
    void browse();
    void browseApiKey();
    void queryActionChanged(int val);
    void queryColIDProcessText();
    void queryKeyValueProcessText();
    void queryCollectionValuesChanged(int val);
    void downloadCollectionIDProcessText();
    void downloadTypeChanged(int val);
    void downloadNameProcessText();
    void directoryChanged(const QString & path);
    void quickSharingChanged(int val);
    void quickCollectionTitleProcessText();
    void quickCollectionEmailsProcessText();
    void quickFrameRateProcessText();
    void quickDownloadTypeChanged(int val);
    void quickDownload();


    void ResetForm();
    void ClearLog();
private:
    QWidget *CreateQuickUploadTab();
    QWidget *CreateUploadTab();
    QWidget *CreateQueryTab();
    QWidget *CreateDownloadTab();
    QWidget *CreateSettingsTab();

    QTabWidget   *tabs;

    // Upload tab
    QWidget      *collectionMode;
    QButtonGroup *collectionModeButtonGroup;
    QLineEdit *collectionID;
    QWidget      *sharing;
    QButtonGroup *sharingButtonGroup;
    QLineEdit *collectionTitle;
    QLineEdit *collectionDescription;
    QLineEdit *collectionCredits;
    QLineEdit *collectionLicense;
    QLineEdit *keyValue;
    QLineEdit *collectionEmails;
    QCheckBox *overwriteFiles;
    QCheckBox *uploadSequenceFile;
    QLineEdit *sequenceTitle;
    QLineEdit *sequenceDescription;
    QCheckBox *createVideo;
    QLineEdit *frameRate;
    QWidget      *queryAction;
    QButtonGroup *queryActionButtonGroup;
    QLineEdit *queryColID;
    QLineEdit *queryKeyValue;
    QWidget      *queryCollectionValues;
    QButtonGroup *queryCollectionValuesButtonGroup;
    QLineEdit *downloadCollectionID;
    QWidget      *downloadType;
    QButtonGroup *downloadTypeButtonGroup;
    QLineEdit *downloadName;
    QWidget      *quickSharing;
    QButtonGroup *quickSharingButtonGroup;
    QLineEdit *quickCollectionTitle;
    QLineEdit *quickCollectionEmails;
    QLineEdit *quickFrameRate;
    QComboBox      *quickDownloadType;
    QLabel *collectionModeLabel;
    QLabel *collectionIDLabel;
    QLabel *sharingLabel;
    QLabel *collectionTitleLabel;
    QLabel *collectionDescriptionLabel;
    QLabel *keyValueLabel;
    QLabel *collectionCreditsLabel;
    QLabel *collectionLicenseLabel;
    QLabel *collectionEmailsLabel;
    QLabel *currentTitleLabel;
    QLabel *currentDescriptionLabel;
    QLabel *sequenceTitleLabel;
    QLabel *sequenceDescriptionLabel;
    QLabel *frameRateLabel;
    QLabel *queryActionLabel;
    QLabel *queryColIDLabel;
    QLabel *queryKeyValueLabel;
    QLabel *queryCollectionValuesLabel;
    QLabel *downloadCollectionIDLabel;
    QLabel *downloadTypeLabel;
    QLabel *downloadNameLabel;

    QLabel *collectionsLink, *quickCollectionsLink;

    QLabel *helpLabel;
    QLabel *helpLabelWarning;

    QLabel *confLocationLabel;

    QFileSystemWatcher *seedmeWatcher;

    QTextBrowser *statusLabel;

    SeedMeAttributes *atts;

    QStringList uploadFiles;

    QPushButton *resetFormButton, *clearLogButton, *quickDownloadButton, *submitButton;

    void updateStatus(QString str);

    QString apikeyFile;

    const char* configFile;
};


#endif
