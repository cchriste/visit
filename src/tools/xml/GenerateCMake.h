/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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

#ifndef GENERATE_CMAKE_H
#define GENERATE_CMAKE_H

#include <QTextStream>
#include "Field.h"
#include <visit-config.h> // for the plugin extension.
#include "Plugin.h"

// ****************************************************************************
//  File:  GenerateCMake
//
//  Purpose:
//    Contains a set of classes which override the default implementation
//    to create cmake input for the plugin.
//
//  Note: This file overrides --
//    Plugin
//
//  Programmer:  Brad Whitlock, 
//  Creation:    Thu Jan 29 13:44:46 PST 2009
//
//  Modifications:
//    Brad Whitlock, Fri Nov  6 11:15:11 PST 2009
//    Handle serial and parallel engine libs.
//
//    Brad Whitlock, Mon Nov 23 15:19:10 PST 2009
//    I added server components and engine only builds.
//
//    David Camp, Thu Jan 14 17:56:29 PST 2010
//    Added the ADD_TARGET_DEFINITIONS function to define ENGINE for plots.
//    
//    Kathleen Bonnell, Tue Jan 26 20:32:55 MST 2010
//    Remove setting of LIBRARY_OUTPUT_PATH, (set by parent instead). Add
//    call to VISIT_PLUGIN_TARGET_PREFIX macro.
//
//    Brad Whitlock, Wed Feb 10 16:36:00 PST 2010
//    I made all of the database plugins use the ADD_TARGET_DEFINITIONS
//    function.
//
//    Eric Brugger, Wed Feb 24 13:00:54 PST 2010
//    I modified the database plugins to list the include paths specified
//    in the xml file before the VTK include paths.
//
//    Eric Brugger, Fri Feb 26 09:47:00 PST 2010
//    I modified the database plugins to list the include paths specified
//    in the xml file before any of the VisIt include paths.  I also modified
//    the database plugins to also treat all flags in CXXFLAGS that start
//    with "-I" as include paths.
//
//    Kathleen Bonnell, Fri May 21 14:15:23 MST 2010 
//    Add DLL_NETCDF, _CGNSDLL EXODUSII_BUILD_SHARED_LIBS defines for 
//    windows projects linking with NETCDF, CGNS or EXODUSII.
//
//    Kathleen Bonnell, Thu May 27 14:59:13 MST 2010 
//    Add some more defines for HDF4, discovered as necessary when compiling
//    with Visual Studio 9.
//
//    Kathleen Bonnell, Fri Sep 24 11:25:32 MST 2010 
//    Add ENGINE target definition for operators if they contain 
//    engine-specific code.
//
//    Kathleen Bonnell, Tue Nov 16 16:26:47 PST 2010
//    Remove logic for mesa.  Add newline after each extraInclude for 
//    legibility in the CMakeLists.txt files.
//
//    David Camp, Wed Nov 17 14:54:02 PST 2010
//    Added the LIBS libraries to the Plot and Operators, did the samething
//    the database code was doing. Also added the link dirs from the ldflags.
//
//    Kathleen Bonnell, Fri Sep 24 11:25:32 MST 2010 
//    Fix windows issues with viewer and gui libs building against an 
//    installed version of VisIt.  Convert Windows paths to CMake paths 
//    since we are creating a CMake file.
//
//    Kathleen Bonnell, Tue Jan  4 08:38:03 PST 2011
//    Fix CGNS dll define, due to update of cgns library.
//    Add call to VISIT_PLUGIN_TARGET_FOLDER for project grouping in VS.
//
//    Eric Brugger, Fri Jan  7 13:38:59 PST 2011
//    I replaced the BOXLIB2D and BOXLIB3D variables with just BOXLIB.
//
//    Kathleen Bonnell, Tue Jan 11 17:06:21 MST 2011 
//    Removed setting EXODUSII_BUILD_SHARED_LIBS definition.
//
//    Kathleen Bonnell, Thu Jan 13 17:54:38 MST 2011
//    Only use VISIT_PLUGIN_TARGET_FOLDER if building from dev.
//
//    Brad Whitlock, Wed Feb 23 15:24:48 PST 2011
//    Enable Fortran language compilation if the user added Fortran code to the
//    list of files.
//
//    Kathleen Biagas, Fri Nov 18 10:09:26 MST 2011
//    Add plugin name to VISIT_PLUGIN_TARGET_FOLDER args. Eases building/
//    debugging individual plugins with Visual Studio when grouped by name.
//
//    Kathleen Biagas, Tue Nov 22 14:39:51 PST 2011
//    Remove VISIT_PLUGIN_TARGET_PREFIX in favor of VISIT_PLUGIN_TARGET_RUNTIME.
//
//    Kathleen Biagas, Mon Jun 18 10:49:07 MST 2012
//    Set VISIT_ARCHIVE_DIR on windows to be /lib. Change minimum CMake
//    version to 2.8.8.
//
//    Kathleen Biagas, Mon Jul 30 15:40:10 MST 2012
//    No longer add definition _HDF5USEDLL_ for hdf5 based plugins, as this
//    is now predefined in an hdf5 header.
//
//    Kathleen Biagas, Wed Oct  9 10:01:15 PDT 2013
//    Added handling of 'Code' and 'Condition' keywords in codefile. 
//    'Condition' allows for conditional includes, definitions and links.
//
//    Kathleen Biagas, Tue Oct 29 16:04:19 MST 2013
//    For extraIncludes specified in CXXFLAGS, check for use of 
//    ${VISIT_INCLUDE_DIR} and correct it if building against public VisIt.
//
//    Eric Brugger, Wed May 21 14:48:11 PDT 2014
//    I added support for EAVL.
//
// ****************************************************************************

class CMakeGeneratorPlugin : public Plugin
{
  public:
    CMakeGeneratorPlugin(const QString &n,const QString &l,const QString &t,
        const QString &vt,const QString &dt, const QString &v, const QString &ifile,
        bool hw, bool ho, bool onlyengine, bool noengine) : 
        Plugin(n,l,t,vt,dt,v,ifile,hw,ho,onlyengine,noengine)
    {
        defaultgfiles.clear();
        defaultsfiles.clear();
        defaultvfiles.clear();
        defaultmfiles.clear();
        defaultefiles.clear();
        defaultwfiles.clear();
        if (type == "database")
        {
            QString filter = QString("avt") + name + "FileFormat.C";
            defaultmfiles.push_back(filter);
            defaultefiles.push_back(filter);
            if (haswriter)
                defaultefiles.push_back(QString("avt") + name + "Writer.C");
            if (hasoptions)
            {
                QString options = QString("avt") + name + QString("Options.C");
                defaultmfiles.push_back(options);
                defaultefiles.push_back(options);
            }
        }
        else if (type == "plot")
        {
            QString filter = QString("avt") + name + "Filter.C";
            defaultvfiles.push_back(filter);
            defaultefiles.push_back(filter);
            QString widgets = QString("Qvis") + name + "PlotWindow.h";
            defaultwfiles.push_back(widgets);
        }
        else if (type == "operator")
        {
            QString filter = QString("avt") + name + "Filter.C";
            defaultefiles.push_back(filter);
        }
    }

    virtual ~CMakeGeneratorPlugin()
    {
    }

    void
    GetFilesWith(const QString &name, const std::vector<QString> &input, 
                 std::set<QString> &output)
    {
         for(size_t i = 0; i < input.size(); ++i)
         {
             if(input[i].indexOf(name) != -1)
                 output.insert(input[i]);
         }
    }

    QString
    ConvertDollarParenthesis(const QString &s) const
    {
        QString retval(s);
        retval = retval.replace("$(", "${");
        retval = retval.replace(")", "}");
        return retval;
    }

    QString
    ToString(const std::vector<QString> &vec, bool withNewline=false) const
    {
        QString s;
        if (withNewline)
        {
            for(size_t i = 0; i < vec.size(); ++i)
                s += (ConvertDollarParenthesis(vec[i]) + "\n");
        }
        else 
        {
            for(size_t i = 0; i < vec.size(); ++i)
                s += (ConvertDollarParenthesis(vec[i]) + " ");
        }
        return s;
    }

    QString
    VisItIncludeDir() const
    {
        return using_dev ? "${VISIT_INCLUDE_DIR}" : "${VISIT_INCLUDE_DIR}/visit";
    }

    QString
    ConvertToProperVisItIncludeDir(const QString &s) const
    {
        QString VID = VisItIncludeDir();
        QString retval(s);
        if (!s.startsWith(VID))
        {
            if (!using_dev && s.startsWith("${VISIT_INCLUDE_DIR}"))
                retval = VID + s.right(s.length()-20);
        }
        return retval;
    }

#ifdef _WIN32
    QString
    ToCMakePath(const QString &s) const
    {
        char exppath[MAX_PATH];
        ExpandEnvironmentStrings(s.toStdString().c_str(), exppath, MAX_PATH);
        QString retval(exppath);
        retval = retval.replace("\\", "/");
        return retval; 
    }
#endif

    bool
    HasCondition(const QString &c) const
    {
        bool retval = false;
        if (atts != NULL && atts->codeFile != NULL)
        {
            retval =  atts->codeFile->HasCondition(c);
        }
        return retval;
    }

    bool
    GetCondition(const QString &c, QString &cond, QString &val) const
    {
        bool retval = false;
        if (HasCondition(c))
        {
            QStringList targets, first, second;
            atts->codeFile->GetCondition(c, targets, first, second);
            for (int i = 0; i < targets.size(); ++i)
            {
                if (targets[i] == "xml2cmake")
                {
                    cond = first[i];
                    val = second[i];
                    retval = true;
                    break;
                }
            }
        }
        return retval;
    }

    void WriteCMake_ConditionalIncludes(QTextStream &out)
    {
        QString condition, incs;
        if(GetCondition("Includes:", condition, incs))
        {
            out << endl;
            out << "IF(" << condition << ")" << endl;
            out << "    INCLUDE_DIRECTORIES(";
            out << incs;
            out << ")" << endl;
            out << "ENDIF(" << condition << ")" << endl;
        }
    }

    void WriteCMake_ConditionalDefinitions(QTextStream &out)
    {
        QString condition, defs;
        if(GetCondition("Definitions:", condition, defs))
        {
            out << "IF(" << condition << ")" << endl;
            out << "    ADD_DEFINITIONS(";
            out << defs;
            out << ")" << endl;
            out << "ENDIF(" << condition << ")" << endl;
        }
    }

    void WriteCMake_ConditionalLinkDirs(QTextStream &out, std::vector<QString> &linkDirs)
    {
        QString condition, dirs;
        if (GetCondition("LinkDirectories:", condition, dirs))
        {
            out << "IF(" << condition << ")" << endl;
            out << "    LINK_DIRECTORIES(" << ToString(linkDirs) << dirs << ")" << endl;
            out << "ELSE(" << condition << ")" << endl;
            out << "    LINK_DIRECTORIES(" << ToString(linkDirs) << ")" << endl;
            out << "ENDIF(" << condition << ")" << endl;
        }
    }

    void WriteCMake_ConditionalTargetLinks(QTextStream &out, const QString &target, const char *libType, const char *plugType, const char *indent)
    {
        QString c(libType);
        c += "LinkLibraries:";
        QString condition, link;
        if (GetCondition(c, condition, link))
        {
            out << indent << "IF(" << condition << ")" << endl;
            out << indent << "    TARGET_LINK_LIBRARIES(" << libType << target << plugType << " " << link << ")" << endl;
            out << indent << "ENDIF(" << condition << ")" << endl;
        }
    }

    void WriteCMake_ConditionalSources(QTextStream &out, const char *libType, const char *indent)
    {
        QString c(libType);
        c += "Sources:";
        QString condition, src;
        if (GetCondition(c, condition, src))
        {
            out << indent << "IF(" << condition << ")" << endl;
            out << indent << "    SET(LIB" << libType << "_SOURCES ${LIB" << libType << "_SOURCES} " << src << ")" << endl;
            out << indent << "ENDIF(" << condition << ")" << endl;
            out << endl;
        }
    }

    void
    WriteCMake_AdditionalCode(QTextStream &out)
    {
        if (atts != NULL && atts->codeFile != NULL)
        {
            QStringList targets, names, first, second;
            atts->codeFile->GetAllCodes(targets, names, first, second);
            for (int i = 0; i < targets.size(); ++i)
            {
                if (targets[i] == "xml2cmake")
                {
                    if (!first[i].isEmpty())
                    {
                        out << first[i] << endl;
                    }
                    if (!second[i].isEmpty())
                    {
                        out << second[i] << endl;
                    }
                }
            }
        }
    }

    void WriteCMake_PlotOperator_Includes(QTextStream &out, bool isOperator)
    {
        // take any ${} from the CXXFLAGS to mean a variable that contains 
        // include directories.
        std::vector<QString> extraIncludes;
        for (size_t i=0; i<cxxflags.size(); i++)
        {
            if(cxxflags[i].startsWith("${"))
                 extraIncludes.push_back(ConvertToProperVisItIncludeDir(cxxflags[i]));
            else if(cxxflags[i].startsWith("$("))
                 extraIncludes.push_back(ConvertToProperVisItIncludeDir(ConvertDollarParenthesis(cxxflags[i])));
        }

        out << endl
            << "IF(VISIT_PYTHON_SCRIPTING)" << endl;
        out << "    SET(PYINCLUDES ${PYTHON_INCLUDE_PATH} " << VisItIncludeDir() << "/visitpy/visitpy)" << endl;
        out << "ENDIF(VISIT_PYTHON_SCRIPTING)" << endl << endl;

        // Includes
        out << "INCLUDE_DIRECTORIES(" << endl;
        out << "${CMAKE_CURRENT_SOURCE_DIR}" << endl;
        out << "${VISIT_COMMON_INCLUDES}" << endl;
        out << VisItIncludeDir() << "/avt/DBAtts/MetaData" << endl;
        out << VisItIncludeDir() << "/avt/DBAtts/SIL" << endl;
        out << VisItIncludeDir() << "/avt/Database/Database" << endl;
        if(isOperator)
        {
            out << VisItIncludeDir() << "/avt/Expressions/Abstract" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/CMFE" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/Conditional" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/Derivations" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/General" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/ImageProcessing" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/Management" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/Math" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/MeshQuality" << endl;
            out << VisItIncludeDir() << "/avt/Expressions/TimeIterators" << endl;
        }
        out << VisItIncludeDir() << "/avt/FileWriter" << endl;
        out << VisItIncludeDir() << "/avt/Filters" << endl;
        out << VisItIncludeDir() << "/avt/IVP" << endl;
        out << VisItIncludeDir() << "/avt/Math" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/AbstractFilters" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Data" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Pipeline" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Sinks" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Sources" << endl;
        out << VisItIncludeDir() << "/avt/Plotter" << endl;
        out << VisItIncludeDir() << "/avt/QtVisWindow" << endl;
        out << VisItIncludeDir() << "/avt/View" << endl;
        out << VisItIncludeDir() << "/avt/VisWindow/Colleagues" << endl;
        out << VisItIncludeDir() << "/avt/VisWindow/Interactors" << endl;
        out << VisItIncludeDir() << "/avt/VisWindow/Proxies" << endl;
        out << VisItIncludeDir() << "/avt/VisWindow/Tools" << endl;
        out << VisItIncludeDir() << "/avt/VisWindow/VisWindow" << endl;
        out << VisItIncludeDir() << "/gui" << endl;
        if(isOperator)
        {
            out << VisItIncludeDir() << "/mdserver/proxy" << endl;
            out << VisItIncludeDir() << "/mdserver/rpc" << endl;
        }
        out << VisItIncludeDir() << "/viewer/main" << endl;
        out << VisItIncludeDir() << "/viewer/proxy" << endl;
        out << VisItIncludeDir() << "/viewer/rpc" << endl;
        out << VisItIncludeDir() << "/winutil" << endl;
        out << VisItIncludeDir() << "/visit_vtk/full" << endl;
        out << VisItIncludeDir() << "/visit_vtk/lightweight" << endl;
        out << "${QT_INCLUDE_DIR}" << endl;
        out << "${QT_QTCORE_INCLUDE_DIR}" << endl;
        out << "${QT_QTGUI_INCLUDE_DIR}" << endl;
        out << "${EAVL_INCLUDE_DIR} " << endl;
        out << "${VTK_INCLUDE_DIRS} " << endl;
        out << "${PYINCLUDES}" << endl;
        if(extraIncludes.size() > 0)
            out << ToString(extraIncludes, true);
        out << ")" << endl;
    }

    bool CustomFilesUseFortran(const std::vector<QString> &files) const
    {
        const char *ext[] = {".f", ".f77", ".f90", ".f95", ".for", 
                             ".F", ".F77", ".F90", ".F95", ".FOR"};
        for(size_t i = 0; i < files.size(); ++i)
        {
            for(int j = 0; j < 10; ++j)
            {
                if(files[i].endsWith(ext[j]))
                    return true;
            }
        }
        return false;
    }

    void WriteCMake_Plot(QTextStream &out, 
                         const QString &guilibname, 
                         const QString &viewerlibname)
    {
        bool useFortran = false;

        out << "PROJECT(" << name<< ")" << endl;
        out << endl;
        if (using_dev)
        {
        out << "INCLUDE(${VISIT_SOURCE_DIR}/CMake/PluginMacros.cmake)" <<endl;
        out << endl;
        }
        out << "SET(COMMON_SOURCES" << endl;
        out << name << "PluginInfo.C" << endl;
        out << name << "CommonPluginInfo.C" << endl;
        out << atts->name << ".C" << endl;
        out << ")" << endl;
        out << endl;
        out << "SET(LIBI_SOURCES " << endl;
        out << name << "PluginInfo.C" << endl;
        out << ")" << endl;
        out << endl;
        WriteCMake_ConditionalSources(out, "I", "");

        // libG sources
        out << "SET(LIBG_SOURCES" << endl;
        out << name << "GUIPluginInfo.C" << endl;
        out << "Qvis" << name << "PlotWindow.C" << endl;
        out << "${COMMON_SOURCES}" << endl;
        if (customgfiles)
        {
            useFortran |= CustomFilesUseFortran(gfiles);
            for (size_t i=0; i<gfiles.size(); i++)
                out << gfiles[i] << endl;
        }
        else
            for (size_t i=0; i<defaultgfiles.size(); i++)
                out << defaultgfiles[i] << endl;
        out << ")" << endl;
        out << "SET(LIBG_MOC_SOURCES" << endl;
        out << "Qvis" << name << "PlotWindow.h" << endl;
        if (customwfiles)
            for (size_t i=0; i<wfiles.size(); i++)
                out << wfiles[i] << endl;
        out << ")" << endl;
        out << endl;
        WriteCMake_ConditionalSources(out, "G", "");

        // libV sources
        out << "SET(LIBV_SOURCES" << endl;
        out << name<<"ViewerPluginInfo.C" << endl;
        out << "avt"<<name<<"Plot.C" << endl;
        if (customvfiles)
        {
            useFortran |= CustomFilesUseFortran(vfiles);
            for (size_t i=0; i<vfiles.size(); i++)
                out << vfiles[i] << endl;
        }
        else
            for (size_t i=0; i<defaultvfiles.size(); i++)
                out << defaultvfiles[i] << endl;
        out << "${COMMON_SOURCES}" << endl;
        out << ")" << endl;
        if (customvwfiles)
        {
            out << "SET(LIBV_MOC_SOURCES" << endl;
            for (size_t i=0; i<vwfiles.size(); i++)
                out << vwfiles[i] << endl;
            out << ")" << endl;
        }
        out << endl;
        WriteCMake_ConditionalSources(out, "V", "");

        // libE sources
        out << "SET(LIBE_SOURCES" << endl;
        out << name<<"EnginePluginInfo.C" << endl;
        out << "avt"<<name<<"Plot.C" << endl;
        if (customefiles)
        {
            useFortran |= CustomFilesUseFortran(efiles);
            for (size_t i=0; i<efiles.size(); i++)
                out << efiles[i] << endl;
        }
        else
            for (size_t i=0; i<defaultefiles.size(); i++)
                out << defaultefiles[i] << endl;
        out << "${COMMON_SOURCES}" << endl;
        out << ")" << endl;
        out << endl;
        WriteCMake_ConditionalSources(out, "E", "");

        if(useFortran)
        {
            out << "ENABLE_LANGUAGE(Fortran)" << endl;
        }

        // Special rules for OpenGL sources.
        std::set<QString> openglFiles;
        GetFilesWith("OpenGL", customvfiles ? vfiles : defaultvfiles, openglFiles);
        GetFilesWith("OpenGL", customefiles ? efiles : defaultefiles, openglFiles);
        if(openglFiles.size() > 0)
        {
            out << "IF (NOT WIN32)" << endl;
            out << "    SET_SOURCE_FILES_PROPERTIES(";
            for(std::set<QString>::iterator it = openglFiles.begin();
                it != openglFiles.end(); ++it)
            {
                 out << *it << " ";
            }
            out << "\n        PROPERTIES" << endl;
            out << "        COMPILE_FLAGS \"-I ${OPENGL_INCLUDE_DIR}\"" << endl;
            out << "    )" << endl;
            out << "ENDIF (NOT WIN32)" << endl;
        }

        WriteCMake_PlotOperator_Includes(out, false);
        WriteCMake_ConditionalIncludes(out);

        // Pass other CXXFLAGS
        for (size_t i=0; i<cxxflags.size(); i++)
        {
            if(!cxxflags[i].startsWith("${") && !cxxflags[i].startsWith("$("))
                 out << "ADD_DEFINITIONS(\"" << cxxflags[i] << "\")" << endl;
        }
        out << endl;

        WriteCMake_ConditionalDefinitions(out);

#if 0
        if (installpublic)
            out << "SET(LIBRARY_OUTPUT_PATH " << visitplugdirpub << ")" << endl;
        else if (installprivate)
            out << "SET(LIBRARY_OUTPUT_PATH " << visitplugdirpri << ")" << endl;
        else
#endif

        out << endl;
        std::vector<QString> linkDirs;
        linkDirs.push_back("${VISIT_LIBRARY_DIR}");
        linkDirs.push_back("${QT_LIBRARY_DIR}");
        linkDirs.push_back("${GLEW_LIBRARY_DIR}");
        linkDirs.push_back("${EAVL_LIBRARY_DIR}");
        linkDirs.push_back("${VTK_LIBRARY_DIRS}");
        // Extract extra link directories from LDFLAGS if they have ${},$(),-L
        for (size_t i=0; i<ldflags.size(); i++)
        {
            if(ldflags[i].startsWith("${") || ldflags[i].startsWith("$("))
                 linkDirs.push_back(ldflags[i]);
            else if(ldflags[i].startsWith("-L"))
                 linkDirs.push_back(ldflags[i].right(ldflags[i].size()-2));
        }

        if (HasCondition("LinkDirectories:"))
        {
            WriteCMake_ConditionalLinkDirs(out, linkDirs);
        }
        else
        {
            out << "LINK_DIRECTORIES(" << ToString(linkDirs) << ")" << endl;
        }
        out << endl;
        out << "ADD_LIBRARY(I"<<name<<"Plot ${LIBI_SOURCES})" << endl;
        out << "TARGET_LINK_LIBRARIES(I"<<name<<"Plot visitcommon)" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "I", "Plot", "");
        out << "SET(INSTALLTARGETS I"<<name<<"Plot)" << endl;
        out << endl;

        out << "IF(NOT VISIT_SERVER_COMPONENTS_ONLY AND NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)" << endl;
        out << "    QT_WRAP_CPP(G" << name << "Plot LIBG_SOURCES ${LIBG_MOC_SOURCES})" << endl;
        out << "    ADD_LIBRARY(G"<<name<<"Plot ${LIBG_SOURCES})" << endl;
        out << "    TARGET_LINK_LIBRARIES(G" << name << "Plot visitcommon "
            << guilibname << " " << ToString(libs) << ToString(glibs) 
            << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "G", "Plot", "    ");
        out << endl;

        if (customvwfiles)
            out << "    QT_WRAP_CPP(V" << name << "Plot LIBV_SOURCES ${LIBV_MOC_SOURCES})" << endl;
        out << "    ADD_LIBRARY(V"<<name<<"Plot ${LIBV_SOURCES})" << endl;
        out << "    TARGET_LINK_LIBRARIES(V" << name << "Plot visitcommon "
            << viewerlibname << " " << ToString(libs) << ToString(vlibs) 
            << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "V", "Plot", "    ");
        out << endl;
        out << "    SET(INSTALLTARGETS ${INSTALLTARGETS} G"<<name<<"Plot V"<<name<<"Plot)" << endl;
        out << endl;
        // libS sources
        out << "    IF(VISIT_PYTHON_SCRIPTING)" << endl;
        out << "        SET(LIBS_SOURCES" << endl;
        out << "            " << name<<"ScriptingPluginInfo.C" << endl;
        out << "            Py"<<atts->name<<".C" << endl;
        if (customsfiles)
            for (size_t i=0; i<sfiles.size(); i++)
                out << "            " << sfiles[i] << endl;
        else
            for (size_t i=0; i<defaultsfiles.size(); i++)
                out << "            " << defaultsfiles[i] << endl;
        out << "            ${COMMON_SOURCES}" << endl;
        out << "        )" << endl;
        WriteCMake_ConditionalSources(out, "S", "        ");
        out << "        ADD_LIBRARY(S"<<name<<"Plot ${LIBS_SOURCES})" << endl;
        out << "        TARGET_LINK_LIBRARIES(S" << name
            << "Plot visitcommon visitpy ${PYTHON_LIBRARY})" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "S", "Plot", "        ");
        out << "        SET(INSTALLTARGETS ${INSTALLTARGETS} S" << name
            << "Plot)" << endl;
        out << "    ENDIF(VISIT_PYTHON_SCRIPTING)" << endl;
        out << endl;
        // Java sources
        out << "    IF(VISIT_JAVA)" << endl;
        out << "        FILE(COPY " << atts->name<<".java " << "DESTINATION ${JavaClient_BINARY_DIR}/src/plots)" << endl;
        out << "        ADD_CUSTOM_TARGET(Java"<<name<<" ALL ${Java_JAVAC_EXECUTABLE} ${VISIT_Java_FLAGS} -d ${JavaClient_BINARY_DIR} -classpath ${JavaClient_BINARY_DIR} -sourcepath ${JavaClient_BINARY_DIR} ";
        if(customjfiles)
        {
            for(size_t i = 0; i < jfiles.size(); ++i)
                out << jfiles[i] << " ";
        }
        out << atts->name<<".java" << endl;
        out << "            DEPENDS JavaClient" << endl;
        out << "            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})" << endl;
        out << "    ENDIF(VISIT_JAVA)" << endl;

        out << "ENDIF(NOT VISIT_SERVER_COMPONENTS_ONLY AND NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)" << endl;
        out << endl;

        out << "ADD_LIBRARY(E"<<name<<"Plot_ser ${LIBE_SOURCES})" << endl;
        out << "TARGET_LINK_LIBRARIES(E"<<name<<"Plot_ser visitcommon avtplotter_ser avtpipeline_ser " << ToString(libs) << ToString(elibsSer) << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "E", "Plot_ser", "");
        out << "SET(INSTALLTARGETS ${INSTALLTARGETS} E"<<name<<"Plot_ser)" << endl;
        out << "ADD_TARGET_DEFINITIONS(E"<<name<<"Plot_ser ENGINE)" << endl;
        out << endl;
        out << "IF(VISIT_PARALLEL)" << endl;
        out << "    ADD_PARALLEL_LIBRARY(E"<<name<<"Plot_par ${LIBE_SOURCES})" << endl;
        out << "    TARGET_LINK_LIBRARIES(E"<<name<<"Plot_par visitcommon avtplotter_par avtpipeline_par " << ToString(libs) << ToString(elibsPar) << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "E", "Plot_par", "    ");
        out << "    SET(INSTALLTARGETS ${INSTALLTARGETS} E"<<name<<"Plot_par)" << endl;
        out << "    ADD_TARGET_DEFINITIONS(E"<<name<<"Plot_par ENGINE)" << endl;
        out << "ENDIF(VISIT_PARALLEL)" << endl;
        out << endl;
        out << "VISIT_INSTALL_PLOT_PLUGINS(${INSTALLTARGETS})" << endl;
        out << "VISIT_PLUGIN_TARGET_RTOD(plots ${INSTALLTARGETS})" << endl;
        if (using_dev)
          out << "VISIT_PLUGIN_TARGET_FOLDER(plots " << name  
              << " ${INSTALLTARGETS})" << endl;
        out << endl;
    }

    void WriteCMake_Operator(QTextStream &out, 
                             const QString guilibname, 
                             const QString viewerlibname)
    {
        bool useFortran = false;

        out << "PROJECT(" << name<< ")" << endl;
        out << endl;
        if (using_dev)
        {
        out << "INCLUDE(${VISIT_SOURCE_DIR}/CMake/PluginMacros.cmake)" <<endl;
        out << endl;
        }
        out << "SET(COMMON_SOURCES" << endl;
        out << name << "PluginInfo.C" << endl;
        out << name << "CommonPluginInfo.C" << endl;
        out << atts->name << ".C" << endl;
        out << ")" << endl;
        out << endl;
        out << "SET(LIBI_SOURCES " << endl;
        out << name << "PluginInfo.C" << endl;
        out << ")" << endl;
        out << endl;
        WriteCMake_ConditionalSources(out, "I", "");

        // libG sources
        out << "SET(LIBG_SOURCES" << endl;
        out << name << "GUIPluginInfo.C" << endl;
        out << "Qvis" << name << "Window.C" << endl;
        out << "${COMMON_SOURCES}" << endl;
        if (customgfiles)
        {
            useFortran |= CustomFilesUseFortran(gfiles);
            for (size_t i=0; i<gfiles.size(); i++)
                out << gfiles[i] << endl;
        }
        else
            for (size_t i=0; i<defaultgfiles.size(); i++)
                out << defaultgfiles[i] << endl;
        out << ")" << endl;
        WriteCMake_ConditionalSources(out, "G", "");
        out << "SET(LIBG_MOC_SOURCES" << endl;
        out << "Qvis" << name << "Window.h" << endl;
        if (customwfiles)
            for (size_t i=0; i<wfiles.size(); i++)
                out << wfiles[i] << endl;
        out << ")" << endl;
        out << endl;

        // libV sources
        out << "SET(LIBV_SOURCES" << endl;
        out << name<<"ViewerPluginInfo.C" << endl;
        if (customvfiles)
        {
            useFortran |= CustomFilesUseFortran(vfiles);
            for (size_t i=0; i<vfiles.size(); i++)
                out << vfiles[i] << endl;
        }
        else
            for (size_t i=0; i<defaultvfiles.size(); i++)
                out << defaultvfiles[i] << endl;
        out << "${COMMON_SOURCES}" << endl;
        out << ")" << endl;
        WriteCMake_ConditionalSources(out, "V", "");
        if (customvwfiles)
        {
            out << "SET(LIBV_MOC_SOURCES" << endl;
            for (size_t i=0; i<vwfiles.size(); i++)
                out << vwfiles[i] << endl;
            out << ")" << endl;
        }
        out << endl;

        // libE sources
        out << "SET(LIBE_SOURCES" << endl;
        out << name<<"EnginePluginInfo.C" << endl;
        if (customefiles)
        {
            useFortran |= CustomFilesUseFortran(efiles);
            for (size_t i=0; i<efiles.size(); i++)
                out << efiles[i] << endl;
        }
        else
            for (size_t i=0; i<defaultefiles.size(); i++)
                out << defaultefiles[i] << endl;
        out << "${COMMON_SOURCES}" << endl;
        out << ")" << endl;
        out << endl;
        WriteCMake_ConditionalSources(out, "E", "");

        if(useFortran)
        {
            out << "ENABLE_LANGUAGE(Fortran)" << endl;
        }

        WriteCMake_PlotOperator_Includes(out, true);
        WriteCMake_ConditionalIncludes(out);

        // Pass other CXXFLAGS
        for (size_t i=0; i<cxxflags.size(); i++)
        {
            if(!cxxflags[i].startsWith("${") && !cxxflags[i].startsWith("$("))
                 out << "ADD_DEFINITIONS(\"" << cxxflags[i] << "\")" << endl;
        }
        out << endl;

        WriteCMake_ConditionalDefinitions(out);

#if 0
        if (installpublic)
            out << "SET(LIBRARY_OUTPUT_PATH " << visitplugdirpub << ")" << endl;
        else if (installprivate)
            out << "SET(LIBRARY_OUTPUT_PATH " << visitplugdirpri << ")" << endl;
        else
#endif

        out << endl;
        // Extract extra link directories from LDFLAGS if they have ${},$(),-L
        std::vector<QString> linkDirs;
        linkDirs.push_back("${VISIT_LIBRARY_DIR}");
        linkDirs.push_back("${QT_LIBRARY_DIR}");
        linkDirs.push_back("${GLEW_LIBRARY_DIR}");
        linkDirs.push_back("${EAVL_LIBRARY_DIR}");
        linkDirs.push_back("${VTK_LIBRARY_DIRS}");
        for (size_t i=0; i<ldflags.size(); i++)
        {
            if(ldflags[i].startsWith("${") || ldflags[i].startsWith("$("))
                 linkDirs.push_back(ldflags[i]);
            else if(ldflags[i].startsWith("-L"))
                 linkDirs.push_back(ldflags[i].right(ldflags[i].size()-2));
        }

        if (HasCondition("LinkDirectories:"))
        {
            WriteCMake_ConditionalLinkDirs(out, linkDirs);
        }
        else
        {
            out << "LINK_DIRECTORIES(" << ToString(linkDirs) << ")" << endl;
        }
        out << endl;
        out << "ADD_LIBRARY(I"<<name<<"Operator ${LIBI_SOURCES})" << endl;
        out << "TARGET_LINK_LIBRARIES(I"<<name<<"Operator visitcommon)" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "I", "Operator", "");
        out << "SET(INSTALLTARGETS I"<<name<<"Operator)" << endl;
        out << endl;

        out << "IF(NOT VISIT_SERVER_COMPONENTS_ONLY AND NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)" << endl;
        out << "    QT_WRAP_CPP(G"<<name<<"Operator LIBG_SOURCES ${LIBG_MOC_SOURCES})" << endl;
        out << "    ADD_LIBRARY(G"<<name<<"Operator ${LIBG_SOURCES})" << endl;
        out << "    TARGET_LINK_LIBRARIES(G" << name << "Operator visitcommon "
            << guilibname << " " << ToString(libs) << ToString(glibs) 
            << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "G", "Operator", "    ");
        out << endl;
        if (customvwfiles)
            out << "    QT_WRAP_CPP(V"<<name<<"Operator LIBV_SOURCES ${LIBV_MOC_SOURCES})" << endl;
        out << "    ADD_LIBRARY(V"<<name<<"Operator ${LIBV_SOURCES})" << endl;
        out << "    TARGET_LINK_LIBRARIES(V" << name << "Operator visitcommon "
            << viewerlibname << " " << ToString(libs) << ToString(vlibs) 
            << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "V", "Operator", "    ");
        out << "    SET(INSTALLTARGETS ${INSTALLTARGETS} G"<<name<<"Operator V"<<name<<"Operator)" << endl;
        out << endl;
        // libS sources
        out << "    IF(VISIT_PYTHON_SCRIPTING)" << endl;
        out << "        SET(LIBS_SOURCES" << endl;
        out << "            " << name<<"ScriptingPluginInfo.C" << endl;
        out << "            Py"<<atts->name<<".C" << endl;
        if (customsfiles)
            for (size_t i=0; i<sfiles.size(); i++)
                out << "            " << sfiles[i] << endl;
        else
            for (size_t i=0; i<defaultsfiles.size(); i++)
                out << "            " << defaultsfiles[i] << endl;
        out << "            ${COMMON_SOURCES}" << endl;
        out << "        )" << endl;
        WriteCMake_ConditionalSources(out, "S", "        ");
        out << "        ADD_LIBRARY(S"<<name<<"Operator ${LIBS_SOURCES})" << endl;
        out << "        TARGET_LINK_LIBRARIES(S"<<name<<"Operator visitcommon visitpy ${PYTHON_LIBRARY})" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "S", "Operator", "        ");
        out << "        SET(INSTALLTARGETS ${INSTALLTARGETS} S"<<name<<"Operator)" << endl;
        out << "    ENDIF(VISIT_PYTHON_SCRIPTING)" << endl;
        out << endl;
        // Java sources
        out << "    IF(VISIT_JAVA)" << endl;
        out << "        FILE(COPY " << atts->name<<".java DESTINATION ${JavaClient_BINARY_DIR}/src/operators)" << endl;
        out << "        ADD_CUSTOM_TARGET(Java"<<name<<" ALL ${Java_JAVAC_EXECUTABLE} ${VISIT_Java_FLAGS} -d ${JavaClient_BINARY_DIR} -classpath ${JavaClient_BINARY_DIR} -sourcepath ${JavaClient_BINARY_DIR} ";
        if(customjfiles)
        {
            for(size_t i = 0; i < jfiles.size(); ++i)
                out << jfiles[i] << " ";
        }
        out << atts->name<<".java" << endl;
        out << "            DEPENDS JavaClient" << endl;
        out << "            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})" << endl;
        out << "    ENDIF(VISIT_JAVA)" << endl;

        out << "ENDIF(NOT VISIT_SERVER_COMPONENTS_ONLY AND NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)" << endl;
        out << endl;

        out << "ADD_LIBRARY(E"<<name<<"Operator_ser ${LIBE_SOURCES})" << endl;
        out << "TARGET_LINK_LIBRARIES(E"<<name<<"Operator_ser visitcommon avtexpressions_ser avtfilters_ser avtpipeline_ser " << ToString(libs) << ToString(elibsSer) << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "E", "Operator_ser", "");
        out << "SET(INSTALLTARGETS ${INSTALLTARGETS} E"<<name<<"Operator_ser)" << endl;
        if (hasEngineSpecificCode)
            out << "ADD_TARGET_DEFINITIONS(E"<<name<<"Operator_ser ENGINE)" << endl;
        out << endl;
        out << "IF(VISIT_PARALLEL)" << endl;
        out << "    ADD_PARALLEL_LIBRARY(E"<<name<<"Operator_par ${LIBE_SOURCES})" << endl;
        out << "    TARGET_LINK_LIBRARIES(E"<<name<<"Operator_par visitcommon avtexpressions_par avtfilters_par avtpipeline_par " << ToString(libs) << ToString(elibsPar) << ")" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "E", "Operator_par", "    ");
        out << "    SET(INSTALLTARGETS ${INSTALLTARGETS} E"<<name<<"Operator_par)" << endl;
        if (hasEngineSpecificCode)
            out << "    ADD_TARGET_DEFINITIONS(E"<<name<<"Operator_par ENGINE)" << endl;
        out << "ENDIF(VISIT_PARALLEL)" << endl;
        out << endl;
        out << "VISIT_INSTALL_OPERATOR_PLUGINS(${INSTALLTARGETS})" << endl;
        out << "VISIT_PLUGIN_TARGET_RTOD(operators ${INSTALLTARGETS})" << endl;
        if (using_dev)
          out << "VISIT_PLUGIN_TARGET_FOLDER(operators " << name 
              << " ${INSTALLTARGETS})" << endl;
        out << endl;
    }

    void WriteCMake_Database(QTextStream &out)
    {
        bool useFortran = false;

        out << "PROJECT("<<name<<")" << endl;
        out << endl;
        if (using_dev)
        {
        out << "INCLUDE(${VISIT_SOURCE_DIR}/CMake/PluginMacros.cmake)" <<endl;
        out << endl;
        }
        out << "SET(COMMON_SOURCES" << endl;
        out << ""<<name<<"PluginInfo.C" << endl;
        out << ""<<name<<"CommonPluginInfo.C" << endl;
        out << ")" << endl;
        out << endl;
        out << "SET(LIBI_SOURCES " << endl;
        out << ""<<name<<"PluginInfo.C" << endl;
        out << ")" << endl;
        out << endl;
        WriteCMake_ConditionalSources(out, "I", "");
        if(!onlyEnginePlugin)
        {
            out << "SET(LIBM_SOURCES" << endl;
            out << ""<<name<<"MDServerPluginInfo.C" << endl;
            out << "${COMMON_SOURCES}" << endl;
            if (custommfiles)
            {
                useFortran |= CustomFilesUseFortran(mfiles);
                for (size_t i=0; i<mfiles.size(); i++)
                    out << mfiles[i] << endl;
            }
            else
                for (size_t i=0; i<defaultmfiles.size(); i++)
                    out << defaultmfiles[i] << endl;
            out << ")" << endl;
            WriteCMake_ConditionalSources(out, "M", "");
            if (customwmfiles)
            {
                useFortran |= CustomFilesUseFortran(wmfiles);
                out << "IF(WIN32)" << endl;
                out << "    SET(LIBM_WIN32_SOURCES" << endl;
                for (size_t i=0; i<wmfiles.size(); i++)
                    out << "    " << wmfiles[i] << endl;
                out << "    )" << endl;
                for (size_t i=0; i<wmfiles.size(); i++)
                {
                    if(wmfiles[i].endsWith(".c"))
                    {
                        out << "    SET_SOURCE_FILES_PROPERTIES("
                            << wmfiles[i] << endl;
                        out << "        PROPERTIES LANGUAGE CXX)" << endl;
                    }
                }
                out << "ENDIF(WIN32)" << endl;
            }
            out << endl;
        }
        if(!noEnginePlugin)
        {
            out << "SET(LIBE_SOURCES" << endl;
            out <<name<<"EnginePluginInfo.C" << endl;
            out << "${COMMON_SOURCES}" << endl;
            if (customefiles)
            {
                useFortran |= CustomFilesUseFortran(efiles);
                for (size_t i=0; i<efiles.size(); i++)
                    out << efiles[i] << endl;
            }
            else
                for (size_t i=0; i<defaultefiles.size(); i++)
                    out << defaultefiles[i] << endl;
            out << ")" << endl;
            WriteCMake_ConditionalSources(out, "E", "");
            if (customwefiles)
            {
                useFortran |= CustomFilesUseFortran(wefiles);
                out << "IF(WIN32)" << endl;
                out << "    SET(LIBE_WIN32_SOURCES" << endl;
                for (size_t i=0; i<wefiles.size(); i++)
                    out << "    " << wefiles[i] << endl;
                out << "    )" << endl;
                for (size_t i=0; i<wefiles.size(); i++)
                {
                    if(wefiles[i].endsWith(".c"))
                    {
                        out << "    SET_SOURCE_FILES_PROPERTIES("
                            << wefiles[i] << endl;
                        out << "        PROPERTIES LANGUAGE CXX)" << endl;
                    }
                }
                out << "ENDIF(WIN32)" << endl;
            }
            out << endl;
        }

        // take any ${} from the CXXFLAGS to mean a variable that contains 
        // include directories.
        std::vector<QString> extraIncludes;
        for (size_t i=0; i<cxxflags.size(); i++)
        {
            if(cxxflags[i].startsWith("${"))
                 extraIncludes.push_back(ConvertToProperVisItIncludeDir(cxxflags[i]));
            else if(cxxflags[i].startsWith("$("))
                 extraIncludes.push_back(ConvertToProperVisItIncludeDir(ConvertDollarParenthesis(cxxflags[i])));
            else if(cxxflags[i].startsWith("-I"))
                 extraIncludes.push_back(ConvertToProperVisItIncludeDir(cxxflags[i].right(cxxflags[i].size()-2)));
        }
        out << "INCLUDE_DIRECTORIES(" << endl;
        out << "${CMAKE_CURRENT_SOURCE_DIR}" << endl;
        if(extraIncludes.size() > 0)
            out << ToString(extraIncludes, true) ;
        out << "${VISIT_COMMON_INCLUDES}" << endl;
        out << VisItIncludeDir() << "/avt/DBAtts/MetaData" << endl;
        out << VisItIncludeDir() << "/avt/DBAtts/SIL" << endl;
        out << VisItIncludeDir() << "/avt/Database/Database" << endl;
        out << VisItIncludeDir() << "/avt/Database/Formats" << endl;
        out << VisItIncludeDir() << "/avt/Database/Ghost" << endl;
        out << VisItIncludeDir() << "/avt/FileWriter" << endl;
        out << VisItIncludeDir() << "/avt/Filters" << endl;
        out << VisItIncludeDir() << "/avt/MIR/Base" << endl;
        out << VisItIncludeDir() << "/avt/MIR/Tet" << endl;
        out << VisItIncludeDir() << "/avt/MIR/Zoo" << endl;
        out << VisItIncludeDir() << "/avt/Math" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/AbstractFilters" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Data" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Pipeline" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Sinks" << endl;
        out << VisItIncludeDir() << "/avt/Pipeline/Sources" << endl;
        out << VisItIncludeDir() << "/avt/VisWindow/VisWindow" << endl;
        out << VisItIncludeDir() << "/visit_vtk/full" << endl;
        out << VisItIncludeDir() << "/visit_vtk/lightweight" << endl;
        out << "${EAVL_INCLUDE_DIR} " << endl;
        out << "${VTK_INCLUDE_DIRS} " << endl;
        out << ")" << endl;
        out << endl;

        WriteCMake_ConditionalIncludes(out);

        // Pass other CXXFLAGS
        for (size_t i=0; i<cxxflags.size(); i++)
        {
            if(!cxxflags[i].startsWith("${") &&
               !cxxflags[i].startsWith("$(") &&
               !cxxflags[i].startsWith("-I"))
                 out << "ADD_DEFINITIONS(\"" << cxxflags[i] << "\")" << endl;
        }
        bool needWindowsDefines = false;
        for (size_t i=0; i<libs.size() && !needWindowsDefines; i++)
        {
            if(libs[i].contains("BOXLIB"))
                 needWindowsDefines = true;
            else if(libs[i].contains("HDF4"))
                 needWindowsDefines = true;
            else if(libs[i].contains("FITS"))
                 needWindowsDefines = true;
            else if(libs[i].contains("NETCDF"))
                 needWindowsDefines = true;
            else if(libs[i].contains("CGNS"))
                 needWindowsDefines = true;
        }
        if (needWindowsDefines)
        {
            out << "IF(WIN32)" << endl;
            bool netcdfAdded = false;
            for (size_t i=0; i<libs.size(); i++)
            {
                if(libs[i].contains("BOXLIB"))
                     out << "  ADD_DEFINITIONS(-DBL_FORT_USE_UPPERCASE)" << endl;
                else if(libs[i].contains("HDF4"))
                     out << "  ADD_DEFINITIONS(-D_HDFDLL_ -D_MFHDFLIB_ -D_HDFLIB_ -DINTEL86)" << endl;
                else if(libs[i].contains("FITS"))
                     out << "  ADD_DEFINITIONS(-D_HDF5USEDLL_)" << endl;
                else if(libs[i].contains("NETCDF")&& !netcdfAdded)
                {
                     out << "  ADD_DEFINITIONS(-DDLL_NETCDF)" << endl;
                     netcdfAdded = true;
                }
                else if(libs[i].contains("CGNS"))
                     out << "  ADD_DEFINITIONS(-DUSE_DLL)" << endl;
            }
            out << "ENDIF(WIN32)" << endl;
        }
        WriteCMake_ConditionalDefinitions(out);

        if(useFortran)
        {
            out << "ENABLE_LANGUAGE(Fortran)" << endl;
        }

        out << endl;
        // Extract extra link directories from LDFLAGS if they have ${},$(),-L
        std::vector<QString> linkDirs;
        linkDirs.push_back("${VISIT_LIBRARY_DIR}");
        linkDirs.push_back("${EAVL_LIBRARY_DIR}");
        linkDirs.push_back("${VTK_LIBRARY_DIRS}");
        for (size_t i=0; i<ldflags.size(); i++)
        {
            if(ldflags[i].startsWith("${") || ldflags[i].startsWith("$("))
                 linkDirs.push_back(ldflags[i]);
            else if(ldflags[i].startsWith("-L"))
                 linkDirs.push_back(ldflags[i].right(ldflags[i].size()-2));
        }
        if (HasCondition("LinkDirectories:"))
        {
            WriteCMake_ConditionalLinkDirs(out, linkDirs);
        }
        else
        {
            out << "LINK_DIRECTORIES(" << ToString(linkDirs) << ")" << endl;
        }
        out << endl;
        out << "ADD_LIBRARY(I"<<name<<"Database ${LIBI_SOURCES})" << endl;
        out << "TARGET_LINK_LIBRARIES(I"<<name<<"Database visitcommon)" << endl;
        WriteCMake_ConditionalTargetLinks(out, name, "I", "Database", "");
        out << "SET(INSTALLTARGETS I"<<name<<"Database)" << endl;
        out << endl;
        if(!onlyEnginePlugin)
        {
            out << "IF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)" << endl;

            out << "    ADD_LIBRARY(M"<<name<<"Database ${LIBM_SOURCES}";
            if (customwmfiles)
                out << "     ${LIBM_WIN32_SOURCES}";
            out << "    )" << endl;
            out << "    TARGET_LINK_LIBRARIES(M"<<name<<"Database visitcommon avtdbatts avtdatabase_ser " << ToString(libs) << ToString(mlibs) << ")" << endl;
            WriteCMake_ConditionalTargetLinks(out, name, "M", "Database", "    ");
            out << "    ADD_TARGET_DEFINITIONS(M"<<name<<"Database MDSERVER)" << endl;
            out << "    SET(INSTALLTARGETS ${INSTALLTARGETS} M"<<name<<"Database)" << endl;
            out << "ENDIF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)" << endl;
            out << endl;
        }
        if(!noEnginePlugin)
        {
            out << "ADD_LIBRARY(E"<<name<<"Database_ser ${LIBE_SOURCES}";
            if (customwefiles)
                out << " ${LIBE_WIN32_SOURCES}";
            out << ")" << endl;
            out << "TARGET_LINK_LIBRARIES(E"<<name<<"Database_ser visitcommon avtdatabase_ser avtpipeline_ser " << ToString(libs) << ToString(elibsSer) << ")" << endl;
            WriteCMake_ConditionalTargetLinks(out, name, "E", "Database_ser", "");
            out << "ADD_TARGET_DEFINITIONS(E"<<name<<"Database_ser ENGINE)" << endl;
            out << "SET(INSTALLTARGETS ${INSTALLTARGETS} E"<<name<<"Database_ser)" << endl;
            out << endl;
            out << "IF(VISIT_PARALLEL)" << endl;
            out << "    ADD_PARALLEL_LIBRARY(E"<<name<<"Database_par ${LIBE_SOURCES})" << endl;
            out << "    TARGET_LINK_LIBRARIES(E"<<name<<"Database_par visitcommon avtdatabase_par avtpipeline_par " << ToString(libs) << ToString(elibsPar) << ")" << endl;
            WriteCMake_ConditionalTargetLinks(out, name, "E", "Database_par", "    ");
            out << "    ADD_TARGET_DEFINITIONS(E"<<name<<"Database_par ENGINE)" << endl;
            out << "    SET(INSTALLTARGETS ${INSTALLTARGETS} E"<<name<<"Database_par)" << endl;
            out << "ENDIF(VISIT_PARALLEL)" << endl;
            out << endl;
        }
        out << "VISIT_INSTALL_DATABASE_PLUGINS(${INSTALLTARGETS})" << endl;
        out << "VISIT_PLUGIN_TARGET_RTOD(databases ${INSTALLTARGETS})" << endl;
        if (using_dev)
          out << "VISIT_PLUGIN_TARGET_FOLDER(databases " << name 
              << " ${INSTALLTARGETS})" << endl;
        out << endl;
    }

    void WriteCMake(QTextStream &out)
    {
        const char *visithome = getenv("VISITARCHHOME");
        if (!visithome && !using_dev)
            throw QString().sprintf("Please set the VISITARCHHOME "
                                    "environment variable.\n"
                                    "You may have it set automatically "
                                    "using 'visit -xml2cmake'.");

        const char *visitplugdirpub = getenv("VISITPLUGININSTPUB");
        if (!visitplugdirpub && installpublic)
            throw QString().sprintf("Please set the VISITPLUGININSTPUB "
                                    "environment variable.\n"
                                    "You may have it set automatically "
                                    "using 'visit -xml2cmake'.");

        const char *visitplugdirpri = getenv("VISITPLUGININSTPRI");
        if (!visitplugdirpri)
        {
           if ((using_dev && installprivate) || !using_dev)
            throw QString().sprintf("Please set the VISITPLUGININSTPRI "
                                    "environment variable.\n"
                                    "You may have it set automatically "
                                    "using 'visit -xml2cmake'.");
        }

        out << "# DO NOT EDIT THIS FILE! THIS FILE IS AUTOMATICALLY GENERATED "
            << "BY xml2cmake" << endl;

        QString qvisithome(visithome);
        QString qvisitplugdirpub(visitplugdirpub);
        QString qvisitplugdirpri(visitplugdirpri);
#ifdef _WIN32
        qvisithome       = ToCMakePath(qvisithome);
        qvisitplugdirpub = ToCMakePath(qvisitplugdirpub);
        qvisitplugdirpri = ToCMakePath(qvisitplugdirpri);
#endif
        // If we're not using a development version then we need to always 
        // include something in the generated output.
        if(!using_dev)
        {
            out << "CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8 FATAL_ERROR)" << endl;
            out << "SET(VISIT_INCLUDE_DIR \"" << qvisithome 
                << "/include\")" << endl;
            out << "SET(VISIT_LIBRARY_DIR \"" << qvisithome 
                << "/lib\")" << endl;
#ifdef _WIN32
            // There is no 'bin' dir for installed VisIt on Windows
            out << "SET(VISIT_BINARY_DIR \""  << qvisithome 
                << "\")" << endl;
            out << "SET(VISIT_ARCHIVE_DIR \"" << qvisithome 
                << "/lib\")" << endl;
#else
            out << "SET(VISIT_BINARY_DIR \""  << qvisithome 
                << "/bin\")" << endl;
            out << "SET(VISIT_ARCHIVE_DIR \"" << qvisithome 
                << "/archives\")" << endl;
#endif
            if(installpublic)
            {
                out << "SET(VISIT_PLUGIN_DIR \"" << qvisitplugdirpub 
                    << "\")" << endl;
            }
            else // installprivate or default
            {
                out << "SET(VISIT_PLUGIN_DIR \"" << qvisitplugdirpri 
                    << "\")" << endl;
            }

            out << "INCLUDE(\"" << qvisithome 
                << "/include/PluginVsInstall.cmake\")" << endl;
            out << "INCLUDE(\"" << qvisithome 
                << "/include/VisItLibraryDependencies.cmake\")" << endl;
            out << endl;
        }
        else
        {
            // We're using a development version but we're installing public 
            // or private.
            if(installpublic)
            {
               out << "SET(VISIT_PLUGIN_DIR " << qvisitplugdirpub << ")" << endl;
            }

            if(installprivate)
            {
               out << "SET(VISIT_PLUGIN_DIR " << qvisitplugdirpri << ")" << endl;
            }
        }

        QString guilibname("gui");
        QString viewerlibname("viewer");
#ifdef WIN32
        if (! using_dev)
        {
            // when calling from an installed version, cmake doesn't know that
            // the gui and viewer lib targets have been renamed to guilib and
            // viewer lib (to prevent conflicts with the exe targets), so they 
            // must be explictily listed by the name of the actual lib created.
            guilibname    = "guilib";
            viewerlibname = "viewerlib";
        }
#endif
        if(type == "plot")
            WriteCMake_Plot(out, guilibname, viewerlibname);
        else if(type == "operator")
            WriteCMake_Operator(out, guilibname, viewerlibname);
        else if(type == "database")
            WriteCMake_Database(out);

        WriteCMake_AdditionalCode(out);
    }
};


// ----------------------------------------------------------------------------
//                           Override default types
// ----------------------------------------------------------------------------
#define Plugin       CMakeGeneratorPlugin

#endif
