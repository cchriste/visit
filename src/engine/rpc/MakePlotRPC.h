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

#ifndef MAKE_PLOT_RPC_H
#define MAKE_PLOT_RPC_H
#include <engine_rpc_exports.h>

#include <VisItRPC.h>
#include <string>
#include <vector>

class MakePlotRPC;

// ****************************************************************************
//  Class:  PreparePlotRPC
//
//  Purpose:
//    Signals the name of the plot about to be created so that
//    the MakePlotRPC has space to store the correct attributes.
//
//  Programmer:  Jeremy Meredith
//  Creation:    March  2, 2001
//
//  Modifications:
//    Brad Whitlock, Wed Mar 21 22:48:12 PST 2007
//    Added name.
//
// ****************************************************************************
class ENGINE_RPC_API PreparePlotRPC : public BlockingRPC
{
  public:
    PreparePlotRPC();
    ~PreparePlotRPC();

    virtual const std::string TypeName() const { return "PreparePlotRPC"; }

    void SetMakePlotRPC(MakePlotRPC*);
    MakePlotRPC *GetMakePlotRPC();

    void operator()(const std::string &plotName, const std::string &pluginID);
    void SelectAll();
    std::string GetID();
    std::string GetName();
  private:
    MakePlotRPC *makePlotRPC;
    std::string id;
    std::string name;
};


// ****************************************************************************
//  Class:  MakePlotRPC
//
//  Purpose:
//    Make a plot.
//
//  Programmer:  Jeremy Meredith
//  Creation:    March  2, 2001
//
//  Modifications:
//    Jeremy Meredith, Fri Nov  9 10:16:51 PST 2001
//    Added a NetworkID as a return type of the rpc.
//
//    Brad Whitlock, Mon Mar 25 09:56:24 PDT 2002
//    Removed SetSocket.
//
//    Jeremy Meredith, Fri Mar 19 14:56:55 PST 2004
//    Modified the rpc to pass data extents.
//
//    Mark C. Miller, Wed Dec 15 21:56:26 PST 2004
//    Added windowID
//
//    Brad Whitlock, Wed Mar 21 22:46:44 PST 2007
//    Added plotName arg to operator().
//
// ****************************************************************************

class ENGINE_RPC_API MakePlotRPC : public BlockingRPC
{
  public:
    struct ENGINE_RPC_API NetworkID : public AttributeSubject
    {
        int id;
      public:
        NetworkID();
        NetworkID(int id_);
        virtual ~NetworkID();
        virtual void SelectAll();

        virtual const std::string TypeName() const
            { return "MakePlotRPC::NetworkID"; }
    };

  public:
    MakePlotRPC();
    virtual ~MakePlotRPC();

    virtual const std::string TypeName() const { return "MakePlotRPC"; }

    int operator()(const std::string &, const std::string&,
                   const AttributeSubject*,
                   const std::vector<double> &, int);

    virtual void SelectAll();

    std::string GetName();
    std::string GetID();
    AttributeSubject *GetAtts();
    const std::vector<double> &GetDataExtents() const;
    PreparePlotRPC &GetPreparePlotRPC();
    int GetWindowID() const;

    void SetAtts(AttributeSubject*);
    void SetDataExtents(const std::vector<double> &);
    void SetWindowID(int id);

    virtual void SetXfer(Xfer *x);

  private:
    AttributeSubject *atts;
    std::vector<double> dataExtents;
    PreparePlotRPC preparePlotRPC;
    NetworkID networkID;
    int       windowID;
};

#endif
