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

#ifndef SIMULATION_COMMAND_RPC_H
#define SIMULATION_COMMAND_RPC_H

#include <engine_rpc_exports.h>
#include <VisItRPC.h>
#include <string>
#include <vector>

// ****************************************************************************
//  Class:  SimulationCommandRPC
//
//  Purpose:
//    Executes a simulation command
//
//  Programmer:  Jeremy Meredith
//  Creation:    March 18, 2005
//
//  Modifications:
//    Brad Whitlock, Fri Mar 27 11:29:39 PDT 2009
//    I simplified it so it only accepts string arguments. I also changed it
//    to nonblocking.
//
// ****************************************************************************

class ENGINE_RPC_API SimulationCommandRPC : public NonBlockingRPC
{
  public:
    SimulationCommandRPC();
    virtual ~SimulationCommandRPC();

    virtual const std::string TypeName() const { return "SimulationCommandRPC"; }

    // Invokation methods
    void operator()(const std::string &command, const std::string&);

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetCommand(const std::string&);
    void SetStringData(const std::string&);

    // Property getting methods
    std::string GetCommand()  const;
    std::string GetStringData() const;

  private:
    std::string command;
    std::string string_data;
};


#endif
