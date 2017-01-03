c-----------------------------------------------------------------------------
c
c Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-442911
c All rights reserved.
c
c This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
c full copyright notice is contained in the file COPYRIGHT located at the root
c of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
c
c Redistribution  and  use  in  source  and  binary  forms,  with  or  without
c modification, are permitted provided that the following conditions are met:
c
c  - Redistributions of  source code must  retain the above  copyright notice,
c    this list of conditions and the disclaimer below.
c  - Redistributions in binary form must reproduce the above copyright notice,
c    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
c    documentation and/or other materials provided with the distribution.
c  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
c    be used to endorse or promote products derived from this software without
c    specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
c ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
c LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
c DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
c DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
c SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
c CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
c LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
c OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
c DAMAGE.
c
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------
c Program: main
c
c Programmer: Brad Whitlock
c Date:       Thu Mar 11 12:13:06 PST 2010
c
c Modifications:
c
c-----------------------------------------------------------------
      program main
      implicit none
      include "visitfortransimV2interface.inc"
ccc   local variables
      integer err

      call simulationarguments()
      err = visitsetupenv()
      err = visitinitializesim("funstructured", 13,
     . "Demonstrates creating an unstructured mesh", 42,
     . "/no/useful/path", 15,
     . VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN,
     . VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN,
     . VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN)
      call mainloop()
      stop
      end

c-----------------------------------------------------------------
c mainloop
c-----------------------------------------------------------------
      subroutine mainloop()
      implicit none
      include "visitfortransimV2interface.inc"
ccc   local variables
      integer visitstate, result, blocking
ccc   SIMSTATE common block
      integer runflag, simcycle
      double precision simtime
      common /SIMSTATE/ runflag, simcycle, simtime
      save /SIMSTATE/

c     main loop
      runflag = 1
      simcycle = 0
      simtime = 0.
      do 10
          if(runflag.eq.1) then
              blocking = 0 
          else
              blocking = 1
          endif

          visitstate = visitdetectinput(blocking, -1)

          if (visitstate.lt.0) then
              goto 1234
          elseif (visitstate.eq.0) then
              call simulate_one_timestep()
          elseif (visitstate.eq.1) then
              runflag = 0
              result = visitattemptconnection()
              if (result.eq.1) then
                  write (6,*) 'VisIt connected!'
              else
                  write (6,*) 'VisIt did not connect!'
              endif
          elseif (visitstate.eq.2) then
              runflag = 0
              if (visitprocessenginecommand().eq.0) then
                  result = visitdisconnect()
                  runflag = 1
              endif
          endif
10    continue
1234  end

      subroutine simulate_one_timestep()
      implicit none
      include "visitfortransimV2interface.inc" 
ccc   SIMSTATE common block
      integer runFlag, simcycle
      double precision simtime
      common /SIMSTATE/ runflag, simcycle, simtime
ccc   UNSTRUCTURED common block (shared with visitgetmesh)
      integer  NNODES, NZONES, LCONN
      parameter (NNODES = 16)
      parameter (NZONES = 5)
      parameter (LCONN = 36)
      real umx(NNODES), umy(NNODES), umz(NNODES)
      integer connectivity(LCONN)
      common /UNSTRUCTURED/ umx, umy, umz, connectivity
      save /UNSTRUCTURED/
c Data values
      data umx/0.,2.,2.,0.,0.,2.,2.,0.,0.,2.,2.,0.,1.,2.,4.,4./
      data umy/0.,0.,0.,0.,2.,2.,2.,2.,4.,4.,4.,4.,6.,0.,0.,0./
      data umz/2.,2.,0.,0.,2.,2.,0.,0.,2.,2.,0.,0.,1.,4.,2.,0./
      data connectivity/VISIT_CELL_HEX, 0,1,2,3,4,5,6,7,
     . VISIT_CELL_HEX, 4,5,6,7,8,9,10,11,
     . VISIT_CELL_PYR, 8,9,10,11,12,
     . VISIT_CELL_WEDGE, 1,14,5,2,15,6,
     . VISIT_CELL_TET, 1,14,13,5/

c Simulate one time step
      simcycle = simcycle + 1
      simtime = simtime + 0.0134
      write (6,*) 'Simulating time step: cycle=',simcycle, 
     .            ' time=', simtime
      call sleep(1)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c These functions must be defined to satisfy the visitfortransimV2interface lib.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c---------------------------------------------------------------------------
c visitcommandcallback
c---------------------------------------------------------------------------
      subroutine visitcommandcallback (cmd, lcmd, args, largs) 
      implicit none
      character*8 cmd, args
      integer     lcmd, largs
      include "visitfortransimV2interface.inc"
ccc   SIMSTATE common block
      integer runflag, simcycle
      double precision simtime
      common /SIMSTATE/ runflag, simcycle, simtime

c     Handle the commands that we define in visitgetmetadata.
      if(visitstrcmp(cmd, lcmd, "halt", 4).eq.0) then
          runflag = 0
      elseif(visitstrcmp(cmd, lcmd, "step", 4).eq.0) then
          call simulate_one_timestep()
      elseif(visitstrcmp(cmd, lcmd, "run", 3).eq.0) then
          runflag = 1
      endif
      end

c---------------------------------------------------------------------------
c visitbroadcastintfunction
c---------------------------------------------------------------------------
      integer function visitbroadcastintfunction(value, sender)
      implicit none
      integer value, sender
c     REPLACE WITH MPI COMMUNICATION IF SIMULATION IS PARALLEL
      visitbroadcastintfunction = 0
      end

c---------------------------------------------------------------------------
c visitbroadcaststringfunction
c---------------------------------------------------------------------------
      integer function visitbroadcaststringfunction(str, lstr, sender)
      implicit none
      character*8 str
      integer     lstr, sender
c     REPLACE WITH MPI COMMUNICATION IF SIMULATION IS PARALLEL
      visitbroadcaststringfunction = 0
      end

c---------------------------------------------------------------------------
c visitslaveprocesscallback
c---------------------------------------------------------------------------
      subroutine visitslaveprocesscallback ()
      implicit none
c     REPLACE WITH MPI COMMUNICATION IF SIMULATION IS PARALLEL
      end

c---------------------------------------------------------------------------
c visitactivatetimestep
c---------------------------------------------------------------------------
      integer function visitactivatetimestep()
      implicit none
      include "visitfortransimV2interface.inc"
      visitactivatetimestep = VISIT_OKAY
      end

c---------------------------------------------------------------------------
c visitgetmetadata
c---------------------------------------------------------------------------
      integer function visitgetmetadata()
      implicit none
      include "visitfortransimV2interface.inc"
ccc   SIMSTATE common block
      integer runflag, simcycle
      double precision simtime
      common /SIMSTATE/ runflag, simcycle, simtime
      integer md, mmd, mat, cmd, err

      if(visitmdsimalloc(md).eq.VISIT_OKAY) then
          err = visitmdsimsetcycletime(md, simcycle, simtime)
          if(runflag.eq.1) then
              err = visitmdsimsetmode(md, VISIT_SIMMODE_RUNNING)
          else
              err = visitmdsimsetmode(md, VISIT_SIMMODE_STOPPED)
          endif

c     Add a 3D unstructured mesh
          err = visitmdmeshalloc(mmd);
          if(err.eq.VISIT_OKAY) then
              err = visitmdmeshsetname(mmd, "unstructured3d", 14)
              err = visitmdmeshsetmeshtype(mmd, 
     .            VISIT_MESHTYPE_UNSTRUCTURED)
              err = visitmdmeshsettopologicaldim(mmd, 3)
              err = visitmdmeshsetspatialdim(mmd, 3)
              err = visitmdmeshsetnumdomains(mmd, 1)
              err = visitmdmeshsetdomaintitle(mmd, "Domains", 7)
              err = visitmdmeshsetdomainpiecename(mmd, "domain", 6)
              err = visitmdmeshsetxunits(mmd, "cm", 2)
              err = visitmdmeshsetyunits(mmd, "cm", 2)
              err = visitmdmeshsetzunits(mmd, "cm", 2)
              err = visitmdmeshsetxlabel(mmd, "Width", 5)
              err = visitmdmeshsetylabel(mmd, "Height", 6)
              err = visitmdmeshsetzlabel(mmd, "Depth", 5)
              err = visitmdmeshsetcellorigin(mmd, 1)
              err = visitmdmeshsetnodeorigin(mmd, 1)
              err = visitmdsimaddmesh(md, mmd)
          endif

c     Add simulation commands
          err = visitmdcmdalloc(cmd)
          if(err.eq.VISIT_OKAY) then
              err = visitmdcmdsetname(cmd, "halt", 4)
              err = visitmdsimaddgenericcommand(md, cmd)
          endif
          err = visitmdcmdalloc(cmd)
          if(err.eq.VISIT_OKAY) then
              err = visitmdcmdsetname(cmd, "step", 4)
              err = visitmdsimaddgenericcommand(md, cmd)
          endif
          err = visitmdcmdalloc(cmd)
          if(err.eq.VISIT_OKAY) then
              err = visitmdcmdsetname(cmd, "run", 3)
              err = visitmdsimaddgenericcommand(md, cmd)
          endif
      endif
      visitgetmetadata = md
      end

c---------------------------------------------------------------------------
c visitgetmesh
c---------------------------------------------------------------------------
      integer function visitgetmesh(domain, name, lname)
      implicit none
      character*8 name
      integer     domain, lname
      include "visitfortransimV2interface.inc" 
ccc   UNSTRUCTURED common block (shared with simulate_one_timestep)
      integer NNODES, NZONES, LCONN
      parameter (NNODES = 16)
      parameter (NZONES = 5)
      parameter (LCONN = 36)
      real umx(NNODES), umy(NNODES), umz(NNODES)
      integer connectivity(LCONN)
      common /UNSTRUCTURED/ umx, umy, umz, connectivity
ccc   local variables
      integer h, x, y, z, conn, err

      h = VISIT_INVALID_HANDLE
      if(visitstrcmp(name, lname, "unstructured3d", 14).eq.0) then          
c Create an unstructured mesh here
          if(visitucdmeshalloc(h).eq.VISIT_OKAY) then
              err = visitvardataalloc(x)
              err = visitvardataalloc(y)
              err = visitvardataalloc(z)
              err = visitvardatasetf(x,VISIT_OWNER_SIM,1,NNODES,umx)
              err = visitvardatasetf(y,VISIT_OWNER_SIM,1,NNODES,umy)
              err = visitvardatasetf(z,VISIT_OWNER_SIM,1,NNODES,umz)

              err = visitvardataalloc(conn)
              err = visitvardataseti(conn,VISIT_OWNER_SIM,1,LCONN,
     .                               connectivity)

              err = visitucdmeshsetcoordsxyz(h, x, y, z)
              err = visitucdmeshsetconnectivity(h, NZONES, conn)
          endif
      endif
      visitgetmesh = h
      end

c---------------------------------------------------------------------------
c visitgetvariable
c---------------------------------------------------------------------------
      integer function visitgetvariable(domain, name, lname)
      implicit none
      character*8 name
      integer     domain, lname
      include "visitfortransimV2interface.inc"
      visitgetvariable = VISIT_INVALID_HANDLE
      end

c---------------------------------------------------------------------------
c visitgetmixedvariable
c---------------------------------------------------------------------------
      integer function visitgetmixedvariable(domain, name, lname)
      implicit none
      character*8 name
      integer     domain, lname
      include "visitfortransimV2interface.inc"
      visitgetmixedvariable = VISIT_INVALID_HANDLE
      end

c---------------------------------------------------------------------------
c visitgetcurve
c---------------------------------------------------------------------------
      integer function visitgetcurve(name, lname)
      implicit none
      character*8 name
      integer     lname
      include "visitfortransimV2interface.inc"
      visitgetcurve = VISIT_INVALID_HANDLE
      end

c---------------------------------------------------------------------------
c visitgetdomainlist
c---------------------------------------------------------------------------
      integer function visitgetdomainlist(name, lname)
      implicit none
      character*8 name
      integer     lname
      include "visitfortransimV2interface.inc"
      visitgetdomainlist = VISIT_INVALID_HANDLE
      end

c---------------------------------------------------------------------------
c visitgetdomainbounds
c---------------------------------------------------------------------------
      integer function visitgetdomainbounds(name, lname)
      implicit none
      character*8 name
      integer     lname
      include "visitfortransimV2interface.inc"
      visitgetdomainbounds = VISIT_INVALID_HANDLE
      end

c---------------------------------------------------------------------------
c visitgetdomainnesting
c---------------------------------------------------------------------------
      integer function visitgetdomainnesting(name, lname)
      implicit none
      character*8 name
      integer     lname
      include "visitfortransimV2interface.inc"
      visitgetdomainnesting = VISIT_INVALID_HANDLE
      end

c---------------------------------------------------------------------------
c visitgetmaterial
c---------------------------------------------------------------------------
      integer function visitgetmaterial(domain, name, lname)
      implicit none
      character*8 name
      integer     domain, lname
      include "visitfortransimV2interface.inc"
      visitgetmaterial = VISIT_INVALID_HANDLE
      end





























