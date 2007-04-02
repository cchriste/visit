c-----------------------------------------------------------------
c Program: main
c
c-----------------------------------------------------------------
      program main
      implicit none
      include "visitfortransiminterface.inc"
ccc   local variables
      integer err

      err = visitsetupenv()
      err = visitinitializesim("funstructured", 13,
     . "Demonstrates creating an unstructured mesh", 42,
     . "/no/useful/path", 15,
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
      include "visitfortransiminterface.inc"
ccc   local variables
      integer visitstate, result, blocking
ccc   SIMSTATE common block
      integer runflag, simcycle
      real simtime
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
      include "visitfortransiminterface.inc" 
ccc   SIMSTATE common block
      integer runFlag, simcycle
      real simtime
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

c Advance time
      simcycle = simcycle + 1
      simtime = simtime + 0.0134
      write (6,*) 'Simulating time step: cycle=',simcycle,
     . ' time=', simtime
      call sleep(1)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c These functions must be defined to satisfy the visitfortransiminterface lib.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c---------------------------------------------------------------------------
c visitcommandcallback
c---------------------------------------------------------------------------
      subroutine visitcommandcallback (cmd, lcmd, intdata, 
     .                                 floatdata, stringdata, 
     .                                 lstringdata)
      implicit none
      character*8 cmd, stringdata
      integer     lcmd, lstringdata, intdata
      real        floatdata
      include "visitfortransiminterface.inc"
ccc   SIMSTATE common block
      integer runflag, simcycle
      real simtime
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
c visitgetmetadata
c---------------------------------------------------------------------------
      integer function visitgetmetadata(handle)
      implicit none
      integer handle
      include "visitfortransiminterface.inc"
ccc   SIMSTATE common block
      integer runflag, simcycle
      real simtime
      common /SIMSTATE/ runflag, simcycle, simtime
      integer err, tdim, sdim, mesh, mt, scalar, curve, mat, e

      err = visitmdsetcycletime(handle, simcycle, simtime)
      if(runflag.eq.1) then
          err = visitmdsetrunning(handle, VISIT_SIMMODE_RUNNING)
      else
          err = visitmdsetrunning(handle, VISIT_SIMMODE_STOPPED)
      endif

c     Add a 3D unstructured mesh
      mt = VISIT_MESHTYPE_UNSTRUCTURED
      tdim = 3
      sdim = 3
      mesh = visitmdmeshcreate(handle, "unstructured3d", 14, mt, tdim,
     . sdim, 1)
      if(mesh.ne.VISIT_INVALID_HANDLE) then
          err = visitmdmeshsetunits(handle, mesh, "cm", 2)
          err = visitmdmeshsetlabels(handle, mesh, "Width", 5,
     .    "Height", 6, "Depth", 5)
          err = visitmdmeshsetblocktitle(handle, mesh, "Domains", 7)
          err = visitmdmeshsetblockpiecename(handle, mesh, "domain", 6)
      endif

c     Add simulation commands
      err = visitmdaddsimcommand(handle, "halt", 4, VISIT_CMDARG_NONE,
     .                           1)
      err = visitmdaddsimcommand(handle, "step", 4, VISIT_CMDARG_NONE,
     .                           1)
      err = visitmdaddsimcommand(handle, "run", 3, VISIT_CMDARG_NONE,
     .                           1)

      visitgetmetadata = VISIT_OKAY
      end

c---------------------------------------------------------------------------
c visitgetmesh
c---------------------------------------------------------------------------
      integer function visitgetmesh(handle, domain, name, lname)
      implicit none
      character*8 name
      integer     handle, domain, lname
      include "visitfortransiminterface.inc" 
ccc   UNSTRUCTURED common block (shared with simulate_one_timestep)
      integer NNODES, NZONES, LCONN
      parameter (NNODES = 16)
      parameter (NZONES = 5)
      parameter (LCONN = 36)
      real umx(NNODES), umy(NNODES), umz(NNODES)
      integer connectivity(LCONN)
      common /UNSTRUCTURED/ umx, umy, umz, connectivity
ccc   local variables
      integer m

      m = VISIT_ERROR
      if(visitstrcmp(name, lname, "unstructured3d", 14).eq.0) then
c Create an unstructured mesh here
          m = visitmeshunstructured(handle, 3, NNODES, NZONES, 0,
     . NZONES-1, umx, umy, umz, LCONN, connectivity)
      endif
      visitgetmesh = m
      end

c---------------------------------------------------------------------------
c visitgetscalar
c---------------------------------------------------------------------------
      integer function visitgetscalar(handle, domain, name, lname)
      implicit none
      character*8 name
      integer     handle, domain, lname
      include "visitfortransiminterface.inc"
      visitgetscalar = VISIT_ERROR
      end


c---------------------------------------------------------------------------
c visitgetcurve
c---------------------------------------------------------------------------
      integer function visitgetcurve(handle, name, lname)
      implicit none
      character*8 name
      integer     handle, lname
      include "visitfortransiminterface.inc"
      visitgetcurve = VISIT_ERROR
      end

c---------------------------------------------------------------------------
c visitgetdomainlist
c---------------------------------------------------------------------------
      integer function visitgetdomainlist(handle)
      implicit none
      integer handle
      include "visitfortransiminterface.inc"
      visitgetdomainlist = VISIT_OKAY
      end

c---------------------------------------------------------------------------
c visitgetmaterial
c---------------------------------------------------------------------------
      integer function visitgetmaterial(handle, domain, name, lname)
      implicit none
      character*8 name
      integer     handle, domain, lname
      include "visitfortransiminterface.inc"
      visitgetmaterial = VISIT_ERROR
      end





























