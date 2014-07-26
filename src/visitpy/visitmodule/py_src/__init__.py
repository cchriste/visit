###############################################################################
#
# Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
# Produced at the Lawrence Livermore National Laboratory
# LLNL-CODE-442911
# All rights reserved.
#
# This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
# full copyright notice is contained in the file COPYRIGHT located at the root
# of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
#
# Redistribution  and  use  in  source  and  binary  forms,  with  or  without
# modification, are permitted provided that the following conditions are met:
#
#  - Redistributions of  source code must  retain the above  copyright notice,
#    this list of conditions and the disclaimer below.
#  - Redistributions in binary form must reproduce the above copyright notice,
#    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
#    documentation and/or other materials provided with the distribution.
#  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
#    be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
# ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
# LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
# DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
# SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
# CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
# LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
# OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
###############################################################################
# file: __init__.py
# Purpose: Main init for the pyside_visithook module.
#
# Programmer: Cyrus Harrison
# Creation: Tue Apr  3 13:40:18 PDT
#
#
# Modifications:
#   Brad Whitlock, Tue May  1 16:49:40 PDT 2012
#   Add writescript.
#
#   Cyrus Harrison, Wed May  2 14:07:07 PDT 2012
#   Avoid creating pyc files when importing the visit module
#
###############################################################################
import sys
__system_bytecode_setting = sys.dont_write_bytecode
sys.dont_write_bytecode = True

#
# note: the frontend is responsible for loading the actual visit module
#

from frontend import *
from evalfuncs import *
from writescript import *

try:
    import PySide
    import pyside_hook
    import pyside_gui
    import pyside_support

    class KeyPressEater(PySide.QtCore.QObject):
        def eventFilter(self, obj, event):
            if event.type() == PySide.QtCore.QEvent.KeyPress:
                return True
            elif event.type() == PySide.QtCore.QEvent.MouseButtonPress:
                return True
            elif event.type() == PySide.QtCore.QEvent.MouseButtonRelease:
                return True
            elif event.type() == PySide.QtCore.QEvent.MouseButtonDblClick:
                return True
            elif event.type() == PySide.QtCore.QEvent.MouseMove:
                return True
            else:
                return PySide.QtCore.QObject.eventFilter(self, obj, event)

    def __VisIt_PySide_Idle_Hook__():
      a = KeyPressEater()
      app = PySide.QtCore.QEventLoop();
      PySide.QtCore.QCoreApplication.instance().installEventFilter(a)
      app.processEvents(PySide.QtCore.QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents);
      PySide.QtCore.QCoreApplication.instance().removeEventFilter(a)

    def IsPySideViewerEnabled(): return True
except ImportError:
    def IsPySideViewerEnabled(): return False
    pass

sys.dont_write_bytecode = __system_bytecode_setting










