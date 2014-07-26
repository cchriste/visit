#*****************************************************************************
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
#*****************************************************************************
"""
 file: property_tree.py
 author: Cyrus Harrison <cyrush@llnl.gov>
 created: 4/15/2010
 description:
    (almost too) flexible property tree class.

"""


class PropertyTree(object):
    def __init__(self,ptype ="tree",init = None):
        self._type = ptype
        self._locked = False
        if self._type == "tree":
            self._value = {}
            if not init is None:
                self.update(init)
        else:
            self._value = init
    def update(self,pval):
        if isinstance(pval,dict):
            for path,value in  pval.items():
                self.add_property(path,value)
        else:
            for path,value in pval.properties().items():
                self.add_property(path,value)
    def clear(self):
        self._locked = False
        if self._type == "tree":
            self._value = {}
        else:
            self._value = None
    def properties(self):
        res = {}
        if self._type == "tree":
            keys = self._value.keys()
            for k in keys:
                curr = self._value[k]
                if curr._type=="tree":
                    for p,v in curr.properties().items():
                        res[k + "/" + p] = v
                else:
                    res[k] = curr._value
        return res
    def children(self):
        res = {}
        if self._type == "tree":
            keys = self._value.keys()
            keys.sort()
            for k in keys:
                curr = self._value[k]
                res[k] = curr
        return res
    def add_property(self,path,value=None):
        idx = path.find("/")
        if idx > 0:
            lpath = path[:idx]
            rpath = path[idx+1:]
            if not lpath in self._value.keys():
                tree = PropertyTree()
                self._value[lpath] = tree
            else:
                tree = self._value[lpath]
            tree.add_property(rpath,value)
        else:
            if value is None:
                self._value[path] = PropertyTree()
            else:
                self._value[path] = PropertyTree("node",value)
    def has_property(self,path):
        node = self.fetch_property(path)
        return not node is None
    def remove_property(self,path):
        # find the proper node in the tree and remove it
        idx = path.find("/")
        if idx > 0:
            lpath = path[:idx]
            rpath = path[idx+1:]
            tree = self._value[lpath]
            tree.remove_property(rpath)
        elif path in self._value.keys():
            del self._value[path]
    def lock(self):
        self._locked = True
        if self._type == "tree":
            for v in self._value.values():
                v.lock()
    def unlock(self):
        self._locked = False
        if self._type == "tree":
            for v in self._value.values():
                v.unlock()
    def __getitem__(self,path):
        node = self.fetch_property(path)
        if node is None:
            if self._locked:
                raise AttributeError(path)
            self.add_property(path)
        node = self.fetch_property(path)
        if node._type == "tree":
            return node
        else:
            return node._value
    def __setitem__(self,path,obj):
        node = self.fetch_property(path)
        if node is None:
            if self._locked:
                raise AttributeError(path)
            self.add_property(path,obj)
        else:
            node._type  = "node"
            node._value = obj
    def fetch_property(self,path):
        idx = path.find("/")
        if idx > 0:
            lpath = path[:idx]
            if lpath in self._value.keys():
                rpath = path[idx+1:]
                tree = self._value[lpath]
                return tree.fetch_property(rpath)
            return None
        elif path in self._value.keys():
            return self._value[path]
        else:
            return None
    def __str__(self):
        return self.__gen_string("")
    def __gen_string(self,path):
        res = ""
        if self._type == "tree":
            for k in self._value.keys():
               npath = path + k + "/"
               res +=  self._value[k].__gen_string(npath)
        else:
            res = path + "%s:%s\n" % (self._type,str(self._value))
        return res
    def __getattr__(self, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        try:
            return self.__getitem__(name)
        except KeyError:
                raise AttributeError(name)
    def __setattr__(self, name,obj):
        if name == "_value" or name == "_type" or name == "_locked":
            self.__dict__[name] = obj
        else:
            self.__setitem__(name, obj)



