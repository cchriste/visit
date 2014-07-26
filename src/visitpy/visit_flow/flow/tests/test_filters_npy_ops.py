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
 file: test_filters_npy_ops.py
 author: Cyrus Harrison <cyrush@llnl.gov>
 created: 3/24/2012
 description:
    unittest test cases for Filters in the npy_ops module.

"""

import unittest
try:
    import numpy as npy
except:
    pass

from flow import *
from flow.filters import npy_ops

from decorators import numpy_test

# uncomment for detailed exe info
#import logging
#logging.basicConfig(level=logging.INFO)


class TestNPyOps(unittest.TestCase):
    def setUp(self):
        print ""
    @numpy_test
    def test_01_simple_workspace(self):
        w = Workspace()
        w.register_filters(npy_ops)
        v_a = npy.array(range(10),dtype=npy.double)
        v_b = npy.array(range(10),dtype=npy.double)
        print ""
        w.registry_add(":src_a",v_a)
        w.registry_add(":src_b",v_b)
        w.add_filter("add","f1")
        w.add_filter("sub","f2")
        w.add_filter("mult","f3")
        w.add_filter("mult","f4")
        w.add_filter("add","f5")
        # f1 = src_a + src_b
        w.connect(":src_a","f1:in_a")
        w.connect(":src_b","f1:in_b")
        # f2 = src_b - src_a
        w.connect(":src_b","f2:in_a")
        w.connect(":src_a","f2:in_b")
        # f3 = f1^2
        w.connect("f1","f3:in_a")
        w.connect("f1","f3:in_b")
        # f4 = f2^2
        w.connect("f2","f4:in_a")
        w.connect("f2","f4:in_b")
        # f5 = f4 + f3
        w.connect("f3","f5:in_a")
        w.connect("f4","f5:in_b")
        print ""
        print w.graph
        print w.execution_plan()
        act_res = w.execute()
        # get output and test
        test_res = npy.power((v_a + v_b),2.0)+ npy.power((v_a - v_b),2.0)
        dsum = npy.sum(act_res - test_res)
        print "Filter Graph Result: %s" % str(act_res)
        print "Test Result:         %s" % str(test_res)
        print "Difference:          %s" % str(dsum)
        self.assertTrue(dsum < 1e-6)
    @numpy_test
    def test_02_more_ops(self):
        w = Workspace()
        w.register_filters(npy_ops)
        v_a = npy.array(range(10),dtype=npy.double)
        v_b = npy.array(range(10),dtype=npy.double)
        v_p = npy.array([2]*10,dtype=npy.double)
        w.registry_add(":src_a",v_a)
        w.registry_add(":src_b",v_b)
        w.registry_add(":src_p",v_p)
        w.add_filter("add","f1")
        w.add_filter("sub","f2")
        w.add_filter("pow","f3")
        w.add_filter("pow","f4")
        w.add_filter("add","f5")
        # f1 = src_a + src_b
        w.connect(":src_a","f1:in_a")
        w.connect(":src_b","f1:in_b")
        # f2 = src_b - src_a
        w.connect(":src_b","f2:in_a")
        w.connect(":src_a","f2:in_b")
        # f3 = f1^2
        w.connect("f1","f3:in_a")
        w.connect(":src_p","f3:in_b")
        # f4 = f2^2
        w.connect("f2","f4:in_a")
        w.connect(":src_p","f4:in_b")
        # f5 = f4 + f3
        w.connect("f3","f5:in_a")
        w.connect("f4","f5:in_b")
        print ""
        act_res = w.execute()
        # get output and test
        test_res = npy.power((v_a + v_b),2.0)+ npy.power((v_a - v_b),2.0)
        dsum = npy.sum(act_res - test_res)
        print "Filter Graph Result: %s" % str(act_res)
        print "Test Result:         %s" % str(test_res)
        print "Difference:          %s" % str(dsum)
        self.assertTrue(dsum < 1e-6)


if __name__ == '__main__':
    unittest.main()

