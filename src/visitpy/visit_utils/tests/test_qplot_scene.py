#!/usr/bin/env python
#
# file: test_basic.py
# author: cdh
#
# Tests for qannote module.
#
#

import unittest
from visit_utils import *

import unittest
import os
from os.path import join as pjoin

from visit_test import *
from visit_utils.qplot import *

try:
    import PySide.QtCore
except:
    pass

output_dir  = pjoin(os.path.dirname(__file__),"_output")
data_dir    = pjoin(os.path.dirname(__file__),"_data")

def patch_scene_input(in_fname,ult_fname):
    #crv_file = os.path.abspath(pjoin(data_dir,"pattern.ult"))
    #qi = open(pjoin(data_dir,"qplot.example.in")).read().replace("$SOURCE_FILE",crv_file)
    #qi_fname = pjoin(output_dir,"qplot.example.in"
    #open(qi_fname,"w").write(qi)
    crv_file = os.path.abspath(pjoin(data_dir,ult_fname))
    qi = open(pjoin(data_dir,in_fname)).read().replace("$SOURCE_FILE",crv_file)
    qi_fname = pjoin(output_dir,in_fname)
    open(qi_fname,"w").write(qi)
    return qi_fname
        

class TestScene(unittest.TestCase):
    def setUp(self):
        self.curves = ult.Curve.load(pjoin(data_dir,"sine.ult"))
    def test_01_basic(self):
        crv = ult.Curve.load(pjoin(data_dir,"pattern.ult"))[0]
        self.assertEqual(crv.interp(7.5),7.5)
        self.assertEqual(crv.interp(5),5)
        self.assertEqual(crv.interp(0),5)
        self.assertEqual(crv.interp(20),20)
        self.assertEqual(crv.interp(30),20)
    @pyside_test
    def test_02_basic(self):
        in_fname = patch_scene_input("qplot.example.in","pattern.ult")
        params = common.load_params(in_fname)
        scene = CurveScene(params.scene)
        opng = pjoin(output_dir,"test.scene.00.png")
        scene.render(opng)
        self.assertTrue(os.path.isfile(opng))
    @pyside_test
    def test_03_basic(self):
        in_fname = patch_scene_input("qplot.example.log.in","pattern.log.ult")
        params = common.load_params(in_fname)
        scene = CurveScene(params.scene)
        opng = pjoin(output_dir,"test.scene.01.png")
        scene.render(opng)
        self.assertTrue(os.path.isfile(opng))
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()

