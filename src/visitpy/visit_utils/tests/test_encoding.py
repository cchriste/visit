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
 file: test_encoding.py
 author: Cyrus Harrison (cyrush@llnl.gov)
 created: 4/09/2010
 description:
    Unit tests for movie encoding helpers.

"""

import unittest
import os
import sys
import glob

from visit_utils import encoding
from visit_utils.common import VisItException

from os.path import join as pjoin

iframes_dir = pjoin(os.path.split(__file__)[0],"_data")
iframes = pjoin(iframes_dir,"wave.movie.%04d.png")
iframes_short_a = pjoin(iframes_dir,"wave.movie.%03d.png")
iframes_short_b = pjoin(iframes_dir,"wave.movie.%d.png")
iframes_stereo  = pjoin(iframes_dir,"wave.movie.stereo.%04d.png")
output_dir  = pjoin(os.path.split(__file__)[0],"_output")

def lst_slnks():
    return glob.glob(pjoin(iframes_dir,"_encode.lnk.*"))

def clean_slnks():
    slnks = lst_slnks()
    for slnk in slnks:
        os.remove(slnk)

def check_encoded_file(path):
    if os.path.isfile(path):
        # make sure the file isn't empty
        st = os.stat(path)
        return st.st_size > 0
    return False

class TestEncoding(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
    def test_encoders(self):
        encoders = encoding.encoders()
        if len(encoders) > 0:
            self.assertTrue("mpg" in encoders)
            self.assertTrue("wmv" in encoders)
    def test_avi(self):
        if "avi" in encoding.encoders():
            ofile = pjoin(output_dir,"wave.movie.avi")
            encoding.encode(iframes,ofile)
            self.assertTrue(check_encoded_file(ofile))
    def test_mpg(self):
        if "mpg" in encoding.encoders():
            ofile = pjoin(output_dir,"wave.movie.mpg")
            encoding.encode(iframes,ofile)
            self.assertTrue(check_encoded_file(ofile))
    def test_mpg_slow(self):
        if "mpg" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.slow.mpg")
            encoding.encode(iframes,ofile,2)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
    def test_mov(self):
        if "mov" in encoding.encoders():
            ofile = pjoin(output_dir,"wave.movie.mov")
            encoding.encode(iframes,ofile)
            self.assertTrue(check_encoded_file(ofile))
    def test_mov_slow(self):
        if "mov" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.slow.mov")
            encoding.encode(iframes,ofile,2)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
    def test_wmv(self):
        if "wmv" in encoding.encoders():
            ofile = pjoin(output_dir,"wave.movie.wmv")
            encoding.encode(iframes,ofile)
            self.assertTrue(check_encoded_file(ofile))
    def test_wmv_slow(self):
        if "wmv" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.slow.wmv")
            encoding.encode(iframes,ofile,2)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
    def test_sm(self):
        if "sm" in encoding.encoders():
            ofile = pjoin(output_dir,"wave.movie.sm")
            encoding.encode(iframes,ofile)
            self.assertTrue(check_encoded_file(ofile))
    def test_sm_slow(self):
        if "sm" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.slow.sm")
            encoding.encode(iframes,ofile,2)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
    def test_swf(self):
        if "swf" in encoding.encoders():
            ofile = pjoin(output_dir,"wave.movie.swf")
            encoding.encode(iframes,ofile )
            self.assertTrue(check_encoded_file(ofile))
    def test_unsupported(self):
            self.assertRaises(VisItException, encoding.encode, iframes,"wave.movie.bad_ext")
    def test_sm_stereo(self):
        if "sm" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.stereo.sm")
            encoding.encode(iframes_stereo,ofile,stereo=True)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
    def test_sm_stereo_slow(self):
        if "sm" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.stereo.slow.sm")
            encoding.encode(iframes_stereo,ofile,2,stereo=True)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
    def test_stereo_uneven_frames_error(self):
            self.assertRaises(VisItException, encoding.encode, iframes,
                                                                pjoin(output_dir,"wave.movie.stereo.bad.sm"),
                                                                stereo=True)
    def test_extract(self):
        if "mpg" in encoding.encoders():
            eframes = pjoin(output_dir,"extract_out_%0m4d.png")
            encoding.encode(iframes,pjoin(output_dir,"wave.movie.mpg"))
            encoding.extract(pjoin(output_dir,"wave.movie.mpg"),eframes)
            ofile = pjoin(output_dir,"wave.movie.extract.and.reencode.mpg")
            encoding.encode(eframes,ofile)
            self.assertTrue(check_encoded_file(ofile))
    def test_pre_lr_stereo(self):
        if "divx" in encoding.encoders():
            iframes = pjoin(iframes_dir,"noise.stereo.left.right.1080p.%04d.png")
            ofile = pjoin(output_dir,"noise.movie.stereo.pre.left.right.avi")
            encoding.encode(iframes,ofile,etype="divx")
            self.assertTrue(check_encoded_file(ofile))
    def test_short_symlinks(self):
        if "mpg" in encoding.encoders():
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.test.seq.pattern.03d.mpg")
            encoding.encode(iframes_short_a,ofile,3)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))
            clean_slnks()
            ofile = pjoin(output_dir,"wave.movie.test.seq.pattern.d.mpg")
            encoding.encode(iframes_short_b,ofile,5)
            self.assertEqual(0,len(lst_slnks()))
            self.assertTrue(check_encoded_file(ofile))

if __name__ == '__main__':
    unittest.main()


