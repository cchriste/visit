#!/usr/bin/env python
#
# file: setup_tests.py
# author: cdh
#
# Provides 'ExecuteTests' a distutils command class that automatically
# collects and launches unittest test cases from scripts in the 'tests'
# subdir.
#
# I created this simple one file module b/c nose & setuptools are not always
# available/installable on all of the systems I develop python modules on.
#
# To use this class in a distutils setup project:
#  A) Create test scripts in the subdir 'tests'
#  B) In your setup.py:
#     #import the test module
#     import setup_tests
#     # Add the custom command to the 'cmdclass' dictonary in your setup call:
#     setup( ...
#            cmdclass = { 'test': setup_tests.ExecuteTests})
#
#  C) Invoke 'python setup.py test' to run tests
#

from distutils.core import Command
from unittest import TextTestRunner, TestLoader
import os, glob, sys

class ExecuteTests(Command):
    description = "locate and execute scripts containing unit tests."
    user_options = [  ('verbose', None, "Verbose output"),
                      ('tests-dir=', None, "Directory containg tests")]
    boolean_options = ['verbose']
    negative_opt = {'quiet' : 'verbose'}
    def initialize_options(self):
        self.verbose = 1
        self.tests_dir = "tests"
    def finalize_options(self):
        self.tests_dir = os.path.join(os.getcwd(),self.tests_dir)
    def run(self):
        """
        Finds and executes unit tests in the 'tests' subdir.
        Because TestLoader imports the tests as a module this method
        automatically creates/updates the 'tests/__init__.py' to
        import all python scripts in the 'tests' subdir.
        """
        self.run_command('build')
        sys.path.insert(0,os.path.join(os.getcwd(),"build","lib"))
        self.tests  = []
        # make sure the 'tests' subdir actually exists.
        if not os.path.isdir(self.tests_dir):
            print "ExecuteTests: <Error> 'tests' subdir not found!"
        else:
            self.find_tests()
            self.gen_tests_init()
            # create a test suite.
            tests = TestLoader().loadTestsFromNames([t[0] for t in self.tests])
            # run the test suite if it actually contains test cases.
            run_verbosity = 2
            if self.verbose == 0:
                run_verbosity = 0
            if tests.countTestCases() > 0:
                runner = TextTestRunner(verbosity=run_verbosity)
                runner.run(tests)
            else:
                print "ExecuteTests: <Warning> No test cases found!"
        sys.path.pop(0)
    def find_tests(self):
        """
        Helper called by 'run' to find python scripts in the 'tests' subdir.
        """
        for f in glob.glob(os.path.join(self.tests_dir,"*.py")):
            if not f.endswith('__init__.py'):
                test = os.path.splitext(os.path.basename(f))[0]
                test = '.'.join(['tests',test])
                self.tests.append( (test,f))
        if len(self.tests) > 0 and self.verbose:
            print "Detected Test Scripts:"
            for t in self.tests:
                print " %s @ %s)" % (t[0],t[1])
    def gen_tests_init(self):
        """
        Helper called by 'run' to update the 'tests/__init__.py' file.

        Note: If no python scripts exist in the 'tests' subdir
              'tests/__init__.py' will be removed.
        """
        tests_init = os.path.join(self.tests_dir,"__init__.py")
        if len(self.tests) > 0:
            f = open(tests_init,"w")
            f.write("# module header auto generated by "
                    "setup_tests.ExecuteTests\n")
            # import each test script as part of the 'tests' module.
            for t in self.tests:
                f.write("import %s\n" % t[0])
        elif os.path.exists(tests_init):
            # remove the module init file if no test scripts are found.
            os.remove(tests_init)