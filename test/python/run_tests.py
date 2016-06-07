"""
Quick Python unit testing
"""

import unittest
import HTMLTestRunner

loader = unittest.TestLoader()

tests = loader.discover('.', pattern="Test*.py")

report_file = file("python_test_report.html", 'wb')
runner = HTMLTestRunner.HTMLTestRunner(stream=report_file, title = "Python Unit Tests")
runner.run(tests)