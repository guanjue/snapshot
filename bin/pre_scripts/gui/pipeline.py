from sys import argv
import sys
import time

from PyQt5.QtCore import (QCoreApplication, QObject, pyqtSignal, QThread)
from snapshot_gui import easy_snapshot

# Here we use pipeline, instead of directly thread
# it is easier to control the process.
class Pipeline(QObject):

    finished = pyqtSignal()

    def __init__(self):
        super(Pipeline, self).__init__()
        self.count_threshold = 1
        self.siglog2 = 'F'
        self.sigscale = 'F'
        self.sigsmallnum = 0.01
        self.function_method = 'mostfreq'
        self.input_dir = ''
        self.output_dir = ''
        self.bin_dir = ''
        print('finish init')

    def set_parameters(self, _count_threshold, _sigsmallnum, _siglog2, _sigscale, _function_method):
        self.count_threshold = _count_threshold
        self.siglog2 = _siglog2
        self.sigsmallnum = _sigsmallnum
        self.sigscale = _sigscale
        self.function_method = _function_method
        print('parameters')

    def set_dir(self, _input_dir, _output_dir, _bin_dir):
        self.input_dir = _input_dir
        self.output_dir = _output_dir
        self.bin_dir = _bin_dir
        print('dir set:')

    def long_running(self):
        print('start long running')
        easy_snapshot(self.count_threshold, self.siglog2, self.sigscale, self.sigsmallnum, self.function_method, self.input_dir, self.output_dir, self.bin_dir)
        #count = 0
        #while count < 5:
        #    time.sleep(1)
        #    print("B Increasing")
        #    count += 1
        #print('my name is', argv[1])
        
        self.finished.emit()

class BackgroundThread(QThread):

    def run(self):
        count = 0
        while count < 10:
            time.sleep(1)
            print("A Increasing")
            count += 1

