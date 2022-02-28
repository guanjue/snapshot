# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, pyqtSignal, QThread
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from color import ColorWindow
from para import ParameterWindow
from visual import VisualWindow
from pipeline import Pipeline
import os
import time

class thread(QThread):
    change_progress = pyqtSignal(int)
 
    def __init__(self):
        super(thread, self).__init__()
        self.start_value = 0
        self.end_value = 100
        self.progress = self.start_value
 
    def run(self):
        time.sleep(0.2)
        for i in range(0, self.end_value):
            self.change_progress.emit(i)
            #time.sleep(0.01)
            
            if i < self.end_value/2:
                time.sleep(0.2)
            elif i >= self.end_value/2 and i < self.end_value * 2 / 3:
                time.sleep(1)
            else:
                time.sleep(5)
            

class faster_thread(QThread):
    change_progress = pyqtSignal(int)
    can_show_pict = pyqtSignal()
 
    def __init__(self):
        super(faster_thread, self).__init__()
        self.start_value = 0
        self.end_value = 101
        self.progress = self.start_value
    
    def set_start_value(self, _start_value):
        self.start_value = _start_value

    def run(self):
        for i in range(self.start_value, self.end_value, 1):
            self.change_progress.emit(i)
            time.sleep(0.05)
        self.can_show_pict.emit()
        


class StartWindow(QtWidgets.QMainWindow):

    # add signal
    changeColor = pyqtSignal()
    changeParameter = pyqtSignal()
    startVisualization = pyqtSignal()

    finish_progress_bar = pyqtSignal()


    def __init__(self):
        super(QtWidgets.QMainWindow, self).__init__()
        # define member
        self.file1 = ""
        self.file2 = ""
        self.file3 = ""
        self.file4 = ""
        self.file5 = ""
        self.para_win = None
        self.color_win = None
        self.visual_win = None

        self.pipeline_thread = None
        self.pipeline = None

        self.bg_thread = None

        # pipeline parameters
        self.count_threshold = 1
        self.siglog2 = 'F'
        self.sigscale = 'F'
        self.sigsmallnum = 0.01
        self.function_method = 'mostfreq'

        self.completed = 0

        #self.input_dir = '/Users/sun/git/snapshot/test_data/input_data/'
        #self.output_dir = '/Users/sun/data/output/'
        self.input_dir = ''
        self.output_dir = ''
        
        self.bin_dir = (os.path.dirname(os.path.realpath(__file__))) + '/../'

        self.setupUi()
        self.setup_signals()

        self.progress_thread = thread()

        self.faster_thread = faster_thread()
        
        self.progress_thread.change_progress.connect(self.set_progress)
        self.faster_thread.change_progress.connect(self.set_progress)
        self.faster_thread.can_show_pict.connect(self.show_picture)


        

    def setupUi(self):
        self.setObjectName("MainWindow")
        self.resize(800, 440)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.button_change_color = QtWidgets.QPushButton(self.centralwidget)
        self.button_change_color.setObjectName("button_change_color")
        self.horizontalLayout.addWidget(self.button_change_color)
        self.button_change_parameter = QtWidgets.QPushButton(self.centralwidget)
        self.button_change_parameter.setObjectName("button_change_parameter")
        self.horizontalLayout.addWidget(self.button_change_parameter)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.line_edit_output_dir = QtWidgets.QLineEdit(self.centralwidget)
        self.line_edit_output_dir.setObjectName("line_edit_output_dir")
        self.gridLayout.addWidget(self.line_edit_output_dir, 4, 2, 1, 1)
        self.button_output_dir = QtWidgets.QPushButton(self.centralwidget)
        self.button_output_dir.setObjectName("button_output_dir")
        self.gridLayout.addWidget(self.button_output_dir, 4, 3, 1, 1)
        self.button_input_dir = QtWidgets.QPushButton(self.centralwidget)
        self.button_input_dir.setObjectName("button_input_dir")
        self.gridLayout.addWidget(self.button_input_dir, 3, 3, 1, 1)
        self.line_edit_input_dir = QtWidgets.QLineEdit(self.centralwidget)
        self.line_edit_input_dir.setObjectName("line_edit_input_dir")
        self.gridLayout.addWidget(self.line_edit_input_dir, 3, 2, 1, 1)
        self.label_file_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_file_2.setObjectName("label_file_2")
        self.gridLayout.addWidget(self.label_file_2, 4, 0, 1, 1)
        self.label_file_1 = QtWidgets.QLabel(self.centralwidget)
        self.label_file_1.setObjectName("label_file_1")
        self.gridLayout.addWidget(self.label_file_1, 3, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem3)
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem4)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem5)
        spacerItem6 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem6)
        spacerItem7 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem7)
        spacerItem8 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem8)
        self.progress_bar = QtWidgets.QProgressBar(self.centralwidget)
        #self.progress_bar.setProperty("Processing", 0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setObjectName("progress_bar")
        self.verticalLayout.addWidget(self.progress_bar)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem9)
        self.button_preprocess = QtWidgets.QPushButton(self.centralwidget)
        self.button_preprocess.setObjectName("button_preprocess")
        self.horizontalLayout_2.addWidget(self.button_preprocess)
        self.button_ok = QtWidgets.QPushButton(self.centralwidget)
        self.button_ok.setObjectName("button_ok")
        self.horizontalLayout_2.addWidget(self.button_ok)
        self.button_cancel = QtWidgets.QPushButton(self.centralwidget)
        self.button_cancel.setObjectName("button_cancel")
        self.horizontalLayout_2.addWidget(self.button_cancel)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        spacerItem10 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem10)
        self.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        self.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(self)
        self.statusbar.setObjectName("statusbar")
        self.setStatusBar(self.statusbar)

        self.retranslateUi()
        #QtCore.QMetaObject.connectSlotsByName(self)

        self.line_edit_input_dir.setEnabled(False)
        self.line_edit_output_dir.setEnabled(False)
        self.progress_bar.setVisible(False)
        #self.progress_bar.setEnabled(False)
        self.button_change_color.setEnabled(False)
        self.button_change_parameter.setEnabled(False)
        self.button_ok.setEnabled(False)

        #self.set_enabled_all_buttons(True)

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.button_change_color.setText(_translate("MainWindow", "Change Display Color"))
        self.button_change_parameter.setText(_translate("MainWindow", "Change Parameters"))
        self.button_output_dir.setText(_translate("MainWindow", "Choose Directory"))
        self.button_input_dir.setText(_translate("MainWindow", "Choose Directory"))
        self.label_file_2.setText(_translate("MainWindow", "Output Directory"))
        self.label_file_1.setText(_translate("MainWindow", "Input Directory"))
        self.button_preprocess.setText(_translate("MainWindow", "Preprocessing"))
        self.button_ok.setText(_translate("MainWindow", "OK"))
        self.button_cancel.setText(_translate("MainWindow", "Exit"))

    def set_enabled_all_buttons(self, true_false):
        self.button_input_dir.setEnabled(true_false)
        self.button_output_dir.setEnabled(true_false)
        self.button_change_color.setEnabled(true_false)
        self.button_change_parameter.setEnabled(true_false)
        self.button_preprocess.setEnabled(true_false)
        self.button_ok.setEnabled(true_false)

    def setup_signals(self):
        #pass
        self.button_change_color.clicked.connect(self.on_button_change_color_clicked)
        self.button_change_parameter.clicked.connect(self.on_button_change_parameter_clicked)
        self.button_input_dir.clicked.connect(self.on_button_input_dir_clicked)
        self.button_output_dir.clicked.connect(self.on_button_output_dir_clicked)
        self.button_preprocess.clicked.connect(self.on_button_preprocess_clicked)
        self.button_ok.clicked.connect(self.on_button_ok_clicked)
        self.button_cancel.clicked.connect(self.on_button_cancel_clicked)
        
        #self.progress_bar.setVisible(False)

    def on_button_change_color_clicked(self):
        # self.changeColor.emit()
        if not self.color_win:
            self.color_win = ColorWindow()
            self.color_win.setupUi(self.input_dir + '/function_color_list.txt')
            self.color_win.setup_signals()
            self.color_win.color_change_signal.connect(self.on_color_changed)
            self.color_win.exec_()
        else:
            self.color_win.exec_()

    def on_button_change_parameter_clicked(self):
        # self.changeParameter.emit()
        if not self.para_win:
            self.para_win = ParameterWindow()
            self.para_win.para_change_signal.connect(self.on_parameter_changed)
            #self.para_win.show()
            self.para_win.exec_()
        else:
            self.para_win.exec_()

    def on_button_input_dir_clicked(self):
        chosen_dir = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if chosen_dir == '':
            return
        self.input_dir = chosen_dir + '/'
        self.line_edit_input_dir.setText(self.input_dir)
        self.set_enabled_all_buttons(False)
        self.button_input_dir.setEnabled(True)
        self.button_preprocess.setEnabled(True)
        self.button_ok.setEnabled(False)
        self.button_change_color.setEnabled(False)
        self.button_change_parameter.setEnabled(False)
        self.button_output_dir.setEnabled(True)

    def on_button_output_dir_clicked(self):
        chosen_dir = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if chosen_dir == '':
            return
        self.output_dir = chosen_dir + '/'
        self.line_edit_output_dir.setText(self.output_dir)

    def open_filename(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            return fileName
        else:
            return "can not find file"

    def on_button_preprocess_clicked(self):
        if self.input_dir == '' and self.output_dir == '':
            QMessageBox.warning(self, 'Warning', 'Please choose input and output directory!')
        elif self.input_dir == '':
            QMessageBox.warning(self, 'Warning', 'Please choose input directory!')
        elif self.output_dir == '':
            QMessageBox.warning(self, 'Warning', 'Please choose output directory!')
        else:
            #self.progress_bar.setVisible(True)
            self.set_enabled_all_buttons(True)
            self.button_preprocess.setEnabled(False)

    def on_button_ok_clicked(self):

        self.button_ok.setEnabled(False)
        
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)

        self.progress_thread.start()
        

        self.pipeline_thread = QThread()
        self.pipeline = Pipeline()
        
        self.pipeline.set_parameters(self.count_threshold, self.sigsmallnum, self.siglog2, self.sigscale, self.function_method)
        self.pipeline.set_dir(self.input_dir, self.output_dir, self.bin_dir)

        self.pipeline.moveToThread(self.pipeline_thread)
        self.pipeline.finished.connect(self.pipeline_thread.quit)
        self.pipeline_thread.started.connect(self.pipeline.long_running)
        self.pipeline_thread.finished.connect(self.on_picture_ready)
        print('thread init finished.')
        self.pipeline_thread.start()
        
        #self.bg_thread = BackgroundThread()
        #self.bg_thread.finished.connect(self.on_picture_ready)
        #self.bg_thread.start()
        

    def on_button_cancel_clicked(self):
        buttonReply = QMessageBox.question(self, 'Exit Application', "Do you want to exit current application?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if buttonReply == QMessageBox.Yes:
            self.close()
        else:
            print('No clicked.')

        if self.pipeline_thread is not None:
            self.pipeline_thread.quit()

        #if self.bg_thread is not None:
        #    self.bg_thread.quit()
            

    def on_picture_ready(self):
        #self.startVisualization.emit()
        if self.pipeline_thread is not None:
            self.pipeline_thread = None

        if self.progress_thread.isRunning():
            self.progress_thread.terminate()
            self.progress_thread.wait()

        self.faster_thread.set_start_value(int(self.progress_bar.value())+1)
        self.faster_thread.start()

    def show_picture(self):

        if not self.visual_win:
            self.visual_win = VisualWindow(self.output_dir)
            #self.visual_win.para_change_signal.connect(self.on_parameter_changed)
            self.visual_win.show()
            #self.setFocus(True)
            #self.visual_win.setFocus(True)
        else:
            self.visual_win.show()

        self.button_ok.setEnabled(True)

        self.progress_bar.setVisible(False)
        image_size = QtCore.QSize(800, 600)
        leftPixelMap = QtGui.QPixmap(self.output_dir + "/snapshot.meansig.png")
        myScaledPixmap = leftPixelMap.scaled(image_size, QtCore.Qt.KeepAspectRatio)
        self.visual_win.label_image.setPixmap(myScaledPixmap)
        self.visual_win.label_image.setAlignment(QtCore.Qt.AlignCenter)

    def on_parameter_changed(self):
        if self.para_win is not None:
            self.count_threshold = int(self.para_win.line_edit_count_threshold.text())
            self.sigsmallnum = float(self.para_win.line_edit_small_number.text())
            if not self.para_win.checkBox_log2.isChecked():
                self.siglog2 = 'F'
            else:
                self.siglog2 = 'T'
            if not self.para_win.checkBox_scale.isChecked():
                self.sigscale = 'F'
            else:
                self.sigscale = 'T'

            if self.para_win.radioButton_frequent.isChecked():
                self.function_method = 'mostfreq'
            elif self.para_win.radioButton_mean.isChecked():
                self.function_method = 'mean'

            #print(self.count_threshold, self.sigsmallnum, self.siglog2, self.sigscale, self.function_method)
            self.para_win.para_change_signal.disconnect()
            self.para_win.close()
            self.para_win.deleteLater()
            self.para_win = None
            #self.updateLabelAnswer(value)

    def on_color_changed(self):
        if self.color_win is not None:
            print('color changed')
            self.color_win.color_change_signal.disconnect()
            self.color_win.close()
            self.color_win.deleteLater()
            self.color_win = None
 
    def set_progress(self,pnumber):
        self.progress_bar.setValue(pnumber)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    ui = StartWindow()
    #ui.setupUi(MainWindow)
    #MainWindow.show()
    ui.show()
    sys.exit(app.exec_())

"""
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
"""

