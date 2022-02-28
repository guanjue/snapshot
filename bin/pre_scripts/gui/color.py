# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'color.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot, pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QMessageBox
import functools
import random

class ColorWindow(QtWidgets.QDialog):
    color_change_signal = pyqtSignal()

    def __init__(self):
        super(QtWidgets.QDialog, self).__init__()
        
        self.labellist = []
        self.colorlist = []
        self.buttonlist = []
        self.color_str_list = []
        self.linelist = []
        self.color_profile_filename = ''
        self.color_id = 0

        # for debug
        #self.setupUi('/Users/sun/git/snapshot/test_data/input_data/function_color_list.txt')

        #self.setup_signals()

    def setupUi(self, _color_profile_filename):
        self.setObjectName("MainWindow")
        self.resize(300, 500)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")

        self.scrollArea = QtWidgets.QScrollArea(self.centralwidget)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 780, 508))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        """
        self.gridLayoutWidget = QtWidgets.QWidget(self.scrollAreaWidgetContents)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 10, 751, 461))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.pushButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 0, 2, 1, 1)
        """
        #self.gridLayoutWidget = QtWidgets.QWidget(self.scrollAreaWidgetContents)
        self.mygroupbox = QtWidgets.QGroupBox('')
        self.mygrid = QtWidgets.QGridLayout()

        self.color_profile_filename = _color_profile_filename
        self.color_id = 0
        with open(self.color_profile_filename) as color_profile:
            for line in color_profile:
                line = line.strip()
                self.linelist.append(line)
                cols = line.split()
                colors = cols[2].split(',')
                self.labellist.append(QtWidgets.QLabel( cols[0] ))
                self.colorlist.append(QtWidgets.QLabel())
                self.buttonlist.append(QtWidgets.QPushButton('Change Color'))
                self.mygrid.addWidget(self.labellist[self.color_id], self.color_id, 0)
                self.mygrid.addWidget(self.colorlist[self.color_id], self.color_id, 1)
                self.mygrid.addWidget(self.buttonlist[self.color_id], self.color_id, 2)
                self.buttonlist[self.color_id].clicked.connect(functools.partial(self.on_button, self.color_id))
                color_name = '#%02X%02X%02X' % (int(colors[0]),int(colors[1]),int(colors[2]))
                self.color_str_list.append(color_name)
                #assert(self.hex2rgb(color_name) == cols[2])
                self.colorlist[self.color_id].setStyleSheet("QWidget { background-color: %s}" % color_name)
                self.color_id += 1

        self.mygroupbox.setLayout(self.mygrid)
        #self.mygroupbox.resize(400,600)

        self.scrollArea.setWidget(self.mygroupbox)
        self.scrollArea.setWidgetResizable(True)
        #self.scrollArea.resize(400, 600)
        #self.scrollArea.setFixedHeight(400)

        self.verticalLayout.addWidget(self.scrollArea)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.button_ok = QtWidgets.QPushButton(self.centralwidget)
        self.button_ok.setObjectName("button_ok")
        self.horizontalLayout.addWidget(self.button_ok)
        self.button_reset = QtWidgets.QPushButton(self.centralwidget)
        self.button_reset.setObjectName("button_reset")
        self.horizontalLayout.addWidget(self.button_reset)
        self.verticalLayout.addLayout(self.horizontalLayout)
        
        #self.setCentralWidget(self.centralwidget)
        #self.menubar = QtWidgets.QMenuBar(self)
        #self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        #self.menubar.setObjectName("menubar")
        #self.setMenuBar(self.menubar)
        
        #self.statusbar = QtWidgets.QStatusBar()
        #self.statusbar.setObjectName("statusbar")
        #self.setStatusBar(self.statusbar)

        self.retranslateUi()
        

    def setup_signals(self):
        self.button_ok.clicked.connect(self.on_button_ok_clicked)
        self.button_reset.clicked.connect(self.on_button_reset_clicked)

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Color Panel"))
        #self.label.setText(_translate("MainWindow", "Color Label"))
        #self.label_2.setText(_translate("MainWindow", "Name Label"))
        #self.pushButton.setText(_translate("MainWindow", "Change Color"))
        self.button_ok.setText(_translate("MainWindow", "Save Colors"))
        self.button_reset.setText(_translate("MainWindow", "Random Colors"))

    def connect_main(self, slider_object):
        slider_object.changeColor.connect(self.show_itself)

    def on_button(self, n):
        #print('Button {0} clicked'.format(n))
        color = QtWidgets.QColorDialog.getColor()
        print color.name()
        self.colorlist[n].setStyleSheet("QWidget { background-color: %s}" % color.name())

    def on_button_ok_clicked(self):
        #print('[TODO] write into color file')
        output_filename, _ = QFileDialog.getSaveFileName(self, 'Save File', self.color_profile_filename, filter ="txt (*.txt *.)")
        if output_filename == '':
            return
        with open(str(output_filename), 'w') as output_file:
            for i in range(len(self.colorlist)):
                cols = self.linelist[i].split()
                cols[2] = self.hex2rgb(self.color_str_list[i])
                new_line = '\t'.join(cols)
                new_line += '\n'
                output_file.write(new_line)

        self.color_change_signal.emit()

    def on_button_reset_clicked(self):
        buttonReply = QMessageBox.question(self, 'Random Colors', "This operation assigns a random generated color to each functional state, and overwrites existing color, do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if buttonReply == QMessageBox.Yes:
            r = lambda: random.randint(0,255)
            for i in range(self.color_id):
                color_name = '#%02X%02X%02X' % (r(),r(),r())
                self.colorlist[i].setStyleSheet("QWidget { background-color: %s}" % color_name)
                self.color_str_list[i] = color_name

    def hex2rgb(self, hex_str):
        rgb_color = ''
        r = int(hex_str[1:3], 16)
        g = int(hex_str[3:5], 16)
        b = int(hex_str[5:], 16)
        rgb_color = str(r) + ',' + str(g) + ',' + str(b)
        return rgb_color


        

    @pyqtSlot()
    def show_itself(self):
        self.show()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    ui = ColorWindow()
    #ui.setupUi(MainWindow)
    ui.show()
    sys.exit(app.exec_())

