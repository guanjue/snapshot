# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'parameter.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot, pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QMessageBox

class ParameterWindow(QtWidgets.QDialog):
    para_change_signal = pyqtSignal()

    def __init__(self):
        super(QtWidgets.QDialog, self).__init__()
        #self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.setupUi()

    def setupUi(self):
        self.setObjectName("MainWindow")
        self.resize(438, 268)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.line_edit_count_threshold = QtWidgets.QLineEdit(self.centralwidget)
        self.line_edit_count_threshold.setObjectName("line_edit_count_threshold")
        self.gridLayout.addWidget(self.line_edit_count_threshold, 0, 1, 1, 1)
        self.checkBox_log2 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_log2.setObjectName("checkBox_log2")
        self.gridLayout.addWidget(self.checkBox_log2, 3, 1, 1, 1)
        self.label_1 = QtWidgets.QLabel(self.centralwidget)
        self.label_1.setObjectName("label_1")
        self.gridLayout.addWidget(self.label_1, 0, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.line_edit_small_number = QtWidgets.QLineEdit(self.centralwidget)
        self.line_edit_small_number.setObjectName("line_edit_small_number")
        self.gridLayout.addWidget(self.line_edit_small_number, 1, 1, 1, 1)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 2, 0, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.radioButton_frequent = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_frequent.setObjectName("radioButton_frequent")
        self.horizontalLayout_2.addWidget(self.radioButton_frequent)
        self.radioButton_mean = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_mean.setObjectName("radioButton_mean")
        self.horizontalLayout_2.addWidget(self.radioButton_mean)
        self.gridLayout.addLayout(self.horizontalLayout_2, 2, 1, 1, 1)
        self.checkBox_scale = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_scale.setObjectName("checkBox_scale")
        self.gridLayout.addWidget(self.checkBox_scale, 4, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.button_ok = QtWidgets.QPushButton(self.centralwidget)
        self.button_ok.setObjectName("button_ok")
        self.horizontalLayout.addWidget(self.button_ok)
        self.button_reset = QtWidgets.QPushButton(self.centralwidget)
        self.button_reset.setObjectName("button_reset")
        self.horizontalLayout.addWidget(self.button_reset)
        self.verticalLayout.addLayout(self.horizontalLayout)
        #self.setCentralWidget(self.centralwidget)
        #self.statusbar = QtWidgets.QStatusBar(self)
        #self.statusbar.setObjectName("statusbar")
        #self.setStatusBar(self.statusbar)

        self.retranslateUi()
        # QtCore.QMetaObject.connectSlotsByName(MainWindow)
        self.radioButton_frequent.setChecked(True)
        self.checkBox_log2.setChecked(False)
        self.checkBox_scale.setChecked(False)        
        
        self.setup_signals()

    def setup_signals(self):
        self.button_ok.clicked.connect(self.on_button_ok_clicked)
        self.button_reset.clicked.connect(self.on_button_reset_clicked)

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Parameter Panel"))
        self.label_3.setText(_translate("MainWindow", "Signal Log2"))
        self.label_4.setText(_translate("MainWindow", "Scale Signal"))
        self.line_edit_count_threshold.setText(_translate("MainWindow", "200"))
        self.checkBox_log2.setText(_translate("MainWindow", "Apply Log2 to Signal"))
        self.label_1.setText(_translate("MainWindow", "Count Threshold (Interger)"))
        self.label_2.setText(_translate("MainWindow", "Add Small Number (Float)"))
        self.line_edit_small_number.setText(_translate("MainWindow", "0.01"))
        self.label.setText(_translate("MainWindow", "Function Method"))
        self.radioButton_frequent.setText(_translate("MainWindow", "Most Frequent"))
        self.radioButton_mean.setText(_translate("MainWindow", "Mean"))
        self.checkBox_scale.setText(_translate("MainWindow", "Scale Signal"))
        self.button_ok.setText(_translate("MainWindow", "OK"))
        self.button_reset.setText(_translate("MainWindow", "Reset to Default"))

    def on_button_ok_clicked(self):
        buttonReply = QMessageBox.question(self, 'Change Parameters', "This operation changes parameters, do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if buttonReply == QMessageBox.Yes:
            self.para_change_signal.emit()

    def on_button_reset_clicked(self):
        self.line_edit_count_threshold.setText('200')
        self.line_edit_small_number.setText('0.01')
        self.radioButton_frequent.setChecked(True)
        self.checkBox_log2.setChecked(False)
        self.checkBox_scale.setChecked(False)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    ui = ParameterWindow()
    #ui.setupUi(MainWindow)
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
