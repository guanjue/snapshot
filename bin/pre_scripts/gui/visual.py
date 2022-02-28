# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'visual2.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class VisualWindow(QtWidgets.QMainWindow):
    def __init__(self, _output_dir):
        super(QtWidgets.QMainWindow, self).__init__()
        self.setupUi()
        self.init_signals()
        self.output_dir = _output_dir

        #image_size = QtCore.QSize(800, 600)
        #leftPixelMap = QtGui.QPixmap(self.output_dir + "/snapshot.meansig.png")
        #myScaledPixmap = leftPixelMap.scaled(image_size, QtCore.Qt.KeepAspectRatio)
        #self.label_image.setPixmap(myScaledPixmap)
        #self.label_image.setAlignment(QtCore.Qt.AlignCenter)

    def setupUi(self):
        self.setObjectName("MainWindow")
        self.resize(810, 606)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        #self.graphicsView = QtWidgets.QGraphicsView(self.centralwidget)
        #self.graphicsView.setObjectName("graphicsView")
        #self.gridLayout.addWidget(self.graphicsView, 0, 0, 1, 1)
        self.label_image = QtWidgets.QLabel(self.centralwidget)
        self.label_image.setObjectName("label_image")
        self.gridLayout.addWidget(self.label_image, 1, 0, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.radioButton_mean = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_mean.setChecked(True)
        self.radioButton_mean.setObjectName("radioButton_mean")
        self.horizontalLayout_2.addWidget(self.radioButton_mean)
        self.radioButton_freq = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_freq.setObjectName("radioButton_freq")
        self.horizontalLayout_2.addWidget(self.radioButton_freq)
        self.gridLayout.addLayout(self.horizontalLayout_2, 0, 0, 1, 1)
        self.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(self)
        self.statusbar.setObjectName("statusbar")
        self.setStatusBar(self.statusbar)
        self.actionResult_1 = QtWidgets.QAction(self)
        self.actionResult_1.setObjectName("actionResult_1")
        self.actionResult_2 = QtWidgets.QAction(self)
        self.actionResult_2.setObjectName("actionResult_2")
        self.actionResult_3 = QtWidgets.QAction(self)
        self.actionResult_3.setObjectName("actionResult_3")

        self.retranslateUi()
        #QtCore.QMetaObject.connectSlotsByName(MainWindow)


    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Visualization"))
        self.radioButton_mean.setText(_translate("MainWindow", "Index-Set Mean Signal"))
        self.radioButton_freq.setText(_translate("MainWindow", "Index-Set Most Frequent Functional Annotation"))
        self.actionResult_1.setText(_translate("MainWindow", "Result 1"))
        self.actionResult_2.setText(_translate("MainWindow", "Result 2"))
        self.actionResult_3.setText(_translate("MainWindow", "Result 3"))

    def init_signals(self):
        self.radioButton_mean.toggled.connect(self.on_radioButton_mean_toggled)
        self.radioButton_freq.toggled.connect(self.on_radioButton_freq_toggled)

    def on_radioButton_mean_toggled(self):
        image_size = QtCore.QSize(800, 600)
        
        leftPixelMap = QtGui.QPixmap(self.output_dir + "/snapshot.meansig.png")
        myScaledPixmap = leftPixelMap.scaled(image_size, QtCore.Qt.KeepAspectRatio)
        self.label_image.setPixmap(myScaledPixmap)
        self.label_image.setAlignment(QtCore.Qt.AlignCenter)
        #self.resize(800, 600)

    def on_radioButton_freq_toggled(self):
        image_size = QtCore.QSize(800, 600)
        
        leftPixelMap = QtGui.QPixmap(self.output_dir + "/snapshot.indexset_fun.png")
        myScaledPixmap = leftPixelMap.scaled(image_size, QtCore.Qt.KeepAspectRatio)
        self.label_image.setPixmap(myScaledPixmap)
        self.label_image.setAlignment(QtCore.Qt.AlignCenter)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    ui = VisualWindow('/Users/sun/data/output/')
    #ui.setupUi(MainWindow)
    ui.show()
    #ui.clearFocus()
    #ui.setFocus(True)
    sys.exit(app.exec_())

