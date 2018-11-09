###                       BIFS GUI                 ###
###                                                ###
###        This file contains the GUI interface    ###
###        for the BIFS package                    ###
###                                                ###

from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import NoNorm 
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from scipy import misc
import bifs
import jsonpickle

from pset_dialogs import Param_Fourier_Space_Dialog,Prior_Dialog
from pset_dialogs import Likelihood_Dialog,Slice3D_Dialog
from pset_dialogs import AddBump_Dialog,DeleteBump_Dialog

class MainWindow(QtGui.QMainWindow):
    """

    Class that generates the main BIFS GUI window
 
    """
    send_fig = QtCore.pyqtSignal(str)

    def __init__(self):
        super(MainWindow, self).__init__()
        self.main_widget = QtGui.QWidget(self)
        # Initialize BIFS object
        self.mybifs = bifs.bifs()
        self.filename = None
        self.didMAP = False
        self.setWindowTitle("Bayesian Imaging in Fourier Space")
        self.fig = Figure()
        plt.rcParams["axes.grid"] = False # turn off grid lines for images
        plt.rcParams["xtick.color"] = (1,1,1,0)
        plt.rcParams["ytick.color"] = (1,1,1,0)
        self.ax1 = self.fig.add_subplot(221,adjustable='box')
        self.ax2 = self.fig.add_subplot(222,adjustable='box')
        self.ax3 = self.fig.add_subplot(223,adjustable='box')
        self.ax4 = self.fig.add_subplot(224,adjustable='box')
        self.axes=[self.ax1, self.ax2, self.ax3, self.ax4]
        self.canvas = FigureCanvas(self.fig)

        self.canvas.setSizePolicy(QtGui.QSizePolicy.Expanding, 
                                   QtGui.QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        self.layout = QtGui.QGridLayout(self.main_widget)
        self.layout.addWidget(self.canvas,0,0,5,5)
        self.createActions()
        self.createMenus()

        # left, top, width, height
        self.setGeometry(30,30,850,850)
        self.setCentralWidget(self.main_widget)
        self.show()
        self.update()

    def getImage(self):
        """
    
        function that calls the BIFS function load_image_file()
        to load an image file.  

        """
        # self.__init__()
        self.fileName = QtGui.QFileDialog.getOpenFileName(self, "Open File",
                                                     QtCore.QDir.currentPath())

        try:
            self.mybifs.load_image_file(self.fileName)
            self.mybifs.image_file_loaded = True
            self.show_initial_image()
            return
        except:
            QtGui.QMessageBox.information(self, "Image Viewer","Cannot load %s % self.fileName")
            return

    def doMAP(self):
        """
    
        function that calls the BIFS function BIFS_MAP()
        to perform a MAP analysis.  

        """
        # Do MAP estimation
        if np.isscalar(self.mybifs.init_image):
            QtGui.QMessageBox.information(self,"MAP Estimator", "Can't perform MAP without an image - load image using Load Initial Image... option from BIFS menu")
            return
        else:
            #try:
            print("Performing k-space MAP estimation")
            # Reinitialize bifs object but keep current parameter setttings
            pft = self.mybifs.param_func_type
            dec = self.mybifs.decay
            pri = self.mybifs.prior
            prs = self.mybifs.prior_scale
            lik = self.mybifs.likelihood
            lis = self.mybifs.likelihood_scale
            bps = self.mybifs.bumps
            v3d = self.mybifs.view3Dslice
            self.mybifs = bifs.bifs()
            self.mybifs.load_image_file(self.fileName)
            self.mybifs.image_file_loaded = True
            self.mybifs.param_func_type = pft
            self.mybifs.decay = dec
            self.mybifs.prior = pri
            self.mybifs.prior_scale = prs
            self.mybifs.likelihood = lik
            self.mybifs.likelihood_scale = lis
            self.mybifs.bumps = bps
            self.mybifs.view3Dslice = v3d
            self.mybifs.BIFS_MAP()
            self.didMAP = True
            self.show_post_proc_images()
            #except:
            #    QtGui.QMessageBox.information(self,"MAP Estimator", "MAP estimate failed") 
            return
            
    def show_initial_image(self):
        """
    
        function that displays the initial image.  

        """
        try:
            if self.mybifs.image_file_loaded == True:
                self.ax1.clear()
                self.canvas.draw()
                self.ax1.set_title("Initial Image")
                if self.mybifs.imdim == 1:
                    self.ax1.plot(self.mybifs.init_image) 
                elif self.mybifs.imdim == 2:
                    self.ax1.imshow(self.mybifs.init_image,cmap = cm.Greys_r,norm=NoNorm())
                else: # assume for now that the only other possibility is 3D
                      # can change later
                    slice_index = np.int(np.round(self.mybifs.view3Dslice[1]*self.mybifs.init_image.shape[self.mybifs.view3Dslice[0]]))
                    if self.mybifs.view3Dslice[0] == 0:
                        init_im_slice = self.mybifs.init_image[slice_index,:,:]
                    elif self.mybifs.view3Dslice[0] == 1:
                        init_im_slice = self.mybifs.init_image[:,slice_index,:]
                    elif self.mybifs.view3Dslice[0] == 2:
                        init_im_slice = self.mybifs.init_image[:,:,slice_index]
                    else:
                        print("Sorry slice index needs to be one of 0,1,2")
                    self.ax1.imshow(init_im_slice, cmap = cm.Greys_r)
                self.canvas.draw()
            return
        except:
            QtGui.QMessageBox.information(self, "Image Viewer","No image loaded %s.")
            return
        
    def show_post_proc_images(self):
        """
    
        function that displays the processed images.  

        """
        # Show initial and post BIFS k-spage images and final image
        try:
            if self.mybifs.image_file_loaded == True:

                self.ax1.set_title("Initial Image") # For some reason plt.cla()
                                                    # clears plot 1 title
                self.ax2.clear()
                self.canvas.draw()
                self.ax2.set_title("Reconstructed Image")
                if self.mybifs.imdim == 1:
                    self.ax2.plot(self.mybifs.final_image) 
                elif self.mybifs.imdim == 2:
                    self.ax2.imshow(self.mybifs.final_image,cmap = cm.Greys_r)
                else: # assume for now that the only other possibility is 3D
                      # can change later
                    slice_index = np.int(np.round(self.mybifs.view3Dslice[1]*self.mybifs.final_image.shape[self.mybifs.view3Dslice[0]]))
                    if self.mybifs.view3Dslice[0] == 0:
                        final_im_slice = self.mybifs.final_image[slice_index,:,:]
                    elif self.mybifs.view3Dslice[0] == 1:
                        final_im_slice = self.mybifs.final_image[:,slice_index,:]
                    elif self.mybifs.view3Dslice[0] == 2:
                        final_im_slice = self.mybifs.final_image[:,:,slice_index]
                    else:
                        print("Sorry slice index needs to be one of 0,1,2")
                    self.ax2.imshow(final_im_slice, cmap = cm.Greys_r)
                self.canvas.draw()
                
                self.ax3.clear()
                self.canvas.draw()
                self.ax3.set_title("Initial K-Space Modulus Image")
                if self.mybifs.imdim == 1:
                    showim1k = np.roll(np.roll(self.mybifs.mod_image,self.mybifs.mod_image.shape[0]//2,0),1)
                    self.ax3.plot(np.log(showim1k)) 
                elif self.mybifs.imdim == 2:
                    showim1k = np.roll(np.roll(self.mybifs.mod_image,self.mybifs.mod_image.shape[0]//2,0),self.mybifs.mod_image.shape[1]//2,1)
                    self.ax3.imshow(np.log(showim1k),cmap = cm.Greys_r)
                else: # assume for now that the only other possibility is 3D
                      # can change later
                    slice_index = np.int(np.round(self.mybifs.view3Dslice[1]*self.mybifs.mod_image.shape[self.mybifs.view3Dslice[0]]))
                    if self.mybifs.view3Dslice[0] == 0:
                        init_mod_im_slice = self.mybifs.mod_image[slice_index,:,:]
                    elif self.mybifs.view3Dslice[0] == 1:
                        init_mod_im_slice = self.mybifs.mod_image[:,slice_index,:]
                    elif self.mybifs.view3Dslice[0] == 2:
                        init_mod_im_slice = self.mybifs.mod_image[:,:,slice_index]
                    else:
                        print("Sorry slice index needs to be one of 0,1,2")

                    showim1k = np.roll(np.roll(init_mod_im_slice,init_mod_im_slice.shape[0]//2,0),init_mod_im_slice.shape[1]//2,1)
                    self.ax3.imshow(np.log(showim1k), cmap = cm.Greys_r)
                self.canvas.draw()

                self.ax4.clear()
                self.canvas.draw()
                self.ax4.set_title("BIFS K-Space Image")
                if self.mybifs.imdim == 1:
                    showim2k = np.roll(np.roll(self.mybifs.bifsk_image,self.mybifs.bifsk_image.shape[0]//2,0),1)
                    self.ax4.plot(np.log(showim2k)) 
                elif self.mybifs.imdim == 2:
                    showim2k = np.roll(np.roll(self.mybifs.bifsk_image,self.mybifs.bifsk_image.shape[0]//2,0),self.mybifs.bifsk_image.shape[1]//2,1)
                    self.ax4.imshow(np.log(showim2k),cmap = cm.Greys_r)
                else: # assume for now that the only other possibility is 3D
                      # can change later
                    slice_index = np.int(np.round(self.mybifs.view3Dslice[1]*self.mybifs.bifsk_image.shape[self.mybifs.view3Dslice[0]]))
                    if self.mybifs.view3Dslice[0] == 0:
                        bifs_mod_im_slice = self.mybifs.bifsk_image[slice_index,:,:]
                    elif self.mybifs.view3Dslice[0] == 1:
                        bifs_mod_im_slice = self.mybifs.bifsk_image[:,slice_index,:]
                    elif self.mybifs.view3Dslice[0] == 2:
                        bifs_mod_im_slice = self.mybifs.bifsk_image[:,:,slice_index]
                    else:
                        print("Sorry slice index needs to be one of 0,1,2")

                    showim2k = np.roll(np.roll(bifs_mod_im_slice,bifs_mod_im_slice.shape[0]//2,0),bifs_mod_im_slice.shape[1]//2,1)
                    self.ax4.imshow(np.log(showim2k), cmap = cm.Greys_r)
                self.canvas.draw()
        except:
            QtGui.QMessageBox.information(self, "Image Viewer","No image loaded %s.")
            return
        return
    
    def saveCurrent(self):
        """
    
        function that calls the BIFS function save_results()
        to save the results of analysis.  

        """
        if np.isscalar(self.mybifs.final_image):
            QtGui.QMessageBox.information(self,"Save Current Results","No Results to Output Yet; Probably need to run - Get MAP Estimate Image - first")
        else:
            print("Saving Current Results")
            self.mybifs.save_results()

    def loadPrevious(self):
        """
    
        function that loads the results of a previous analysis
        saved via the BIFS function save_results().  

        """        
        fileName = QtGui.QFileDialog.getOpenFileName(self, "Load Existing Paramter File",QtCore.QDir.currentPath())
        # Check if it's a bifs parameter file
        if fileName.split('.')[1] == 'bifspout':
            bifs_handle = open(fileName,"r")
            bifsj = bifs_handle.read()
            bifs_handle.close()
            self.mybifs = jsonpickle.decode(bifsj)
            # Show intial image
            self.show_initial_image()
            # Show initial and post BIFS k-spage images and final image
            self.show_post_proc_images()                    
        else:
            QtGui.QMessageBox.information(self, "Parameter File Loader","Cannot load %s., doesn't look like a bifs parameter file" % fileName)
            return
        return
    
    def close(self):
        sys.exit(app.exec_())
        
    def zoomIn(self):
        self.scaleImage(1.25)

    def zoomOut(self):
        self.scaleImage(0.8)

    def normalSize(self):
        self.imageLabel.adjustSize()
        self.scaleFactor = 1.0

    def fitToWindow(self):
        fitToWindow = self.fitToWindowAct.isChecked()
        self.scrollArea.setWidgetResizable(fitToWindow)
        if not fitToWindow:
            self.normalSize()

        self.updateActions()

    def createActions(self):

        self.setParamFuncAct = QtGui.QAction("&Parameter Space Function", self, triggered=self.setParamFunc)

        self.setPriorAct = QtGui.QAction("&Prior Distribution", self, triggered=self.setPrior)

        self.setLikelihoodAct = QtGui.QAction("&Likelihood Distribution", self, triggered=self.setLikelihood)

        self.setAddBumpAct =  QtGui.QAction("&Add K-Space Bump", self, triggered=self.setAddBump)

        self.setDeleteBumpAct =  QtGui.QAction("&Delete K-Space Bump", self, triggered=self.setDeleteBump)
        
        self.set3DSliceAct = QtGui.QAction("&3D Slice to View", self, triggered=self.set3DSlice)
        
        self.getImageAct = QtGui.QAction("&Load Initial Image...", self, triggered=self.getImage)

        self.doMapAct = QtGui.QAction("&Get MAP Estimate Image...", self,triggered=self.doMAP)

        self.saveCurrentAct = QtGui.QAction("&Save Current Results...", self,triggered=self.saveCurrent)

        self.loadPreviousAct = QtGui.QAction("&Load Previous Results...", self,triggered=self.loadPrevious)
        
        self.exitAct = QtGui.QAction("&Exit", self, shortcut="Ctrl+Q",
                triggered=self.close)

        self.zoomInAct = QtGui.QAction("Zoom &In (25%)", self,
                shortcut="Ctrl++", enabled=False, triggered=self.zoomIn)

        self.zoomOutAct = QtGui.QAction("Zoom &Out (25%)", self,
                shortcut="Ctrl+-", enabled=False, triggered=self.zoomOut)

        self.normalSizeAct = QtGui.QAction("&Normal Size", self,
                shortcut="Ctrl+S", enabled=False, triggered=self.normalSize)

        self.fitToWindowAct = QtGui.QAction("&Fit to Window", self,
                enabled=False, checkable=True, shortcut="Ctrl+F",
                triggered=self.fitToWindow)

        # self.aboutAct = QtGui.QAction("&About", self, triggered=self.about)

        # self.aboutQtAct = QtGui.QAction("&About Qt", self,
        #        triggered=QtGui.qApp.aboutQt)

    def createMenus(self):
        self.bifsMenu = QtGui.QMenu("&BIFS Operations", self)
        self.bifsMenu.addAction(self.getImageAct)
        self.bifsMenu.addAction(self.doMapAct)
        self.bifsMenu.addAction(self.saveCurrentAct)
        self.bifsMenu.addAction(self.loadPreviousAct)
        self.bifsMenu.addSeparator()
        self.bifsMenu.addAction(self.exitAct)

        self.paramMenu = QtGui.QMenu("&Parameter Set", self)
        self.paramMenu.addAction(self.setParamFuncAct)
        self.paramMenu.addAction(self.setPriorAct)
        self.paramMenu.addAction(self.setLikelihoodAct)
        self.paramMenu.addAction(self.setAddBumpAct)
        self.paramMenu.addAction(self.setDeleteBumpAct)
        self.paramMenu.addAction(self.set3DSliceAct)
        self.paramMenu.addSeparator()
        self.paramMenu.addAction(self.exitAct)
        
        self.viewMenu = QtGui.QMenu("&View", self)
        self.viewMenu.addAction(self.zoomInAct)
        self.viewMenu.addAction(self.zoomOutAct)
        self.viewMenu.addAction(self.normalSizeAct)
        self.viewMenu.addSeparator()
        self.viewMenu.addAction(self.fitToWindowAct)

        # self.helpMenu = QtGui.QMenu("&Help", self)
        # self.helpMenu.addAction(self.aboutAct)
        # self.helpMenu.addAction(self.aboutQtAct)
        
        self.menuBar().addMenu(self.bifsMenu)
        self.menuBar().addMenu(self.paramMenu)
        self.menuBar().addMenu(self.viewMenu)
        # self.menuBar().addMenu(self.helpMenu)

    def updateActions(self):
        self.zoomInAct.setEnabled(not self.fitToWindowAct.isChecked())
        self.zoomOutAct.setEnabled(not self.fitToWindowAct.isChecked())
        self.normalSizeAct.setEnabled(not self.fitToWindowAct.isChecked())
                    
    def adjustScrollBar(self, scrollBar, factor):
        scrollBar.setValue(int(factor * scrollBar.value()
                                + ((factor - 1) * scrollBar.pageStep()/2)))

    def setParamFunc(self):
        if self.mybifs.basis == "Fourier": # add else clauses as bases are added
            dialog = StartParam_Fourier_Space_Dialog(self)
            if dialog.exec_() == 1: # If OK button pressed
                self.mybifs.param_func_type = dialog.getFunc()
                self.mybifs.decay = dialog.getDecay()
        else:
            print("Sorry, don't know how to handle",self.mybifs.basis,"basis yet.")

    def setPrior(self):
        dialog = StartPrior_Dialog(self)
        if dialog.exec_() == 1: # If OK button pressed
            self.mybifs.prior = dialog.getDist()
            self.mybifs.prior_scale = dialog.getScale()

    def setLikelihood(self):
        dialog = StartLikelihood_Dialog(self)
        if dialog.exec_() == 1: # If OK button pressed
            self.mybifs.likelihood = dialog.getDist()
            self.mybifs.likelihood_scale = dialog.getScale()

    def setAddBump(self):
        if np.isscalar(self.mybifs.init_image):
            QtGui.QMessageBox.information(self,"Add Bump", "Can't add bump without k-space structure, need to first load image using Load Initial Image... option from BIFS menu")
            return
        else:
            dialog = StartAddBump_Dialog(self)
            if dialog.exec_() == 1: # If OK button pressed
                get_key = dialog.getKey()
                get_position = dialog.getPosition()
                get_amplitude = dialog.getAmplitude()
                get_width = dialog.getWidth()
                self.mybifs.add_bump(get_key,get_position,get_amplitude,get_width)
                return

#                QtGui.QMessageBox.information(self, "Add Bump","Failed")
#                return
    def setDeleteBump(self):
        if self.mybifs.bumps:
            self.bumplist = []
            for keys in self.mybifs.bumps:
                self.bumplist.append(str(keys)+"  "+str(self.mybifs.bumps[keys][0]) + "  " + str(self.mybifs.bumps[keys][1]) + "  " + str(self.mybifs.bumps[keys][2]))
            dialog = StartDeleteBump_Dialog(self)
            if dialog.exec_() == 1: # If OK button pressed
                delstr = dialog.getDeleteString()
                delkey = str(delstr).split(" ")[0].strip("'[")
                del self.mybifs.bumps[delkey]
            return
        else:
            print("No bumps to delete yet")
            return
        
    def set3DSlice(self):
        if self.mybifs.imdim != 3:
            print("Sorry can only view slices for 3D data")
            return
        else:
            dialog = Start3DSlice_Dialog(self)
            if dialog.exec_() == 1: # If OK button pressed
                self.mybifs.view3Dslice[0] = dialog.getSliceIndex()
                self.mybifs.view3Dslice[1] = dialog.getSlicePercent()
                if self.didMAP:
                    self.show_initial_image()
                    self.show_post_proc_images()
                else:
                    if self.mybifs.image_file_loaded:
                        self.show_initial_image()
                return

            
class StartParam_Fourier_Space_Dialog(QtGui.QDialog, Param_Fourier_Space_Dialog):
    def __init__(self,parent=MainWindow):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

    def getDecay(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.decay))
        except:
            print("Couldn't get decay value")
            
    def getFunc(self):
        try:
            return str(self.param_spaceBox.currentText())
        except:
            print("Couldn't get parameter function")

class StartPrior_Dialog(QtGui.QDialog, Prior_Dialog):
    def __init__(self,parent=MainWindow):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

    def getScale(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.scale))
        except:
            print("Couldn't get scale value")
            
    def getDist(self):
        try:
            return str(self.priorBox.currentText())
        except:
            print("Couldn't get prior distribution")
            
class StartLikelihood_Dialog(QtGui.QDialog, Likelihood_Dialog):
    def __init__(self,parent=MainWindow):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

    def getScale(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.scale))
        except:
            print("Couldn't get scale value")
            
    def getDist(self):
        try:
            return str(self.likelihoodBox.currentText())
        except:
            print("Couldn't get likelihood distribution")

class Start3DSlice_Dialog(QtGui.QDialog, Slice3D_Dialog):
    def __init__(self,parent=MainWindow):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

    def getSliceIndex(self):
        try:
            return np.int(QtGui.QLineEdit.text(self.slice_index))
        except:
            print("Couldn't get slice index")
            
    def getSlicePercent(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.slice_percent))
        except:
            print("Couldn't get slice percent")

class StartAddBump_Dialog(QtGui.QDialog, AddBump_Dialog):
    def __init__(self,parent=MainWindow):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

    def getKey(self):
        try:
            return str(self.addBumpKeyBox.currentText())
        except:
            print("Couldn't get bump key")
            
    def getPosition(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.Position))
        except:
            print("Couldn't get bump position")
            
    def getAmplitude(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.Amplitude))
        except:
            print("Couldn't get bump position")
            
    def getWidth(self):
        try:
            return np.float(QtGui.QLineEdit.text(self.Width))
        except:
            print("Couldn't get bump width")
            
class StartDeleteBump_Dialog(QtGui.QDialog, DeleteBump_Dialog):
    def __init__(self,parent=MainWindow):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

    def getDeleteString(self):
        try:
            return str(self.deleteBumpKeyBox.currentText())
        except:
            print("Couldn't get bump function information")

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    win = MainWindow()
    sys.exit(app.exec_())
