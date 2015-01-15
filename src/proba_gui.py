#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import Tkinter as Tk
import tkFileDialog

import sys
from copy import deepcopy
from angryboys import AngryBoys
from radiatingboys import RadiatingBoys
from exponentialboys import ExponentialBoys
from exactsolver import ExactSolver
from mpssolver import MpsSolver
from exactmeasurement import ExactMeasurement
from mpsmeasurement import MpsMeasurement
from projectionboys import ProjectionBoys
import pylab

class gui_main(Tk.Tk):
    def __init__(self):
        Tk.Tk.__init__(self)
        
        self.wm_title("MPS on Probability")
        self.setPlotFrame()
        self.setControlFrame()
        self.setMeasureFrame()
        self.setMeasureControlFrame()
        self.initVariables()
        
    def initVariables(self):
        pass
           
    def setPlotFrame(self):
        self.plot_frame = Tk.Frame(self)
        self.plot_frame.grid(row=0, column=0)
        self.fig = Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.plot_frame )
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
                
    def setControlFrame(self):
        self.control_frame = Tk.Frame(self)
        self.control_frame.grid(row=1, column=0)
        Tk.Label(master=self.control_frame, text="Total Run Time").grid(row=0, column=0)
        
        # total time entry
        self.v_total_time = Tk.IntVar()
        self.total_time_entry = Tk.Entry(master=self.control_frame, textvariable = self.v_total_time)
        self.total_time_entry.grid(row=0, column=1)
        self.v_total_time.set(100)
        Tk.Label(master=self.control_frame, text=" steps").grid(row=0, column=2)

    
        # bound dimension entry
        self.v_bound_dim = Tk.IntVar()
        Tk.Label(master=self.control_frame, text="Bound dimension").grid(row=1, column=0)
        self.bound_dim_entry = Tk.Entry(master=self.control_frame, textvariable = self.v_bound_dim)
        self.bound_dim_entry.grid(row=1, column=1)
        self.v_bound_dim.set(10)
    
    
        # chain length entry
        self.v_chain_len = Tk.IntVar()
        Tk.Label(master=self.control_frame, text="Chain length").grid(row=2, column=0)
        self.chain_len_entry = Tk.Entry(master=self.control_frame, textvariable = self.v_chain_len)
        self.chain_len_entry.grid(row=2, column=1)
        self.v_chain_len.set(10)
        
        self.model_list = ["AngryBoys", "RadiatingBoys", "ExponentialBoys", "ProjectionBoys"]
        
        self.v_model = ""
        for i, text in enumerate(self.model_list):
            Tk.Radiobutton(self.control_frame, text=text, variable=self.v_model, value=text).grid(row=3, column=i)

        self.compare_with_exact = Tk.IntVar()
        Tk.Checkbutton(self.control_frame, text = "Compare with exact solution", variable = self.compare_with_exact).grid(row=4,column=1)

        self.quit_button = Tk.Button(master=self.control_frame, text='Quit', command=self._quit)
        self.quit_button.grid(row=10, column=0)
    
        self.start_button = Tk.Button(master=self.control_frame, text='Start!', command=self.start)
        self.start_button.grid(row=10, column=1)
    
        self.reset_button = Tk.Button(master=self.control_frame, text='Reset', command=self.reset)
        self.reset_button.grid(row=10, column=2)
        
    def start(self):
        pass
        
    def reset(self):
        pass
        
    def setMeasureFrame(self):
        self.measure_frame = Tk.Frame(self)
        self.measure_frame.grid(row=0, column=1)
        self.measure_list = Tk.Listbox(self.measure_frame, height=25, width=25)
        self.measure_list.pack()
        
    def setMeasureControlFrame(self):
        self.measure_control_frame = Tk.Frame(self)
        self.measure_control_frame.grid(row=1, column=1)
        
        self.add_proba = Tk.Button(self.measure_control_frame, text="Add Probability Measurement", width=28, command=self.addProba)
        self.add_proba.pack()
        self.add_mean = Tk.Button(self.measure_control_frame, text="Add Mean Measurement", width=28, command=self.addMean)
        self.add_mean.pack()
        self.add_variance = Tk.Button(self.measure_control_frame, text="Add Variance Measurement", width=28, command=self.addVariance)
        self.add_variance.pack()
    
        self.delete_measure = Tk.Button(self.measure_control_frame, text="Delete A Measurement", width=28, command=lambda lb=self.measure_list: lb.delete(Tk.ANCHOR))
        self.delete_measure.pack()
        
    def addProba(self):
        pass
    def addMean(self):
        pass
    def addVariance(self):
        pass
    
    def _quit(self):
        self.quit()     # stops mainloop
        self.destroy()  # this is necessary on Windows to prevent

    def getModelFile(self):
        self.filename=tkFileDialog.askopenfilename()
        print self.filename

    



if __name__=="__main__":
    gui = gui_main()
    gui.mainloop()
