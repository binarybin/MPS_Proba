#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import Tkinter as Tk
import tkFileDialog


def _quit(tkobj):
    tkobj.quit()     # stops mainloop
    tkobj.destroy()  # this is necessary on Windows to prevent

def getModelFile(filename):
    filename[0]=tkFileDialog.askopenfilename()

def initialization():
    root = Tk.Tk()
    root.wm_title("MPS on Probability")
    
    # the frame for the plot
    plot_frame = Tk.Frame(root) 
    plot_frame.grid(row=0, column=0)
    f = Figure(figsize=(5,4), dpi=100)

    canvas = FigureCanvasTkAgg(f, master=plot_frame)
    canvas.show()
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg( canvas, plot_frame )
    toolbar.update()
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    # the frame for the control
    control_frame = Tk.Frame(root)
    control_frame.grid(row=1, column=0)
    
    Tk.Label(master=control_frame, text="Total Run Time").grid(row=0, column=0)
    total_time_entry = Tk.Entry(master=control_frame)
    total_time_entry.grid(row=0, column=1)
    Tk.Label(master=control_frame, text=" steps").grid(row=0, column=2)
    
    Tk.Label(master=control_frame, text="Bound dimension").grid(row=1, column=0)
    bound_dim_entry = Tk.Entry(master=control_frame)
    bound_dim_entry.grid(row=1, column=1)
    
    Tk.Label(master=control_frame, text="Chain length").grid(row=2, column=0)
    chain_len_entry = Tk.Entry(master=control_frame)
    chain_len_entry.grid(row=2, column=1)
    
    filename = [""]
    model_file_button = Tk.Button(master=control_frame, text="Choose model file", command=lambda: getModelFile(filename))
    model_file_button.grid(row=3, column=1)

    compare_with_exact = Tk.IntVar()
    Tk.Checkbutton(control_frame, text = "Compare with exact solution", variable = compare_with_exact).grid(row=4,column=1)

    quit_button = Tk.Button(master=control_frame, text='Quit', command=lambda: _quit(root))
    quit_button.grid(row=10, column=0)
    
    start_button = Tk.Button(master=control_frame, text='Start!', command=lambda: _quit(root))
    start_button.grid(row=10, column=1)
    
    start_button = Tk.Button(master=control_frame, text='Reset', command=lambda: _quit(root))
    start_button.grid(row=10, column=2)
    
    # the frame for the measurement list
    measure_frame = Tk.Frame(root)
    measure_frame.grid(row=0, column=1)
    measure_list = Tk.Listbox(measure_frame, height=25, width=25)
    measure_list.pack()
    
    measure_list.insert(Tk.END, "a list entry")

    for item in ["one", "two", "three", "four"]:
        measure_list.insert(Tk.END, item)
    
    # the frame for the measurement control
    measure_control_frame = Tk.Frame(root)
    measure_control_frame.grid(row=1, column=1)
    add_proba = Tk.Button(measure_control_frame, text="Add Probability Measurement", width=28, command=lambda measure_list=measure_list: measure_list.delete(Tk.ANCHOR))
    add_proba.pack()
    add_mean = Tk.Button(measure_control_frame, text="Add Mean Measurement", width=28, command=lambda measure_list=measure_list: measure_list.delete(Tk.ANCHOR))
    add_mean.pack()
    add_variance = Tk.Button(measure_control_frame, text="Add Variance Measurement", width=28, command=lambda measure_list=measure_list: measure_list.delete(Tk.ANCHOR))
    add_variance.pack()
    
    delete_measure = Tk.Button(measure_control_frame, text="Delete A Measurement", width=28, command=lambda measure_list=measure_list: measure_list.delete(Tk.ANCHOR))
    delete_measure.pack()
    
    Tk.mainloop()


if __name__=="__main__":
    initialization()
