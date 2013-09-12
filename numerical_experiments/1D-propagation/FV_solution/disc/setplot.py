
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    

    # Figure for strain
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Stress', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Stress'
    plotaxes.scaled = False   
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = stress
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    

    plotfigure = plotdata.new_plotfigure(name='x-velocity', figno=1)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'x-velocity'
    plotaxes.scaled = False    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = xvel
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?

    plotfigure = plotdata.new_plotfigure(name='y-velocity', figno=2)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'y-velocity'
    plotaxes.scaled = False     
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = yvel
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?

    plotfigure = plotdata.new_plotfigure(name='left-going', figno=4)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'left-going part'
    plotaxes.scaled = False      
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = lwave
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?

    plotfigure = plotdata.new_plotfigure(name='right-going', figno=5)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'right-going part'
    plotaxes.scaled = False      
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = rwave
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?

    plotfigure = plotdata.new_plotfigure(name='sound speed', figno=6)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'sound speed'
    plotaxes.scaled = False      
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sound_speed
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?

    plotfigure = plotdata.new_plotfigure(name='velocity', figno=3)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'velocity'
    plotaxes.scaled = False      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
    plotitem.quiver_var_x = xvel
    plotitem.quiver_var_y = yvel
    #plotitem.quiver_coarsening_x = 4
    #plotitem.quiver_coarsening_y = 1
    plotitem.show = True       # show on plot?
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'pdf'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 1           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

def stress(current_data):    
    import numpy as np
    from psystem import setaux
    aux = setaux(current_data.x[:,0],current_data.y[0,:])
    q = current_data.q
    #return np.exp(aux[1,...]*q[0,...])-1.
    return aux[1,...]*q[0,...]

def xvel(current_data):    
    from psystem import setaux
    aux = setaux(current_data.x[:,0],current_data.y[0,:])
    q = current_data.q
    return q[1,...]/aux[0,...]

def yvel(current_data):    
    from psystem import setaux
    aux = setaux(current_data.x[:,0],current_data.y[0,:])
    q = current_data.q
    return q[2,...]/aux[0,...]

def lwave(current_data):
    from psystem import setaux
    import numpy as np
    aux = setaux(current_data.x[:,0],current_data.y[0,:])
    z = np.sqrt(aux[0,...]*aux[1,...])
    sigma = stress(current_data)
    u = xvel(current_data)
    return -(z*u+sigma)/(2.*z)

def rwave(current_data):
    from psystem import setaux
    import numpy as np
    aux = setaux(current_data.x[:,0],current_data.y[0,:])
    z = np.sqrt(aux[0,...]*aux[1,...])
    sigma = stress(current_data)
    u = xvel(current_data)
    return -(z*u-sigma)/(2.*z)

def sound_speed(current_data):
    from psystem import setaux
    import numpy as np
    aux = setaux(current_data.x[:,0],current_data.y[0,:])
    c = np.sqrt(aux[1,...]/aux[0,...])
    return c

