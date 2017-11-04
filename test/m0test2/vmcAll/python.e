Traceback (most recent call last):
  File "plan2.py", line 11, in <module>
    plt.savefig("alpha.eps")
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/pyplot.py", line 687, in savefig
    fig = gcf()
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/pyplot.py", line 578, in gcf
    return figure()
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/pyplot.py", line 527, in figure
    **kwargs)
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/backends/backend_qt4agg.py", line 46, in new_figure_manager
    return new_figure_manager_given_figure(num, thisFig)
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/backends/backend_qt4agg.py", line 53, in new_figure_manager_given_figure
    canvas = FigureCanvasQTAgg(figure)
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/backends/backend_qt4agg.py", line 76, in __init__
    FigureCanvasQT.__init__(self, figure)
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/backends/backend_qt4.py", line 68, in __init__
    _create_qApp()
  File "/home/luca.parisi-1/anaconda2/lib/python2.7/site-packages/matplotlib/backends/backend_qt5.py", line 138, in _create_qApp
    raise RuntimeError('Invalid DISPLAY variable')
RuntimeError: Invalid DISPLAY variable
