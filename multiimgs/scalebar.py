import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from tools import *

"""
**Update log**
    
    1. add automatic barlen and outpu_unit
"""

def scalebar(**kwargs):
    """Add scalebar (:class:`~matplotlib.patches.Patch`) on current axis 
    (:class:`~matplotlib.axes.Axes`).
    
    Keyword arguments:
        *axis*: scalars (left, right, bottom, top)
            Frame limits on current axis. If given, the function will take the 
            region delineated by axis as the region to display. If not given, 
            the default value is the current axis given by plt.gca()
        
        *loc*: ['lowerleft' | 'lowerright' | 'upperleft' | 'upperright']
            Set the location of the scalebar on the frame. Default='lowerleft'.
        
        *color*: any matplotlib color
            Set the color of scalebar and text on scalebar. Default='black'
        
        *textcolor*: any matplotlib color
            Set the color of text on scalebar. Default=color
        
        *barcolor*: any matplotlib color
            Set the color of scalebar. Default=color
        
        *unit*: string
            Unit used in current axis. Default='mm'
        
        *barlen*: scaler
            Length of the scalebar in unit designated in 'unit'. If not given,
            automatically choose the most appropriate length.
            
        *direction*: ['horizontal' | 'vertical']
            Direction of the scalebar. Default='horizontal' 
            
        *text_direction*: ['horizontal' | 'vertical']
            Direction of the text on scalebar. Default='horizontal' 
            
        *texton*: [True | False]
            Show text on scalebar or not. Default=True
            
        *output_unit*: [*None* | string]
            The unit to be shown on the scalebar. If output_unit='auto', it will
            automatically choose the most appropriate unit.
        
    **Example**
        import numpy as np
        import matplotlib.pyplot as plt
        x=np.arange(100)
        y=np.arange(100)
        x,y=np.meshgrid(x,y)
        img=x*y
        axis=[-5.2,-15.1,-1.3,-6.4]
        plt.imshow(img,origin='lower',extent=axis)
        scalebar(color='white',unit='cm',output_unit='auto')
        plt.show()
        """
    
    plt.hold(True)
    axis=kwargs.get('axis',plt.axis())
    loc=kwargs.get('loc','lowerleft')
    color=kwargs.get('color','black')
    textcolor=kwargs.get('textcolor',color)
    barcolor=kwargs.get('barcolor',color)
    unit=kwargs.get('unit','mm')
    direction=kwargs.get('direction','horizontal')
    text_direction=kwargs.get('text_direction',direction)
    barlen=kwargs.get('barlen',search(abs(axis[1]-axis[0])*0.15))
    texton=kwargs.get('texton',True)
    output_unit=kwargs.get('output_unit')
    
    if output_unit=='auto':
        barlen_text,unit=auto_unit(barlen,unit)
    elif output_unit:
        barlen_text=int_out(float(str(convert_unit(barlen,unit,new_unit=output_unit))))
        unit=output_unit
    else:
        barlen_text=barlen


    #scalebar location
    if loc=='lowerleft':
        xloc=axis[0]+(axis[1]-axis[0])*0.07
        yloc=axis[2]+(axis[3]-axis[2])*0.07
        islower=1
        isleft=1
    elif loc=='upperleft':
        xloc=axis[0]+(axis[1]-axis[0])*0.07
        yloc=axis[3]-(axis[3]-axis[2])*0.07
        islower=-1
        isleft=1
    elif loc=='lowerright':
        xloc=axis[1]-(axis[1]-axis[0])*0.07
        yloc=axis[2]+(axis[3]-axis[2])*0.07
        islower=1
        isleft=-1
    elif loc=='upperright':
        xloc=axis[1]-(axis[1]-axis[0])*0.07
        yloc=axis[3]-(axis[3]-axis[2])*0.07
        islower=-1
        isleft=-1
    elif type(loc) is not str:
        xloc=loc[0]
        yloc=loc[1]
        islower=1
        isleft=1
    #scalebar length
    x_axis_direction=np.sign(axis[1]-axis[0])
    y_axis_direction=np.sign(axis[3]-axis[2])
    if direction=='vertical':
        ylen=barlen*y_axis_direction
        xlen=(axis[1]-axis[0])*0.02
        if texton:
            x_text=xloc+isleft*xlen*2.5
            y_text=yloc+islower*ylen/2. 
            plt.text(x_text,y_text,str(barlen_text)+' '+unit,rotation=text_direction,color=textcolor,verticalalignment='center',withdash=True,horizontalalignment='center')     
    else:
        xlen=barlen*x_axis_direction
        ylen=(axis[3]-axis[2])*0.02
        if texton:
            x_text=xloc+isleft*xlen/2.
            y_text=yloc+islower*ylen*2.5 
            plt.text(x_text,y_text,str(barlen_text)+' '+unit,rotation=text_direction,color=textcolor,verticalalignment='center',withdash=True,horizontalalignment='center')   
    current_axis=plt.gca()
    current_axis.add_patch(Rectangle([xloc,yloc],isleft*xlen,islower*ylen,color=barcolor))
    
"""
#*Example*
    
x=np.arange(100)

y=np.arange(100)

x,y=np.meshgrid(x,y)

img=x*y
axis=np.array([-5.2,-15.1,-1.3,-6.4])*0.01
plt.imshow(img,origin='lower',extent=axis)
#cl=plt.axis()
#ca=plt.gca()
#ca.add_patch(Rectangle([-6,-4],1,1,facecolor='red'))
#scalebar(direction='horizontal',barcolor='white',textcolor='red',loc='upperleft',barlen=2,unit='mm')
scalebar(color='white',output_unit='auto')
plt.show()
"""