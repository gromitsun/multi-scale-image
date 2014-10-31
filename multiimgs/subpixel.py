import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage

from tools import *
    
#test pattern    
a=np.arange(5)
a=np.append(a,a[::-1])
b=np.copy(a)
x,y=np.meshgrid(a,b)
im=x*y


#image
im=plt.imread('t.jpg')
im=np.mean(im,2)

#im_up=ndimage.zoom(im,1)
#im_up1=crop(im_up,new_axis=[50,200,30,160])
im_up=fft_resample(im,5)
#im_shift=shiftimg(im_up,[5.6,0])
#plt.subplot(121)
plt.imshow(im_up)
#plt.subplot(122)
#plt.imshow(im_up1)
plt.show()