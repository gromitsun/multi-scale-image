#coding: utf-8
"""
**Update**

1. add mapsimg.out() method
2. add output matplotlib.figure.Figure in MultiImgs.show()
3. fix attribute transfer from mapsimg to MultiImg. Relative citation changed to
absolute citation.
4. reverse what's done in 3, add set() method in MultiImg
5. use inheritance in MultiImg: MultiImg inherits mapsim;
change order of classes
6. remove 'return fig' in MultiImgs.MultiImgs.show(); add a kwarg 'display' = True/False.
7. added back 'return fig' in MultiImgs.MultiImgs.show(), if 'display' == False.
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py

from scalebar import *
from image_reg import *
from load_hdf import *

class mapsimg(object):
	def __init__(self,image=None,name=None,unit=None,xaxis=None,yaxis=None,xaxis_unit='mm',yaxis_unit='mm'):
		self.image=image
		self.name=name
		self.unit=unit
		self.xaxis=xaxis
		self.yaxis=yaxis
		self.xaxis_unit=xaxis_unit
		self.yaxis_unit=yaxis_unit
		if np.ndim(self.image) == 2:
			self.pixels = self.image.shape[::-1]
		elif np.ndim(self.image) == 3:
			self.pixels = self.image[0].shape[::-1]
		try:
			_x_dir = np.sign(xaxis[-1] - xaxis[0])
			_y_dir = np.sign(yaxis[-1] - yaxis[0])
			self.pixel_size = float(str(abs(self.xaxis[-1]-self.xaxis[0])*1./(self.pixels[0]-1)))###
			self.axis=[xaxis[0] - _x_dir*0.5*self.pixel_size,
						xaxis[-1] + _x_dir*0.5*self.pixel_size,
						yaxis[0] - _y_dir*0.5*self.pixel_size,
						yaxis[-1] + _y_dir*0.5*self.pixel_size]
		except TypeError:
			self.axis=None
	def loadh5(self,filename,element=None):
		self.__init__(readh5(filename,element))
	def show(self,**kwargs):
		element=kwargs.get('element',None)
		axis=kwargs.get('axis','on')
		scale_bar=kwargs.get('scalebar','on')
		colorbar=kwargs.get('colorbar','on')
		hold=kwargs.get('hold',False)
		vmin=kwargs.get('vmin')
		vmax=kwargs.get('vmax')
		extent=kwargs.get('extent',self.axis)
		if element!=None:
			n=self.search(element)
			image=self.image[n]
		else:
			image=self.image
		plt.imshow(image,extent=self.axis,aspect='equal',vmin=vmin,vmax=vmax, interpolate = 'none')
		if axis=='off':
			plt.axis('off')
		if colorbar=='on':
			plt.colorbar()
		if scale_bar=='on':
			scalebar(axis=extent,unit=self.xaxis_unit,**kwargs)
		plt.xlabel(self.xaxis_unit)
		plt.ylabel(self.yaxis_unit)
		if hold:
			plt.hold(True)
		else:
			plt.xlim(extent[:2])
			plt.ylim(extent[-2:])
			plt.show()
	def search(self,element):
		n=np.where(self.name==element)[0][0]
		return n
	def out(self,element=None):
		#x=object.__new__(type('mapsimg',(object,),dict(mapsimg.__dict__)))
		x=mapsimg()
		x.__dict__=self.__dict__.copy()
		if element:
			x.image=x.image[self.search(element)]
			x.name=x.name[self.search(element)]
		return x

class MultiImg(mapsimg):
	def __init__(self,image=None,layer=None,**kwargs):
		citation=kwargs.get('citation','relative')
		if isinstance(image,mapsimg):
			if citation != 'absolute':
				for attr in dir(image):
					try:
						exec('self.'+attr+'=image.'+attr)
					except AttributeError: #attribute '__weakref__' is not writable
						continue
				self.data=image
			else:
				self.__dict__=image.__dict__
		else:
			self.image=image
			self.axis=kwargs.get('axis', [0,len(image[0]),0,len(image)])
			self.name=kwargs.get('name')
			self.unit=kwargs.get('unit', 'a.u.')
			self.xaxis_unit=kwargs.get('xaxis_unit', 'mm')
			self.yaxis_unit=kwargs.get('yaxis_unit', 'mm')
		self.layer=layer
	def set(self,attr,new_value):
		exec('self.data.'+attr+'='+str(new_value))

class MultiImgs(object):
	def __init__(self,image_set=[],**kwargs):
		self.image_set=image_set
		self.image_sorted=self.sort()
		self.layers=self.image_sorted.keys()
		self.nlayers=len(self.layers)
	def list_layers(self):
		n=0
		layers=[]
		for x in self.image_set:
			if x.layer not in layers:
				layers.append(x.layer)
				n+=1
		return layers
	def sort(self):
		image_sorted={}
		for x in self.image_set:
			if x.layer in image_sorted.keys():
				image_sorted[x.layer].append(x)
			else:
				image_sorted[x.layer]=[x]
		return image_sorted
	def add_image(self,new_image):
		self.image_set.append(new_image)
		"""
		if new_image.layer in self.image_sorted.keys():
			self.image_sorted[new_image.layer].append(new_image)
		else:
			self.image_sorted[new_image.layer]=[new_image]
		"""
		self.__init__(self.image_set)
	def show(self,**kwargs):
		display = kwargs.get('display', True)
		show_layers = kwargs.get('show_layers',self.layers)
		try:
			show_layers=sorted(show_layers)
		except TypeError:
			show_layers=[show_layers]
		extent=kwargs.get('extent', 
						max_axis(*tuple(_image.axis for _image in self.image_sorted[self.layers[0]])))
		vmin=kwargs.get('vmin')
		vmax=kwargs.get('vmax')
		fig = plt.figure(figsize=(8, 8*abs((extent[3]-extent[2])*1./(extent[1]-extent[0]))))
		for layer in show_layers:
			for image in self.image_sorted[layer]:
				if layer==show_layers[0] and image==self.image_sorted[layer][0]:
					if not vmin:
						kwargs['vmin']=np.nanmin(image.image)
						vmin=np.nanmin(image.image)
					if not vmax:
						kwargs['vmax']=np.nanmax(image.image)
						vmax=np.nanmax(image.image)
					image.show(hold=True,**kwargs)
				else:
					image.show(hold=True,vmin=vmin,vmax=vmax,scalebar='off',colorbar='off')
		plt.xlim(extent[:2])
		plt.ylim(extent[-2:])
		if display:
			plt.show()
		else:
			return fig
	def imreg(self,MultiImg1,MultiImg2,**kwargs):
		kwargs['axis1']=kwargs.get('axis1',MultiImg1.axis)
		kwargs['axis2']=kwargs.get('axis2',MultiImg2.axis)
		if kwargs.get('overwrite'):
			img2.axis=reg(MultiImg1.image,MultiImg2.image,axis1=MultiImg1.axis,axis2=MultiImg2.axis)
		return reg(MultiImg1.image,MultiImg2.image,**kwargs)
		
	def savefig(self,**kwargs):
		format = kwargs.get('format', 'eps')
		filename = kwargs.get('filename', 'multiscaleimage.'+format)
		fig = self.show(display = False)
		fig.savefig(filename, format = format)

