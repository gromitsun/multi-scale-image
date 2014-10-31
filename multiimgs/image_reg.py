"""
**update**
1. move zeropadding to xcorr
2. corrected typo 'tuckey' -> 'tukey'
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.patches import Rectangle
from tools import *

def max_axis(axis1, axis2 = None, *args):
	out = []
	if axis2 == None:
		return axis1
	if len(args) > 0:
		axis2 = max_axis(axis2, *args)
	x_min = np.min((axis1[:2], axis2[:2]))
	x_max = np.max((axis1[:2], axis2[:2]))
	y_min = np.min((axis1[-2:], axis2[-2:]))
	y_max = np.max((axis1[-2:], axis2[-2:]))
	if axis1[0] < axis1[1]:
		out += [x_min, x_max]
	else:
		out += [x_max, x_min]
	if axis1[2] < axis1[3]:
		out += [y_min, y_max]
	else:
		out += [y_max, y_min]
	return out

def min_axis(axis1, axis2, *args):
	if axis2 == None:
		return axis1
	if len(args) > 0:
		axis2 = min_axis(axis2, *args)
	if axis1[0] < axis1[1]:
		x_left = max(axis1[0], axis2[0])
		x_right = min(axis1[1], axis2[1])
	else:
		x_left = min(axis1[0], axis2[0])
		x_right = max(axis1[1], axis2[1])
	if axis1[2] < axis1[3]:
		y_bottom = max(axis1[2], axis2[2])
		y_top = min(axis1[3], axis2[3])
	else:
		y_bottom = min(axis1[2], axis2[2])
		y_top = max(axis1[3], axis2[3])
	return [x_left, x_right, y_bottom, y_top]
	

def reg(img1,img2,**kwargs):
	
	#get kwargs
	t=kwargs.get('tukey')
	axis1=kwargs.get('axis1',[0,len(img1[0]),0,len(img1)])
	axis2=kwargs.get('axis2',[0,len(img2[0]),0,len(img2)])
	precision=kwargs.get('precision',abs(axis1[1]-axis1[0])*1./len(img1[0]))
	method=kwargs.get('method','fft')
	method_us=kwargs.get('method_us',method)
	method_ds=kwargs.get('method_ds',method)
	remove_dc=kwargs.get('remove_dc',True)
	region=kwargs.get('region')
	output=kwargs.get('output','new_axis2')

	x1_dir=np.sign(axis1[1]-axis1[0])
	y1_dir=np.sign(axis1[3]-axis1[2])
	x2_dir=np.sign(axis2[1]-axis2[0])
	y2_dir=np.sign(axis2[3]-axis2[2])
	
	#crop img1 according to 'region' and calculate axis1_compute
	img1_copy=img1[:]
	if not region:
		region = axis1
	axis_compute = [region[0], region[1]+x1_dir*abs(axis2[1]-axis2[0]), region[2]-y1_dir*abs(axis2[3]-axis2[2]), region[3]]
	axis1_crop = min_axis(axis1, axis_compute)
	if axis1_crop != axis1:
		img1, axis1_crop = crop(img1, axis = axis1, new_axis = axis1_crop, actual_axis = True)
	
	#process NaNs
	if np.isnan(img1).any():
		img1=nan_to_value(img1)
	if np.isnan(img2).any():
		img2=nan_to_value(img2)

	#choose up- and down- sampling method
	if method_us==('scipy' or 'interpolation'):
		fun_us=ndimage.zoom
	elif method_us=='fft':
		fun_us=fft_resample
	else:
		raise ValueError("Method %s is not recognized" % method_us)
	if method_ds==('scipy' or 'interpolation'):
		fun_ds=ndimage.zoom
	elif method_ds=='fft':
		fun_ds=fft_resample
	else:
		raise ValueError("Method %s is not recognized" % method_ds)	
	
	#rescaling
	res1=abs(region[3]-region[2])*1./len(img1)
	res2=abs(axis2[3]-axis2[2])*1./len(img2)
	if precision<res1:
		img1=fun_us(img1,res1/precision, **kwargs)
	elif precision>res1:
		img1=fun_ds(img1,res1/precision, **kwargs)
	if precision<res2:
		img2=fun_us(img2,res2/precision, **kwargs)
	elif precision>res2:
		img2=fun_ds(img2,res2/precision, **kwargs)
		
	#remove dc
	if remove_dc:
		img1 = img1 - np.average(img1)
		img2 = img2 - np.average(img2)
	
	#tukey window
	if t:
		try:
			t1=t[0]
			t2=t[1]
		except TypeError:
			t1=t
			t2=t
		img1=img1*tukey2d(t1,len(img1),len(img1[0]))
		img2=img2*tukey2d(t2,len(img2),len(img2[0]))
	
	#zeropad img1 to axis1_compute
	if axis1_crop != axis_compute:
		img1, axis_compute = zeropadding(img1, axis = axis1_crop, new_axis = axis_compute, actual_axis = True)
	
	
	#xcorr
	cor=xcorr(img2,img1)
	
	#crop cor according to region
	cor, region = crop(cor, axis = axis_compute, new_axis = region, actual_axis = True)
	
	#best location
	y,x=np.where(cor==np.max(cor))
	x=x[0]
	y=y[0]
	xcor_max = cor[y, x]
	tempx=(axis_compute[1]-axis_compute[0])*x*1./len(img1[0])+axis_compute[0]
	tempy=(axis_compute[2]-axis_compute[3])*y*1./len(img1)+axis_compute[3]
	new_axis2=[tempx,axis2[1]-axis2[0]+tempx,tempy-axis2[3]+axis2[2],tempy]
	
	#output
	def out_options(x,**kwargs):
		box_color=kwargs.get('box_color')
		if x=='xcor_map':
			out=cor
		elif x=='region':
			out=region
		elif x=='new_axis2':
			out=new_axis2
		elif x=='loc_map':
			loc_map=plt.figure()
			plt.imshow(img1_copy,extent=axis1)
			current_axis=plt.gca()
			current_axis.add_patch(Rectangle([new_axis2[0],new_axis2[2]],new_axis2[1]-new_axis2[0],new_axis2[3]-new_axis2[2],fill=False,color=box_color))
			out=loc_map
		elif x=='resampled_img1':
			out=img1
		elif x=='resampled_img2':
			out=img2
		elif x=='xcor_max':
			out=xcor_max
		elif x=='translation':
			out=new_axis2[0]-axis2[0], new_axis2[2]-axis2[2]
		return out
	if type(output) is str:
		return out_options(output,**kwargs)
	else:
		out=[]
		for i in output:
			out.append(out_options(i,**kwargs))
		return tuple(out)

def show_image(img1, axis1=None, img2=None, axis2=None, **kwargs):
	if axis1 is None and img1 is not None:
		axis1 = [0, len(img1[0]), 0, len(img1)]
	if axis2 is None and img2 is not None:
		axis2 = [0, len(img2[0]), 0, len(img2)]
	plt.figure(figsize=(8, 8))
	if img1 is not None:
		plt.imshow(img1,
					interpolation='none',
					extent=axis1,
					**kwargs)
	plt.imshow(img2,
				interpolation='none',
				extent=axis2,
					**kwargs)
	plt.xlabel('X-axis')
	plt.ylabel('Y-axis')
	plt.axis('equal')
	plt.axis('image')
	plt.colorbar()
	plt.grid()
	plt.show()