"""
**update**
1. update zeropadding()
2. add automatic zeropadding into xcorr and xcor
3. add infmax and infmin
4. add method keyword in zeropadding, can do corner/center zeropadding.
5. zeropadding can be done by providing axis and new_axis.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import gaussian
from scipy.ndimage.filters import _ni_support


##math tools

def int_out(x):
	"""If x is an integer, convert x to int."""
	if x.is_integer():
		x=int(x)
	return x

def roundn(x,n=1,**kwargs): 
	"""Round x to the nth effective digit."""
	out=round(x,-int(np.floor(np.log10(x)))+n-1)
	if kwargs.get('int_output')==True:
		out=int_out(out)
	return out
	


def find_nearest(x,lib=[1,2,5]): 
	"""round x to the nearest number in lib"""
	return lib[np.abs(np.array(lib)-x).argmin()]

def search(x,lib=[1,2,5],period=10):
	"""find the nearest number in lib*period**n for all integer n"""
	p=int(np.floor(np.log(x)/np.log(period)))
	return find_nearest(x*period**(-p),lib)*period**p
	
def convert_unit(value,unit='mm',new_unit='nm'):
	"""value (in unit) -> new_value (in new_unit)"""
	unitlib={'km':1000.,'m':1.,'dm':0.1,'cm':0.01,'mm':0.001,'um':10**-6,'nm':10**-9,'$\AA$':10**-10}
	return value*unitlib[unit]/unitlib[new_unit]

def inverse_dict(dic):
	return {value:key for key,value in dic.items()}

def auto_unit(value,unit='mm',**kwargs):
	"""value, unit -> appropriate new_value, new_unit"""
	unitlib={'km':1000.,'m':1.,'dm':0.1,'cm':0.01,'mm':0.001,'um':10**-6,'nm':10**-9,'$\AA$':10**-10}
	ivd=inverse_dict(unitlib)
	new_unit=ivd[search(unitlib[unit]*value,lib=[1],period=1000)]
	new_value=int_out(float(str(convert_unit(value,unit,new_unit))))
	return new_value,new_unit


def nanmean(x):
	return np.nansum(x)*1./np.sum(np.isfinite(x))

def nan_to_value(x,value='average'):
	if value=='average':
		value=nanmean(x)
	x=np.array(x)
	x[np.isnan(x)]=value
	return x

def infmax(arr):
	x, y = np.where(arr == np.inf)
	arr[x, y] = -np.inf
	return np.max(arr)

def infmin(arr):
	x, y = np.where(arr == np.inf)
	arr[x, y] = np.inf
	return np.min(arr)

def infnanmax(arr):
	x, y = np.where(arr == np.inf)
	arr[x, y] = -np.inf
	return np.nanmax(arr)

def infnanmin(arr):
	x, y = np.where(arr == np.inf)
	arr[x, y] = np.inf
	return np.nanmin(arr)

def simplify_array(arr):
	for x in range(len(arr)):
		arr[x] = float(str(arr[x]))
	return arr

def gaussian_func(x, sigma):
	return np.exp(-0.5*(float(x)/sigma)**2)

def gaussian_kernel1d(sigma, axis=-1, order=0, output=None,
					  mode="reflect", cval=0.0):
	"""One-dimensional Gaussian filter.

	Parameters
	----------
	%(input)s
	sigma : scalar
		standard deviation for Gaussian kernel
	%(axis)s
	order : {0, 1, 2, 3}, optional
		An order of 0 corresponds to convolution with a Gaussian
		kernel. An order of 1, 2, or 3 corresponds to convolution with
		the first, second or third derivatives of a Gaussian. Higher
		order derivatives are not implemented
	%(output)s
	%(mode)s
	%(cval)s

	Returns
	-------
	gaussian_filter1d : ndarray

	"""
	if order not in range(4):
		raise ValueError('Order outside 0..3 not implemented')
	sd = float(sigma)
	# make the length of the filter equal to 4 times the standard
	# deviations:
	lw = int(4.0 * sd + 0.5)
	weights = [0.0] * (2 * lw + 1)
	weights[lw] = 1.0
	sum = 1.0
	sd = sd * sd
	# calculate the kernel:
	for ii in range(1, lw + 1):
		tmp = np.exp(-0.5 * float(ii * ii) / sd)
		weights[lw + ii] = tmp
		weights[lw - ii] = tmp
		sum += 2.0 * tmp
	for ii in range(2 * lw + 1):
		weights[ii] /= sum
	# implement first, second and third order derivatives:
	if order == 1:  # first derivative
		weights[lw] = 0.0
		for ii in range(1, lw + 1):
			x = float(ii)
			tmp = -x / sd * weights[lw + ii]
			weights[lw + ii] = -tmp
			weights[lw - ii] = tmp
	elif order == 2:  # second derivative
		weights[lw] *= -1.0 / sd
		for ii in range(1, lw + 1):
			x = float(ii)
			tmp = (x * x / sd - 1.0) * weights[lw + ii] / sd
			weights[lw + ii] = tmp
			weights[lw - ii] = tmp
	elif order == 3:  # third derivative
		weights[lw] = 0.0
		sd2 = sd * sd
		for ii in range(1, lw + 1):
			x = float(ii)
			tmp = (3.0 - x * x / sd) * x * weights[lw + ii] / sd2
			weights[lw + ii] = -tmp
			weights[lw - ii] = tmp
	return weights


	
	
#image tools

def tukey(alpha,N):
	w=np.empty(N)
	for n in range(N):
		if 0<=n<alpha*(N-1)/2.:
			w[n]=1/2.*(1+np.cos(np.pi*(2*n/alpha/(N-1.)-1)))
		elif alpha*(N-1)/2.<=n<=(N-1)*(1-alpha/2.):
			w[n]=1
		else:
			w[n]=1/2.*(1+np.cos(np.pi*(2*n/alpha/(N-1.)-2./alpha+1)))
	return w


def tukey2d(alpha,M,N):
	x=tukey(alpha,N)
	y=tukey(alpha,M)
	X,Y=np.meshgrid(x,y)
	return X*Y

def xcorr(x,y,**kwargs):
	"""cross correlation by rfft"""
	x = np.asarray(x)
	y = np.asarray(y)
	if np.ndim(x) == np.ndim(y):
		shape=kwargs.get('shape',np.max((x.shape, y.shape), axis = 0))
		return np.fft.irfftn(np.conjugate(np.fft.rfftn(x,s=shape))*np.fft.rfftn(y,s=shape))
	elif np.ndim(y) == 1:
		axis = kwargs.get('axis', 0)
		shape=kwargs.get('shape', max(x.shape[axis], len(y)))
		shape+=shape%2
		outshape = np.array(x.shape[:])
		outshape[axis] = shape
		out = np.zeros(outshape)
		y = np.fft.ifftshift(np.pad(y, pad_width = (int((shape-len(y)+1)/2), int((shape-len(y))/2)), mode = 'constant'))
		y_fft = np.fft.rfft(y, n=shape)
		x_fft = np.fft.rfft(x, n=shape, axis=axis)
		if axis == 0:
			for ii in range(len(x_fft[0])):
				out[:,ii] = np.fft.irfft(x_fft[:,ii]*np.conjugate(y_fft))
		else:
			for ii in range(len(x_fft)):
				out[ii] = np.fft.irfft(x_fft[ii]*np.conjugate(y_fft))
		return out
	else:
		raise ValueError('Only inputs with dimensions of 1 or 2 can be processed.')

def xcor(x,y,**kwargs):
	"""cross correlation by fft"""
	if np.ndim(x) == 1 and np.ndim(y) == 1:
		shape=kwargs.get('shape',max(len(x),len(y)))
		return np.fft.ifft(np.conjugate(np.fft.fft(x,s=shape))*np.fft.fft(y,s=shape))
	elif np.ndim(x) == 2 and np.ndim(y) == 2:
		shape=kwargs.get('shape',[max(len(x),len(y)),max(len(x[0]),len(y[0]))])
		return np.fft.ifft2(np.conjugate(np.fft.fft2(x,s=shape))*np.fft.fft2(y,s=shape))
	else:
		raise ValueError('Only inputs with dimensions of 1 or 2 can be processed.')


def fftfreq2d(m,n=None,d=1,zero='center',axis=0):
	x=np.fft.fftfreq(m,d)
	if n==None: n=m
	y=np.fft.fftfreq(n,d)
	X,Y=np.meshgrid(x,y)
	if axis==0: result=X
	else: result=Y
	if zero=='center':
		result=np.fft.fftshift(result,1-axis)
	return result

def shiftimg_fft(im_fft,shifts,zero='front'):
	"""fft of image -> shifted fft of image by x=shifts[0], y=shifts[1]"""
	im_fft_shift=im_fft[:]
	for n in range(len(shifts)):
		im_fft_shift=im_fft_shift*np.exp(1j*fftfreq2d(np.shape(im_fft)[1],np.shape(im_fft)[0],axis=n,zero=zero)*2*np.pi*(-shifts[n]))
	return im_fft_shift
	
def shiftimg(im,shifts,**kwargs):
	"""image -> shifted image by x=shifts[0], y=shifts[1]"""
	im_fft=np.fft.fft2(im)
	im_fft_shift=shiftimg_fft(im_fft,shifts,**kwargs)
	return np.real(np.fft.ifft2(im_fft_shift))

def zeropadding(im,new_shape=None,ratio=None, axis=None, new_axis=None, **kwargs):
	"""Pad ndarray im with zeros arround the original ndarray.
	Shape of padded image indicated in 'new_shape'.
	Input the ratio of new image to original (new_shape/np.shape(im)) through ratio is also acceptable."""
	actual_axis = kwargs.get('actual_axis')
	if ratio or new_shape:
		method = kwargs.get('method', 'center')
		if ratio:
			new_shape=np.round(np.array(np.shape(im))*np.array(ratio))
		im_padded=np.zeros(new_shape)
		if np.iscomplex(im).any():
			im_padded=im_padded*complex(1)
		if method == 'center':
			x1 = -(new_shape[1]-len(im[0])+(new_shape[1]-len(im[0]))%2)/2
			x2 = x1 + new_shape[1]
			y1 = -(new_shape[0]-len(im)+(new_shape[0]-len(im))%2)/2
			y2 = y1 + new_shape[0]	
		else:
			x1 = 0
			x2 = new_shape[1]
			y1 = 0
			y2 = new_shape[0]
	elif axis and new_axis:
		x1=round((new_axis[0]-axis[0])*1./(axis[1]-axis[0])*len(im[0]))
		x2=round((new_axis[1]-axis[0])*1./(axis[1]-axis[0])*len(im[0]))
		y1=round((new_axis[3]-axis[3])*1./(axis[2]-axis[3])*len(im))
		y2=round((new_axis[2]-axis[3])*1./(axis[2]-axis[3])*len(im))
		im_padded = np.zeros([y2 - y1, x2 - x1])
	else:
		raise TypeError("Either 'new_shape' or 'ratio' or ('axis' and 'new_axis') must be given")
	
	try:
		im_padded[-y1:-y1+len(im),-x1:-x1+len(im[0])] = im[:]
	except IndexError or ValueError:
		raise ValueError('axis should be included in new_axis')
	
	if actual_axis and (axis is not None):
		actual_x1=(axis[1]-axis[0])*x1*1./len(im[0])+axis[0]
		actual_x2=(axis[1]-axis[0])*x2*1./len(im[0])+axis[0]
		actual_y2=(axis[2]-axis[3])*y1*1./len(im)+axis[3]
		actual_y1=(axis[2]-axis[3])*y2*1./len(im)+axis[3]
		return im_padded, [actual_x1, actual_x2, actual_y1, actual_y2]
	else:
		return im_padded

def crop(image,**kwargs):
	new_shape=kwargs.get('new_shape')
	scale=kwargs.get('scale')
	axis=kwargs.get('axis',[0,len(image[0]),0,len(image)])
	new_axis=kwargs.get('new_axis')
	slice_pixel = kwargs.get('slice_pixel')
	actual_axis = kwargs.get('actual_axis')
	if scale or new_shape:
		if scale:
			new_shape=np.round(np.array(np.shape(image))*np.array(scale))
		x1 = int((len(image[0])-new_shape[1])/2)
		x2 = x1+new_shape[1]
		y1 = int((len(image)-new_shape[0])/2)
		y2 = y1+new_shape[0]
	elif slice_pixel:
		y1, y2, x1, x2 = tuple(slice_pixel)
	elif axis and new_axis:
		actual_axis = kwargs.get('actual_axis')
		x1=round((new_axis[0]-axis[0])*1./(axis[1]-axis[0])*len(image[0]))
		x2=round((new_axis[1]-axis[0])*1./(axis[1]-axis[0])*len(image[0]))
		y1=round((new_axis[3]-axis[3])*1./(axis[2]-axis[3])*len(image))
		y2=round((new_axis[2]-axis[3])*1./(axis[2]-axis[3])*len(image))
	else:
		raise ValueError("'new_shape' or 'scale' or 'new_axis' must be given")
	
	if actual_axis and (axis is not None):
		actual_x1=(axis[1]-axis[0])*x1*1./len(image[0])+axis[0]
		actual_x2=(axis[1]-axis[0])*x2*1./len(image[0])+axis[0]
		actual_y2=(axis[2]-axis[3])*y1*1./len(image)+axis[3]
		actual_y1=(axis[2]-axis[3])*y2*1./len(image)+axis[3]
		return image[y1:y2,x1:x2], [actual_x1, actual_x2, actual_y1, actual_y2]
	else:
		return image[y1:y2,x1:x2]
		

def fft_resample(image, scale, **kwargs):
	"""resample 'image' to 'scale' using fft"""
	conserve_min = kwargs.get('conserve_min')
	conserve_range = kwargs.get('conserve_range')
	rescale_range = kwargs.get('rescale_range')
	zoomin_filter = kwargs.get('zoomin_filter')
	sigma = kwargs.get('gaussian_sigma', 1./scale)
	
	im_fft=np.fft.fft2(image)
	if zoomin_filter != 'gaussian':
		im_fft=np.fft.fftshift(im_fft,axes=[0,1])
	
	if scale>1:
		im_fft=zeropadding(im_fft,ratio=scale)
	elif scale<1:
		if zoomin_filter == 'gaussian':
			#_x = np.array(gaussian(len(im_fft[0]), sigma), ndmin = 2)
			#_y = np.array(gaussian(len(im_fft), sigma), ndmin = 2)
			#_gaussian = _x*_y.transpose()
			_temp = np.array(gaussian_kernel1d(sigma), ndmin = 2)
			_gaussian = _temp.transpose()*_temp
			_gaussian = zeropadding(_gaussian, new_shape = im_fft.shape)
			_gaussian = np.fft.fft2(np.fft.ifftshift(_gaussian))
			im_fft = np.fft.fftshift(im_fft*_gaussian)
		im_fft=crop(im_fft,scale=scale)
	im_fft=np.fft.ifftshift(im_fft,axes=[0,1])
	out = np.real(np.fft.ifft2(im_fft))
	if conserve_min:
		out -= np.min(out)
	if conserve_range:
		out *= (np.max(image) - np.min(image))*1./(np.max(out)-np.min(out))
	if rescale_range:
		out *= scale
	return out

def overlay_display(img1, img2, axis1 = None, axis2 = None, translation = None, remove_background = True, normalize = True, **kwargs):
	if remove_background:
		img1 = img1 - np.min(img1)
		img2 = img2 - np.min(img2)
	if normalize:
		img1 = img1/np.max(img1)
		img2 = img2/np.max(img2)
	if (axis1 is not None) and (axis2 is not None):
		translation = [0, 0]
		translation[0] = round((axis2[0]-axis1[0])*1./(axis1[1]-axis1[0])*len(img1[0]))
		translation[1] = round((axis2[3]-axis1[3])*1./(axis1[2]-axis1[3])*len(img1))
	x_start = min(0, translation[0])
	y_start = min(0, translation[1])
	x_end = max(np.shape(img1)[1], np.shape(img2)[1] + translation[0])
	y_end = max(np.shape(img1)[0], np.shape(img2)[0] + translation[1])
	disp = np.zeros([y_end - y_start, x_end - x_start, 3])
	disp[-y_start:-y_start+np.shape(img1)[0], -x_start:-x_start+np.shape(img1)[1], 0] = img1[:]
	disp[-y_start + translation[1]:-y_start + translation[1] +np.shape(img2)[0], -x_start + translation[0]:-x_start + translation[0]+np.shape(img2)[1], 1] = img2[:]
	plt.figure()
	plt.imshow(disp, **kwargs)
	plt.show()

def gaussian_filter1d(input, sigma, axis=-1, order=0, output=None,
					  mode="reflect", cval=0.0):
	"""One-dimensional Gaussian filter.

	Parameters
	----------
	%(input)s
	sigma : scalar
		standard deviation for Gaussian kernel
	%(axis)s
	order : {0, 1, 2, 3}, optional
		An order of 0 corresponds to convolution with a Gaussian
		kernel. An order of 1, 2, or 3 corresponds to convolution with
		the first, second or third derivatives of a Gaussian. Higher
		order derivatives are not implemented
	%(output)s
	%(mode)s
	%(cval)s

	Returns
	-------
	gaussian_filter1d : ndarray

	"""
	#weights = gaussian(sigma)
	weights = gaussian_kernel1d(sigma, axis=-1, order=0, output=None,
					  							mode="reflect", cval=0.0)
	return xcorr(input, weights, axis=axis)

def gaussian_filter(input, sigma, order=0, output=None,
				  mode="reflect", cval=0.0):
	"""Multidimensional Gaussian filter.

	Parameters
	----------
	%(input)s
	sigma : scalar or sequence of scalars
		Standard deviation for Gaussian kernel. The standard
		deviations of the Gaussian filter are given for each axis as a
		sequence, or as a single number, in which case it is equal for
		all axes.
	order : {0, 1, 2, 3} or sequence from same set, optional
		The order of the filter along each axis is given as a sequence
		of integers, or as a single number.  An order of 0 corresponds
		to convolution with a Gaussian kernel. An order of 1, 2, or 3
		corresponds to convolution with the first, second or third
		derivatives of a Gaussian. Higher order derivatives are not
		implemented
	%(output)s
	%(mode)s
	%(cval)s

	Returns
	-------
	gaussian_filter : ndarray
		Returned array of same shape as `input`.

	Notes
	-----
	The multidimensional filter is implemented as a sequence of
	one-dimensional convolution filters. The intermediate arrays are
	stored in the same data type as the output. Therefore, for output
	types with a limited precision, the results may be imprecise
	because intermediate results may be stored with insufficient
	precision.

	"""
	input = np.asarray(input)
	"""
	output, return_value = _ni_support._get_output(output, input)
	orders = _ni_support._normalize_sequence(order, input.ndim)
	if not set(orders).issubset(set(range(4))):
		raise ValueError('Order outside 0..4 not implemented')
	sigmas = _ni_support._normalize_sequence(sigma, input.ndim)
	"""
	
	axes = list(range(input.ndim))
	if type(sigma) is int:
		sigmas = np.ones(len(axes))*sigma
	else:
		sigmas = sigma[:]
	if type(order) is int:
		orders = np.ones(len(axes))*order
	else:
		orders = order[:]
	axes = [(axes[ii], sigmas[ii], orders[ii])
						for ii in range(len(axes)) if sigmas[ii] > 1e-15]
	if len(axes) > 0:
		for axis, sigma, order in axes:
			output = gaussian_filter1d(input, sigma, axis, order, output,
							  mode, cval)
			input = output
	else:
		output[...] = input[...]
	return output

def zoom(img, zoom):
	return fft_resample(img, scale = zoom, zoomin_filter = 'gaussian', conserve_min=True, rescale_range = True)