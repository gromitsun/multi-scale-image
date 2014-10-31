import h5py
import numpy as np
import glob
from sys import argv
from tools import *
try:
	script, output_filename, path= argv
except ValueError:
	path='./'
	output_filename='msi.h5'

def link_exchange_files(output_filename='linked_ex.h5', path='./'):
	"""Link all exchange HDF5 files. Output a single file with 'HDF SoftLink's.
	
	Keyword arguments:
	output_filename -- filename of the output file. (Default 'linked_ex2_msi.h5')
	path -- path to files to be linked. (Default './')
	"""
	filelist = glob.glob(path+'2xfm*.h*')
	f=h5py.File(output_filename,'w')
	for sub_file in filelist:
		f[sub_file+'/images']=h5py.ExternalLink(sub_file,'/exchange/images')
		f[sub_file+'/x_axis']=h5py.ExternalLink(sub_file,'/exchange/x_axis')
		f[sub_file+'/y_axis']=h5py.ExternalLink(sub_file,'/exchange/y_axis')
		f[sub_file+'/names']=h5py.ExternalLink(sub_file,'/exchange/images_names')
		f[sub_file+'/units']=h5py.ExternalLink(sub_file,'/exchange/images_units')
	f.close()

def exchange_to_ex2(output_filename='_ex2.h5', path='./', separate_files=True):
	"""Convert xfm hdf5 files with 'exchange' group to 'exchange2' hdf5 files.
	
	Keyword arguments:
	*path*:
		Path of the hdf5 files to be converted.
	*output_filename*:
		Filename of the output file. 
		If output are separate files, the value in this argument will be 
		appended to the end of the original filename.
	*separate_files*:
		If True, create one file for each original file.
		If False, create one single file for all files.
	"""
	filelist = glob.glob(path+'2xfm_*.h*')
	if not separate_files:
		f_msi=h5py.File(output_filename,'w')
	for sub_file in filelist:
		print 'Processing file \''+sub_file+'\'...'
		f_exchange=h5py.File(sub_file,'r')
		try:
			images=f_exchange['/exchange/images']
			names=f_exchange['/exchange/images_names']
			units=f_exchange['/exchange/images_units']
			x_axis=f_exchange['/exchange/x_axis']
			y_axis=f_exchange['/exchange/y_axis']
			n_xcoord=np.where(np.array(names)=='x_coord')[0][0]
			n_ycoord=np.where(np.array(names)=='y_coord')[0][0]
		except KeyError:
			print 'Error: Not a standard exchange file! File skipped.'
			continue
		if separate_files:
			sub_file=sub_file[:sub_file.find('.h')]
			f_msi=h5py.File(sub_file+output_filename,'w')
			group=f_msi.create_group(sub_file)
		else:
			group=f_msi.create_group('image_set/'+sub_file)
		group.attrs['axis']=[x_axis[0],x_axis[-1],y_axis[-1],y_axis[0]]
		group.attrs['axis_unit']=[units[n_xcoord],units[n_ycoord]]
		for n in range(len(images)):
			dset=group.create_dataset(names[n],(len(images[n]),len(images[n][0]),),dtype='float')
			dset[...]=images[n]
			dset.attrs['name']=names[n]
			dset.attrs['unit']=units[n]
		f_exchange.close()
		if separate_files:
			f_msi.close()
		print 'Done!'
	if not separate_files:
		f_msi.close()
	print 'All done!'

def link_ex2_files(output_filename='linked_ex2_msi.h5', path='./'):
	"""Link exchange2 HDF5 files. Output a single file with 'HDF SoftLink's.
	
	ex2 files -> msi files
	Keyword arguments:
	output_filename -- filename of the output file. (Default 'linked_ex2_msi.h5')
	path -- path to files to be linked. (Default './')
	"""
	filelist = glob.glob(path+'2xfm*ex2.h*')
	f=h5py.File(output_filename,'w')
	for sub_file in filelist:
		sub_file=sub_file[:sub_file.find('_ex2.h5')]
		f['image_set/'+sub_file]=h5py.ExternalLink(sub_file+'_ex2.h5',sub_file)
	f.close()
	
def msi_add_sorted(filename='linked_ex2_msi.h5'):
	"""Add 'image_sorted' group to linked ex2 files (msi files).
	
	Put sort images in group 'image_set' into different 'layers' according to 
	their pixel size. Layers numbered starting from 0, and from largest pixel 
	size to lowest.

	Keyword arguments:
	filename -- filename of the input file. (Default 'linked_ex2_msi.h5')
	
	Notes:
	Pixel size of each image in group 'image_set' is calculated by dividing 
	the extent of 'axis' in the x direction (first 2 values of 'axis') by 
	the number of pixels in the same direction;
	Pixel sizes are rounded to the first effective digit before sorting.
	"""
	f=h5py.File(filename,'r+')
	dic={}
	i=0
	for image in f['image_set']:
		i+=1
		pixel_size=abs(f['image_set'][image].attrs['axis'][1]-f['image_set'][image].attrs['axis'][0])*1./len(f['image_set'][image]['x_coord'][0])
		if i==1:
			unit=f['image_set'][image].attrs['axis_unit'][0]
		else:
			unit_temp=f['image_set'][image].attrs['axis_unit'][0]
			if unit != unit_temp:
				pixel_size=convert_unit(pixel_size, unit_temp, unit)
		pixel_size=roundn(pixel_size)
		if pixel_size in dic.keys():
			dic[pixel_size].append(image)
		else:
			dic[pixel_size]=[image]
	dic_sorted=sorted(dic.items())[::-1]
	for n in range(len(dic)):
		layer=f.create_group('image_sorted/layer'+str(n))
		layer.attrs['pixel_size']=dic_sorted[n][0]
		layer.attrs['pixel_size_unit']=unit
		for image in dic_sorted[n][1]:
			layer[image]=h5py.SoftLink('/image_set/'+image)
			layer[image].attrs['layer']=n
	f.close()
		

#link_exchange_files(output_filename)
#link_msi_files(output_filename,path)
#exchange_to_msi(path=path,separate_files=False)