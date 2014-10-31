import glob
from MultiImgs import *

def _exchange_to_mapsimg(exchange_group, element = None, filename = None):
	if filename == None:
		filename = exchange_group.name
	name=np.array(exchange_group['images_names'])
	xaxis=list(exchange_group['x_axis'])[::1]
	yaxis=list(exchange_group['y_axis'])[::-1]
	n=np.where(name=='x_coord')[0][0]
	xaxis_unit=exchange_group['images_units'][n]
	n=np.where(name=='y_coord')[0][0]
	yaxis_unit=exchange_group['images_units'][n]
	if element!=None:
		n=np.where(name==element)[0][0]
		image=exchange_group['images'][n]
		name = filename + ' - ' + name[n]
		unit=exchange_group['images_units'][n]
	else:
		image=list(exchange_group['images'])
		unit=list(exchange_group['images_units'])
		name = filename
	return mapsimg(image,name,unit,xaxis,yaxis,xaxis_unit,yaxis_unit)

def readh5(filename,element=None):
	f=h5py.File(filename,'r')
	name=np.array(f['exchange']['images_names'])
	xaxis=list(f['exchange']['x_axis'])[::1]
	yaxis=list(f['exchange']['y_axis'])[::-1]
	n=np.where(name=='x_coord')[0][0]
	xaxis_unit=f['exchange']['images_units'][n]
	n=np.where(name=='y_coord')[0][0]
	yaxis_unit=f['exchange']['images_units'][n]
	if element!=None:
		n=np.where(name==element)[0][0]
		image=f['exchange']['images'][n]
		name=name[n]
		unit=f['exchange']['images_units'][n]
	else:
		image=list(f['exchange']['images'])
		unit=list(f['exchange']['images_units'])
	f.close()
	return image,name,unit,xaxis,yaxis,xaxis_unit,yaxis_unit

def loadh5(filename, element=None, path = './'):
	f=h5py.File(path + filename,'r')
	name=np.array(f['exchange']['images_names'])
	xaxis=list(f['exchange']['x_axis'])[::1]
	yaxis=list(f['exchange']['y_axis'])[::-1]
	n=np.where(name=='x_coord')[0][0]
	xaxis_unit=f['exchange']['images_units'][n]
	n=np.where(name=='y_coord')[0][0]
	yaxis_unit=f['exchange']['images_units'][n]
	if element!=None:
		n=np.where(name==element)[0][0]
		image=f['exchange']['images'][n]
		name = filename + ' - ' + name[n]
		unit=f['exchange']['images_units'][n]
	else:
		image=list(f['exchange']['images'])
		unit=list(f['exchange']['images_units'])
		name = filename
	f.close()
	return mapsimg(image,name,unit,xaxis,yaxis,xaxis_unit,yaxis_unit)
   
def load_all(element, path = './', simplify_filename = True):
	filelist = glob.glob(path + '*.h*')
	datadict={}
	for f in filelist:
		if simplify_filename:
			n = f.find('_', len(path))
			num = f[n+1:n+5]
		else:
			num = f[:]
		print 'Loading %s ...' % f
		f = f[len(path):]
		try:
			datadict[num] = loadh5(f, element=element, path = path)
		except KeyError:
			print 'File type error. File skipped.'
	print 'Successfully loaded %d files.' % len(datadict)
	return datadict

def sort(datadict, return_pixel_size_unit = False):
	dict_sorted = {}
	i = 0
	for num, data in datadict.items():
		i += 1
		_pixel_size = data.pixel_size
		if i == 1:
			unit = data.xaxis_unit
		else:
			unit_temp = data.xaxis_unit
			if unit != unit_temp:
				_pixel_size = convert_unit(_pixel_size, unit_temp, unit)
			_pixel_size = roundn(_pixel_size, 2)
			if _pixel_size in dict_sorted.keys():
				dict_sorted[_pixel_size].append(data)
			else:
				dict_sorted[_pixel_size] = [data]
	if return_pixel_size_unit:
		return dict_sorted, unit
	else:
		return dict_sorted

def write_link_file(dict_sorted, pixel_size_unit = 'mm', datapath = './', filename = 'link_file.h5', element = None):
	f=h5py.File(filename, 'w')
	n = 0
	try:
		for pixel_size, dataset in sorted(dict_sorted.items())[::-1]:
			level = f.create_group('Level_' + str(n))
			level.attrs['pixel_size'] = int_out(float(str(pixel_size)))
			level.attrs['pixel_size_unit'] = pixel_size_unit
			for data in dataset:
				_name = data.name
				_n = _name.find(' - ')
				if _n != -1:
					_name = _name[:_n]
				level[_name] = h5py.ExternalLink(datapath + _name, '/exchange')
				level[_name].attrs['name'] = _name
				level[_name].attrs['axis'] = data.axis
				level[_name].attrs['axis_unit'] = data.xaxis_unit ###
			n += 1
	finally:
		f.close()

def load_link_file(filename, path = './', element = None, return_pixel_size_unit = False):
	f = h5py.File(path + filename, 'r')
	dict_sorted = {}
	try:
		for level in f.values():
			if 'Level_' in level.name:
				level_num = int(level.name[7:])
				dict_sorted[level_num] = [_exchange_to_mapsimg(level[filename], element, filename) for filename in level.keys()]
				pixel_size_unit = level.attrs['pixel_size_unit']
		if return_pixel_size_unit:
			return dict_sorted, pixel_size_unit
		else:
			return dict_sorted
	finally:
		f.close()
	
def load(element = None, datapath = './', linkfilepath = './', link_file_name = 'link_file.h5', from_link_file = True, create_link_file = True):
	
	if from_link_file:
		try:
			dict_sorted, pixel_size_unit = load_link_file(link_file_name, path = linkfilepath, element = element, return_pixel_size_unit = True)
		except IOError:
			from_link_file = False
	
	if not from_link_file:
		dict_sorted, pixel_size_unit = sort(load_all(element, path = datapath), return_pixel_size_unit = True)
	
	m = MultiImgs()
	
	for level_num, dataset in dict_sorted.items():
		for data in dataset:
			m.add_image(MultiImg(data, level_num))
	if create_link_file and not from_link_file:
		write_link_file(dict_sorted, pixel_size_unit = pixel_size_unit, datapath = datapath, filename = linkfilepath + link_file_name, element = element)
	return m

def print_axis(datadict):
	for num, data in sorted(datadict.items()):
		print num, data.axis

def print_pixels(datadict):
	for num, data in sorted(datadict.items()):
		print num, data.image.shape[::-1]

def print_pixel_size(datadict):
	for num, data in sorted(datadict.items()):
		print num, data.pixel_size
