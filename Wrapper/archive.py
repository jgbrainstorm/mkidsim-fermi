# Copyright 2012 ETHZ.ch Lukas Gamper, gamperl@gmail.com
import numpy as np
import tables as tbl
from collections import Set

class nodeList(Set):

	def __init__(self, nodes):
		self.elements = [node._v_name for node in nodes]

	def __str__(self):
		return str(self.elements)

	def __iter__(self):
		return iter(self.elements)

	def __contains__(self, data):
         return data in self.elements

	def __len__(self):
		return len(self.elements)

class archive:

	def __init__(self, filename, mode = 'r'):
		self.file = tbl.openFile(filename, mode)

	@property
	def isopen(self):
		return self.file.isopen == 1

	def close(self):
		self.file.close()

	def allocate(self, path, dtype, shape):
		segments = self.__splitPath(path)
		group = self.file.root
		for segment in segments[:-1]:
			try:
				group = self.file.getNode(group, name=segment, classname='Group')
			except tbl.exceptions.NoSuchNodeError:
				group = self.file.createGroup(group, segment)
		return self.file.createCArray(group, segments[-1], atom=tbl.Atom.from_dtype(dtype), shape=shape)

	def __enter__(self):
		if not self.isopen:
			raise tbl.exceptions.ClosedFileError('HDF5 file is closed')
		return self

	def __exit__(self, type, data, traceback):
		self.close()

	def __setitem__(self, arg, data):
		segments = self.__splitPath(arg[0] if isinstance(arg, tuple) else arg)
		group = self.file.root
		for segment in segments[:-1]:
			try:
				group = self.file.getNode(group, name=segment, classname='Group')
			except tbl.exceptions.NoSuchNodeError:
				group = self.file.createGroup(group, segment)
		if isinstance(arg, tuple):
			self.file.getNode(group, name=segments[-1], classname='CArray')[arg[1:]] = data
		else:
			try:
				self.file.createArray(group, segments[-1], data)
			except tbl.exceptions.NodeError:
				self.file.removeNode(group, name=segments[-1], recursive=True)
				self.file.createArray(group, segments[-1], data)

	def __getitem__(self, arg):
		path = arg[0] if isinstance(arg, tuple) else arg
		if path == '/':
			return nodeList(self.file.listNodes(self.file.root))
		segments = self.__splitPath(path)
		group = self.file.root
		for segment in segments[:-1]:
			group = self.file.getNode(group, name=segment, classname='Group')
		try:
			if isinstance(arg, tuple):
				return self.file.getNode(group, name=segments[-1], classname='CArray')[arg[1:]]
			else:
				return self.file.getNode(group, name=segments[-1], classname='Array').read()
		except tbl.exceptions.NoSuchNodeError:
			return nodeList(self.file.listNodes(self.file.getNode(group, name=segments[-1], classname='Group')))

	def __splitPath(self, path):
		segments = path.split('/')
		if len(segments[0]) > 0 or min([len(segment) for segment in segments[1:]]) == 0:
			raise tbl.exceptions.NodeError('Invalid path ' + path)
		return segments[1:]

# with archive('test.h5', 'w') as ar:
# 	ar['/path/to/numpy/data'] = np.random.random((2,2,2))
# 	ar['/path/to/i'] = 1
# 	ar['/path/to/c'] = 5 + 1j
# 	ar['/path/to/s'] = 'test'

# 	i = ar['/path/to/i']
# 	s = ar['/path/to/s']
# 	c = ar['/path/to/c']
# 	n = ar['/path/to/numpy/data']
# 	p = ar['/path']

# 	print i, type(i)
# 	print s, type(s)
# 	print c, type(c)
# 	print n, type(n)
# 	print p, type(p)

# 	print ar['/']
# 	print ar['/path']
# 	print ar['/path/to']

# 	#check if path is group:
# 	print isinstance(ar['/'], nodeList)
# 	print isinstance(ar['/path/to/c'], nodeList)

# 	print ar.isopen

# 	ar['/path/to/i'] = 1.5
# 	i = ar['/path/to/i']
# 	print i, type(i)

# 	#write hyperslab
# 	data = np.ones((2, 2));
# 	ar.allocate('/hyperslab', dtype=data.dtype, shape=(3,3))
# 	ar['/hyperslab', 1:3, 1:3] = data
# 	ar['/hyperslab', 0:3, 1] = np.ones((1, 3))

# 	#read hyperslab
# 	print ar['/hyperslab', 0:2, 1:3]
# 	print ar['/hyperslab']
