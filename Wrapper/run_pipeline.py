import os
import sys
import time
from archive import nodeList
import archive
import tests
import numpy as np
import string
import re
import argparse
import ConfigParser
import shutil

import psutil
import subprocess

def copyHDF5Path(source, path, dest):
	for file in source[path]:
		if (path != "/"):
			file = "/" + str(file)
		if isinstance(source[path + file], nodeList) == False:
			dest[path + file] = source[path + file]
		else:
			copyHDF5Path(source, path + file, dest)

class ProcessTimer:
  def __init__(self,command):
    self.command = command
    self.execution_state = False

  def execute(self):
    self.max_vms_memory = 0
    self.max_rss_memory = 0

    self.t1 = None
    self.t0 = time.time()
    self.p = subprocess.Popen(self.command,shell=True)
    self.execution_state = True

  def poll(self):
    if not self.check_execution_state():
      return False

    self.t1 = time.time()

    try:
      pp = psutil.Process(self.p.pid)

      #obtain a list of the subprocess and all its descendants
      descendants = list(pp.get_children(recursive=True))
      descendants = descendants + [pp]

      rss_memory = 0
      vms_memory = 0

      #calculate and sum up the memory of the subprocess and all its descendants 
      for descendant in descendants:
        try:
          mem_info = descendant.get_memory_info()

          rss_memory += mem_info[0]
          vms_memory += mem_info[1]
        except psutil.error.NoSuchProcess:
          #sometimes a subprocess descendant will have terminated between the time
          # we obtain a list of descendants, and the time we actually poll this
          # descendant's memory usage.
          pass
      self.max_vms_memory = max(self.max_vms_memory,vms_memory)
      self.max_rss_memory = max(self.max_rss_memory,rss_memory)

    except psutil.error.NoSuchProcess:
      return self.check_execution_state()

    return self.check_execution_state()

  def is_running(self):
    return psutil.pid_exists(self.p.pid) and self.p.poll() == None
  def check_execution_state(self):
    if not self.execution_state:
      return False
    if self.is_running():
      return True
    self.executation_state = False
    self.t1 = time.time()
    return False

  def close(self,kill=False):

    try:
      pp = psutil.Process(self.p.pid)
      if kill:
        pp.kill()
      else:
        pp.terminate()
    except psutil.error.NoSuchProcess:
      pass

class SpokesConfig:
	def __init__ (self, config, commandline, cwd):
		self.config = config
		self.cwd = cwd
		self.commandline = commandline		
	
	def mergeDatafile(self, datafile):
		#merge datafiles into the databank
		with archive.archive(self.config.get("options", "data_dir") + self.config.get("options", "data_bank"), 'a') as dest:
			with archive.archive(datafile, 'r') as src:
				print "... copying " + self.config.get("data_bank", datafile) + " into " + self.config.get("options", "data_bank")
				copyHDF5Path(src, "/", dest)
	
	def mergeParamfile(self, file):
		param = ConfigParser.RawConfigParser()
		param.read(file)
		with archive.archive('../../data/data_bank.h5', 'a') as ar:
			print "Importing sections:", param.sections(), "from parameter file:", file
			for section in param.sections():
				for option in param.options(section):
					ar["/" + section + "/" + option] = eval(param.get(section, option))
	
	def mergeGlobalDatafiles(self):
		#merge datafiles into the databank
		for datafile in self.config.options("data_bank"):
			self.mergeDatafile(datafile)

	def mergeGlobalParamfile(self):
		if self.config.has_option('options','param'):
			if os.path.exists(self.config.get("options", "base_dir") + self.config.get('options','param')) == False:
				raise Exception('Global param file not found: ' + self.config.get("options", "base_dir") + self.config.get('options','param'))
				
			self.mergeParamfile(self.config.get("options", "base_dir") + self.config.get('options','param'))

	def mergeCmdLineParam(self):
		if (self.commandline.param != None):
			print "Merging command line param file"
			os.chdir(self.cwd)
			if (os.path.exists(self.commandline.param)):
				print "... Reading param file " + self.commandline.param
				self.mergeParamfile(self.commandline.param)
			else:
				raise Exception('Command line parameter file not found')

	def mergeModParam(self, module, paramfile, datafile):
		os.chdir(self.cwd)
	
		#h5 datafile
		if datafile != None:
			if os.path.exists(datafile) == False:
				raise Exception('Module data file not found: ' + datafile + ' (' + module + ')')
			print "... Reading data file for " + module
			self.mergeDatafile(datafile)
		
		#parameterfile
		if paramfile != None:
			if os.path.exists(paramfile) == False:
				raise Exception('Module param file not found: ' + paramfile + ' (' + module + ')')
			print "... Reading param file for " + module
			self.mergeParamfile(paramfile)
	

	def mergeModParams(self, prerun, modules, postrun):
		for module in prerun:
			self.mergeModParam(module['name'], module['param'], module['datafile'])
			
		for module in modules:
			self.mergeModParam(module['name'], module['param'], module['datafile'])
			
		for module in postrun:
			self.mergeModParam(module['name'], module['param'], module['datafile'])						
	
	def read(self, spokesModules):
		#merge Parameters and data files
		self.mergeGlobalDatafiles()
		self.mergeModParams(spokesModules.prerun, spokesModules.modules, spokesModules.postrun)
		self.mergeGlobalParamfile()
		self.mergeCmdLineParam()


class SpokesModules:
	#constructor
	def __init__ (self, config, commandline, cwd):
		self.config = config
		self.cwd = cwd
		self.commandline = commandline		
		self.nextModule = 0
		self.prerun = []
		self.modules = []
		self.postrun = []
		self.tests = tests.spokes_tests()
	
		#read module List
		if self.commandline.noConvert != True:
			module_list = eval(self.config.get('prerun','modules'))
			print module_list, type(module_list)
			for module_name in module_list:
				module = self.readModConfig(module_name)
				self.prerun.append(module)
		
		module_list = eval(self.config.get('run','modules'))
		print module_list, type(module_list)
		for module_name in module_list:
			module = self.readModConfig(module_name)
			self.modules.append(module)
			
		module_list = eval(self.config.get('postrun','modules'))
		print module_list, type(module_list)
		for module_name in module_list:
			module = self.readModConfig(module_name)
			self.postrun.append(module)

	#read config data for a module
	def readModConfig(self, mod):
		modDir = self.config.get('options', 'base_dir') + self.config.get(mod, 'directory') + '/'
		return {
			'name': self.config.get(mod, 'name') if self.config.has_option(mod, 'name') else None,
			'directory': modDir,
			'param': modDir + self.config.get(mod, 'param') if self.config.has_option(mod, 'param') else None,
			'datafile': modDir + self.config.get(mod, 'datafile') if self.config.has_option(mod, 'datafile') else None,
			'command': self.config.get(mod, 'command') if self.config.has_option(mod, 'command') else None,
			'datafile': self.config.get(mod, 'datafile') if self.config.has_option(mod, 'datafile') else None,
			'make': self.config.get(mod, 'make') if self.config.has_option(mod, 'make') else None,
			'profile': self.config.get(mod, 'profile') if self.config.has_option(mod, 'profile') else None,
			'test': self.config.get(mod, 'test') if self.config.has_option(mod, 'test') else None}	
			
	def list(self):
		print "\nSpokes pipeline module list \n"
		print "\nPre-Run \n"
		for i in range(len(self.prerun)):
			print self.prerun[i]['name']
		print "\nSpokes Pipeline"
		for i in range(len(self.modules)):
			print '[' + str(i) + ']: ' + self.modules[i]['name']
		print "\n"
		print "\nPost-Run \n"
		for i in range(len(self.postrun)):
			print self.postrun[i]['name']
		print "\n"
		
	def timings(self, runId):
		sumTime = 0
		for i in range(len(self.modules)):
			try:
				sumTime +=  self.modules[i]['endtime'] - self.modules[i]['starttime']
			except KeyError:
				pass
		timelist=[]
		modulelist=[]

		os.chdir(self.cwd)

		timingdata = []

		for i in range(len(self.modules)):
			try:
				print "[" + str(i) + "] %s\t%.2fs\t%.1f%%\t%.0fMB" % (self.modules[i]['name'].ljust(25, " "), self.modules[i]['endtime'] - self.modules[i]['starttime'], (self.modules[i]['endtime'] - self.modules[i]['starttime']) / sumTime *100, self.modules[i]['memory_rss'] / (1024*1024))
				modulelist.append(self.modules[i]['name'])
				timelist.append(self.modules[i]['endtime'] - self.modules[i]['starttime'])
				timingdata.append([self.modules[i]['name'].encode("ascii","ignore"), str(self.modules[i]['endtime'] - self.modules[i]['starttime']), str(self.modules[i]['starttime']), str(self.modules[i]['endtime']), str(self.modules[i]['memory_rss'])])
			except KeyError:
				print "[" + str(i) + "] " + (self.modules[i]['name'].ljust(25, " "))
				timingdata.append(['', '', '', '', ''])

		print "%s\t%.2fs\t%.1f%%" % ("Total".ljust(25, " "), sumTime, 100)

		os.chdir(self.cwd)
		with archive.archive('../../data/data_bank.h5', 'a') as ar:
			ar['/Glue/run/timings'] = timingdata
			ar['/Timing/timelist'] = timelist
			ar['/Timing/modulelist'] = modulelist
			if (runId == None):
				ar['/Glue/run/runId'] = ''
			else:
				ar['/Glue/run/runId'] = runId		
	
	def execProfile(self, module):
		match = re.match("^(i*)python (.*)", module['command'])
		if match:
			os.system("PYTHONPATH=$PYTHONPATH:. python " + self.cwd + "/lsprofcalltree.py -o " + self.cwd + "/profile." + module['name'] + ".log " + match.group(2))
		module['memory_rss'] = 0
	
	def execMemProfile(self, module):
		ptimer = ProcessTimer(module['command'])
		try:
			ptimer.execute()
			while ptimer.poll():
				time.sleep(.2)
		finally:
		#make sure that we don't leave the process dangling?
			ptimer.close()
		if (ptimer.p.returncode != 0):
			#wont happen if you use ipython
			print '... error in ' + module['name']
			sys.exit(1)
		module['memory_rss'] = ptimer.max_rss_memory
	
	def execMod(self, module):	
		if subprocess.check_call(module['command'], shell=True) != 0:
			raise Exception('Error in module: ' + module['name'])
		module['memory_rss'] = 0
	
	def testModule(self, module):
		if (module['test'] != None):
			os.chdir(self.cwd)
			print "\n************************** testing " + module['name'] + " **************************\n"
			if (getattr(self.tests, module['test'])() == False):
				raise Exception('Error in Module result: ' + module['name'])
			else:
				print "test succeeded for: " + module['name'] + " \n\n"
	
	def runModule(self, module):
		module['starttime'] = time.time()
		print "\n\n************************** starting " + module['name'] + " **************************\n\n"
		os.chdir(self.cwd)
		os.chdir(module['directory'])
		
		if module['profile'] != None and module['profile'] == 'True':
			self.execProfile(module)
		elif self.config.has_option("options", "memprofile") and (self.config.get("options", "memprofile") == 'True'):
			self.execMemProfile(module)
		else:
			self.execMod(module)

		module['endtime'] = time.time()
		print "\n\n************************** finished " + module['name'] + " **************************"
		print module['name'] + ' runtime: ' + str(round(module['endtime']-module['starttime'], 2)) + " seconds\n"

		#TODO: insert switch for testcases
		self.testModule(module)
		
	def setup(self, runFrom, runId):
		spokesConfig = SpokesConfig(self.config, self.commandline, self.cwd)
		if (runFrom > 0):
			if runId != None:
				filename = 'data_bank' + str(runFrom - 1) + '_run' + str(runId) + '.h5'
			else:
				filename = 'data_bank' + str(runFrom - 1) + '.h5'
			if os.path.exists(self.config.get("options", "data_dir") + filename):
				print '... copying snapshot to data bank'
				shutil.copy(self.config.get("options", "data_dir") + filename, self.config.get("options", "data_dir") + 'data_bank.h5')
			else:
				raise Exception('Data bank snapshot ' + filename + ' not found, run the pipeline from Module 0')
			#reread config after copy
			spokesConfig.read(self)
		else:
			#read config for convert
			spokesConfig.read(self)
			#run convert etc
			print "\n\n************************** prerun modules **************************\n"
			for i in range(len(self.prerun)):
				print '[' + str(i) + ']: ' + self.prerun[i]['name']
			print " "				
			for module in self.prerun:
				self.runModule(module)
			#reread config after convert
			spokesConfig.read(self)				

	def finish(self, runUntil, runId):
		if runUntil != len(self.modules) and self.commandline.noCopy != True:
			print '... copying data bank snapshot, disable with --noCopy'
			os.chdir(self.cwd)
			if runId != None:
				dest = self.config.get("options", "data_dir") + 'data_bank' + str(runUntil-1) + '_run' + str(runId) + '.h5';
			else:
				dest = self.config.get("options", "data_dir") + 'data_bank' + str(runUntil-1) + '.h5';
			if self.config.has_option("options", "ptrepack"):
				os.system(self.config.get("options", "ptrepack") + ' -o ' + self.config.get("options", "data_dir") + 'data_bank.h5:/ ' + dest + ':/')
			else:
				os.system('ptrepack -o ' + self.config.get("options", "data_dir") + 'data_bank.h5:/ ' + dest + ':/')				
			
		#saving run data for plotting
		if runUntil == len(self.modules):
			print '... saving data for runId ' + str(runId)
			os.chdir(self.cwd)
			with archive.archive('../../data/data_bank.h5') as src:
				with archive.archive('../../data/run' + str(runId) + '.h5', 'w') as dest:
					for f in src['/']:
						#Only Copy the config (first character uppercase)
						if (f[0].isupper()):
							print '/' , f
							copyHDF5Path(src, '/' + f, dest)
					#TODO: save more interesting data (nselectedgal, etc)
					dest['/Glue/run/runId'] = runId
					dest['/Glue/run/ngal'] = len(src['/gal/galaxy_index'])

	
	def run(self, runFrom, runUntil, runId):
		
		self.setup(runFrom, runId)
		
		print "\n\n************************** running modules **************************\n"
		for i in range(runFrom, runUntil):
			print '[' + str(i) + ']: ' + self.modules[i]['name']
		print " "	
	
		#run pipeline
		for i in range(runFrom, runUntil):
			self.runModule(self.modules[i])
		
		print "\n\n************************** postrun modules **************************\n"
		for i in range(len(self.postrun)):
			print '[' + str(i) + ']: ' + self.postrun[i]['name']
		print " "
						
		#run reporting etc
		if self.commandline.report == True or runUntil == len(self.modules): 		
			self.timings(runId)	
			for module in self.postrun:
				self.runModule(module)
			
		#show timings 
		self.timings(runId)
		self.finish(runUntil, runId)
		
				
#Commandline Arguments
parser = argparse.ArgumentParser(description="SPOKES pipeline glue")
parser.add_argument('--config', help="Config file for the glue", required=False)
parser.add_argument('--param', help="Param file overriding other param files", required=False)
parser.add_argument('--list', help="Print a list of the modules with the ids and exit", required=False, action='store_true', default=False)
parser.add_argument('--runModule', type=int, help="Run module X using a saved databank and creates a snapshot afterwards. Reads the params of the module", required=False)
parser.add_argument('--runUntil', type=int, help="Run until module X and create a snapshot of the data bank", required=False)
parser.add_argument('--runFrom', type=int, help="Start with saved databank at module X and continue to run the pipeline", required=False)
parser.add_argument('--report', help="Generate a report ", required=False, action='store_true', default=False)
parser.add_argument('--noCopy', help="Prevent the pipeline from copying the databank after the run", required=False, action='store_true', default=False)
parser.add_argument('--noConvert', help="Don't run convert", required=False, action='store_true', default=False)
parser.add_argument('--runId', help="A run ID that identifies the pipline run", required=False)
parser.add_argument('--getParams', help="Read all the Parameter Files and dump them into the given file and exit", required=False)
parser.add_argument('--compareParams', help="Read all the Parameter Files, compare them to the given file and exit", required=False)

class glue:
	#constructor
	def __init__ (self):
		self.cwd = os.getcwd()
		self.config = ConfigParser.RawConfigParser()
		self.commandline = parser.parse_args()

		if (self.commandline.config != None):
			if (os.path.exists(self.commandline.config)):
				self.config.read(self.commandline.config)
			else:
				raise Exception("Config File " + self.commandline.config + " not found")
		elif (os.path.exists("glue.mine.ini")):
			self.config.read("glue.mine.ini")
		else:
			self.config.read('glue.ini')

		self.runId = None
		if self.config.has_option('options', 'run_id'):
			self.runId = self.config.get('options', 'run_id')
		if self.commandline.runId != None:
			self.runId = self.commandline.runId
		
		#loading modules
		self.modules = SpokesModules(self.config, self.commandline, self.cwd)
		self.modules.list()

	def runPipeline(self, runFrom, runUntil):
		
		self.modules.run(runFrom, runUntil, self.runId)

	def dumpParamFile(self, file):
		paramFile = ConfigParser.RawConfigParser()

		os.chdir(self.cwd)

		for module in range(len(self.modules.modules)):
			if (self.modules.modules[module]['param'] != None and os.path.exists(self.modules.modules[module]['param'])):
				print "... Reading Param File for " + self.modules.modules[module]['name']
				param = ConfigParser.RawConfigParser()
				param.read(self.modules.modules[module]['param']);
				for section in param.sections():
					paramFile.add_section(section)
					for option in param.options(section):
						paramFile.set(section, option, param.get(section, option))
		with open(file, 'wb') as f:
			paramFile.write(f)

	def compareParamFile(self, file):
		paramFile = ConfigParser.RawConfigParser()
		paramFile.read(file)

		os.chdir(self.cwd)
		data = {}

		for module in range(len(self.modules.modules)):
			if (self.modules.modules[module]['param'] != None and os.path.exists(self.modules.modules[module]['param'])):
				print "... Comparing Param File for " + self.modules.modules[module]['name']
				param = ConfigParser.RawConfigParser()
				param.read(self.modules.modules[module]['param']);
				for section in param.sections():
					data[section] = param.options(section)
					for option in param.options(section):
						if paramFile.has_option(section, option) == False:
							print "Missing option in section " + section + ": " + option
						elif paramFile.get(section, option) != param.get(section, option):
							print "Different value in section " + section + ": " + option

		#compare the other way
		for section in paramFile.sections():
			if section not in data:
				print "Unnecessary section " + section
			else:
				for option in paramFile.options(section):
					if option not in data[section]:
						print "Unnecessary option in section " + section + ": " + option

	def glue(self):
		#dispatching
		
		if self.commandline.getParams != None:
			self.dumpParamFile(self.commandline.getParams)
			sys.exit(0)
		elif g.commandline.compareParams != None:
			self.compareParamFile(self.commandline.compareParams)
			sys.exit(0)
		elif self.commandline.list == True:
			self.modules.list()
			sys.exit(0)
		
		if (self.commandline.runModule != None):
			self.runPipeline(self.commandline.runModule, self.commandline.runModule + 1)
		else:
			if (self.commandline.runFrom == None):
				runFrom = 0
			else:
				runFrom = self.commandline.runFrom

			if (self.commandline.runUntil == None):
				runUntil = len(self.modules.modules)
			else:
				runUntil = self.commandline.runUntil+1

			self.runPipeline(runFrom, runUntil)

if __name__ == '__main__':

	print " "
	print "************************** Starting the spokes pipeline **************************"
	print " "
	g = glue()
	g.glue()

	sys.exit(0)

