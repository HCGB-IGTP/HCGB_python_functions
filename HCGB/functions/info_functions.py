#!/usr/bin/env python3
############################################################
## Jose F. Sanchez										##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################

import json
import os
from HCGB import functions
import hashlib
from hmac import compare_digest


######################################
def dump_info_run(out_folder, module, userinput, runInfo, debug):
	"""
	Prints information of the job run in json format
	
	:param out_folder: Absolute path to folder to dump results
	:param module: Module or script name executed
	:param runInfo: Contains configuration detail information to dump
	
	:type out_folder: string 
	:type module: string 
	:type runInfo: dict
	
	:returns: Prints information in json format in the output folder provided.
	
	Example of information to include:
	
	userinput = {"filename":infiles, 
				"database":path_to_database,
				"file_format":file_format}
				
	runInfo = { "module":module, 
				"analysis":example,
				"date":date, 
				"time":time}
	
	Original idea extracted from https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/kmerfinder.py
	
	"""
	
	# Collect all data of run-analysis
	
	## convert to dictionary
	userinput_dict = vars(userinput)
	del userinput_dict['func']
	
	## merge dictionaries:: data = {**runInfo, **userinput_dict}
	
	## add sub dictionary
	data = runInfo.copy()
	data['input_options'] = userinput_dict
	
	## debug messages
	if debug:
		functions.aesthetics_functions.debug_message("Dump information:", 'yellow')
		print()
		functions.aesthetics_functions.debug_message("runInfo:", 'yellow')
		print(runInfo)
		print(type(runInfo))
		print()
		
		functions.aesthetics_functions.debug_message("userinput_dict:", 'yellow')
		print(userinput_dict)
		print(type(userinput_dict))
		print()
		
		functions.aesthetics_functions.debug_message("data:", 'yellow')
		print(data)
		print(type(data))
		print()
	
	# Save json output
	result_json_file = os.path.join(out_folder, module + ".json") 
	with open(result_json_file, "w") as outfile:  
		json.dump(data, outfile, indent=3)

######################################
def compare_info_dict(dict1, dict2, listItems, debug):
	"""Given a pair of dictionaries and a list of items, compare values."""
	##
	print()

######################################
def dump_info_conda(folder, module, package_name="BacterialTyper", debug=False):
	"""
	Prints information of the conda environment in yml format
	
	:param folder: Absolute path to folder to dump conda environment file 
	:param module: Module or script name executed
	
	:type folder: string 
	:type module: string 
	
	:returns: Prints information in yml format in the output folder provided.
	
	Example of information to include:
	
	...
	name: Bacterialtyper_dev
	channels:
	  - bioconda
	  - defaults
	dependencies:
	  - _libgcc_mutex=0.1=conda_forge
	  - _openmp_mutex=4.5=2_gnu
	  - _sysroot_linux-64_curr_repodata_hack=3=h69a702a_13
	  - alsa-lib=1.2.10=hd590300_0
	  ...

	"""
	
	## Check if file exist and contains information
	file_yml = os.path.join(folder, package_name + '.yml')
	tmp_file_yml = os.path.join(folder, 'tmp_' + package_name + ".yml")
	## if file exists and contains information
	if functions.files_functions.is_non_zero_file(file_yml):
		## check filehash 
		file_hash = read_filehash(file_yml)
		
		## Create a new one
		conda_export_call_tmp = 'conda env export -f ' + tmp_file_yml
		functions.system_call_functions.system_call(conda_export_call_tmp,  True, message=False)
		
		file_hash2 = read_filehash(tmp_file_yml)
		
		## if different, create a new file containing module name provided
		if debug:
			functions.aesthetics_functions.debug_message("Filehash:")
			print("file_yml: " + file_yml)
			print(file_hash.hexdigest())
			
			print("file_yml2: " + tmp_file_yml)
			print(file_hash2.hexdigest())
			
		if file_hash.hexdigest() == file_hash2.hexdigest():
			if debug:
				functions.aesthetics_functions.debug_message("Equal filehashes > Return")
				
			## remove tmp_file generated
			os.remove(tmp_file_yml)
			
			return()
		else:
			data_string = functions.time_functions.create_human_timestamp()
			new_tag = data_string + '_' + module
			new_file_yml = os.path.join(folder, new_tag + '_' + package_name + '.yml')
	
			## rename file_yml + '_tmp' to new_file_yml
			os.rename(tmp_file_yml, new_file_yml)
	else:
		## Create a new one
		conda_export_call_tmp = 'conda env export -f ' + file_yml
		functions.system_call_functions.system_call(conda_export_call_tmp,  True, message=False)
		


######################################
def read_filehash(file_given):
	"""
	Obtains file hash for the given file using python hashlib
	
	:param file_given: Absolute path to file to check
	:type file_given: string 
	
	:returns: file_hash generated (hashlib.sha256())
	
	Example of information generated:
	
	Parameters
	----------
	file_given : /path/fo/file/file1.py

	Returns
	-------
	f1adc236059bf58ea6cf2f6838cafb1daaf0353c3e1d4fcac5cff1ff6a20d2ab

	"""

	## original code from: https://nitratine.net/blog/post/how-to-hash-files-in-python/
	BLOCK_SIZE = 65536 # The size of each read from the file

	file_hash = hashlib.sha256() # Create the hash object, can use something other than `.sha256()` if you wish
	with open(file_given, 'rb') as f: # Open the file to read it's bytes
		fb = f.read(BLOCK_SIZE) # Read from the file. Take in the amount declared above
		while len(fb) > 0: # While there is still data being read from the file
			file_hash.update(fb) # Update the hash
			fb = f.read(BLOCK_SIZE) # Read the next block from the file
	
	return(file_hash)


