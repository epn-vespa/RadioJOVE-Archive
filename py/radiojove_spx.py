#! /usr/bin/python
import numpy 
import os
import struct
import pprint as pp
import datetime
from astropy import time as astime
from spacepy import pycdf
import json

################################################################################
# Loading local config file 
################################################################################
def load_local_config(config_file,debug=False):
	
	if debug:
		print "### [load_local_config]"
	
	with open(config_file) as f:
		config = json.load(f)
	
	if debug:
		print config
	
	return config

################################################################################
# Decoding the primary header section of RSP files
################################################################################
def load_radiojove_spx_header(hdr_raw,debug=False):

	if debug:
		print "### [load_radiojove_spx_header]"
	
	hdr_fmt = '<10s6d1h10s20s20s40s1h1i'
	hdr_values = struct.unpack(hdr_fmt,hdr_raw[0:156])
	
	header = {}
	header['sft_version']  = hdr_values[0]
	header['start_jdtime'] = hdr_values[1]-0.5+2415020 # + julday(1,0,1900)
	header['start_time']   = astime.Time(header['start_jdtime'],format='jd').datetime
	header['stop_jdtime']  = hdr_values[2]-0.5+2415020 # + julday(1,0,1900)
	header['stop_time']    = astime.Time(header['stop_jdtime'],format='jd').datetime
	header['latitude']     = hdr_values[3]
	header['longitude']    = hdr_values[4]
	header['chartmax']     = hdr_values[5]
	header['chartmin']     = hdr_values[6]
	header['timezone']     = hdr_values[7]
	header['source']       = (hdr_values[8].strip('\x00')).strip(' ')
	header['author']       = (hdr_values[9].strip('\x00')).strip(' ')
	header['obsname']      = (hdr_values[10].strip('\x00')).strip(' ')
	header['obsloc']       = (hdr_values[11].strip('\x00')).strip(' ')
	header['nchannels']    = hdr_values[12]
	header['note_length']  = hdr_values[13]

	return header

################################################################################
# decoding the SPS notes header section
################################################################################
def extract_radiojove_sps_notes(raw_notes,debug=False):

	if debug:
		print "### [extract_radiojove_sps_notes]"
	
	notes = {}
	
	# list of valid metadata keys
	key_list = ['SWEEPS','LOWF','HIF','STEPS','RCVR','DUALSPECFILE','COLORRES','BANNER','ANTENNATYPE','ANTENNAORIENTATION','COLORFILE','COLOROFFSET','COLORGAIN','CORRECTIONFILENAME','CAXF','CAX1','CAX2','CLOCKMSG']
	
	# list of metadata keys with multiple values
	key_list_multi = ['BANNER','COLOROFFSET','COLORGAIN','CAXF','CAX1','CAX2','CLOCKMSG']

	# list of metadata keys with interger values
	key_list_int = ['SWEEPS','LOWF','HIF','STEPS','RCVR','COLORRES','COLOROFFSET','COLORGAIN']
	
	# stripping raw note text stream from "*[[*" and "*]]*" delimiters, and splitting with '\xff'
	start_index = raw_notes.find('*[[*')
	stop_index  = raw_notes.find('*]]*')
	notes['free_text'] = raw_notes[0:start_index]
	note_list = raw_notes[start_index+4:stop_index].strip('\xff').split('\xff')

	# Looping on note items
	for note_item in note_list:
		if debug: 
			print 'Current Item = %s' % note_item
		
		# looping on valid key items
		for key_item in key_list:
		
			# getting length of key name
			key_len = len(key_item)

			# checking if current note item contains current key item
			if note_item[0:key_len] == key_item:
			
				if debug: 
					print 'Detected Key = %s' % key_item
				# if current key item has multiple values, do this
				if key_item in key_list_multi:
					
					# if current key item has multiple values, initializing the list
					if key_item not in notes.keys():
						notes[key_item] = [] 

					# checking specific cases 
					if key_item[0:3] == 'CAX':
						note_item = note_item.split('|')
						note_index = int(note_item[0][key_len:])
						note_value = note_item[1]
					elif key_item[0:8] == 'CLOCKMSG':
						note_item = note_item.split(' ')
						note_index = int(note_item[0][key_len:])
						note_value = ' '.join(note_item[1:])
					else:
						note_index = int(note_item[key_len:key_len+1])
						note_value = note_item[key_len+1:]
					if debug: 
						print 'Index = %s' % note_index
						print 'Value = %s' % note_value

					# setting value to note item 
					notes[key_item].append(note_value)

				else:
					note_value = note_item[key_len:]
					if key_item in key_list_int:
						if (key_item == 'RCVR') & (note_value == ''):
							note_value = '-1'
						if debug: 
							print 'Value = %s' % note_value
						notes[key_item] = int(note_value)
					else:
						if debug: 
							print 'Value = %s' % note_value
						notes[key_item] = note_value

	if 'DUALSPECFILE' not in notes.keys():
		notes['DUALSPECFILE'] = False
	
	return notes

################################################################################
# decoding the SPD notes header section
################################################################################
def extract_radiojove_spd_notes(raw_notes,debug=False):

	if debug:
		print "### [extract_radiojove_spd_notes]"

	notes = {}
	notes['CHL'] = {}
	notes['CHO'] = {}
	notes['MetaData'] = {}

	start_index = raw_notes.find('*[[*')
	stop_index  = raw_notes.find('*]]*')
	notes['free_text'] = raw_notes[0:start_index]
	note_list = raw_notes[start_index+4:stop_index].split('\xff')
	for note_item in note_list:
		
		if note_item == 'Logged Using UT':
			notes['Logged Using UT'] = True
		else:
			notes['Logged Using UT'] = False
		
		if note_item == 'No Time Stamps':
			notes['No Time Stamps'] = True
		else:
			notes['No Time Stamps'] = False

		if note_item[0:3] == 'CHL':
			notes['CHL'][int(note_item[3])] = note_item[4:]

		if note_item[0:3] == 'CHO':
			notes['CHO'][int(note_item[3])] = note_item[4:]

		if note_item == 'Integer Save':
			notes['Integer Save'] = True
		else:
			notes['Integer Save'] = False

		if note_item[0:7] == 'XALABEL':
			notes['XALABEL'] = note_item[7:]

		if note_item[0:7] == 'YALABEL':
			notes['YALABEL'] = note_item[7:]

		if note_item[0:9] == 'MetaData_':
			item_metadata = note_item.split('\xc8')
			notes['MetaData'][item_metadata[0][9:].strip(' ').strip(':').strip(' ')] = item_metadata[1]

	return notes

################################################################################
# function to display header information
################################################################################
def display_header(file_spx,debug=False):

	if debug:
		print "### [display_header]"
		
	# Opening file:
	prim_hdr_length = 156
	lun = open(file_spx,'rb')
	# Reading header:
	prim_hdr_raw = lun.read(prim_hdr_length)
	header = load_radiojove_spx_header(prim_hdr_raw)
	pp.pprint(header)

	header['file_type'] = file_spx[-3:].upper()
	# Reading notes:
	notes_raw = lun.read(header['note_length'])
	print notes_raw
	if header['file_type'] == 'SPS':
		notes = extract_radiojove_sps_notes(notes_raw,debug)
	if header['file_type'] == 'SPD':
		notes = extract_radiojove_spd_notes(notes_raw)
	pp.pprint(notes)
	return
	
################################################################################
# Open SPx file and return header,notes
################################################################################
def open_radiojove_spx(file_info,debug=False):

	if debug:
		print "### [open_radiojove_spx]"
	
	file_info['size'] = os.path.getsize(file_info['name'])
	
	# Opening file:
	file_info['prim_hdr_length'] = 156
	file_info['lun'] = open(file_info['name'],'rb')

	# Reading header:
	file_info['prim_hdr_raw'] = file_info['lun'].read(file_info['prim_hdr_length'])
	header = load_radiojove_spx_header(file_info['prim_hdr_raw'])
	header['file_name'] = file_info['name']
	header['file_type'] = file_info['name'][-3:].upper()

	# Reading notes:
	file_info['notes_raw'] = file_info['lun'].read(header['note_length'])
	if header['file_type'] == 'SPS':
		notes = extract_radiojove_sps_notes(file_info['notes_raw'],debug)
		
	if header['file_type'] == 'SPD':
		notes = extract_radiojove_spd_notes(file_info['notes_raw'],debug)
		header['nfreq'] = 1

	if header['obsname'] == 'AJ4CO DPS':
		header['obsty_id'] = 'AJ4CO'
		header['instr_id'] = 'DPS'
		header['gain0'] = 1.95
		header['gain1'] = 1.95
		header['offset0'] = 1975
		header['offset1'] = 1975
	else:
		header['obsty_id'] = 'ABCDE'
		header['instr_id'] = 'XXX'
		header['gain0'] = 2.00
		header['gain1'] = 2.00
		header['offset0'] = 2000
		header['offset1'] = 2000

	if header['file_type'] == 'SPS':
		header['level'] = 'EDR'
	if header['file_type'] == 'SPD':
		header['level'] = 'DDR'

	if debug:
		print header
		print notes

	# Reading data:

	file_info['data_length'] = file_info['size'] - file_info['prim_hdr_length'] - header['note_length']

	# nfeed = number of observation feeds 
	# nfreq = number of frequency step (1 for SPD)  
	# nstep = number of sweep (SPS) or time steps (SPD)
	
	# SPS files
	header['feeds'] = {}

	if header['file_type'] == 'SPS':
		header['nfreq'] = header['nchannels']
		if notes['DUALSPECFILE']:
			header['nfeed'] = 2
			header['feeds'][0] = {}
			header['feeds'][0]['FIELDNAM'] = 'RR'
			header['feeds'][0]['CATDESC']  = 'RCP Flux Density'
			header['feeds'][0]['LABLAXIS'] = 'RCP Power Spectral Density'
			header['feeds'][1] = {}
			header['feeds'][1]['FIELDNAM'] = 'LL'
			header['feeds'][1]['CATDESC']  = 'LCP Flux Density'
			header['feeds'][1]['LABLAXIS'] = 'LCP Power Spectral Density'
		else:
			header['nfeed'] = 1
			header['feeds'][0] = {}
			header['feeds'][0]['FIELDNAM'] = 'RR'
			header['feeds'][0]['CATDESC']  = 'RCP Flux Density'
			header['feeds'][0]['LABLAXIS'] = 'RCP Power Spectral Density'
		
		file_info['bytes_per_step'] = (header['nfreq'] * header['nfeed'] + 1) * 2
		file_info['data_format']    = '>%sH' % (file_info['bytes_per_step']/2)

		header['fmin'] = notes['LOWF']/1.E6   # MHz
		header['fmax'] = notes['HIF']/1.E6    # MHz
		frequency = header['fmax']-(numpy.arange(header['nfreq'])/float(header['nfreq'])*(header['fmax']-header['fmin']))

	# SPD files
		
	if header['file_type'] == 'SPD':
		header['nfreq'] = 1
		header['nfeed'] = header['nchannels']
		for i in range(header['nchannels']):
			header['feeds'][i] = {}
			header['feeds'][i]['FIELDNAM'] = 'CH{:02d}'.format(i)
			header['feeds'][i]['CATDESC']  = 'CH{:02d} Flux Density'.format(i)
			header['feeds'][i]['LABLAXIS'] = 'CH{:02d} Flux Density'.format(i)
		
		if notes['INTEGER_SAVE_FLAG']:
			file_info['bytes_per_step'] = 2 
			file_info['data_format']    = '%sh' % (header['nfeed'])
		else:
			file_info['nbytes_per_sample'] = 8
			file_info['data_format']    = '%sd' % (header['nfeed'])
		
		if notes['NO_TIME_STAMPS_FLAG']:
			file_info['bytes_per_step'] = header['nfeed']*file_info['nbytes_per_sample']
			file_info['data_format'] = '<%s' % (file_info['data_format'])
		else:
			file_info['bytes_per_step'] = header['nfeed']*file_info['nbytes_per_sample'] + 8
			file_info['data_format'] = '<1d%s' % (file_info['data_format'])
		
		frequency = 20.1
		header['fmin'] = frequency   # MHz
		header['fmax'] = frequency   # MHz
		
	if header['file_type'] == 'SPS':
		header['product_type'] = ('sp{}_{}'.format(header['nfeed'],header['nfreq']))
		file_info['record_data_offset'] = 0
	if header['file_type'] == 'SPD':
		header['product_type'] = ('ts{}'.format(header['nfeed']))
		if notes['NO_TIME_STAMPS_FLAG']:
			file_info['record_data_offset']=0
		else:
			file_info['record_data_offset']=1
	
	header['nstep'] = file_info['data_length'] / file_info['bytes_per_step']

	if header['file_type'] == 'SPS':
		time_step = (header['stop_jdtime']-header['start_jdtime']) / float(header['nstep'])
		time = numpy.arange(header['nstep'])*time_step+header['start_jdtime']
	if header['file_type'] == 'SPD':
		if notes['NO_TIME_STAMPS_FLAG']:
			time_step = (header['stop_jdtime']-header['start_jdtime']) / float(header['nstep'])
			time = numpy.arange(header['nstep'])*time_step+header['start_jdtime']
		else:
			time = numpy.array()
			for i in range(header['nstep']):
				time.append(data_raw[i][0])
			time_step = numpy.median(time[1:nstep]-time[0:nstep-1])
	
	# transforming times from JD to datetime
	time = astime.Time(time,format='jd').datetime

	# time sampling step in seconds
	header['time_step'] = time_step*86400.
	header['time_integ'] = header['time_step'] # this will have to be checked at some point
	
	if debug:
		print "nfeed : {}".format(header['nfeed'])
		print "nfreq : {} ({})".format(header['nfreq'],len(frequency))
		print "nstep : {} ({})".format(header['nstep'],len(time))

	return header,notes,time,frequency
	
################################################################################
# Initialize raw data object
################################################################################
def init_radiojove_data(header,notes,debug=False):
	if debug:
		print "### [init_radiojove_data]"

	data = {}
	for i in range(header['nfeed']):
		data[header['feeds'][i]['FIELDNAM']] = {}
	
	return data
	
################################################################################
# Read SPx sweep
################################################################################
def read_radiojove_spx_sweep(file_info,packet_size,debug=False):

	if debug:
		print "### [read_radiojove_spx_sweep]"
		print "loading packet of {} step(s), with format `{}`.".format(packet_size,file_info['data_format'])

	raw = []
	for i in range(packet_size):
		raw.append(struct.unpack(file_info['data_format'],file_info['lun'].read(file_info['bytes_per_step'])))
		if raw[i][-1] != 65278:
			print "WARNING ! wrong end of sweep delimiter. (Got 0x{:04X} instead of 0x{:04X})".format(raw[i][-1],65278) 
			

	if debug:
		print "Size of loaded data: {}".format(len(raw))
	
	return raw

################################################################################
# Close SPx file
################################################################################
def close_radiojove_spx(file_info,debug=False):
	if debug:
		print "### [close_radiojove_spx]"
	
	file_info['lun'].close()
	return

################################################################################
# Check CDF file with PDS script
################################################################################
def check_radiojove_cdf(file_info,debug=False):
	if debug:
		print "### [check_radiojove_cdf]"
	
	os.system(config)
	return

################################################################################
# Init CDF output file
################################################################################
def init_radiojove_cdf(file_info,header,start_time,config,debug=False):
	if debug:
		print "### [init_radiojove_cdf]"

#	Setting CDF output name 
	if file_info['daily']:
		file_info['cdfout_file'] = "radiojove_{}_{}_{}_{}_{:%Y%m%d}_V{}.cdf".format(header['obsty_id'], header['instr_id'], header['level'], header['product_type'], start_time.date() ,config['vers']['cdf_version']).lower()
	else:
		file_info['cdfout_file'] = "radiojove_{}_{}_{}_{}_{:%Y%m%d%H%M}_V{}.cdf".format(header['obsty_id'], header['instr_id'], header['level'], header['product_type'], start_time ,config['vers']['cdf_version']).lower()
	if os.path.exists(config['path']['cdfout_path']+file_info['cdfout_file']):
		os.remove(config['path']['cdfout_path']+file_info['cdfout_file'])
	
	print "CDF file output: {}".format(config['path']['cdfout_path']+file_info['cdfout_file'])
		
#	Opening CDF object 
	pycdf.lib.set_backward(False)
	cdfout = pycdf.CDF(config['path']['cdfout_path']+file_info['cdfout_file'],'')
	cdfout.col_major(True)
	cdfout.compress(pycdf.const.NO_COMPRESSION)

	return cdfout
	
################################################################################
# Close CDF output file
################################################################################
def close_radiojove_cdf(cdfout,debug=False):
	if debug:
		print "### [close_radiojove_cdf]"

	cdfout.close()
	return
	
################################################################################
# Global Attributes for CDF 
################################################################################
def write_gattr_radiojove_cdf(cdfout,header,time,freq,config,debug=False):
	if debug:
		print "### [write_gattr_radiojove_cdf]"

	# Creating Time and Frequency Axes 
	ndata = len(time)
	jul_date = astime.Time(time,format="datetime",scale="utc").jd.tolist()
	
    
	# SETTING ISTP GLOBAL ATTRIBUTES
	cdfout.attrs['Project']         = ["PDS>Planetary Data System","PADC>Paris Astronomical Data Centre"]
	cdfout.attrs['Discipline']      = "Space Physics>Magnetospheric Science"
	cdfout.attrs['Data_type']       = "{}_{}".format(header['level'],header['product_type']).upper()
	cdfout.attrs['Descriptor']      = "{}_{}".format(header['obsty_id'],header['instr_id']).upper()
	cdfout.attrs['Data_version']    = config['vers']['cdf_version']
	cdfout.attrs['Instrument_type'] = "Radio Telescope"
	cdfout.attrs['Logical_source']  = "radiojove_{}_{}".format(cdfout.attrs['Descriptor'],cdfout.attrs['Data_type']).lower()
	cdfout.attrs['Logical_file_id'] = "{}_00000000_v00".format(cdfout.attrs['Logical_source'])
	cdfout.attrs['Logical_source_description'] = obs_description(header['obsty_id'], header['instr_id'])
	cdfout.attrs['File_naming_convention'] ="source_descriptor_datatype_yyyyMMdd_vVV"	
	cdfout.attrs['Mission_group']   = "RadioJOVE"
	cdfout.attrs['PI_name']         = header['author']
	cdfout.attrs['PI_affiliation']  = "RadioJOVE"
	cdfout.attrs['Source_name']     = "RadioJOVE"
	cdfout.attrs['TEXT']            = "RadioJOVE Project data. More info at http://radiojove.org and http://radiojove.gsfc.nasa.gov" 
	cdfout.attrs['Generated_by']    = ["SkyPipe","RadioJOVE","PADC"]	
	cdfout.attrs['Generation_date'] = "{:%Y%m%d}".format(datetime.datetime.now())
	cdfout.attrs['LINK_TEXT']       = ["Radio-SkyPipe Software available on ","More info on RadioJOVE at ","More info on Europlanet at " ]
	cdfout.attrs['LINK_TITLE']      = ["Radio-SkyPipe website","NASA/GSFC web page","Paris Astronomical Data Centre"]
	cdfout.attrs['HTTP_LINK']       = ["http://www.radiosky.com/skypipeishere.html","http://radiojove.gsfc.nasa.gov","http://www.europlanet-vespa.eu"]
	cdfout.attrs['MODS']            = ""
	cdfout.attrs['Rules_of_use']    = "RadioJOVE Data are provided for scientific use. As part of a amateur community project, the RadioJOVE data should be used with careful attention. The RadioJOVE observer of this particular file must be cited or added as a coauthor if the data is central to the study. The RadioJOVE team (radiojove-data@lists.nasa.gov) should also be contacted for any details about publication of studies using this data."
	cdfout.attrs['Skeleton_version'] = config['vers']['cdf_version']
	cdfout.attrs['Sotfware_version'] = config['vers']['sft_version']
	cdfout.attrs['Time_resolution'] = "{} Seconds".format(str(header['time_step']))
	cdfout.attrs['Acknowledgement'] = "This study is using data from RadioJOVE project data, that are distributed by NASA/PDS/PPI and PADC at Observatoire de Paris (France)."
	cdfout.attrs['ADID_ref']        = ""
	cdfout.attrs['Validate']        = ""
	cdfout.attrs['Parent']          = os.path.basename(header['file_name'])
	cdfout.attrs['Software_language'] = 'python'

	# SETTING PDS GLOBAL ATTRIBUTES
	cdfout.attrs['PDS_Start_time']  = time[0].isoformat()+'Z'
	cdfout.attrs['PDS_Stop_time']   = time[ndata-1].isoformat()+'Z'
	cdfout.attrs['PDS_Observation_target'] = 'Jupiter'
	cdfout.attrs['PDS_Observation_type'] = 'Radio'
	
	# SETTING VESPA GLOBAL ATTRIBUTES
	cdfout.attrs['VESPA_dataproduct_type'] = "DS>Dynamic Spectra"
	cdfout.attrs['VESPA_target_class'] = "planet"
	cdfout.attrs['VESPA_target_region'] = "Magnetosphere"
	cdfout.attrs['VESPA_feature_name'] = "Radio Emissions#Aurora"

	cdfout.attrs['VESPA_time_min']  = jul_date[0]
	cdfout.attrs['VESPA_time_max']  = jul_date[ndata-1]
	cdfout.attrs['VESPA_time_sampling_step'] = header['time_step'] #numpy.median([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.
	cdfout.attrs['VESPA_time_exp']  = header['time_integ']           #numpy.median([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.

	cdfout.attrs['VESPA_spectral_range_min']  = numpy.amin(freq)*1e6
	cdfout.attrs['VESPA_spectral_range_max']  = numpy.amax(freq)*1e6
	cdfout.attrs['VESPA_spectral_sampling_step'] = numpy.median([freq[i+1]-freq[i] for i in range(len(freq)-1)])*1e6
	cdfout.attrs['VESPA_spectral_resolution'] = 50.e3

	cdfout.attrs['VESPA_instrument_host_name'] = header['obsty_id']
	cdfout.attrs['VESPA_instrument_name'] = header['instr_id']
	cdfout.attrs['VESPA_measurement_type'] = "phys.flux;em.radio"
	cdfout.attrs['VESPA_access_format'] = "application/x-cdf"
		
	# SETTING RADIOJOVE GLOBAL ATTRIBUTES

	cdfout.attrs['RadioJOVE_observer_name']   = header['author']
	cdfout.attrs['RadioJOVE_observatory_loc'] = header['obsloc']
	cdfout.attrs['RadioJOVE_observatory_lat'] = header['latitude']
	cdfout.attrs['RadioJOVE_observatory_lon'] = header['longitude']
	cdfout.attrs['RadioJOVE_sft_version']     = header['sft_version']
	cdfout.attrs['RadioJOVE_chartmin']  = header['chartmin']
	cdfout.attrs['RadioJOVE_chartmax']  = header['chartmax']
	cdfout.attrs['RadioJOVE_nchannels'] = header['nfeed']
	cdfout.attrs['RadioJOVE_rcvr']   	= -1
	cdfout.attrs['RadioJOVE_banner0']   = ""
	cdfout.attrs['RadioJOVE_banner1']   = ""
	cdfout.attrs['RadioJOVE_antenna_type']    = ""
	cdfout.attrs['RadioJOVE_antenna_orientation'] = ""
	cdfout.attrs['RadioJOVE_color_file']      = ""
	cdfout.attrs['RadioJOVE_color_offset0']   = header['offset0']
	cdfout.attrs['RadioJOVE_color_offset1']   = header['offset1']
	cdfout.attrs['RadioJOVE_color_gain0']     = header['gain0']
	cdfout.attrs['RadioJOVE_color_gain1']     = header['gain1']
	cdfout.attrs['RadioJOVE_correction_filename'] = ""
	cdfout.attrs['RadioJOVE_caxf']   	= ""
	cdfout.attrs['RadioJOVE_cax1']   	= ""
	cdfout.attrs['RadioJOVE_cax2']   	= ""
	cdfout.attrs['RadioJOVE_clockmsg']  = ""
	
	
	if debug:
		print cdfout.attrs
	return
	
################################################################################
# EPOCH variable for CDF
################################################################################
def write_epoch_radiojove_cdf(cdfout,time,debug=False):
	if debug:
		print "### [write_epoch_radiojove_cdf]"

	ndata = len(time)
	date_start = time[0]
	date_stop  = time[ndata-1]
	date_start_round = time[0].replace(minute=0, second=0, microsecond=0)
	date_stop_round  = time[ndata-1].replace(minute=0, second=0, microsecond=0)+datetime.timedelta(hours=1)
	
	# SETTING UP VARIABLES AND VARIABLE ATTRIBUTES
	cdfout.new('EPOCH',data=time,type=pycdf.const.CDF_TIME_TT2000,compress=pycdf.const.NO_COMPRESSION)
	cdfout['EPOCH'].attrs.new('VALIDMIN',data=datetime.datetime(2000,1,1),type=pycdf.const.CDF_TIME_TT2000)
	cdfout['EPOCH'].attrs.new('VALIDMAX',data=datetime.datetime(2100,1,1),type=pycdf.const.CDF_TIME_TT2000)
	cdfout['EPOCH'].attrs.new('SCALEMIN',data=date_start_round,type=pycdf.const.CDF_TIME_TT2000)
	cdfout['EPOCH'].attrs.new('SCALEMAX',data=date_stop_round,type=pycdf.const.CDF_TIME_TT2000)
	cdfout['EPOCH'].attrs['CATDESC']  = "Default time (TT2000)"
	cdfout['EPOCH'].attrs['FIELDNAM'] = "Epoch"
	cdfout['EPOCH'].attrs.new('FILLVAL',data=-9223372036854775808,type=pycdf.const.CDF_TIME_TT2000)
	cdfout['EPOCH'].attrs['LABLAXIS'] = "Epoch"
	cdfout['EPOCH'].attrs['UNITS']    = "ns" 
	cdfout['EPOCH'].attrs['VAR_TYPE'] = "support_data"
	cdfout['EPOCH'].attrs['SCALETYP'] = "linear" 
	cdfout['EPOCH'].attrs['MONOTON']  = "INCREASE"
	cdfout['EPOCH'].attrs['TIME_BASE'] = "J2000" 
	cdfout['EPOCH'].attrs['TIME_SCALE'] = "UTC" 
	cdfout['EPOCH'].attrs['REFERENCE_POSITION'] = "Earth"
	cdfout['EPOCH'].attrs['SI_CONVERSION'] = "1.0e-9>s" 
	cdfout['EPOCH'].attrs['UCD']      = "time.epoch"
	
	if debug:
		print cdfout['EPOCH']
		print cdfout['EPOCH'].attrs
	return 
	
################################################################################
# FREQUENCY variable for CDF
################################################################################
def write_frequency_radiojove_cdf(cdfout,header,freq,debug=False):
	if debug:
		print "### [write_frequency_radiojove_cdf]"

	cdfout.new('FREQUENCY',data=freq, type=pycdf.const.CDF_FLOAT,compress=pycdf.const.NO_COMPRESSION, recVary=False)
	cdfout['FREQUENCY'].attrs['CATDESC'] = "Frequency"
	cdfout['FREQUENCY'].attrs['DICT_KEY'] = "electric_field>power"
	cdfout['FREQUENCY'].attrs['FIELDNAM'] = "FREQUENCY" 
	cdfout['FREQUENCY'].attrs.new('FILLVAL',data=-1.0e+31,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs['FORMAT'] = "F6.3"
	cdfout['FREQUENCY'].attrs['LABLAXIS'] = "Frequency" 
	cdfout['FREQUENCY'].attrs['UNITS'] = "MHz" 
	cdfout['FREQUENCY'].attrs.new('VALIDMIN',data=0.,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs.new('VALIDMAX',data=40.,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs['VAR_TYPE'] = "support_data"
	cdfout['FREQUENCY'].attrs['SCALETYP'] = "linear"
	cdfout['FREQUENCY'].attrs.new('SCALEMIN',data=header['fmin'],type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs.new('SCALEMAX',data=header['fmax'],type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs['SI_CONVERSION'] = "1.0e6>Hz" 
	cdfout['FREQUENCY'].attrs['UCD'] = "em.freq"

	if debug:
		print cdfout['FREQUENCY']
		print cdfout['FREQUENCY'].attrs
	return

################################################################################
# Data variables for CDF
################################################################################
def write_data_radiojove_cdf(cdfout,header,time,freq,file_info,packet_size,debug=False):
	if debug:
		print "### [write_data_radiojove_cdf]"

	nt = header['nstep'] # len(time)
	nf = header['nfreq'] # len(freq)
	
	# defining variables
	for i in header['feeds'].keys():
		var_name = header['feeds'][i]['FIELDNAM']
		if debug:
			print "Creating {} variable".format(var_name)
		cdfout.new(var_name,data=numpy.zeros((nt,nf)),type=pycdf.const.CDF_UINT2,compress=pycdf.const.NO_COMPRESSION)
		cdfout[var_name].attrs['CATDESC'] = header['feeds'][i]['CATDESC']
		cdfout[var_name].attrs['DEPEND_0'] = "EPOCH"
		cdfout[var_name].attrs['DEPEND_1'] = "FREQUENCY"
		cdfout[var_name].attrs['DICT_KEY'] = "electric_field>power"
		cdfout[var_name].attrs['DISPLAY_TYPE' ] = "spectrogram" 
		cdfout[var_name].attrs['FIELDNAM'] = var_name
		cdfout[var_name].attrs.new('FILLVAL',data=65535,type=pycdf.const.CDF_UINT2)
		cdfout[var_name].attrs['FORMAT'] = "E12.2"
		cdfout[var_name].attrs['LABLAXIS'] = header['feeds'][i]['LABLAXIS']
		cdfout[var_name].attrs['UNITS'] = "ADU" 
		cdfout[var_name].attrs.new('VALIDMIN',data=0,type=pycdf.const.CDF_UINT2)
		cdfout[var_name].attrs.new('VALIDMAX',data=4096,type=pycdf.const.CDF_UINT2)
		cdfout[var_name].attrs['VAR_TYPE'] = "data"
		cdfout[var_name].attrs['SCALETYP'] = "linear"
		cdfout[var_name].attrs.new('SCALEMIN',data=2050,type=pycdf.const.CDF_UINT2)
		cdfout[var_name].attrs.new('SCALEMAX',data=2300,type=pycdf.const.CDF_UINT2)
		cdfout[var_name].attrs['SI_CONVERSION'] = " " 
		cdfout[var_name].attrs['UCD'] = "phys.flux;em.radio"

	# reading sweeps structure
	for j in range(0,header['nstep'],packet_size):
		j1 = j
		j2 = j+packet_size
		if j2 > nt:
			j2=nt

		if debug:
			if packet_size == 1:
				print "Loading record #{}".format(j)		
			else: 
				print "Loading records #{} to #{}".format(j1,j2)
		
		data_raw = numpy.array(read_radiojove_spx_sweep(file_info,j2-j1,debug))[:,file_info['record_data_offset']:file_info['record_data_offset']+header['nfreq']*header['nfeed']].reshape(j2-j1,header['nfreq'],header['nfeed'])
				
		for i in range(header['nfeed']):
			cdfout[header['feeds'][i]['FIELDNAM']][j1:j2,:] = data_raw[:,:,i]

	return

################################################################################
# Observatory Descriptions
################################################################################

def obs_description(obsty,instr,debug=False):
	if debug:
		print "### [obs_description]"

	desc = "RadioJOVE {}".format(obsty.upper())
	if instr.upper() == 'DPS':
		desc = "{} Dual Polarization Spectrograph".format(desc)
	return desc
	
################################################################################
# Main SPX to CDF 
################################################################################
def spx_to_cdf(file_spx,config_file='local_config_bc.json',daily=False,debug=False):
	if debug:
		print "### [spx_to_cdf]"

#	setting up local paths and versions 
	config = load_local_config(config_file,debug)

	#file_sps='/Users/baptiste/Projets/VOParis/RadioJove/data/CDF/data/dat/V01/spectrogram/AJ4CO_DPS_150101071000_corrected_using_CA_2014_12_18_B.sps'
	#file_spd='/Users/baptiste/Projets/VOParis/RadioJove/data/CDF/data/dat/V01/timeseries/AJ4CO_RSP_UT150101000009.spd'

	file_info = {}
	file_info['name'] = file_spx
	file_info['daily'] = daily

	# Opening file, initializing file info and loading header + notes
	header,notes,time,frequency = open_radiojove_spx(file_info,debug)
	
	# initializing data structure
	#data = init_radiojove_data(header,notes)

	#edr = {}
	#edr['header'] = header
	#edr['notes'] = notes
	#edr['data'] = data
	#edr['freq'] = frequency
	#edr['time'] = time

	# initializing CDF file
	cdfout = init_radiojove_cdf(file_info,header,time[0],config,debug)

	write_gattr_radiojove_cdf(cdfout,header,time,frequency,config,debug)
	write_epoch_radiojove_cdf(cdfout,time,debug)
	write_data_radiojove_cdf(cdfout,header,time,frequency,file_info,config['proc']['packet_size'],debug)
	write_frequency_radiojove_cdf(cdfout,header,frequency,debug)

	close_radiojove_cdf(cdfout,debug)
	
	check_radiojove_cdf(file_info,config)
	
	return
	
		
	
