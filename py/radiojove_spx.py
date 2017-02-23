#! /usr/bin/python
import numpy 
import os
import struct
import pprint as pp
import datetime
from astropy import Time as astime
from spacepy import pycdf

# Decoding the primary header section of RSP files
def load_radiojove_spx_header(hdr_raw):
	
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

# decoding the SPS notes header section
def extract_radiojove_sps_notes(raw_notes,debug):
	
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

# decoding the SPD notes header section
def extract_radiojove_spd_notes(raw_notes):
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

# function to display header information
def display_header(file_spx,debug=False):
		
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
	
# function to read radiojove files (.spd or.sps)
def read_radiojove_spx(file_spx,debug=False):
	#file_sps='/Users/baptiste/Projets/VOParis/RadioJove/data/CDF/data/dat/V01/spectrogram/AJ4CO_DPS_150101071000_corrected_using_CA_2014_12_18_B.sps'
	#file_spd='/Users/baptiste/Projets/VOParis/RadioJove/data/CDF/data/dat/V01/timeseries/AJ4CO_RSP_UT150101000009.spd'
	file_size = os.path.getsize(file_spx)

	# Opening file:
	prim_hdr_length = 156
	lun = open(file_spx,'rb')

	# Reading header:
	prim_hdr_raw = lun.read(prim_hdr_length)
	header = load_radiojove_spx_header(prim_hdr_raw)
	header['file_name'] = os.path.basename(file_spx)
	header['file_type'] = file_spx[-3:].upper()

	# Reading notes:
	notes_raw = lun.read(header['note_length'])
	if header['file_type'] == 'SPS':
		notes = extract_radiojove_sps_notes(notes_raw,debug)
		
	if header['file_type'] == 'SPD':
		notes = extract_radiojove_spd_notes(notes_raw)
		header['nfreq'] = 1

	header['obsty_id'] = 'AJ4CO'
	header['instr_id'] = 'DPS'
	if header['file_type'] == 'SPS':
		header['level'] = 'EDR'
	if header['file_type'] == 'SPD':
		header['level'] = 'DDR'
	
	print header
	print notes
	# Reading data:

	data_length = file_size - prim_hdr_length - header['note_length']

	# nfeed = number of observation feeds 
	# nfreq = number of frequency step (1 for SPD)  
	# nstep = number of sweep (SPS) or time steps (SPD)
	
	# SPS files

	if header['file_type'] == 'SPS':
		header['nfreq'] = header['nchannels']
		if notes['DUALSPECFILE']:
			header['nfeed'] = 2
		else:
			header['nfeed'] = 1
		
		nbytes_per_step = (header['nfreq'] * header['nfeed'] + 1) * 2
		data_fmt = '<%sH' % (nbytes_per_step/2)

		header['fmin'] = notes['LOWF']/1.E6   # MHz
		header['fmax'] = notes['HIF']/1.E6    # MHz
		frequency = header['fmax']-(numpy.arange(header['nfreq'])/float(header['nfreq'])*(header['fmax']-header['fmin']))

	# SPD files
		
	if header['file_type'] == 'SPD':
		header['nfreq'] = 1
		header['nfeed'] = header['nchannels']
		
		if notes['INTEGER_SAVE_FLAG']:
			nbytes_per_sample = 2 
			data_fmt = '%sh' % (header['nfeed'])
		else:
			nbytes_per_sample = 8
			data_fmt = '%sd' % (header['nfeed'])
		
		if notes['NO_TIME_STAMPS_FLAG']:
			nbytes_per_step = header['nfeed']*nbytes_per_sample
			data_fmt = '<%s' % (data_fmt)
		else:
			nbytes_per_step = header['nfeed']*nbytes_per_sample + 8
			data_fmt = '<1d%s' % (data_fmt)
		
		frequency = 20.1
		header['fmin'] = frequency   # MHz
		header['fmax'] = frequency   # MHz
		
	if header['file_type'] == 'SPS':
		header['product_type'] = ('sp%s_%s' % (header['nfeed'],header['nfreq']))
	if header['file_type'] == 'SPD':
		header['product_type'] = ('ts%s' % header['nfeed'])
	
	header['nstep'] = data_length / nbytes_per_step
	data_raw = []
	
	for i in range(header['nstep']):
		data_raw.append(struct.unpack(data_fmt,lun.read(nbytes_per_step)))
	
	lun.close()
	
	data = {}
	if header['file_type'] == 'SPS':
		for i in range(header['nfeed']):
			data['CH%s' % (i+1)] = {}
			for j in range(header['nstep']):
				data['CH%s' % (i+1)][j] = [data_raw[j][i+2*k] for k in range(header['nfreq'])]

	if header['file_type'] == 'SPD':
		if notes['NO_TIME_STAMPS_FLAG']:
			jj=0
		else:
			jj=1
		
		for i in range(header['nfeed']):
			data['CH%s' % (i+1)] = {}
			for j in range(header['nstep']):
				data['CH%s' % (i+1)][j] = data_raw[j][i+jj]
	
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
	header['time_sampling_step'] = time_step*86400.
	
	output = {}
	output['header'] = header
	output['notes'] = notes
	output['data'] = data
	output['frequency'] = frequency
	output['time'] = time
	
	return output
	
#def master_cdf_name(obsty_id, instr_id, level, cdf_version):
#	return 'radiojove_{}_{}_{}_000000000000_000000000000_V{}.cdf'.format(obsty_id, instr_id, level, cdf_version)
#
#def master_skt_name(obsty_id, instr_id, level, cdf_version):
#	return 'radiojove_{}_{}_{}_000000000000_000000000000_V{}.skt'.format(obsty_id, instr_id, level, cdf_version)

def obs_decription(obsty,instr):
	desc = "RadioJOVE {}".format(obsty.upper())
	if instr.upper() == 'DPS':
		desc = "{} Dual Polarization Spectrograph".format(desc)
	return desc
	
def edr_to_cdf(edr,conf,build_cdf_master):

#	setting up variables
#	config = load_local_config(conf)
#	if config['data_path']['relative_path']:
#		path_seed = 
	master_path = 'master/'
	cdfbin_path = '/Applications/cdf/cdf36_0-dist/bin/'
	cdfout_path = '../data/cdf/V10/'
	cdf_version = '10'
	dat_version = '00'
	sft_version = '03'

#	Setting SKT and CDF names 
#	master_cdf = master_cdf_name(edr['header']['obsty_id'], edr['header']['instr_id'], 'edr', cdf_version)
#	skelet_cdf = master_skt_name(edr['header']['obsty_id'], edr['header']['instr_id'], 'edr', cdf_version)
#	
#	Creating CDF Master (removing it if already there)
#	if build_cdf_master:
#		if os.path.exists(master_path+master_cdf):
#			os.remove(master_path+master_cdf)
#		os.system(cdfbin_path+'skeletoncdf -cdf '+master_path+master_cdf+' '+master_path+skelet_cdf)
#
#	print "Master CDF file name:"
#	print master_cdf

#	Creating Time and Frequency Axes 
	ndata = len(edr['time'])
	jul_date = Time(edr['time'],format="datetime",scale="utc").jd.tolist()
	frequency = edr['freq'][0:-1]
	
#	Setting CDF output name 
	cdfout_file = "radiojove_{}_{}_{}_{}_{:%Y%m%d}_V{}.cdf".format(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['level'], edr['header']['product_type'], edr['time'][0].date() ,cdf_version)
	if os.path.exists(cdfout_path+cdfout_file):
		os.remove(cdfout_path+cdfout_file)

#	Opening CDF object 
	pycdf.lib.set_backward(False)
	cdfout = pycdf.CDF(cdfout_path+cdfout_file, '')
	cdfout.col_major(True)
	cdfout.compress(pycdf.const.NO_COMPRESSION)
    
	# SETTING ISTP GLOBAL ATTRIBUTES
	cdfout.attrs['Project']         = ["PDS>Planetary Data System","PADC>Paris Astronomical Data Centre"]
	cdfout.attrs['Discipline']      = "Space Physics>Magnetospheric Science"
	cdfout.attrs['Data_type']       = "{}_{}".format(edr['header']['level'],edr['header']['product_type']).upper()
	cdfout.attrs['Descriptor']      = "{}_{}".format(edr['header']['obsty_id'], edr['header']['instr_id']).upper()

	cdfout.attrs['Data_version']    = cdf_version
	cdfout.attrs['Instrument_type'] = "Radio Telescope"
	cdfout.attrs['Logical_file_id'] = "{}_00000000_v00".format(cdfout.attrs['Logical_source'])
	cdfout.attrs['Logical_source']  = "radiojove_{}_{}".format(cdfout.attrs['Descriptor'].lower(),cdfout.attrs['Data_type'].lower())

	cdfout.attrs['Logical_source_description'] = obs_description(edr['header']['obsty_id'], edr['header']['instr_id'])
	cdfout.attrs['File_naming_convention'] ="source_descriptor_datatype_yyyyMMdd_vVV"	
	cdfout.attrs['Mission_group']   = "RadioJOVE"
	cdfout.attrs['PI_name']         = "Dave Typinski"
	
	cdfout.attrs['PI_affiliation']  = "RadioJOVE"
	cdfout.attrs['Source_name']     = "RadioJOVE"
	cdfout.attrs['TEXT']            = "RadioJOVE Project data. More info at http://radiojove.org and http://radiojove.gsfc.nasa.gov" 
	cdfout.attrs['Generated_by]     = ["SkyPipe","RadioJOVE","PADC"]
	
	cdfout.attrs['Generation_date]  = "{:%Y%m%d}".format(datetime.datetime.now())
	cdfout.attrs['LINK_TEXT']       = ["Radio-SkyPipe Software available on ","More info on RadioJOVE at ","More info on Europlanet at " ]
	cdfout.attrs['LINK_TITLE']      = ["Radio-SkyPipe website","NASA/GSFC web page","Paris Astronomical Data Centre"]
	cdfout.attrs['HTTP_LINK']       = ["http://www.radiosky.com/skypipeishere.html","http://radiojove.gsfc.nasa.gov","http://www.europlanet-vespa.eu" } .
	
	cdfout.attrs['MODS']            = ""
	cdfout.attrs['Rules_of_use']    = "RadioJOVE Data are provided for scientific use. As part of a amateur community project, the RadioJOVE data should be used with careful attention. The RadioJOVE observer of this particular file must be cited or added as a coauthor if the data is central to the study. The RadioJOVE team (radiojove-data@lists.nasa.gov) should also be contacted for any details about publication of studies using this data."
	cdfout.attrs['Skeleton_version'] = cdf_version
	cdfout.attrs['Sotfware_version'] = sft_version
	
	cdfout.attrs['Time_resolution'] = "{} Seconds".format(str(edr['header']['time_sampling_step']))
	cdfout.attrs['Acknowledgement'] = "This study is using data from RadioJOVE project data, that are distributed by NASA/PDS/PPI and PADC at Observatoire de Paris (France)."
	cdfout.attrs['ADID_ref']        = ""
	cdfout.attrs['Validate']        = ""

	cdfout.attrs['Parent']          = edr['header']['file_name']
	cdfout.attrs['Software_language'] = 'python'

	# SETTING PDS GLOBAL ATTRIBUTES
	cdfout.attrs['PDS_Start_time']  = edr['time'][0].isoformat()+'Z'
	cdfout.attrs['PDS_Stop_time']   = edr['time'][ndata-1].isoformat()+'Z'
	cdfout.attrs['PDS_Observation_target'] = 'Jupiter'
	cdfout.attrs['PDS_Observation_type'] = 'Radio'
	
	# SETTING VESPA GLOBAL ATTRIBUTES
	cdfout.attrs['VESPA_dataproduct_type'] = "DS>Dynamic Spectra"
	cdfout.attrs['VESPA_target_class'] = "planet"
	cdfout.attrs['VESPA_target_region'] = "Magnetosphere"
	cdfout.attrs['VESPA_feature_name'] = "Radio Emissions#Aurora"

	cdfout.attrs['VESPA_time_min']  = jul_date[0]
	cdfout.attrs['VESPA_time_max']  = jul_date[ndata-1]
	cdfout.attrs['VESPA_time_sampling_step'] = header['time_step'] #np.median([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.
	cdfout.attrs['VESPA_time_exp']  = header['time_integ']           #np.median([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.

	cdfout.attrs['VESPA_spectral_range_min']  = np.amin(frequency)*1e6
	cdfout.attrs['VESPA_spectral_range_max']  = np.amax(frequency)*1e6
	cdfout.attrs['VESPA_spectral_sampling_step'] = np.median(frequency[i+1]-frequency[i] for i in range(len(frequency)-1))*1e6
	cdfout.attrs['VESPA_spectral_resolution'] = 50.e3

	cdfout.attrs['VESPA_instrument_host_name'] = edr['header']['obsty_id']
	cdfout.attrs['VESPA_instrument_name'] = edr['header']['instr_id']
	cdfout.attrs['VESPA_measurement_type'] = "phys.flux;em.radio
	cdfout.attrs['VESPA_access_format'] = "application/x-cdf"
		
	# SETTING RADIOJOVE GLOBAL ATTRIBUTES

	cdfout.attrs['RadioJOVE_observer_name'] = edr['header']['author']
	cdfout.attrs['RadioJOVE_observatory_location'] = edr['header']['obsloc']
#	cdfout.attrs['RadioJOVE_observatory_latitude'] = edr['header']['']
#	cdfout.attrs['RadioJOVE_observatory_longitude'] = edr['header']['']
#	cdfout.attrs['RadioJOVE_dipole_orientation_angle'] = edr['header']['']
#	cdfout.attrs['RadioJOVE_n_channel'] = edr['header']['']

	date_start = edr['time'][0]
	date_stop = edr['time'][ndata-1]
	date_start_round = edr['time'][0].replace(minute=0, second=0, microsecond=0)
	date_stop_round = edr['time'][ndata-1].replace(minute=0, second=0, microsecond=0)+datetime.timedelta(hours=1)

	# SETTING UP VARIABLES AND VARIABLE ATTRIBUTES
	cdfout.new('EPOCH',data=edr['time'],type=pycdf.const.CDF_TIME_TT2000,compress=pycdf.const.NO_COMPRESSION)
	cdfout['EPOCH'].attrs.new('VALIDMIN',data=date_start,type=pycdf.const.CDF_TIME_TT2000)
	cdfout['EPOCH'].attrs.new('VALIDMAX',data=date_stop,type=pycdf.const.CDF_TIME_TT2000)
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


	cdfout.new('RR',data=edr['data']['spectrum_A'], type=pycdf.const.CDF_REAL4,compress=pycdf.const.NO_COMPRESSION)
    cdfout['RR'].attrs['CATDESC'] = "RCP Flux Density"
    cdfout['RR'].attrs['DEPEND_0'] = "EPOCH"
    cdfout['RR'].attrs['DEPEND_1'] = "FREQUENCY"
	cdfout['RR'].attrs['DICT_KEY'] = "electric_field>power"
	cdfout['RR'].attrs['DISPLAY_TYPE' ] = "spectrogram" 
	cdfout['RR'].attrs['FIELDNAM'] = "RR" 
	cdfout['RR'].attrs.new('FILLVAL',data=-1.0e+31,type=pycdf.const.CDF_REAL4}
	cdfout['RR'].attrs['FORMAT'] = "E12.2"
	cdfout['RR'].attrs['LABLAXIS'] = "RCP Power Spectral Density" 
	cdfout['RR'].attrs['UNITS'] = "ADU" 
	cdfout['RR'].attrs.new('VALIDMIN',data=0.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs.new('VALIDMAX',data=5000.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs['VAR_TYPE'] = "data"
	cdfout['RR'].attrs['SCALETYP'] = "linear"
	cdfout['RR'].attrs.new('SCALEMIN',data=0.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs.new('SCALEMAX',data=5000.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs['SI_CONVERSION'] = " " 
	cdfout['RR'].attrs['UCD'] = "phys.flux;em.radio"

	cdfout.new('LL',data=edr['data']['spectrum_A'], type=pycdf.const.CDF_REAL4,compress=pycdf.const.NO_COMPRESSION)
    cdfout['RR'].attrs['CATDESC'] = "LCP Flux Density"
    cdfout['RR'].attrs['DEPEND_0'] = "EPOCH"
    cdfout['RR'].attrs['DEPEND_1'] = "FREQUENCY"
	cdfout['RR'].attrs['DICT_KEY'] = "electric_field>power"
	cdfout['RR'].attrs['DISPLAY_TYPE' ] = "spectrogram" 
	cdfout['RR'].attrs['FIELDNAM'] = "LL" 
	cdfout['RR'].attrs.new('FILLVAL',data=-1.0e+31,type=pycdf.const.CDF_REAL4}
	cdfout['RR'].attrs['FORMAT'] = "E12.2"
	cdfout['RR'].attrs['LABLAXIS'] = "LCP Power Spectral Density" 
	cdfout['RR'].attrs['UNITS'] = "ADU" 
	cdfout['RR'].attrs.new('VALIDMIN',data=0.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs.new('VALIDMAX',data=5000.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs['VAR_TYPE'] = "data"
	cdfout['RR'].attrs['SCALETYP'] = "linear"
	cdfout['RR'].attrs.new('SCALEMIN',data=0.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs.new('SCALEMAX',data=5000.0,type=pycdf.const.CDF_REAL4)
	cdfout['RR'].attrs['SI_CONVERSION'] = " " 
	cdfout['RR'].attrs['UCD'] = "phys.flux;em.radio"

	cdfout.new('FREQUENCY',data=frequency, type=pycdf.const.CDF_FLOAT,compress=pycdf.const.NO_COMPRESSION, recVary=False)
    cdfout['FREQUENCY'].attrs['CATDESC'] = "Frequency"
	cdfout['FREQUENCY'].attrs['DICT_KEY'] = "electric_field>power"
	cdfout['FREQUENCY'].attrs['FIELDNAM'] = "FREQUENCY" 
	cdfout['FREQUENCY'].attrs.new('FILLVAL',data=-1.0e+31,type=pycdf.const.CDF_REAL4}
	cdfout['FREQUENCY'].attrs['FORMAT'] = "F6.3"
	cdfout['FREQUENCY'].attrs['LABLAXIS'] = "Frequency" 
	cdfout['FREQUENCY'].attrs['UNITS'] = "MHz" 
	cdfout['FREQUENCY'].attrs.new('VALIDMIN',data=1.0,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs.new('VALIDMAX',data=40.,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs['VAR_TYPE'] = "support_data"
	cdfout['FREQUENCY'].attrs['SCALETYP'] = "linear"
	cdfout['FREQUENCY'].attrs.new('SCALEMIN',data=10.,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs.new('SCALEMAX',data=40.,type=pycdf.const.CDF_REAL4)
	cdfout['FREQUENCY'].attrs['SI_CONVERSION'] = "1.0e6>Hz" 
	cdfout['FREQUENCY'].attrs['UCD'] = "em.freq"
	
	cdfout.close()

	
		
	
