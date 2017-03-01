import radiojove_spx as rj
from spacepy import pycdf

file='/Users/baptiste/Projets/VOParis/RadioJove/data/CDF/data/dat/V01/spectrogram/AJ4CO_DPS_150228020000_corrected_using_CA_2014_12_18_B.sps'

rj.spx_to_cdf(file, '../config/local_config_bc.json', True, False)

cdf = pycdf.CDF('/Users/baptiste/Projets/VOParis/RadioJove/data/CDF/data/cdf/tmp/radiojove_aj4co_dps_edr_sp2_300_20150301_v10.cdf')
print cdf['EPOCH'].attrs
print cdf['RR'].attrs
print cdf['LL'].attrs
print cdf['FREQUENCY'].attrs
cdf.close()