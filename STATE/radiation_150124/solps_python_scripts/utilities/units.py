def initialise():

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	# LaTeX conventions

	units = {}

	###########################

	# b2fplasmf

	units['bb'] = 'T'

	units['d']    = 'm^2 \\cdot s^{-1}'
	units['chii'] = 'm^2 \\cdot s^{-1}'
	units['chie'] = 'm^2 \\cdot s^{-1}'

	units['na']    = 'm^{-3}'
	units['na3da'] = 'm^{-3}'
	units['na3di'] = 'm^{-3}'
	units['na3dl'] = 'm^{-3}'
	units['na3dr'] = 'm^{-3}'	
	units['na3dtl'] = 'm^{-3}'
	units['na3dtr'] = 'm^{-3}'

	units['ne']    = 'm^{-3}'
	units['ne3da'] = 'm^{-3}'
	units['ne3di'] = 'm^{-3}'
	units['ne3dl'] = 'm^{-3}'
	units['ne3dr'] = 'm^{-3}'
	units['ne3dtl'] = 'm^{-3}'
	units['ne3dtr'] = 'm^{-3}'

	units['nae'] = '-'
	units['nae3da'] = '-'
	units['nae3di'] = '-'
	units['nae3dl'] = '-'
	units['nae3dr'] = '-'

	units['Zeff'] = '-'

	units['ti']    = 'eV'
	units['ti3da'] = 'eV'
	units['ti3di'] = 'eV'
	units['ti3dl'] = 'eV'
	units['ti3dr'] = 'eV'
	units['ti3dtl'] = 'eV'
	units['ti3dtr'] = 'eV'

	units['te']    = 'eV'
	units['te3da'] = 'eV'
	units['te3di'] = 'eV'
	units['te3dl'] = 'eV'
	units['te3dr'] = 'eV'
	units['te3dtl'] = 'eV'
	units['te3dtr'] = 'eV'

	# uits['tie']    = '-'
	units['tie3da'] = '-'
	units['tie3di'] = '-'
	units['tie3dl'] = '-'
	units['tie3dr'] = '-'
	units['tie3dtl'] = '-'
	units['tie3dtr'] = '-'

	units['po']    = 'eV'
	units['po3da'] = 'eV'
	units['po3di'] = 'eV'
	units['po3dl'] = 'eV'
	units['po3dr'] = 'eV'
	units['po3dtl'] = 'eV'
	units['po3dtr'] = 'eV'

	units['ft3da'] = 'W \\cdot m^{-2}'
	units['ft3di'] = 'W \\cdot m^{-2}'
	units['ft3dl'] = 'W \\cdot m^{-2}'
	units['ft3dr'] = 'W \\cdot m^{-2}'
	units['ft3dtl'] = 'W \\cdot m^{-2}'
	units['ft3dtr'] = 'W \\cdot m^{-2}'

	units['fi3da'] = 'W \\cdot m^{-2}'
	units['fi3di'] = 'W \\cdot m^{-2}'
	units['fi3dl'] = 'W \\cdot m^{-2}'
	units['fi3dr'] = 'W \\cdot m^{-2}'
	units['fi3dtl'] = 'W \\cdot m^{-2}'
	units['fi3dtr'] = 'W \\cdot m^{-2}'

	units['fe3da'] = 'W \\cdot m^{-2}'
	units['fe3di'] = 'W \\cdot m^{-2}'
	units['fe3dl'] = 'W \\cdot m^{-2}'
	units['fe3dr'] = 'W \\cdot m^{-2}'
	units['fe3dtl'] = 'W \\cdot m^{-2}'
	units['fe3dtr'] = 'W \\cdot m^{-2}'

	units['fo3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dr'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dtl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dtr'] = 'm^{-2} \\;\\cdot s^{-1}'

	units['fl3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dr'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dtl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dtr'] = 'm^{-2} \\;\\cdot s^{-1}'

	units['ft3da'] = 'W \\cdot m^{-2}'
	units['ft3di'] = 'W \\cdot m^{-2}'
	units['ft3dlP'] = 'W \\cdot m^{-2}'
	units['ft3drP'] = 'W \\cdot m^{-2}'
	units['ft3dtlP'] = 'W \\cdot m^{-2}'
	units['ft3dtrP'] = 'W \\cdot m^{-2}'

	units['fi3da'] = 'W \\cdot m^{-2}'
	units['fi3di'] = 'W \\cdot m^{-2}'
	units['fi3dlP'] = 'W \\cdot m^{-2}'
	units['fi3drP'] = 'W \\cdot m^{-2}'
	units['fi3dtlP'] = 'W \\cdot m^{-2}'
	units['fi3dtrP'] = 'W \\cdot m^{-2}'

	units['fe3da'] = 'W \\cdot m^{-2}'
	units['fe3di'] = 'W \\cdot m^{-2}'
	units['fe3dlP'] = 'W \\cdot m^{-2}'
	units['fe3drP'] = 'W \\cdot m^{-2}'
	units['fe3dtlP'] = 'W \\cdot m^{-2}'
	units['fe3dtrP'] = 'W \\cdot m^{-2}'

	units['fo3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dlP'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3drP'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dtlP'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fo3dtrP'] = 'm^{-2} \\;\\cdot s^{-1}'

	units['fl3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dlP'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3drP'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dtlP'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['fl3dtrP'] = 'm^{-2} \\;\\cdot s^{-1}'

	units['ga3dr'] = 'degrees'
	units['ga3dl'] = 'degrees'
	units['ga3dtr'] = 'degrees'
	units['ga3dtl'] = 'degrees'

	units['sei3dl'] = '-'
	units['sei3dr'] = '-'
	units['sii3dl'] = '-'
	units['sii3dr'] = '-'
	units['set3dl'] = '-'
	units['set3dr'] = '-'
	units['sit3dl'] = '-'
	units['sit3dr'] = '-'
	units['stt3dl'] = '-'
	units['stt3dr'] = '-'
	units['siik3dl'] = '-'
	units['siik3dr'] = '-'
	units['sitk3dl'] = '-'
	units['sitk3dr'] = '-'
	units['sttk3dl'] = '-'
	units['sttk3dr'] = '-'

	units['dab23da'] = 'm^{-3}'
	units['dab23di'] = 'm^{-3}'
	units['dab23dl'] = 'm^{-3}'
	units['dab23dr'] = 'm^{-3}'
	units['dmb23da'] = 'm^{-3}'
	units['dmb23di'] = 'm^{-3}'
	units['dmb23dl'] = 'm^{-3}'
	units['dmb23dr'] = 'm^{-3}'
	units['dib23da'] = 'm^{-3}'
	units['dib23di'] = 'm^{-3}'
	units['dib23dl'] = 'm^{-3}'
	units['dib23dr'] = 'm^{-3}'

	units['tab23da'] = 'eV'
	units['tab23di'] = 'eV'
	units['tab23dl'] = 'eV'
	units['tab23dr'] = 'eV'
	units['tmb23da'] = 'eV'
	units['tmb23di'] = 'eV'
	units['tmb23dl'] = 'eV'
	units['tmb23dr'] = 'eV'
	units['tib23da'] = 'eV'
	units['tib23di'] = 'eV'
	units['tib23dl'] = 'eV'
	units['tib23dr'] = 'eV'

	units['rfluxa3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['rfluxa3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['rfluxa3dl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['rfluxa3dr'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxa3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxa3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxa3dl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxa3dr'] = 'm^{-2} \\;\\cdot s^{-1}'

	units['rfluxm3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['rfluxm3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['rfluxm3dl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['rfluxm3dr'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxm3da'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxm3di'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxm3dl'] = 'm^{-2} \\;\\cdot s^{-1}'
	units['pfluxm3dr'] = 'm^{-2} \\;\\cdot s^{-1}'

	units['emiss3da'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emiss3di'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emiss3dl'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emiss3dr'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'

	units['emissmol3da'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emissmol3di'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emissmol3dl'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emissmol3dr'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'

	units['eneutrad3da'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['eneutrad3di'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['eneutrad3dl'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['eneutrad3dr'] = 'W \\rightarrow W \\cdot m^{-3}'

	units['emolrad3da'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['emolrad3di'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['emolrad3dl'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['emolrad3dr'] = 'W \\rightarrow W \\cdot m^{-3}'

	units['eionrad3da'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['eionrad3di'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['eionrad3dl'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['eionrad3dr'] = 'W \\rightarrow W \\cdot m^{-3}'

	units['edissml3da'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['edissml3di'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['edissml3dl'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['edissml3dr'] = 'W \\rightarrow W \\cdot m^{-3}'

	units['srcml3da'] = 'A \\rightarrow m^{-3}'
	units['srcml3di'] = 'A \\rightarrow m^{-3}'
	units['srcml3dl'] = 'A \\rightarrow m^{-3}'
	units['srcml3dr'] = 'A \\rightarrow m^{-3}'

	units['rqrad'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['rqbrm'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['rqradtot'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['rqbrmtot'] = 'W \\rightarrow W \\cdot m^{-3}'

	###########################

	# fort.44

	units['dab2'] = 'm^{-3}'
	units['tab2'] = 'eV'
	units['pab2'] = 'Pa'

	units['dmb2'] = 'm^{-3}'
	units['tmb2'] = 'eV'
	units['pmb2'] = 'Pa'

	units['dib2'] = 'm^{-3}'
	units['tib2'] = 'eV'
	units['pib2'] = 'Pa'

	units['dnb2'] = 'm^{-3}'
	units['tnb2'] = 'eV'
	units['pnb2'] = 'Pa'

	units['daeb2'] = '-'
	units['dmeb2'] = '-'
	units['dneb2'] = '-'

	units['taib2'] = '-'
	units['tmib2'] = '-'
	units['tnib2'] = '-'

	units['vrab2'] = 'm \\cdot s^{-1}'
	units['vpab2'] = 'm \\cdot s^{-1}'
	units['vtab2'] = 'm \\cdot s^{-1}'
	
	units['vrmb2'] = 'm \\cdot s^{-1}'
	units['vpmb2'] = 'm \\cdot s^{-1}'
	units['vtmb2'] = 'm \\cdot s^{-1}'

	units['e0ab2'] = 'eV'
	units['e0mb2'] = 'eV'

	units['emiss']    = 'photons \\cdot m^{-3} \\cdot s^{-1}'
	units['emissmol'] = 'photons \\cdot m^{-3} \\cdot s^{-1}'

	units['eneutrad'] = 'W \\rightarrow W \\cdot m^{-3}'
	units['emolrad']  = 'W \\rightarrow W \\cdot m^{-3}'
	units['eionrad']  = 'W \\rightarrow W \\cdot m^{-3}'
	units['etotrad']  = 'W \\rightarrow W \\cdot m^{-3}'

	###########################

	# fort.46

	units['pdena'] = 'm^{-3}'
	units['pdenm'] = 'm^{-3}'
	units['pdeni'] = 'm^{-3}'
	units['pdenn'] = 'm^{-3}'

	units['tdena'] = 'eV'
	units['tdenm'] = 'eV'
	units['tdeni'] = 'eV'
	units['tdenn'] = 'eV'

	units['edena'] = 'Pa'
	units['edenm'] = 'Pa'
	units['edeni'] = 'Pa'
	units['edenn'] = 'Pa'

	###########################

	# databases

	# units['AMJUEL'] = 'm^{-3} \\cdot s^{-1}'
	# units['HYDHEL'] = 'm^{-3} \\cdot s^{-1}'
	# units['AMMONX'] = 'm^{-3} \\cdot s^{-1}'
	# units['ADAS']   = 'm^{-3} \\cdot s^{-1}'

	return units

def units(of = None):

	try:    
		return initialise()[of]
	except:
		if 'AMJUEL' in of or 'HYDHEL' in of or 'AMMONX' in of:
			if 'mfp' in of:
				return 'm'
			else:
				return 'm^{-3} \\cdot s^{-1}'
		# if 'relative variation' in of:
		# 	return '\\%'
		print('Unit of ' + of + ' currently missing! Add it yourself in units.py :)')
		return '-'
