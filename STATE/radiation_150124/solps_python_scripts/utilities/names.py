def initialise():

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	# LaTeX conventions

	names = {}

	###########################

	# b2fplasmf

	names['bb'] = 'Magnetic field'

	names['diff'] = 'Particle diffusivity $D^{AN}$'
	names['chi_i'] = 'Ion heat diffusivity $\\chi_i^{AN}$'
	names['chi_e'] = 'Electron heat diffusivity $\\chi_e^{AN}$'

	names['na']    = 'density'
	names['na3da'] = 'Ion density (OMP)'
	names['na3di'] = 'Ion density (IMP)'
	names['na3dl'] = 'Ion density (IT)'
	names['na3dr'] = 'Ion density (OT)'
	names['na3dtl'] = 'Ion density (IT top)'
	names['na3dtr'] = 'Ion density (OT top)'

	names['ne']    = 'density'
	names['ne3da'] = 'Electron density (OMP)'
	names['ne3di'] = 'Electron density (IMP)'
	names['ne3dl'] = 'Electron density (IT)'
	names['ne3dr'] = 'Electron density (OT)'
	names['ne3dtl'] = 'Electron density (IT top)'
	names['ne3dtr'] = 'Electron density (OT top)'

	names['nae'] = '-to-$e^-$ density ratio'
	names['nae3da'] = '-to-$e^-$ density ratio (OMP)'
	names['nae3di'] = '-to-$e^-$ density ratio (IMP)'
	names['nae3dl'] = '-to-$e^-$ density ratio (IT)'
	names['nae3dr'] = '-to-$e^-$ density ratio (OT)'

	names['Zeff'] = 'Effective charge $Z_{eff}$'

	names['ti']    = 'Ion temperature'
	names['ti3da'] = 'Ion temperature (OMP)'
	names['ti3di'] = 'Ion temperature (IMP)'
	names['ti3dl'] = 'Ion temperature (IT)'
	names['ti3dr'] = 'Ion temperature (OT)'
	names['ti3dtl'] = 'Ion temperature (IT top)'
	names['ti3dtr'] = 'Ion temperature (OT top)'

	names['te']    = 'Electron temperature'
	names['te3da'] = 'Electron temperature (OMP)'
	names['te3di'] = 'Electron temperature (IMP)'
	names['te3dl'] = 'Electron temperature (IT)'
	names['te3dr'] = 'Electron temperature (OT)'
	names['te3dtl'] = 'Electron temperature (IT top)'
	names['te3dtr'] = 'Electron temperature (OT top)'

	names['tie']    = 'Ion-to-electron temperature ratio'
	names['tie3da'] = 'Ion-to-electron temperature ratio (OMP)'
	names['tie3di'] = 'Ion-to-electron temperature ratio (IMP)'
	names['tie3dl'] = 'Ion-to-electron temperature ratio (IT)'
	names['tie3dr'] = 'Ion-to-electron temperature ratio (OT)'
	names['tie3dtl'] = 'Ion-to-electron temperature ratio (IT top)'
	names['tie3dtr'] = 'Ion-to-electron temperature ratio (OT top)'

	names['po']    = 'Electric potential'
	names['po3da'] = 'Electric potential (OMP)'
	names['po3di'] = 'Electric potential (IMP)'
	names['po3dl'] = 'Electric potential (IT)'
	names['po3dr'] = 'Electric potential (OT)'
	names['po3dtl'] = 'Electric potential (IT top)'
	names['po3dtr'] = 'Electric potential (OT top)'

	names['ft3da'] = 'Total perpendicular power flux (OMP)'
	names['ft3di'] = 'Total perpendicular power flux (IMP)'
	names['ft3dl'] = 'Total perpendicular power flux (IT)'
	names['ft3dr'] = 'Total perpendicular power flux (OT)'
	names['ft3dtl'] = 'Total perpendicular power flux (IT top)'
	names['ft3dtr'] = 'Total perpendicular power flux (OT top)'

	names['fi3da'] = 'Ion perpendicular power flux (OMP)'
	names['fi3di'] = 'Ion perpendicular power flux (IMP)'
	names['fi3dl'] = 'Ion perpendicular power flux (IT)'
	names['fi3dr'] = 'Ion perpendicular power flux (OT)'
	names['fi3dtl'] = 'Ion perpendicular power flux (IT top)'
	names['fi3dtr'] = 'Ion perpendicular power flux (OT top)'

	names['fe3da'] = 'Electron perpendicular power flux (OMP)'
	names['fe3di'] = 'Electron perpendicular power flux (IMP)'
	names['fe3dl'] = 'Electron perpendicular power flux (IT)'
	names['fe3dr'] = 'Electron perpendicular power flux (OT)'
	names['fe3dtl'] = 'Electron perpendicular power flux (IT top)'
	names['fe3dtr'] = 'Electron perpendicular power flux (OT top)'

	names['fo3da'] = 'Ion perpendicular particle flux (OMP)'
	names['fo3di'] = 'Ion perpendicular particle flux (IMP)'
	names['fo3dl'] = 'Ion perpendicular particle flux (IT)'
	names['fo3dr'] = 'Ion perpendicular particle flux (OT)'
	names['fo3dtl'] = 'Ion perpendicular particle flux (IT top)'
	names['fo3dtr'] = 'Ion perpendicular particle flux (OT top)'

	names['fl3da'] = 'Electron perpendicular particle flux (OMP)'
	names['fl3di'] = 'Electron perpendicular particle flux (IMP)'
	names['fl3dl'] = 'Electron perpendicular particle flux (IT)'
	names['fl3dr'] = 'Electron perpendicular particle flux (OT)'
	names['fl3dtl'] = 'Electron perpendicular particle flux (IT top)'
	names['fl3dtr'] = 'Electron perpendicular particle flux (OT top)'

	names['ft3da'] = 'Total parallel power flux (OMP)'
	names['ft3di'] = 'Total parallel power flux (IMP)'
	names['ft3dlP'] = 'Total parallel power flux (IT)'
	names['ft3drP'] = 'Total parallel power flux (OT)'
	names['ft3dtlP'] = 'Total parallel power flux (IT top)'
	names['ft3dtrP'] = 'Total parallel power flux (OT top)'

	names['fi3da'] = 'Ion parallel power flux (OMP)'
	names['fi3di'] = 'Ion parallel power flux (IMP)'
	names['fi3dlP'] = 'Ion parallel power flux (IT)'
	names['fi3drP'] = 'Ion parallel power flux (OT)'
	names['fi3dtlP'] = 'Ion parallel power flux (IT top)'
	names['fi3dtrP'] = 'Ion parallel power flux (OT top)'

	names['fe3da'] = 'Electron parallel power flux (OMP)'
	names['fe3di'] = 'Electron parallel power flux (IMP)'
	names['fe3dlP'] = 'Electron parallel power flux (IT)'
	names['fe3drP'] = 'Electron parallel power flux (OT)'
	names['fe3dtlP'] = 'Electron parallel power flux (IT top)'
	names['fe3dtrP'] = 'Electron parallel power flux (OT top)'

	names['fo3da'] = 'Ion parallel particle flux (OMP)'
	names['fo3di'] = 'Ion parallel particle flux (IMP)'
	names['fo3dlP'] = 'Ion parallel particle flux (IT)'
	names['fo3drP'] = 'Ion parallel particle flux (OT)'
	names['fo3dtlP'] = 'Ion parallel particle flux (IT top)'
	names['fo3dtrP'] = 'Ion parallel particle flux (OT top)'

	names['fl3da'] = 'Electron parallel particle flux (OMP)'
	names['fl3di'] = 'Electron parallel particle flux (IMP)'
	names['fl3dlP'] = 'Electron parallel particle flux (IT)'
	names['fl3drP'] = 'Electron parallel particle flux (OT)'
	names['fl3dtlP'] = 'Electron parallel particle flux (IT top)'
	names['fl3dtrP'] = 'Electron parallel particle flux (OT top)'

	names['ga3dl'] = 'Grazing angle (IT)'
	names['ga3dr'] = 'Grazing angle (OT)'
	names['ga3dtl'] = 'Grazing angle (IT top)'
	names['ga3dtr'] = 'Grazing angle (OT top)'

	names['sei3dl'] = 'Electron internal energy transmission coefficient (IT)'
	names['sei3dr'] = 'Electron internal energy transmission coefficient (OT)'
	names['sii3dl'] = 'Ion internal energy transmission coefficient (IT)'
	names['sii3dr'] = 'Ion internal energy transmission coefficient (OT)'
	names['set3dl'] = 'Electron total energy transmission coefficient (IT)'
	names['set3dr'] = 'Electron total energy transmission coefficient (OT)'
	names['sit3dl'] = 'Ion total energy transmission coefficient (IT)'
	names['sit3dr'] = 'Ion total energy transmission coefficient (IT)'
	names['stt3dl'] = 'Plasma total energy transmission coefficient (IT)'
	names['stt3dr'] = 'Plasma total energy transmission coefficient (OT)'
	names['siik3dl'] = 'Kin. ion internal energy transmission coefficient (IT)'
	names['siik3dr'] = 'Kin. ion internal energy transmission coefficient (OT)'
	names['sitk3dl'] = 'Kin. ion total energy transmission coefficient (IT)'
	names['sitk3dr'] = 'Kin. ion total energy transmission coefficient (OT)'
	names['sttk3dl'] = 'Kin. plasma total energy transmission coefficient (IT)'
	names['sttk3dr'] = 'Kin. plasma total energy transmission coefficient (OT)'

	names['dab23da'] = 'Atom density (OMP)'
	names['dab23di'] = 'Atom density (IMP)'
	names['dab23dl'] = 'Atom density (IT)'
	names['dab23dr'] = 'Atom density (OT)'
	names['dmb23da'] = 'Molecule density (OMP)'
	names['dmb23di'] = 'Molecule density (IMP)'
	names['dmb23dl'] = 'Molecule density (IT)'
	names['dmb23dr'] = 'Molecule density (OT)'
	names['dib23da'] = 'Molecular ion density (OMP)'
	names['dib23di'] = 'Molecular ion density (IMP)'
	names['dib23dl'] = 'Molecular ion density (IT)'
	names['dib23dr'] = 'Molecular ion density (OT)'

	names['tab23da'] = 'Atom temperature (OMP)'
	names['tab23di'] = 'Atom temperature (IMP)'
	names['tab23dl'] = 'Atom temperature (IT)'
	names['tab23dr'] = 'Atom temperature (OT)'
	names['tmb23da'] = 'Molecule temperature (OMP)'
	names['tmb23di'] = 'Molecule temperature (IMP)'
	names['tmb23dl'] = 'Molecule temperature (IT)'
	names['tmb23dr'] = 'Molecule temperature (OT)'
	names['tib23da'] = 'Molecular ion temperature (OMP)'
	names['tib23di'] = 'Molecular ion temperature (IMP)'
	names['tib23dl'] = 'Molecular ion temperature (IT)'
	names['tib23dr'] = 'Molecular ion temperature (OT)'

	names['rfluxa3da'] = 'Radial atom flux (OMP)'
	names['rfluxa3di'] = 'Radial atom flux (IMP)'
	names['rfluxa3dl'] = 'Radial atom flux (IT)'
	names['rfluxa3dr'] = 'Radial atom flux (OT)'

	names['rfluxm3da'] = 'Radial molecule flux (OMP)'
	names['rfluxm3di'] = 'Radial molecule flux (IMP)'
	names['rfluxm3dl'] = 'Radial molecule flux (IT)'
	names['rfluxm3dr'] = 'Radial molecule flux (OT)'

	names['pfluxa3da'] = 'Poloidal atom flux (OMP)'
	names['pfluxa3di'] = 'Poloidal atom flux (IMP)'
	names['pfluxa3dl'] = 'Poloidal atom flux (IT)'
	names['pfluxa3dr'] = 'Poloidal atom flux (OT)'

	names['pfluxm3da'] = 'Poloidal molecule flux (OMP)'
	names['pfluxm3di'] = 'Poloidal molecule flux (IMP)'
	names['pfluxm3dl'] = 'Poloidal molecule flux (IT)'
	names['pfluxm3dr'] = 'Poloidal molecule flux (OT)'

	names['emiss3da'] = 'Atomic $H_{\\alpha}$ emission from atomic hydrogen (OMP)'
	names['emiss3di'] = 'Atomic $H_{\\alpha}$ emission from atomic hydrogen (IMP)'
	names['emiss3dl'] = 'Atomic $H_{\\alpha}$ emission from atomic hydrogen (IT)'
	names['emiss3dr'] = 'Atomic $H_{\\alpha}$ emission from atomic hydrogen (OT)'

	names['emissmol3da'] = 'Atomic $H_{\\alpha}$ emission from molecular hydrogen (OMP)'
	names['emissmol3di'] = 'Atomic $H_{\\alpha}$ emission from molecular hydrogen (IMP)'
	names['emissmol3dl'] = 'Atomic $H_{\\alpha}$ emission from molecular hydrogen (IT)'
	names['emissmol3dr'] = 'Atomic $H_{\\alpha}$ emission from molecular hydrogen (OT)'

	names['eneutrad3da'] = 'Cumulative atomic emission (OMP)'
	names['eneutrad3di'] = 'Cumulative atomic emission (IMP)'
	names['eneutrad3dl'] = 'Cumulative atomic emission (IT)'
	names['eneutrad3dr'] = 'Cumulative atomic emission (OT)'

	names['emolrad3da'] = 'Cumulative molecular emission (OMP)'
	names['emolrad3di'] = 'Cumulative molecular emission (IMP)'
	names['emolrad3dl'] = 'Cumulative molecular emission (IT)'
	names['emolrad3dr'] = 'Cumulative molecular emission (OT)'

	names['eionrad3da'] = 'Cumulative molecular ion emission (OMP)'
	names['eionrad3di'] = 'Cumulative molecular ion emission (IMP)'
	names['eionrad3dl'] = 'Cumulative molecular ion emission (IT)'
	names['eionrad3dr'] = 'Cumulative molecular ion emission (OT)'

	names['edissml3da'] = 'Power density to atoms due to molecule DS (OMP)'
	names['edissml3di'] = 'Power density to atoms due to molecule DS (IMP)'
	names['edissml3dl'] = 'Power density to atoms due to molecule DS (IT)'
	names['edissml3dr'] = 'Power density to atoms due to molecule DS (OT)'

	names['srcml3da'] = 'Molecule particle source (OMP)'
	names['srcml3di'] = 'Molecule particle source (IMP)'
	names['srcml3dl'] = 'Molecule particle source (IT)'
	names['srcml3dr'] = 'Molecule particle source (OT)'

	names['rqrad'] = 'Line radiation emission'
	names['rqbrm'] = 'Bremsstrahlung emission'

	###########################

	# fort.44

	names['dab2'] = 'density'
	names['tab2'] = 'temperature'
	names['pab2'] = 'pressure'

	names['dmb2'] = 'density'
	names['tmb2'] = 'temperature'
	names['pmb2'] = 'pressure'

	names['dib2'] = 'ion density'
	names['tib2'] = 'ion temperature'
	names['pib2'] = 'ion pressure'
	
	names['dnb2'] = 'density'
	names['tnb2'] = 'temperature'
	names['pnb2'] = 'pressure'

	names['daeb2'] = 'to $e^-$ density ratio'
	names['dmeb2'] = 'to $e^-$ density ratio'
	names['dneb2'] = 'neutral to $e^-$ density ratio'

	names['taib2'] = 'to $i^+$ temperature ratio'
	names['tmib2'] = 'to $i^+$ temperature ratio'
	names['tnib2'] = 'neutral to $i^+$ temperature ratio'

	names['vrab2'] = 'radial velocity'
	names['vpab2'] = 'poloidal velocity'
	names['vtab2'] = 'total velocity'

	names['vrmb2'] = 'radial velocity'
	names['vpmb2'] = 'poloidal velocity'
	names['vtmb2'] = 'total velocity'

	names['e0ab2'] = 'kinetic energy'
	names['e0mb2'] = 'kinetic energy'

	names['emiss']    = 'driven alpha emission'
	names['emissmol'] = 'driven alpha emission'

	names['eneutrad'] = 'total emission'
	names['emolrad']  = 'total emission'
	names['eionrad']  = 'total emission'

	###########################

	# fort.46

	names['pdena'] = 'density'
	names['pdenm'] = 'density'
	names['pdeni'] = 'density'
	names['pdenn'] = 'density'

	names['tdena'] = 'mean energy'
	names['tdenm'] = 'mean energy'
	names['tdeni'] = 'ion mean energy'
	names['tdenn'] = 'mean energy'

	names['edena'] = 'pressure'
	names['edenm'] = 'pressure'
	names['edeni'] = 'ion pressure'
	names['edenn'] = 'pressure'

	return names

def name(of = None):

	try:    
		return initialise()[of]
	except:
		# print('Name of '+ of + ' currently missing! Add it yourself in names.py :)')
		return None
