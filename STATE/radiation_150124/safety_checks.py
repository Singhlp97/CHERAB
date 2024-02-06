def safety_checks(cfg = None):

    if cfg['plasma']['homemade']['use_homemade_emission'] is True:
        if cfg['plasma']['type_radiation'] != 'halpha_total_radiation':
            raise ValueError('Can NOT use homemade emission with total radiation! Only Dalpha allowed...')
        if cfg['raytracing']['sampling']['padding']['use_padding'] is True:
            raise ValueError('Can NOT use homemade emission with padding!')

    if cfg['plasma']['SOLPS']['weights']['use_weighted_emission'] is True and cfg['plasma']['SOLPS']['use_extra_emission'] is False:
    	raise ValueError('Extra emission needed to use weights!')

    if cfg['plasma']['SOLPS']['use_EIRENE_emission'] is True and cfg['plasma']['type_radiation'] == 'total_radiation':
        raise ValueError('Do NOT combine halpha_*_radiation from EIRENE and total_radiation from B2 ;)') 

    if cfg['baserun'] == 'h_alpha_camera' and cfg['plasma']['type_radiation'] == 'total_radiation':
        raise ValueError('Using total_radiation for h_alpha_camera?? o.O') 

    if cfg['baserun'] == 'radiation_load' and cfg['plasma']['type_radiation'] == 'halpha_total_radiation' or cfg['plasma']['type_radiation'] == 'halpha_mol_radiation':
        raise ValueError('Using halpha_*_radiation for radiation_load?? o.O') 

    if cfg['baserun'] == 'radiation_load' and cfg['plasma']['absorption']['use_absorption_function'] is True:
        raise ValueError('Not yet implemented... maybe will never be... :(') 

    return
