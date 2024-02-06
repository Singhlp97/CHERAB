from import_ASTRA_data import load_astra_data
from import_SOLPS_data import load_solps_data
from plasma_emitter    import make_plasma_emitter
from homemade_emitter  import make_homemade_emitter

######################################################################################################################
######################################################################################################################
######################################################################################################################

def load_plasma_radiation(cfg = None, parent = None):

    if parent is None:
      raise ValueError('A parent node, e.g. World(), must be provided :)')

    use_ASTRA_emission  = cfg['plasma']['ASTRA']['use_ASTRA_emission']
    use_B2_emission     = cfg['plasma']['SOLPS']['use_B2_emission']
    use_EIRENE_emission = cfg['plasma']['SOLPS']['use_EIRENE_emission']

    use_homemade_emission = cfg['plasma']['homemade']['use_homemade_emission']

    ################

    if use_homemade_emission is False:

        astra_simulation = None
        solps_simulation = None

        if use_ASTRA_emission is True: astra_simulation = load_astra_data(cfg = cfg)
        if use_B2_emission    is True: solps_simulation = load_solps_data(cfg = cfg)

        plasmaEmitter = make_plasma_emitter(cfg = cfg,
                                            parent = parent,
                                            astra_simulation = astra_simulation,
                                            solps_simulation = solps_simulation)
    else:
        plasmaEmitter = make_homemade_emitter(cfg = cfg,
                                              parent = parent)

    return