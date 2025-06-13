import numpy as np

from src.components import phytoplankton as phyto
from src.components import dim
from src.components import heterotrophs as het
from src.components import detritus
from src.components import dom
from src.components import flocs
from src.utils import phys
from src.utils import functions as fns
from src.config_model import varinfos
from src.config_model import base_config

# Onur22 = Kerimoglu et al 2022
# =============================
# Model units: [d-1], [mmol m-3] for all tracers except Chl [mg m-3]


"""
Notes on values and units
202501 - theta_max:
0.07 gChl/gC in Onur
Saumya's analysis: at MOW1, 15-55 gC/cChl = 0.018-0.066 gChl/gC
"""

iniC = 3.
MOW_ICs = [{'Phy': {'initialization': {'C': iniC,
                                       'N': iniC / 6.,
                                       'Chl': 3.5,
                                       'P': iniC / 106,
                                       'Si': iniC * 0.3,  }}},
           {'BacF': {'initialization':{'C': 1. * 1e-1,
                                       'N': 1. / 4.9 * 1e-1,
                                       'P': 1. / 77. * 1e-1}}},
           {'BacA': {'initialization':{'C': 1. * 1e-1,
                                       'N': 1. / 4.9 * 1e-1,
                                       'P': 1. / 77. * 1e-1}}},
           {'HF': {'initialization': {'C': 1. * 1e-1,
                                      'N': 1. / 7. * 1e-1,
                                      'P': 1. / 180. * 1e-1}}},
           {'Cil': {'initialization': {'C': 1. * 1e-1,
                                       'N': 1. / 7. * 1e-1,
                                       'P': 1. / 180. * 1e-1}}},
           {'DSi': {'initialization': {'concentration': 20.}}},
           {'DIP': {'initialization': {'concentration': 1.}}},
           {'NH4': {'initialization': {'concentration': 4.5}}},
           {'NO3': {'initialization': {'concentration': 45}}},
           ]

MOW1_some_pars = [{'DetL': {'parameters': {'kleak': 0.}}},
                       {'Phy': {'parameters': {'alpha': 7. * 1e-3 / varinfos.molmass_C,  # 5.8*1e-2/varinfos.molmass_C,
                                               }}}]

MOW1_noFlocs = fns.deep_update(base_config.Onur, *MOW_ICs, *MOW1_some_pars)



# FLOCS
# =====
# Model units: s-1 hence conversion to d-1 is needed and done through
# Np_0=2.25e11 #  2.25e11

prFloc_df_ini = 32.0e-6
"""
CHECK 202501
SPMC are 50-200 mg/l on a climatology basis
Hence 120 mg/l is a good start
--> 0.12 in g/l 
"""

pconc = 0.12  # [g/l] -> = [kg/m3]
pdens = 2500  # [kg/m³ or g/l]
pfrac_ini = 0.9999  # [-] # Initial fraction of "free" microflocs
df_ini = prFloc_df_ini

prFloc_dp = 5e-6  # [m]
prFloc_nfrac = 2.1  # [-]

nc_ini = (df_ini / prFloc_dp) ** prFloc_nfrac  # Initial number of microflocs per macroflocs
pntot = pconc / (
            pdens * 1.0 / 6.0 * np.pi * prFloc_dp ** 3)  # Initial total number concentration of microflocs = free+in_macro

Flocs = {
    'formulation': 'Onur22',
    'Microflocs': {'class': flocs.Flocs,
                   'parameters': {
                       'alpha_PP': 0.02,  # prFloc.alpha
                       'alpha_PF': 0.02,  # prFloc.alpha
                       'alpha_FF': 0.02,  # prFloc.alpha
                       'p_exp': 0.513,  # prFloc.brk_p
                       'q_exp': 1.004,  # prFloc.brk_q
                       'f_frac_floc_break': 0.467,  # prFloc.brk_f
                       'efficiency_break': 1.5e-4,  # /np.sqrt(3600*24), #prFloc.brk_eb
                       'fyflocstrength': 1e-10,  # prFloc.brk_fy
                       'mu_viscosity': 1e-6,  # /3600/24, #prFloc.brk_p
                       'd_p_microflocdiam': 5e-6,  # prFloc.dp
                       'nf_fractal_dim': 2.1  # prFloc.nfrac
                   },
                   'coupling': {
                       'coupled_Nf': 'Macroflocs',
                       'coupled_Nt': 'Micro_in_Macro'
                   },
                   'initialization': {
                       'numconc': pntot * pfrac_ini},
                   'diagnostics': ['volconcentration',
                                   'massconcentration', 'alpha_FF', 'fyflocstrength', 'alpha_PF', 'alpha_PP', 'Ncnum'],
                   },
    'Macroflocs': {'class': flocs.Flocs,
                   'parameters': {
                       'alpha_PP': 0.02,  # prFloc.alpha
                       'alpha_PF': 0.02,  # prFloc.alpha
                       'alpha_FF': 0.02,  # prFloc.alpha
                       'p_exp': 0.513,  # prFloc.brk_p
                       'q_exp': 1.004,  # prFloc.brk_q
                       'f_frac_floc_break': 0.467,  # prFloc.brk_f
                       'efficiency_break': 1.5e-4,  # /np.sqrt(3600*24), #prFloc.brk_eb
                       'fyflocstrength': 1e-10,  # prFloc.brk_fy
                       'mu_viscosity': 1e-6,  # /3600/24, #prFloc.brk_p
                       'd_p_microflocdiam': 5e-6,  # prFloc.dp
                       'nf_fractal_dim': 2.1  # prFloc.nfrac
                   },
                   'coupling': {
                       'coupled_Np': 'Microflocs',
                       'coupled_Nt': 'Micro_in_Macro'
                   },
                   'initialization': {
                       'numconc': pntot * (1.0 - pfrac_ini) / nc_ini},
                   'diagnostics': ['volconcentration', 'diam', 'settling_vel', 'massconcentration',
                                   'alpha_FF', 'fyflocstrength', 'alpha_PF', 'alpha_PP', 'Ncnum'],

                   },
    'Micro_in_Macro': {'class': flocs.Flocs,
                       'parameters': {
                           'alpha_PP': 0.02,  # prFloc.alpha
                           'alpha_PF': 0.02,  # prFloc.alpha
                           'alpha_FF': 0.02,  # prFloc.alpha
                           'p_exp': 0.513,  # prFloc.brk_p
                           'q_exp': 1.004,  # prFloc.brk_q
                           'f_frac_floc_break': 0.467,  # prFloc.brk_f
                           'efficiency_break': 1.5e-4,  # /np.sqrt(3600*24), #prFloc.brk_eb
                           'fyflocstrength': 1e-10,  # prFloc.brk_fy
                           'mu_viscosity': 1e-6,  # /3600/24, #prFloc.brk_p
                           'd_p_microflocdiam': 5e-6,  # prFloc.dp
                           'nf_fractal_dim': 2.1  # prFloc.nfrac
                       },
                       'coupling': {
                           'coupled_Np': 'Microflocs',
                           'coupled_Nf': 'Macroflocs'
                       },
                       'initialization': {
                           'numconc': pntot * (1.0 - pfrac_ini)},
                       'diagnostics': ['volconcentration',
                                       'massconcentration', 'alpha_FF', 'fyflocstrength', 'alpha_PF', 'alpha_PP', 'Ncnum'],

                       },
}

Flocs_TEP_FF_only = fns.deep_update(Flocs, {
    'Microflocs': {
        'parameters': {
            'alpha_FF_ref': 0.05,
            'K_glue': 15.0,
            'deltaFymax': 1e-10
        },
        'coupling': {
            'coupled_glue': 'TEPC'
        },
    },
    'Macroflocs':  {'parameters': {'sinking_leak': 1e-8},
                    'coupling': { 'coupled_glue': 'TEPC'  }},
    'Micro_in_Macro': {'coupling': {'coupled_glue': 'TEPC'}},
})

Flocs_TEP_all_alphas = fns.deep_update(Flocs_TEP_FF_only, {
    'Microflocs': {
        'parameters': {
            'alpha_PF_ref': 0.05,
            'alpha_PP_ref': 0.05
        },
        # 'diagnostics': Flocs_TEP_FF_only['Microflocs']['diagnostics'] + ['alpha_PF', 'alpha_PP']
    },
})

SPM_phyto_light_attenuation = [{'Phy': {'coupling': {'coupled_SPM': ['Microflocs', 'Micro_in_Macro']}}}]




MOW1 = fns.deep_update(base_config.Onur, *MOW_ICs, *MOW1_some_pars, Flocs)

phyto_Flocs_2ways = fns.deep_update(MOW1, Flocs_TEP_all_alphas, *SPM_phyto_light_attenuation)
# fns.print_dict(phyto_Flocs_2ways)

newDOCaggregpars =  [{'DOCS': {'parameters': {'alpha_TEPC': 0.1}}},
                     {'DOCL': {'parameters': {'alpha_TEPC': 0.3}}}]

phyto_Flocs_2ways_modifiedDOCalphas = fns.deep_update(phyto_Flocs_2ways, *newDOCaggregpars)


