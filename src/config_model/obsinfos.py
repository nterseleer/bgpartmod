from . import varinfos

"""
REF UNITS = MODEL UNITS
- trsfrm: factor applied to transform observation data into model units
- oprt: to build a new variable from existing variable (e.g. a stoichiometric ratio)
"""

mow1fname = 'DATA_AllStations.txt'

# UNITS IN DATA_AllStations.txt
# {'SPMC': 'mg/l', 'error': '%', 'POC': 'mg/l', 'error.1': '%', 'PON': 'mg/l', 'error.2': '%',
#  'ChlA': 'µg/l', 'error.3': '%', 'ChlB': 'µg/l', 'error.4': '%', 'PheoA': 'µg/l', 'error.5': '%',
#  'PheoB': 'µg/l', 'error.6': '%', 'TEP': 'mgXG/l', 'error.7': '%', 'Salinit': '°C', 'T@pH': 'mg/l',
#  'DOC': 'mg/l', 'DIC': 'µmol/l', 'TP': 'µmol/l', 'TN': 'µmol/l', 'NOx': 'µmol/l', 'NO2': 'µmol/l',
#  'NH4': 'µmol/l', 'PO4': 'µmol/l', 'Si': 'mg/l', 'Alkal': nan, 'pH': nan}


mow1conv = {'ChlA': {'correspondingVars': ['Phy_Chl'],
                     # 'trsfrm': 1 # µg l-1 = mg m-3 OK
                     },
            'NH4': {'correspondingVars': ['NH4_concentration'],
                    # 'trsfrm': 1 # µmol l-1 = mmol m-3 OK
                    },
            'PO4': {'correspondingVars': ['DIP_concentration'],
                    # 'trsfrm': 1 # µmol l-1 = mmol m-3 OK
                    },
            'Si': {'correspondingVars': ['DSi_concentration'],
                   # 'trsfrm': 1 # µmol l-1 = mmol m-3 OK
                   },
            'TEP': {'correspondingVars': ['TEPC_C'],
                    },
            'TN': {'correspondingVars': ['N_tot'],
                   # 'trsfrm': 1 # µmol l-1 = mmol m-3 OK
                   },
            'TP': {'correspondingVars': ['P_tot'],
                   # 'trsfrm': 1 # µmol l-1 = mmol m-3 OK
                   },
            'SPMC': {'correspondingVars': ['SPMC'],
                     # 'trsfrm': 1e3 # mg m-3 / varinfos.molmass_C  # mg l-1 to mmol m-3
                     },
            'POC': {'correspondingVars': ['POC'],
                    'trsfrm': 1e3 / varinfos.molmass_C  # mg l-1 to mmol m-3
                    },
            'PON': {'correspondingVars': ['PON'],
                    'trsfrm': 1e3 / varinfos.molmass_N  # mg l-1 to mmol m-3
                    },
            'DOC': {'correspondingVars': ['DOC'],
                    'trsfrm': 1e3 / varinfos.molmass_C  # mg l-1 to mmol m-3
                    },
            'DIN': {'correspondingVars': ['DIN_concentration'],
                    # 'oprt': 'NH4+NOx'
                    },
            'NOx': {'correspondingVars': ['NO3_concentration'],
                    # 'oprt': 'NH4+NOx'
                    },



            # 'thetaN': {'oprt': 'mPhy_Chl/mPhy_N/1e3'},
            # 'thetaC': {'oprt': 'mPhy_Chl/mPhy_C/1e3'}
            }

mow1tepname = 'MOW1_TEPf + TEPm estimates.txt'
""" 
UNITS
mg XG eq/l --> cf Fettweis et al 2022 Appendix C:
For fresh TEP, 1 mgC/l ~ 0.55 mgXGeq/l
(or ~0.57 for overall TEP concentrations)
"""
mow1tepconv = {'Tep_fresh': {'correspondingVars': ['TEPC_C'],
                             'trsfrm': 1e3 * 0.55 / varinfos.molmass_C   # mg XG eq/l to mmol C/m3
                             },
               }

# mow1spmcname = 'MOW1_SPMC (2005-2016).txt'
# mow1spmcconv = {'SPMC': {'correspondingVars': ['SPMC'],
#                              },
#                }

mow1d50name = 'MOW1_D50-floc size (2007-2019).txt'
mow1d50conv = {'D50': {'correspondingVars': ['Floc_diam'],
                       },
               }

st330conv = {
    'DA': {'correspondingVars': ['Phy_C'],
           'trsfrm': 1 / varinfos.molmass_C},   # mgC/m3 to mmol C/m3
    'NF': {'correspondingVars': ['HF_C'],
           'trsfrm': 1 / varinfos.molmass_C},   # mgC/m3 to mmol C/m3
    # 'OP': {'correspondingVars': ['Phy_Phaeocystis_C'],
    #            'trsfrm': 1 / varinfos.molmass_C},   # mgC/m3 to mmol C/m3
    'MZ': {'correspondingVars': ['Cil_C'],
           'trsfrm': 1 / varinfos.molmass_C},   # mgC/m3 to mmol C/m3
    # 'CP': {'correspondingVars': ['Copepods_C'],
    #            'trsfrm': 1 / varinfos.molmass_C},   # mgC/m3 to mmol C/m3
    'BC': {'correspondingVars': ['BacF_C'],
           'trsfrm': 1 / varinfos.molmass_C},   # mgC/m3 to mmol C/m3
    'NO3': {'correspondingVars': ['NO3_concentration']}, # OK mmol/m3
    'NH4': {'correspondingVars': ['NH4_concentration']}, # OK mmol/m3
    'PO4': {'correspondingVars': ['DIP_concentration']}, # OK mmol/m3
    'SiO': {'correspondingVars': ['DSi_concentration']}, # OK mmol/m3
}

# Flynn et al 1994 used in Geider et al 1998 Fig. 12
flynnfname = 'Flynn et al 1994.xlsx'
flynnconv = {'Ammonium': {'correspondingVars': ['DIN', 'DIN_conc', 'DIM_DINconc'],
                          'trsfrm': 1e3 / varinfos.molmass_N},
             'Cell C': {'correspondingVars': ['Phy_C', 'totPhy_C'],
                        'trsfrm': 1e3 / varinfos.molmass_C},
             'Cell N': {'correspondingVars': ['Phy_N', 'totPhy_N'],
                        'trsfrm': 1e3 / varinfos.molmass_N},
             'Chlorophyll': {'correspondingVars': ['Phy_Chl', 'totPhy_Chl']},
             'QN': {'oprt': 'mPhy_N/mPhy_C'},
             'thetaN': {'oprt': 'mPhy_Chl/mPhy_N/1e3'},
             'thetaC': {'oprt': 'mPhy_Chl/mPhy_C/1e3'}
             }
