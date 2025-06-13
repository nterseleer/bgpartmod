# Conversion
molmass_Chl = 893.489  # [g mol-1]
molmass_C = 12.0107  # [g mol-1]
molmass_N = 14.0067  # [g mol-1]
molmass_NH4 = 18.03846  # [g mol-1]
CNredfield = 106. / 16.  # [gC gN-1]
NPredfield = 16.  # [gN gP-1]
thetaNSchartau = 1.56 / molmass_N  # [gChla gN-1]
degCtoK = 273.15
boltz = 8.6210e-5


ref_values = {
    'Phy+mu_max': {
        'reference_value': 5.2,
        'symbol': '\\mu_{max}',
        'units': 'd-1',
        'complete_name': 'Maximum growth rate'
    },
    'Phy+mortrate': {
        'reference_value': 0.05,
        'symbol': 'd^l_{Dia}',
        'units': 'd-1',
        'complete_name': 'Specific mortality rate'
    },
    'Phy+lysrate': {
        'reference_value': 0.1,
        'symbol': 'l^l_{Dia}',
        'units': 'd-1',
        'complete_name': 'Specific lysis rate'
    },
    'Phy+zeta_resp_base': {
        'reference_value': 0.01,
        'symbol': '\\zeta_{pb}',
        'units': 'd-1',
        'complete_name': 'Basal respiration rate'
    },
    'Phy+zeta_resp_prod': {
        'reference_value': 0.01,
        'symbol': '\\zeta_{pp}',
        'units': '-',
        'complete_name': 'Production specific respiration rate'
    },
    'Phy+gamma_C_exud_base': {
        'reference_value': 0.02,
        'symbol': '\\gamma_{pb}',
        'units': 'd-1',
        'complete_name': 'Basal exudation rate'
    },
    'Phy+gamma_C_exud_prod': {
        'reference_value': 0.11,
        'symbol': '\\gamma_{pp}',
        'units': '-',
        'complete_name': 'Production specific exudation rate'
    },
    'Phy+QP_min': {
        'reference_value': 0.003,
        'symbol': 'Q^P_{min}',
        'units': 'molP molC-1',
        'complete_name': 'Subsistence P quota'
    },
    'Phy+QN_min': {
        'reference_value': 0.05,
        'symbol': 'Q^N_{min}',
        'units': 'molN molC-1',
        'complete_name': 'Subsistence N quota'
    },
    'Phy+QSi_min': {
        'reference_value': 0.06,
        'symbol': 'Q^{Si}_{min}',
        'units': 'molSi molC-1',
        'complete_name': 'Subsistence Si quota'
    },
    'Phy+QP_max': {
        'reference_value': 0.012,
        'symbol': 'Q^P_{max}',
        'units': 'molP molC-1',
        'complete_name': 'Maximum quota for P'
    },
    'Phy+QN_max': {
        'reference_value': 0.15,
        'symbol': 'Q^N_{max}',
        'units': 'molN molC-1',
        'complete_name': 'Maximum quota for N'
    },
    'Phy+QSi_max': {
        'reference_value': 0.18,
        'symbol': 'Q^{Si}_{max}',
        'units': 'molSi molC-1',
        'complete_name': 'Maximum quota for Si'
    },
    'Phy+alpha': {
        'reference_value': 7e-6,
        'symbol': '\\alpha^{Chl}',
        'units': 'gC gChl-1 / (µE m-2)',
        'complete_name': 'Chl-specific slope of P-I curve'
    },
    'Phy+theta_max': {
        'reference_value': 0.07,
        'symbol': '\\theta_{max}',
        'units': 'gChl gC-1',
        'complete_name': 'Maximum Chl:C ratio'
    },
    'DOCS+alpha_TEPC': {
        'reference_value': 0.85,
        'symbol': '\\alpha^{A1}_{DOC_S-TEP}',
        'units': '-',
        'complete_name': 'DOC_S-TEP attachment probability for A1'
    },
    'DOCS+beta_TEPC': {
        'reference_value': 0.032,
        'symbol': '\\beta^{A1}_{DOC_S-TEP}',
        'units': 'm3 mmolC-1 d-1',
        'complete_name': 'DOC_S-TEP collision kernel for A1'
    },
    'DOCL+alpha_TEPC': {
        'reference_value': 1.0,
        'symbol': '\\alpha^{A1}_{DOC_L-TEP}',
        'units': '-',
        'complete_name': 'DOC_L-TEP attachment probability for A1'
    },
    'DOCL+beta_TEPC': {
        'reference_value': 0.064,
        'symbol': '\\beta^{A1}_{DOC_S-TEP}',
        'units': 'm3 mmolC-1 d-1',
        'complete_name': 'DOC_L-TEP collision kernel for A1'
    },
    'DetL+beta_max': {
        'reference_value': 0.033,
        'symbol': '\\beta^{A2}_{max}',
        'units': 'm3 mmolC-1 d-1',
        'complete_name': 'Maximum collision kernel for A2'
    },
    'TEPC+rho_TEP': {
        'reference_value': 0.1,
        'symbol': 'r^{TEP}',
        'units': 'd-1',
        'complete_name': 'Specific degradation rate of TEP'
    },
    'Microflocs+applyalphafact_PP': {
        'reference_value': 0.,
        'symbol': '\\alpha_{factor}^{PP}',
        'units': '-',
        'complete_name': "Glue's exponent for PP factor"
    },
    'Microflocs+applyalphafact_PF': {
        'reference_value': 0.,
        'symbol': '\\alpha_{factor}^{PF}',
        'units': '-',
        'complete_name': "Glue's exponent for PF factor"
    },
    'Microflocs+applyalphafact_FF': {
        'reference_value': 0.,
        'symbol': '\\alpha_{factor}^{FF}',
        'units': '-',
        'complete_name': "Glue's exponent for FF factor"
    },
    'Microflocs+nf_fractal_dim': {
        'reference_value': 2.1,
        'symbol': 'n{f}',
        'units': '-',
        'complete_name': "Fractal dimension of macroflocs [-]"
    },
    'Microflocs+alpha_PP': {
        'reference_value': 0.02,
        'symbol': '\\alpha_{PP}',
        'units': '-',
        'complete_name': "PP collision efficiency [-]"
    },
    'Microflocs+alpha_PF': {
        'reference_value': 0.02,
        'symbol': '\\alpha_{PF}',
        'units': '-',
        'complete_name': "PF collision efficiency [-]"
    },
    'Microflocs+alpha_FF': {
        'reference_value': 0.02,
        'symbol': '\\alpha_{FF}',
        'units': '-',
        'complete_name': "FF collision efficiency [-]"
    },
    'Macroflocs+sinking_leak': {
        'reference_value': 0.,
        'symbol': 'sinking_{leak}',
        'units': 'm^{¸1}',
        'complete_name': "Sinking leak [m^{¸1}]"
    },
    'Microflocs+f_frac_floc_break': {
        'reference_value': 0.1,
        'symbol': 'f',
        'units': '-',
        'complete_name': "Fraction of microflocs by breakage [-]"
    },
    'Microflocs+efficiency_break': {
        'reference_value': 1.5e-4,
        'symbol': 'E_b',
        'units': 's^{0.5}/m',
        'complete_name': "Efficiency factor for breakage [s^{0.5}/m]"
    },
    'Microflocs+fyflocstrength': {
        'reference_value': 1e-10,
        'symbol': 'F_y',
        'units': 'Pa',
        'complete_name': "Yield strength of flocs [Pa]"
    },
    'Microflocs+K_glue': {
        'reference_value': 15,
        'symbol': 'K_glue',
        'units': 'mmol TEP_C m-3',
        'complete_name': "Half saturation constant for TEP effect on flocculation"
    },
    'Microflocs+deltaFymax': {
        'reference_value': 1e-9,
        'symbol': 'deltaFymax',
        'units': '---',
        'complete_name': "deltaFymax"
    }


}






"""
Definition of output variables: Each output variable refers to a dictionary potentially containing:
    - 'units': units used for output (not necessarily the same as in the model)
    - 'munits': units used in the model while running
    - 'trsfrm': transformation factor to convert from model units to output units
    - 'cleanname': Short name written in a format to be interpreted by bgpartfunctions.cleantext()
    - 'longname': Long name 
    - 'oprt': operation to compute the output variable. To be interpreted by bgpartfunctions.eval_expr()

Notes:
    - When no computation is needed to access the variable (i.e., no 'oprt'), it is not needed to explicitely define
    the model equivalent of a variable. E.g. : "mPhy_C" is directly available from "Phy_C" infos (hence the need to 
    define 'units' and 'munits', but "totmPhy_C" has to be defined entirely here and cannot be estimated automatically
    from "totPhy_C".

"""
doutput = {"Phy_C": {'units': 'mmol C m-3',
                     'munits': 'mmol C m-3',
                     'longname': 'Phytoplankton C'},
           "Phy_Chl": {'units': 'mg Chla m-3',
                       'cleanname': 'Phy_{Chla}',
                       'longname': 'Phytoplankton Chla'},
           "Phy_N": {'units': 'mmol N m-3',
                     'munits': 'mmol N m-3',
                     'longname': 'Phytoplankton N'},
           "Het_C": {'units': 'mmol C m-3',
                     'munits': 'mmol C m-3',
                     'longname': 'Heterotrophs C'},
           "Het_N": {'units': 'mmol N m-3',
                     'munits': 'mmol N m-3',
                     'longname': 'Heterotrophs N'},
           "Det_C": {'units': 'mmol C m-3',
                     'munits': 'mmol C m-3',
                     'longname': 'Detritus C'},
           "Det_N": {'units': 'mmol N m-3',
                     'munits': 'mmol N m-3',
                     'longname': 'Detritus N'},
           "BacA_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'cleanname': 'BAC_{attached}',
                      'longname': 'Attached Bacteria'},
           "BacF_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'cleanname': 'BAC_{free}',
                      'longname': 'Free-living Bacteria'},
           "HF_C": {'units': 'mmol C m-3',
                    'munits': 'mmol C m-3',
                    'cleanname': 'Flagellates',
                    'longname': 'Heterotrophic Flagellates'},
           "Cil_C": {'units': 'mmol C m-3',
                     'munits': 'mmol C m-3',
                     'cleanname': 'Ciliates',
                     'longname': 'Ciliates'},
           "NH_4": {'units': 'mmol N m-3',
                    'munits': 'mmol N m-3',
                    'perPhyto': False,
                    'longname': 'Ammonium'},
           "NO3_concentration": {'units': 'mmol N m-3',
                           'munits': 'mmol N m-3',
                           'perPhyto': False,
                           'cleanname': 'NO_3',
                           'longname': 'NO_3'
                                 },
           "NH4_concentration": {'units': 'mmol N m-3',
                           'munits': 'mmol N m-3',
                           'perPhyto': False,
                           'cleanname': 'NH_4',
                           'longname': 'NH_4'
                                 },
           "DIP_concentration": {'units': 'mmol P m-3',
                           'munits': 'mmol P m-3',
                           'perPhyto': False,
                           'cleanname': 'DIP',
                           'longname': 'Dissolved Inorganic P'
                                 },
           "DIN_concentration": {'units': 'mmol N m-3',
                           'munits': 'mmol N m-3',
                                 'oprt': 'mNO3_concentration+mNH4_concentration',
                                 'perPhyto': False,
                           'cleanname': 'DIN',
                           'longname': 'Dissolved Inorganic N'
                                 },
           "DSi_concentration": {'units': 'mmol Si m-3',
                           'munits': 'mmol Si m-3',
                           'perPhyto': False,
                           'cleanname': 'DSi',
                           'longname': 'Dissolved Inorganic Si'
                                 },
           "DOCS_C": {'units': 'mmol C m-3',
                           'munits': 'mmol C m-3',
                           'cleanname': 'DOC_{small}',
                           'longname': 'Dissolved Organic C (small)'
                                 },
           "DOCL_C": {'units': 'mmol C m-3',
                           'munits': 'mmol C m-3',
                           'cleanname': 'DOC_{large}',
                           'longname': 'Dissolved Organic C (large)'
                                 },
           "Macroflocs_settling_vel": {'units': 'm d-1',
                           'munits': 'm d-1',
                           'cleanname': 'Settling velocity',
                           'longname': 'Flocs settling velocity'
                                 },


           "DIM_DINconc": {'units': 'mmol N m-3',
                           'munits': 'mmol N m-3',
                           'perPhyto': False,
                           'cleanname': 'DIN_{concentration}',
                           'longname': 'Dissolved Inorganic Nitrogen'},
           "DIM_DICconc": {'units': 'mmol C m-3',
                           'munits': 'mmol C m-3',
                           'perPhyto': False,
                           'cleanname': 'DIC_{concentration}',
                           'longname': 'Dissolved Inorganic Carbon'},
           "DOM_DONconc": {'units': 'mmol N m-3',
                           'munits': 'mmol N m-3',
                           'perPhyto': False,
                           'longname': 'Dissolved Organic NItrogen'},
           "DOM_PCHOconc": {'units': 'mmol C m-3',
                            'munits': 'mmol C m-3',
                            'cleanname': 'PCHO',
                            'perPhyto': False,
                            'longname': 'Dissolved Polysaccharides Carbon'},
           "TEPC_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'cleanname': 'TEP',
                      'perPhyto': False,
                      'longname': 'Transparent Exopolymer Particles'},
           "DIM_TAconc": {'units': '[]',
                          'munits': '[]]',
                          'perPhyto': False,
                          'longname': 'Total Alkalinity'},
           "DOM_resDOCconc": {'units': 'mmol C m-3',
                              'munits': 'mmol C m-3',
                              'perPhyto': False,
                              'longname': 'Residual Dissolved Carbon'},

           # "I_t": {'units': 'µmol quanta m-2 d-1',
           #         'oprt': 'setup.PAR /3600./24 * np.exp(-setup.k_att * SUMALLmPhy_Chl * setup.z)',
           #         'longname': 'Ambiant light level (during Light phases)\n'},
           "QN": {'units': 'molN molC-1',
                  'oprt': 'mPhy_N/mPhy_C',
                  'cleanname': 'Q^N',
                  'longname': 'Phytoplankton N:C ratio'},
           "CN": {'units': 'molC molN-1',
                  'oprt': 'mPhy_C/mPhy_N',
                  'cleanname': 'C:N',
                  'longname': 'Phytoplankton C:N ratio'},
           "CP": {'units': 'molC molP-1',
                  'oprt': 'mPhy_C/mPhy_P',
                  'cleanname': 'C:P',
                  'longname': 'Phytoplankton C:P ratio'},
           "CSi": {'units': 'molC molSi-1',
                  'oprt': 'mPhy_C/mPhy_Si',
                  'cleanname': 'C:Si',
                  'longname': 'Phytoplankton C:Si ratio'},
           "QP": {'units': 'molP molC-1',
                  'oprt': 'mPhy_P/mPhy_C',
                  'cleanname': 'Q^P',
                  'longname': 'Phytoplankton P:C ratio'},
           "QSi": {'units': 'molSi molC-1',
                   'oprt': 'mPhy_N/mPhy_C',
                   'cleanname': 'Q^{Si}',
                   'longname': 'Phytoplankton Si:C ratio'},
           "BacA_QN": {'units': 'gN gC-1',
                       'oprt': 'mBacA_N/mBacA_C',
                       'cleanname': 'Q^N',
                       'longname': 'BacA N:C ratio'},
           "BacF_QN": {'units': 'gN gC-1',
                       'oprt': 'mBacF_N/mBacF_C',
                       'cleanname': 'Q^N',
                       'longname': 'BacF N:C ratio'},
           "HF_QN": {'units': 'gN gC-1',
                     'oprt': 'mHF_N/mHF_C',
                     'cleanname': 'Q^N',
                     'longname': 'HF N:C ratio'},
           "Cil_QN": {'units': 'gN gC-1',
                      'oprt': 'mCil_N/mCil_C',
                      'cleanname': 'Q^N',
                      'longname': 'Cil N:C ratio'},

           "thetaN": {'units': 'gChla gN-1',
                      'oprt': 'mPhy_Chl/mPhy_N*14',
                      'cleanname': '\\theta^N',
                      'longname': 'Phytoplankton Chla:N ratio'},
           "thetaC": {'units': 'gChla gC-1',
                      'oprt': 'mPhy_Chl/mPhy_C/12',
                      'cleanname': '\\theta^C',
                      'longname': 'Phytoplankton Chla:C ratio'},
           "CChl": {'units': 'gC gChla-1',
                    'oprt': 'mPhy_C*12/mPhy_Chl',
                    'cleanname': 'C:Chl',
                    'longname': 'Phytoplankton C:Chla ratio'},
           "totPhy_Chl": {'units': 'µg Chla l-1',
                          'oprt': 'SUMALL(Phy_Chl)',
                          'cleanname': 'Total Phy_{Chla}',
                          'longname': 'Total Chla'},
           "totPhy_C": {'units': 'mmol C m-3',
                        'oprt': 'SUMALL(Phy_C)',
                        'cleanname': 'Total Phy_C',
                        'longname': 'Total Phytoplankton Carbon'},
           "totmPhy_C": {'units': 'mmol C m-3',
                         'oprt': 'SUMALL(mPhy_C)',
                         'cleanname': 'Total Phy_C',
                         'longname': 'Total Phytoplankton Carbon'},
           "totPhy_N": {'units': 'mmol N m-3',
                        'oprt': 'SUMALL(Phy_N)',
                        'cleanname': 'Total Phy_N',
                        'longname': 'Total Phytoplankton Nitrogen'},
           "totmPhy_N": {'units': 'mmol N m-3',
                         'oprt': 'SUMALL(mPhy_N)',
                         'cleanname': 'Total mPhy_N',
                         'longname': 'Total Phytoplankton Nitrogen'},
           "C_tot": {'units': 'mmol C m-3',
                     'cleanname': 'C_{tot}^{Sys}',
                     'longname': 'Total C in the system'},
           "N_tot": {'units': 'mmol N m-3',
                     'cleanname': 'N_{tot}^{Sys}',
                     'longname': 'Total N in the system'},
           "PON": {'units': 'mmol N m-3',
                   'longname': 'PON'},
           "POC": {'units': 'mmolC m-3',
                   'longname': 'Particulate Organic C'},
           "POC:PON": {'units': 'molC molN-1',
                       'oprt': 'POC/PON',
                       'longname': 'POC:PON'},
           "Microflocs_numconc": {'units': 'm-3',
                                'cleanname': 'Microflocs',
                                'munits': 'm-3',
                                'longname': 'Microflocs'},

           # "Microflocs_concentration": {'units': 'g m-3',
           #                              'oprt': 'Microflocs_numconc * np.pi / 6. ((Micro_in_Macro_numconc/Macroflocs_numconc)**(1/2.1)*5e-6*1e6)*((Micro_in_Macro_numconc/Macroflocs_numconc)**(1/2.1)*5e-6*1e6)*((Micro_in_Macro_numconc/Macroflocs_numconc)**(1/2.1)*5e-6*1e6)',
           #                          'cleanname': 'Microflocs',
           #                          'munits': 'g m-3',
           #                          'longname': 'Microflocs'},
           "Macroflocs_numconc": {'units': 'm-3',
                                'cleanname': 'Macroflocs',
                                'munits': 'm-3',
                                'longname': 'Macroflocs'},
           "Micro_in_Macro_numconc": {'units': 'm-3',
                                    'cleanname': 'Microflocs in flocs',
                                    'munits': 'm-3',
                                    'longname': 'Microflocs in flocs'},
           "Floc_diam": {'units': 'µm',
                         'cleanname': 'Floc_{diam}',
                         # 'oprt': '(Micro_in_Macro_numconc/Macroflocs_numconc)**(1/2.1)*5e-6*1e6', # should not be
                         # used since Macroflocs_diam is now available (with non-static fractal num etc)
                         'oprt': 'Macroflocs_diam*1e6',
                         'longname': 'Floc diameter'},
           "Microflocs_massconcentration": {'units': 'mg l-1',
                         'cleanname': 'MassConcentration_{free microflocs}',
                         'trsfrm': 1e3,
                         'longname': 'Free Microflocs mass concentration'},
            "Micro_in_Macro_massconcentration": {'units': 'mg l-1',
                         'cleanname': 'MassConcentration_{microflocs in flocs}',
                         'trsfrm': 1e3,
                         'longname': 'Microflocs in flocs mass concentration'},
            "SPMC": {'units': 'mg l-1',
                         'cleanname': 'SPMC_{mineral}',
                         'trsfrm': 1e3,
                         'oprt': 'Micro_in_Macro_massconcentration + Microflocs_massconcentration',
                         'longname': 'Suspended Particulate Matter Concentration'},
            "Macroflocs_fyflocstrength": {'units': '-',
                         'cleanname': 'F_y floc strenght',
                         'longname': 'Yield Floc strength'},
           "Macroflocs_diam": {'units': 'µm',
                         'cleanname': 'Floc_{diam}',
                               'trsfrm': 1e6,
                         'longname': 'Floc diameter'},
            # "Microflocs_TOT": {'units': 'm-3',
            #              'cleanname': 'SPMC_{mineral}',
            #              'oprt': 'Micro_in_Macro_numconc + Microflocs_numconc',
            #              'longname': 'Tot free + flocs Microflocs numeric concentration'},
           # "Ffrac": {'units': '%',
           #           'oprt': 'Floc_Micro_in_Macronumconc/(Floc_Micro_in_Macronumconc+Floc_Microflocsnconc)',
           #           'longname': 'Proportion of microflocs in flocs'},
           }

# OLD 20231010 - Change model units for Sch07 - below is used for GMK98
doutput_MASS = {"Phy_C": {'units': 'µM',
                          'munits': 'gC m-3',
                          'trsfrm': 1e3 / molmass_C,
                          'longname': 'Phytoplankton C'},
                "Phy_Chl": {'units': 'µg Chla l-1',
                            'trsfrm': 1e3,
                            'cleanname': 'Phy_{Chla}',
                            'longname': 'Phytoplankton Chla'},
                "Phy_Chla": {'units': 'µg Chla l-1',
                             'trsfrm': 1e3,
                             'cleanname': 'Phy_{Chla}',
                             'longname': 'Phytoplankton Chla'},
                "Phy_N": {'units': 'µM',
                          'munits': 'gN m-3',
                          'trsfrm': 1e3 / molmass_N,
                          'longname': 'Phytoplankton N'},
                "Het_C": {'units': 'µM',
                          'munits': 'gC m-3',
                          'trsfrm': 1e3 / molmass_C,
                          'longname': 'Heterotrophs C'},
                "Het_N": {'units': 'µM',
                          'munits': 'gN m-3',
                          'trsfrm': 1e3 / molmass_N,
                          'longname': 'Heterotrophs N'},
                "Det_C": {'units': 'µM',
                          'munits': 'gC m-3',
                          'trsfrm': 1e3 / molmass_C,
                          'longname': 'Detritus C'},
                "Det_N": {'units': 'µM',
                          'munits': 'gN m-3',
                          'trsfrm': 1e3 / molmass_N,
                          'longname': 'Detritus N'},
                "NH_4": {'units': 'µM',
                         'munits': 'gN m-3',
                         'trsfrm': 1e3 / molmass_N,
                         'perPhyto': False,
                         'longname': 'Ammonium'},
                "DIM_DINconc": {'units': 'µM',
                                'munits': 'gN m-3',
                                'trsfrm': 1e3 / molmass_N,
                                'perPhyto': False,
                                'longname': 'Dissolved Inorganic Nitrogen'},
                "DIM_DICconc": {'units': 'µM',
                                'munits': 'gC m-3',
                                'trsfrm': 1e3 / molmass_C,
                                'perPhyto': False,
                                'longname': 'Dissolved Inorganic Carbon'},
                "DOM_DONconc": {'units': 'µM',
                                'munits': 'gN m-3',
                                'trsfrm': 1e3 / molmass_N,
                                'perPhyto': False,
                                'longname': 'Dissolved Organic NItrogen'},
                "DOM_PCHOconc": {'units': 'µM',
                                 'munits': 'gC m-3',
                                 'trsfrm': 1e3 / molmass_C,
                                 'perPhyto': False,
                                 'longname': 'Dissolved Polysaccharides Carbon'},
                "DOM_TEPCconc": {'units': 'µM',
                                 'munits': 'gC m-3',
                                 'trsfrm': 1e3 / molmass_C,
                                 'perPhyto': False,
                                 'longname': 'Transparent Exopolymer Particles Carbon'},
                "DIM_TAconc": {'units': '[]',
                               'munits': '[]]',
                               'trsfrm': 1,
                               'perPhyto': False,
                               'longname': 'Total Alkalinity'},
                "DOM_resDOCconc": {'units': 'µM',
                                   'munits': 'gC m-3',
                                   'trsfrm': 1e3 / molmass_C,
                                   'perPhyto': False,
                                   'longname': 'Residual Dissolved Carbon'},

                "I_t": {'units': 'µmol quanta m-2 d-1',
                        'oprt': 'setup.I /3600./24 * np.exp(-setup.k_att * SUMALLmPhy_Chl * setup.z)',
                        'longname': 'Ambiant light level (during Light phases)\n'},
                "QN": {'units': 'gN gC-1',
                       'oprt': 'mPhy_N/mPhy_C',
                       'cleanname': 'Q^N',
                       'longname': 'Phytoplankton N:C ratio'},
                "thetaN": {'units': 'gChla gN-1',
                           'oprt': 'mPhy_Chl/mPhy_N',
                           'cleanname': '\\theta^N',
                           'longname': 'Phytoplankton Chla:N ratio'},
                "thetaC": {'units': 'gChla gC-1',
                           'oprt': 'mPhy_Chl/mPhy_C',
                           'cleanname': '\\theta^C',
                           'longname': 'Phytoplankton Chla:C ratio'},
                "totPhy_Chl": {'units': 'µg Chla l-1',
                               'oprt': 'SUMALL(Phy_Chl)',
                               'cleanname': 'Total Phy_{Chla}',
                               'longname': 'Total Chla'},
                "totPhy_C": {'units': 'µM',
                             'oprt': 'SUMALL(Phy_C)',
                             'cleanname': 'Total Phy_C',
                             'longname': 'Total Phytoplankton Carbon'},
                "totmPhy_C": {'units': 'gC m-3',
                              'oprt': 'SUMALL(mPhy_C)',
                              'cleanname': 'Total Phy_C',
                              'longname': 'Total Phytoplankton Carbon'},
                "totPhy_N": {'units': 'µM',
                             'oprt': 'SUMALL(Phy_N)',
                             'cleanname': 'Total Phy_N',
                             'longname': 'Total Phytoplankton Nitrogen'},
                "totmPhy_N": {'units': 'gN m-3',
                              'oprt': 'SUMALL(mPhy_N)',
                              'cleanname': 'Total mPhy_N',
                              'longname': 'Total Phytoplankton Nitrogen'},
                # "Ntot": {'units': 'gN m-3',
                #          'oprt': 'mNH_4 + SUMALLmPhy_N',
                #          'cleanname': 'N_{tot}^{Sys}',
                #          'longname': 'Total N in the system'},
                "Ctot": {'units': 'gC m-3',
                         'oprt': 'DIM_DICconc + SUMALLmPhy_C',
                         'cleanname': 'C_{tot}^{Sys}',
                         'longname': 'Total C in the system'},
                "C_tot": {'units': 'gC m-3',
                          'cleanname': 'C_{tot}^{Sys}',
                          'longname': 'Total C in the system'},
                "N_tot": {'units': 'gN m-3',
                          'cleanname': 'N_{tot}^{Sys}',
                          'longname': 'Total N in the system'}
                }
