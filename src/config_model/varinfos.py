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
    'Phy+divide_water_depth_ratio': {
        'reference_value': 1.,
        'symbol': 'depth^{water}_{ratio}',
        'units': '[-]',
        'complete_name': 'Water Depth Reduction Ratio for Phytoplankton'
    },
    'Phy+grazing_loss_max': {
        'reference_value': 0.,
        'symbol': 'grazing^{max}_{loss}',
        'units': '[mmolC m-3 d-1]',
        'complete_name': 'Max Grazing loss by zooplankton'
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
    'DetL+kleak': {
        'reference_value': 0.0,
        'symbol': 'k_{leak}^{DetL}',
        'units': 'd-1',
        'complete_name': 'DetL leakage rate'
    },
    'DetL+beta_max': {
        'reference_value': 0.033,
        'symbol': '\\beta^{A2}_{max}',
        'units': 'm3 mmolC-1 d-1',
        'complete_name': 'Maximum collision kernel for A2'
    },
    'TEPC+rho_TEP': {
        'reference_value': 0.1,
        'symbol': '\\rho^{TEP}',
        'units': 'd-1',
        'complete_name': 'Specific degradation rate of TEP'
    },
    'TEPC+kleak': {
        'reference_value': 0.,
        'symbol': 'k_{leak}^{TEP}',
        'units': 'd-1',
        'complete_name': 'TEP leakage rate'
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
    'Microflocs+alpha_PP_base': {
        'reference_value': 0.002,
        'symbol': '\\alpha_{PP}^{base}',
        'units': '-',
        'complete_name': "Base PP collision efficiency (mineral only) [-]"
    },
    'Microflocs+delta_alpha_PP': {
        'reference_value': 0.028,
        'symbol': '\\Delta\\alpha_{PP}',
        'units': '-',
        'complete_name': "TEP increment for PP collision efficiency [-]"
    },
    'Microflocs+alpha_FF_base': {
        'reference_value': 0.002,
        'symbol': '\\alpha_{FF}^{base}',
        'units': '-',
        'complete_name': "Base FF collision efficiency (mineral only) [-]"
    },
    'Microflocs+delta_alpha_FF': {
        'reference_value': 0.028,
        'symbol': '\\Delta\\alpha_{FF}',
        'units': '-',
        'complete_name': "TEP increment for FF collision efficiency [-]"
    },
    'Microflocs+alpha_PF_base': {
        'reference_value': 0.002,
        'symbol': '\\alpha_{PF}^{base}',
        'units': '-',
        'complete_name': "Base PF collision efficiency (mineral only) [-]"
    },
    'Microflocs+delta_alpha_PF': {
        'reference_value': 0.028,
        'symbol': '\\Delta\\alpha_{PF}',
        'units': '-',
        'complete_name': "TEP increment for PF collision efficiency [-]"
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
    'Microflocs+fyflocstrength_base': {
        'reference_value': 1e-10,
        'symbol': 'F_y^{base}',
        'units': 'N',
        'complete_name': "Base yield strength of flocs (mineral only) [N]"
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
    },
    'Microflocs+eps_kd': {
        'reference_value': 0.066 * 1e3,
        'symbol': 'eps_{kd}',
        'units': 'm-1 (mg l-1)-1',
        'complete_name': "SPM extinction coefficient"
    },
    'Microflocs+tau_cr_base': {
        'reference_value': 0.5,
        'symbol': '\\tau_{cr}^{base}',
        'units': 'Pa',
        'complete_name': "Base Critical shear stress"
    },
    'Microflocs+delta_tau_cr': {
        'reference_value': 1.,
        'symbol': 'Δ \\tau_{cr}^{base}',
        'units': 'Pa',
        'complete_name': "Max TEP-increased Critical shear stress"
    },
    'Microflocs+delta_nf_fractal_dim': {
        'reference_value': 0.,
        'symbol': 'Δnf',
        'units': '-',
        'complete_name': "Max TEP-induced decrease in fractal dimension"
    },


    'Macroflocs+resuspension_rate': {
        'reference_value': 1e5,
        'symbol': 'rate_{resusp}',
        'units': '[]',
        'complete_name': "Resuspension rate"
    },
    'NH4+k_remin': {
        'reference_value': 0,
        'symbol': 'K_remin^{NH_4}',
        'units': '--',
        'complete_name': "k_remin_NH4"
    },
    'DIP+k_remin': {
        'reference_value': 0,
        'symbol': 'K_remin^{DIP}',
        'units': '--',
        'complete_name': "k_remin_DIP"
    },
    'DSi+k_remin': {
        'reference_value': 0,
        'symbol': 'K_remin^{DSi}',
        'units': '--',
        'complete_name': "k_remin_DSi"
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
           "Phy_P": {'units': 'mmol P m-3',
                     'munits': 'mmol P m-3',
                     'longname': 'Phytoplankton P'},
           "Phy_Si": {'units': 'mmol NSi m-3',
                      'munits': 'mmol Si m-3',
                      'longname': 'Phytoplankton Si'},
           "Het_C": {'units': 'mmol C m-3',
                     'munits': 'mmol C m-3',
                     'longname': 'Heterotrophs C'},
           "Het_N": {'units': 'mmol N m-3',
                     'munits': 'mmol N m-3',
                     'longname': 'Heterotrophs N'},
           "BacA_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'cleanname': 'BAC_{attached}',
                      'longname': 'Attached Bacteria'},
           "BacA_N": {'units': 'mmol N m-3',
                      'munits': 'mmol N m-3',
                      'longname': 'Attached bacteria nitrogen',
                      'cleanname': 'BacA N'},

           "BacA_P": {'units': 'mmol P m-3',
                      'munits': 'mmol P m-3',
                      'longname': 'Attached bacteria phosphorus',
                      'cleanname': 'BacA P'},
           "BacF_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'cleanname': 'BAC_{free}',
                      'longname': 'Free-living Bacteria'},
           "BacF_N": {'units': 'mmol N m-3',
                      'munits': 'mmol N m-3',
                      'longname': 'Free-living bacteria nitrogen',
                      'cleanname': 'BacF N'},
           "BacF_P": {'units': 'mmol P m-3',
                      'munits': 'mmol P m-3',
                      'longname': 'Free-living bacteria phosphorus',
                      'cleanname': 'BacF P'},
           "HF_C": {'units': 'mmol C m-3',
                    'munits': 'mmol C m-3',
                    'cleanname': 'Flagellates',
                    'longname': 'Heterotrophic Flagellates'},
           "HF_N": {'units': 'mmol N m-3',
                    'munits': 'mmol N m-3',
                    'longname': 'Heterotrophic flagellates nitrogen',
                    'cleanname': 'HF N'},
           "HF_P": {'units': 'mmol P m-3',
                    'munits': 'mmol P m-3',
                    'longname': 'Heterotrophic flagellates phosphorus',
                    'cleanname': 'HF P'},

           "Cil_C": {'units': 'mmol C m-3',
                     'munits': 'mmol C m-3',
                     'cleanname': 'Ciliates',
                     'longname': 'Ciliates'},
           "Cil_N": {'units': 'mmol N m-3',
                     'munits': 'mmol N m-3',
                     'longname': 'Ciliates nitrogen',
                     'cleanname': 'Cil N'},
           "Cil_P": {'units': 'mmol P m-3',
                     'munits': 'mmol P m-3',
                     'longname': 'Ciliates phosphorus',
                     'cleanname': 'Cil P'},

           "DetS_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'longname': 'Small detritus carbon',
                      'cleanname': 'DetS C'},

           "DetS_N": {'units': 'mmol N m-3',
                      'munits': 'mmol N m-3',
                      'longname': 'Small detritus nitrogen',
                      'cleanname': 'DetS N'},

           "DetS_P": {'units': 'mmol P m-3',
                      'munits': 'mmol P m-3',
                      'longname': 'Small detritus phosphorus',
                      'cleanname': 'DetS P'},

           "DetL_C": {'units': 'mmol C m-3',
                      'munits': 'mmol C m-3',
                      'longname': 'Large detritus carbon',
                      'cleanname': 'DetL C'},

           "DetL_N": {'units': 'mmol N m-3',
                      'munits': 'mmol N m-3',
                      'longname': 'Large detritus nitrogen',
                      'cleanname': 'DetL N'},

           "DetL_P": {'units': 'mmol P m-3',
                      'munits': 'mmol P m-3',
                      'longname': 'Large detritus phosphorus',
                      'cleanname': 'DetL P'},

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
           "Macroflocs_settling_vel": {'units': 'm s-1',
                                       'munits': 'm s-1',
                                       'cleanname': 'Settling velocity',
                                       'longname': 'Flocs settling velocity'
                                       },
           "Macroflocs_settling_vel_base": {'units': 'mm s-1',
                                            'munits': 'm s-1',
                                            'trsfrm': 1e3,
                                            'cleanname': 'Settling velocity',
                                            'longname': 'Winterwerp\'s Flocs settling velocity'
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
           #         'oprt': 'setup.PAR /3600./24 * np.exp(-setup.k_att * SUMALLmPhy_Chl * setup.water_depth)',
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

           "Macroflocs_numconc": {'units': 'm-3',
                                  'cleanname': 'Macroflocs',
                                  'munits': 'm-3',
                                  'longname': 'Macroflocs'},
           "Micro_in_Macro_numconc": {'units': 'm-3',
                                      'cleanname': 'Microflocs in flocs',
                                      'munits': 'm-3',
                                      'longname': 'Microflocs in flocs'},

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
            "Macroflocs_nf_fractal_dim": {'units': '-',
                               'cleanname': 'n_f',
                               'longname': 'Fractal dimension'},

           "Phy_limNUT": {'units': '-',
                          'munits': '-',
                          'longname': 'Phytoplankton nutrient limitation',
                          'cleanname': 'Phy_{lim}^{NUT}'},

           "Phy_limT": {'units': '-',
                        'munits': '-',
                        'longname': 'Phytoplankton temperature limitation',
                        'cleanname': 'Phy_{lim}^{T}'},

           "Phy_limI": {'units': '-',
                        'munits': '-',
                        'longname': 'Phytoplankton light limitation',
                        'cleanname': 'Phy_{lim}^{I}'},

           "Phy_kd": {'units': 'm^{-1}',
                      'munits': 'm^{-1}',
                      'longname': 'Light attenuation coefficient',
                      'cleanname': 'k_d'},

           # Contributions to kd from different components
           "kd_contrib_Phy": {'units': 'm^{-1}',
                              'munits': 'm^{-1}',
                              'oprt': "model.components['Phy'].eps_kd * mPhy_C",
                              'longname': 'Light attenuation from Phytoplankton',
                              'cleanname': 'k_d^{Phy}'},

           "kd_contrib_DetL": {'units': 'm^{-1}',
                               'munits': 'm^{-1}',
                               'oprt': "model.components['DetL'].eps_kd * mDetL_C",
                               'longname': 'Light attenuation from Large Detritus',
                               'cleanname': 'k_d^{DetL}'},

           "kd_contrib_DetS": {'units': 'm^{-1}',
                               'munits': 'm^{-1}',
                               'oprt': "model.components['DetS'].eps_kd * mDetS_C",
                               'longname': 'Light attenuation from Small Detritus',
                               'cleanname': 'k_d^{DetS}'},

           "kd_contrib_BacA": {'units': 'm^{-1}',
                               'munits': 'm^{-1}',
                               'oprt': "model.components['BacA'].eps_kd * mBacA_C",
                               'longname': 'Light attenuation from Attached Bacteria',
                               'cleanname': 'k_d^{BacA}'},

           "kd_contrib_BacF": {'units': 'm^{-1}',
                               'munits': 'm^{-1}',
                               'oprt': "model.components['BacF'].eps_kd * mBacF_C",
                               'longname': 'Light attenuation from Free Bacteria',
                               'cleanname': 'k_d^{BacF}'},

           "kd_contrib_HF": {'units': 'm^{-1}',
                             'munits': 'm^{-1}',
                             'oprt': "model.components['HF'].eps_kd * mHF_C",
                             'longname': 'Light attenuation from Heterotrophic Flagellates',
                             'cleanname': 'k_d^{HF}'},

           "kd_contrib_Cil": {'units': 'm^{-1}',
                              'munits': 'm^{-1}',
                              'oprt': "model.components['Cil'].eps_kd * mCil_C",
                              'longname': 'Light attenuation from Ciliates',
                              'cleanname': 'k_d^{Cil}'},

           "kd_contrib_Microflocs": {'units': 'm^{-1}',
                                     'munits': 'm^{-1}',
                                     'oprt': "model.components['Microflocs'].eps_kd * Microflocs_massconcentration",
                                     'longname': 'Light attenuation from free Microflocs (SPM)',
                                     'cleanname': 'k_d^{Microflocs}'},

           "kd_contrib_Micro_in_Macro": {'units': 'm^{-1}',
                                         'munits': 'm^{-1}',
                                         'oprt': "model.components['Micro_in_Macro'].eps_kd * Micro_in_Macro_massconcentration",
                                         'longname': 'Light attenuation from Microflocs in Macroflocs (SPM)',
                                         'cleanname': 'k_d^{Micro in Macro}'},

           "Phy_lim_N": {'units': '-',
                         'munits': '-',
                         'longname': 'Phytoplankton N limitation',
                         'cleanname': 'Phy_{lim}^{N}'},

           "Phy_lim_P": {'units': '-',
                         'munits': '-',
                         'longname': 'Phytoplankton P limitation',
                         'cleanname': 'Phy_{lim}^{P}'},

           "Phy_lim_Si": {'units': '-',
                          'munits': '-',
                          'longname': 'Phytoplankton Si limitation',
                          'cleanname': 'Phy_{lim}^{Si}'},

           "Phy_limQUOTA.N": {'units': '-',
                              'munits': '-',
                              'longname': 'Phytoplankton N quota limitation',
                              'cleanname': 'Phy_{lim}^{QN}'},

           "Phy_limQUOTA.P": {'units': '-',
                              'munits': '-',
                              'longname': 'Phytoplankton P quota limitation',
                              'cleanname': 'Phy_{lim}^{QP}'},

           "Phy_limQUOTA.Si": {'units': '-',
                               'munits': '-',
                               'longname': 'Phytoplankton Si quota limitation',
                               'cleanname': 'Phy_{lim}^{QSi}'},

           "Phy_limQUOTAmin.N": {'units': '-',
                                 'munits': '-',
                                 'longname': 'Phytoplankton N minimum quota',
                                 'cleanname': 'Phy_{lim}^{QN,min}'},

           "Phy_limQUOTAmin.P": {'units': '-',
                                 'munits': '-',
                                 'longname': 'Phytoplankton P minimum quota',
                                 'cleanname': 'Phy_{lim}^{QP,min}'},

           "Phy_limQUOTAmin.Si": {'units': '-',
                                  'munits': '-',
                                  'longname': 'Phytoplankton Si minimum quota',
                                  'cleanname': 'Phy_{lim}^{QSi,min}'},

           "Phy_mmNH4": {'units': '-',
                         'munits': '-',
                         'longname': 'Phytoplankton NH4 external lim',
                         'cleanname': 'Phy_{mm}^{NH4}'},

           "Phy_mmNO3": {'units': '-',
                         'munits': '-',
                         'longname': 'Phytoplankton NO3 external lim',
                         'cleanname': 'Phy_{mm}^{NO3}'},

           "Phy_mmDIP": {'units': '-',
                         'munits': '-',
                         'longname': 'Phytoplankton DIP external lim',
                         'cleanname': 'Phy_{mm}^{DIP}'},

           "Phy_mmDSi": {'units': '-',
                         'munits': '-',
                         'longname': 'Phytoplankton DSi external lim',
                         'cleanname': 'Phy_{mm}^{DSi}'},

           "Phy_PC_max": {'units': 'd^{-1}',
                          'munits': 'd^{-1}',
                          'longname': 'Maximum photosynthetic rate',
                          'cleanname': 'P_C^{max}'},
           "Phy_PAR_t": {'units': 'µE m-2 d^{-1}',
                         'munits': 'µE m-2 d^{-1}',
                         'longname': 'Incident PAR',
                         'cleanname': 'PAR_{incident}'},
           "Phy_PAR_t_water_column": {'units': 'µE m-2 d^{-1}',
                                      'munits': 'µE m-2 d^{-1}',
                                      'longname': 'PAR in the water column',
                                      'cleanname': 'PAR_{water}'},

           "Phy_PC": {'units': 'd^{-1}',
                      'munits': 'd^{-1}',
                      'longname': 'Photosynthetic rate',
                      'cleanname': 'P_C'},

           "Phy_QN": {'units': 'molN:molC',
                      'munits': 'molN:molC',
                      'longname': 'Phytoplankton N:C ratio',
                      'cleanname': 'Q_N'},

           "Phy_QP": {'units': 'molP:molC',
                      'munits': 'molP:molC',
                      'longname': 'Phytoplankton P:C ratio',
                      'cleanname': 'Q_P'},

           "Phy_QSi": {'units': 'molSi:molC',
                       'munits': 'molSi:molC',
                       'longname': 'Phytoplankton Si:C ratio',
                       'cleanname': 'Q_{Si}'},

           "Phy_thetaC": {'units': 'gChl gC^{-1}',
                          'munits': 'mgChl mmolC^{-1}',
                          'longname': 'Phytoplankton Chl:C ratio',
                          'cleanname': 'θ_C',
                          'trsfrm': 1 / molmass_C},

           "Phy_source_PP.C": {'units': 'mmolC m^{-3} d^{-1}',
                               'munits': 'mmolC m^{-3} d^{-1}',
                               'longname': 'Primary production (C_{source})',
                               'cleanname': 'PP source'},

           "Phy_source_Chlprod.Chl": {'units': 'mgChl m^{-3} d^{-1}',
                                      'munits': 'mgChl m^{-3} d^{-1}',
                                      'longname': 'Chlorophyll production',
                                      'cleanname': 'Chl prod'},

           "Phy_rho_Chl": {'units': 'mgChl mmolN^{-1}',
                           'munits': 'mgChl mmolN^{-1}',
                           'longname': 'Chlorophyll synthesis rate (N_{spec})',
                           'cleanname': 'ρ_{Chl}'},

           "Macroflocs_alpha_FF": {'units': '-',
                                   'munits': '-',
                                   'longname': 'Floc-floc collision efficiency',
                                   'cleanname': 'α_{FF}'},

           "Macroflocs_alpha_PF": {'units': '-',
                                   'munits': '-',
                                   'longname': 'Primary-floc collision efficiency',
                                   'cleanname': 'α_{PF}'},

           "Macroflocs_alpha_PP": {'units': '-',
                                   'munits': '-',
                                   'longname': 'Primary-primary collision efficiency',
                                   'cleanname': 'α_{PP}'},

           "Macroflocs_Ncnum": {'units': '-',
                                'munits': '-',
                                'longname': 'Number of microflocs per macrofloc',
                                'cleanname': 'N_c'},

           # Resuspension diagnostics variables
           "Macroflocs_g_shear_rate_at_t": {'units': 's^{-1}',
                                            'munits': 's^{-1}',
                                            'longname': 'Shear rate at time t',
                                            'cleanname': 'g_{shear}'},

           "Macroflocs_tau_cr": {'units': 'Pa',
                                 'munits': 'Pa',
                                 'longname': 'Critical shear stress',
                                 'cleanname': '\\tau_{cr}'},



           "Macroflocs_sedimentation": {'units': ' m^{-3} d^{-1}',
                                        'munits': ' m^{-3} d^{-1}',
                                        'longname': 'Sedimentation rate',
                                        'cleanname': 'Sedimentation'},

           "Macroflocs_resuspension": {'units': ' m^{-3} d^{-1}',
                                       'munits': ' m^{-3} d^{-1}',
                                       'longname': 'Resuspension rate',
                                       'cleanname': 'Resuspension'},

           # Micro_in_Macro resuspension diagnostics
           "Micro_in_Macro_g_shear_rate_at_t": {'units': 's^{-1}',
                                                'munits': 's^{-1}',
                                                'longname': 'Shear rate at time t (micro in macro)',
                                                'cleanname': 'g_{shear}'},

           "Micro_in_Macro_tau_cr": {'units': 'Pa',
                                     'munits': 'Pa',
                                     'longname': 'Critical shear stress (micro in macro)',
                                     'cleanname': '\\tau_{cr}'},

           "Micro_in_Macro_sedimentation": {'units': ' m^{-3} d^{-1}',
                                            'munits': ' m^{-3} d^{-1}',
                                            'longname': 'Sedimentation rate (micro in macro)',
                                            'cleanname': 'Sedimentation'},

           "Micro_in_Macro_resuspension": {'units': ' m^{-3} d^{-1}',
                                           'munits': ' m^{-3} d^{-1}',
                                           'longname': 'Resuspension rate (micro in macro)',
                                           'cleanname': 'Resuspension'},

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
                        'oprt': 'setup.I /3600./24 * np.exp(-setup.k_att * SUMALLmPhy_Chl * setup.water_depth)',
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
