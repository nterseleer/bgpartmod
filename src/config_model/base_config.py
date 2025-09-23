import numpy as np

from src.components import phytoplankton as phyto
from src.components import dim
from src.components import heterotrophs as het
from src.components import detritus
from src.components import dom
from src.config_model import varinfos
from src.utils import functions as fns

# Onur22 = Kerimoglu et al 2022
# =============================
# Model units: [d-1], [mmol m-3] for all tracers except Chl [mg m-3]

# 20240319 - All targets and consumers are given for coupling, but the actual preference p_i,j determines whether
#            there is an effective target-consumer relationship or not. The code should be standard enough to adapt.
heterotrophs_list = ['BacF', 'BacA', 'Cil', 'HF']
potential_targets = ['DOCS', 'DOCL', 'TEPC', 'DetS', 'DetL'] + heterotrophs_list

# ===============================================================================
# PHYTOPLANKTON DIAGNOSTICS CONFIGURATIONS - SEPARATE FOR PERFORMANCE CONTROL
# ===============================================================================

# Full phytoplankton diagnostics
Phy_diagnostics_full = {
    'Phy': {
        'diagnostics': [
            'lim_N', 'lim_P', 'lim_Si', 'QN', 'QP', 'QSi', 'thetaC',
            'limNUT', 'limQUOTA.N', 'limQUOTA.P', 'limQUOTA.Si',
            'limQUOTAmin.N', 'limQUOTAmin.P', 'limQUOTAmin.Si',
            'mmNH4', 'mmNO3', 'mmDIP', 'mmDSi',
            'kd', 'PC_max', 'PC', 'source_PP.C', 'rho_Chl', 'limI', 'limT',
            'PAR_t', 'PAR_t_water_column',
            'source_uptake.NH4', 'source_uptake.NO3', 'source_uptake.N',
            'source_uptake.P', 'source_uptake.Si',
            'source_Chlprod.Chl',
            'sink_lysis.C', 'sink_lysis.N', 'sink_lysis.Chl', 'sink_lysis.P', 'sink_lysis.Si',
            'sink_mortality.C', 'sink_mortality.N', 'sink_mortality.Chl', 'sink_mortality.P',
            'sink_mortality.Si',
            'sink_exudation.C', 'sink_exudation.N', 'sink_exudation.Chl', 'sink_exudation.P',
            'sink_exudation.Si',
            'frac_exud_small',
            'sink_respiration.C', 'sink_respiration.N', 'sink_respiration.Chl', 'sink_respiration.P',
            'sink_respiration.Si',
            'sink_ingestion.C', 'sink_ingestion.N', 'sink_ingestion.Chl', 'sink_ingestion.P',
            'sink_ingestion.Si',
            'sink_aggregation.C', 'sink_aggregation.N', 'sink_aggregation.Chl', 'sink_aggregation.P',
            'sink_aggregation.Si',
        ]
    }
}


# ===============================================================================
# BASE CONFIGURATION (based on KERIMOGLU ET AL 2022)
# ===============================================================================

Onur = {
    'formulation': 'Onur22',
    'BacF': {'class': het.Heterotrophs,
             'parameters':
                 {'g_max': 4.,  # [-] ?!
                  'K_i': 10.,  # [-] ?!
                  'eff_C': 0.6,  # [-]
                  'eff_N': 1.,  # [-]
                  'eff_P': 1.,  # [-]
                  'mortrate_lin': 0.,  # [d-1]
                  'mortrate_quad': 0.,  # [d-1]
                  'lysrate_lin': 0.1,  # [d-1]
                  'lysrate_quad': 0.1,  # [d-1]  ! Difference: 0 in the paper, 0.1 in Onur's code
                  'f_unass_excr': 0.8,  # [-]
                  'f_unass_Si': 0.9,  # [-]
                  'zeta_resp': 0.05,  # [d-1]
                  'A_E': 0.65,  # [-]
                  'T_ref': 283.15,  # [K]
                  'eps_kd': 0.012,  # [m2 mmolC-1]
                  'pref_DOCS': 0.5,  # [-]
                  'pref_DOCL': 0.5,  # [-]
                  'pref_TEPC': 0.,  # [-]
                  'pref_DetS': 0.,  # [-]
                  'pref_DetL': 0.,  # [-]
                  'pref_BacA': 0.,  # [-]
                  'pref_BacF': 0.,  # [-]
                  'pref_HF': 0.,  # [-]
                  'pref_Cil': 0.,  # [-]
                  },
             'coupling': {
                 'coupled_targets': potential_targets,
                 'coupled_consumers': heterotrophs_list,
             },
             'initialization': {
                 'C': 1.,
                 'N': 1. / 4.9,
                 'P': 1. / 77.
             },
             'aggregate':
                 {'N_tot': 'N',
                  'C_tot': 'C',
                  'P_tot': 'P',
                  'POC': 'C',
                  'PON': 'N',
                  'POP': 'P'}
             },
    'BacA': {'class': het.Heterotrophs,
             'parameters':
                 {'g_max': 3.,  # [-] ?!
                  'K_i': 10.,  # [-] ?!
                  'eff_C': 0.6,  # [-]
                  'eff_N': 1.,  # [-]
                  'eff_P': 1.,  # [-]
                  'mortrate_lin': 0.,  # [d-1]
                  'mortrate_quad': 0.,  # [d-1]
                  'lysrate_lin': 0.1,  # [d-1]
                  'lysrate_quad': 0.1,  # [d-1]  ! Difference: 0 in the paper, 0.1 in Onur's code
                  'f_unass_excr': 0.8,  # [-]
                  'f_unass_Si': 0.9,  # [-]
                  'zeta_resp': 0.05,  # [d-1]
                  'A_E': 0.65,  # [-]
                  'T_ref': 283.15,  # [K]
                  'eps_kd': 0.012,  # [m2 mmolC-1]
                  'pref_DOCS': 0.,  # [-]
                  'pref_DOCL': 0.,  # [-]
                  'pref_TEPC': 0.2,  # [-]
                  'pref_DetS': 0.4,  # [-]
                  'pref_DetL': 0.4,  # [-]
                  'pref_BacA': 0.,  # [-]
                  'pref_BacF': 0.,  # [-]
                  'pref_HF': 0.,  # [-]
                  'pref_Cil': 0.,  # [-]
                  },
             'coupling': {
                 'coupled_targets': potential_targets,
                 'coupled_consumers': heterotrophs_list,
             },
             'initialization': {
                 'C': 1.,
                 'N': 1. / 4.9,
                 'P': 1. / 77.
             },
             'aggregate':
                 {'N_tot': 'N',
                  'C_tot': 'C',
                  'P_tot': 'P',
                  'POC': 'C',
                  'PON': 'N',
                  'POP': 'P'}
             },
    'HF': {'class': het.Heterotrophs,
           'parameters':
               {'g_max': 4.,  # [-] ?!
                'K_i': 20.,  # [-] ?!
                'eff_C': 0.4,  # [-]
                'eff_N': 1.,  # [-]
                'eff_P': 1.,  # [-]
                'mortrate_lin': 0.05,  # [d-1]
                'mortrate_quad': 0.0,  # [d-1] ! Difference: 0 in the paper, 0.1 in Onur's code
                'lysrate_lin': 0.,  # [d-1]
                'lysrate_quad': 0.06,  # [d-1]
                'f_unass_excr': 1.,  # [-]
                'f_unass_Si': 0.95,  # [-]
                'zeta_resp': 0.05,  # [d-1]
                'A_E': 0.65,  # [-]
                'T_ref': 283.15,  # [K]
                'eps_kd': 0.012,  # [m2 mmolC-1]
                'pref_DOCS': 0.,  # [-]
                'pref_DOCL': 0.,  # [-]
                'pref_TEPC': 0.1,  # [-]
                'pref_DetS': 0.1,  # [-]
                'pref_DetL': 0.,  # [-]
                'pref_BacA': 0.4,  # [-]
                'pref_BacF': 0.4,  # [-]
                'pref_HF': 0.,  # [-]
                'pref_Cil': 0.,  # [-]
                },
           'coupling': {
               'coupled_targets': potential_targets,
               'coupled_consumers': heterotrophs_list,
           },
           'initialization': {
               'C': .1,
               'N': .1 / 7.,
               'P': .1 / 180.
           },
           'aggregate':
               {'N_tot': 'N',
                'C_tot': 'C',
                'P_tot': 'P',
                'POC': 'C',
                'PON': 'N',
                'POP': 'P'}
           },
    'Cil': {'class': het.Heterotrophs,
            'parameters':
                {'g_max': 2.4,  # [-] ?!
                 'K_i': 20.,  # [-] ?!
                 'eff_C': 0.4,  # [-]
                 'eff_N': 1.,  # [-]
                 'eff_P': 1.,  # [-]
                 'mortrate_lin': 0.05,  # [d-1]
                 'mortrate_quad': 0.0,  # [d-1] ! Difference: 0 in the paper, 0.1 in Onur's code
                 'lysrate_lin': 0.,  # [d-1]
                 'lysrate_quad': 0.02,  # [d-1]
                 'f_unass_excr': 1.,  # [-]
                 'f_unass_Si': 0.95,  # [-]
                 'zeta_resp': 0.05,  # [d-1]
                 'A_E': 0.65,  # [-]
                 'T_ref': 283.15,  # [K]
                 'eps_kd': 0.012,  # [m2 mmolC-1]
                 'pref_DOCS': 0.,  # [-]
                 'pref_DOCL': 0.,  # [-]
                 'pref_TEPC': 0.1,  # [-]
                 'pref_DetS': 0.1,  # [-]
                 'pref_DetL': 0.,  # [-]
                 'pref_BacA': 0.1,  # [-]
                 'pref_BacF': 0.1,  # [-]
                 'pref_HF': 0.6,  # [-]
                 'pref_Cil': 0.,  # [-]
                 },
            'coupling': {
                'coupled_targets': potential_targets,
                'coupled_consumers': heterotrophs_list,
            },
            'initialization': {
                'C': .1,
                'N': .1 / 7,
                'P': .1 / 180.
            },
            'aggregate':
                {'N_tot': 'N',
                 'C_tot': 'C',
                 'P_tot': 'P',
                 'POC': 'C',
                 'PON': 'N',
                 'POP': 'P'}
            },

    'Phy': {'class': phyto.Phyto,
            'parameters':
                {'mu_max': 5.2,  # [d-1] !OK
                 'alpha': 7.e-6,  # converted to [mgC mgChl-1 µE-1 m2] from [mgC mgChl-1 E-1 m2]
                 # 'thetaN_max': 0.07 / 0.15 * varinfos.molmass_C,  # [mgChl mmolN-1] from theta_max/QNmax
                 'theta_max': 0.07 * varinfos.molmass_C,  # Converted to [gChl molC-1] from [gChl gC-1]
                 # Conversion is needed as theta_max/QN_max must give [gChl molN-1]
                 # to have appropriate units for rho_chl [gChl/mol]
                 'QN_max': 0.15,  # [molN:molC] !OK
                 'QP_max': 0.012,  # [molP:molC] !OK
                 'QSi_max': 0.18,  # [molSi:molC] !OK
                 'QN_min': 0.05,  # [molN:molC] !OK
                 'QP_min': 0.003,  # [molP:molC] !OK
                 'QSi_min': 0.06,  # [molSi:molC] !OK
                 'v_max_N': 0.78,  # [molN molC-1 d-1]
                 'v_max_P': 0.075,  # [molP molC-1 d-1]
                 'v_max_Si': 1.2,  # [molSi molC-1 d-1]
                 'KNO3': 2.,  # [mmolN m-3]
                 'KNH4': 0.5,  # [mmolN m-3]
                 'KDIP': 0.05,  # [mmolP m-3]
                 'KDSi': 0.43,  # [mmolSi m-3]
                 'mortrate': 0.05,  # [d-1] !OK
                 'zeta_resp_base': 0.01,  # [d-1] !OK
                 'zeta_resp_prod': 0.01,  # [-] !OK
                 'gamma_C_exud_base': 0.02,  # [d-1] !OK
                 'gamma_C_exud_prod': 0.11,  # [-] !OK
                 'lysrate': 0.1,  # [d-1] !OK
                 'A_E': 0.32,  # [-] !OK
                 'T_ref': 283.15,  # [K] !OK
                 'eps_kd': 0.024,  # [m2 mmolC-1]
                 'kdvar': True,
                 },
            'coupling':
                {
                    'coupled_aggreg_target': 'DetL',
                    'coupled_NH4': 'NH4',
                    'coupled_NO3': 'NO3',
                    'coupled_DIP': 'DIP',
                    'coupled_DSi': 'DSi',
                    'coupled_light_attenuators': ['DetL', 'DetS', 'BacA', 'BacF', 'HF', 'Cil'],
                },
            'initialization':
                {'C': 20.,
                 'N': 3.,
                 'Chl': 5.,
                 'P': 0.24,
                 'Si': 1.8, },
            'aggregate':
                {'C_tot': 'C',
                 'N_tot': 'N',
                 'P_tot': 'P',
                 'Si_tot': 'Si',
                 'POC': 'C',
                 'PON': 'N',
                 'POP': 'P',
                 'Chl_tot': 'Chl',
                 'Cphy_tot': 'C'},
            },

    'DOCS': {'class': dom.DOM,
             'parameters': {
                 'beta_TEPC': 0.032,  # [m3 mmolC-1 d-1]
                 'alpha_TEPC': 0.85,  # [-]
             },
             'coupling': {
                 'coupled_exud_sources_phyto': 'Phy',
                 'coupled_sloppy_feeding_sources': heterotrophs_list,
                 'coupled_lysis_sources': ['Phy'] + heterotrophs_list,
                 'coupled_TEPC': 'TEPC',
                 'coupled_consumers': heterotrophs_list,
             },
             'initialization': {
                 'C': 20.,
                 'N': 0.,
                 'P': 0.
             },
             'aggregate':
                 {'C_tot': 'C',
                  'N_tot': 'N',
                  'P_tot': 'P',
                  'DOC': 'C',
                  }
             },
    'DOCL': {'class': dom.DOM,
             'parameters': {
                 'alpha': 0.001,  # [-]
                 'alpha_TEPC': 1.,  # [-]
                 'beta': 0.86,  # [m3 mmolC-1 d-1]
                 'beta_TEPC': 0.064,  # [m3 mmolC-1 d-1]
             },
             'coupling': {
                 'coupled_exud_sources_phyto': 'Phy',
                 'coupled_breakdown_sources': ['TEPC'],
                 'coupled_lysis_sources': ['Phy'] + heterotrophs_list,
                 'coupled_TEPC': 'TEPC',
                 'coupled_consumers': heterotrophs_list,
             },
             'initialization': {
                 'C': 1.,
             },
             'aggregate':
                 {'C_tot': 'C',
                  'POC': 'C',
                  'DOC': 'C',
                  }
             },

    'DetL': {'class': detritus.Detritus,
             'parameters':
                 {'kleak': 0.25,
                  'eps_kd': 0.012,  # [m2 mmolC-1]
                  'beta_max': 0.033,  # [m3 mmolC-1 d-1] !OK
                  'KA2': 57.48,  # [mmolC m-3] !OK
                  },
             'coupling': {
                 'coupled_Phy': 'Phy',
                 'coupled_TEPC': 'TEPC',
                 'coupled_smaller_Det': 'DetS',
                 'coupled_consumers': heterotrophs_list,
             },
             'initialization': {
                 'C': 10.,
                 'N': 0.9,
                 'P': 0.065,
                 'Si': 1.
             },
             'aggregate':
                 {'N_tot': 'N',
                  'C_tot': 'C',
                  'P_tot': 'P',
                  'Si_tot': 'Si',
                  'POC': 'C',
                  'PON': 'N',
                  'POP': 'P',
                  }
             },

    'TEPC': {'class': dom.DOM,
             'parameters': {
                 'rho_TEP': 0.1,
             },
             'coupling': {
                 'coupled_aggreg_sources': ['DOCS', 'DOCL'],
                 'coupled_aggreg_target': 'DetL',
                 'coupled_consumers': heterotrophs_list
             },
             'initialization': {
                 'C': 0.62,
             },
             'aggregate':
                 {'C_tot': 'C',
                  'POC': 'C',
                  }
             },

    'DetS': {'class': detritus.Detritus,
             'parameters':
                 {'eps_kd': 0.012,  # [m2 mmolC-1]
                  },
             'coupling': {
                 'coupled_mortality_sources': ['Phy'] + heterotrophs_list,
                 'coupled_consumers': heterotrophs_list,
                 'coupled_sloppy_feeding_sources': heterotrophs_list,
                 'coupled_larger_Det': 'DetL',
             },
             'initialization': {
                 'C': 20.,
                 'N': 1.6,
                 'P': 0.1,
                 'Si': 1.8
             },
             'aggregate':
                 {'C_tot': 'C',
                  'N_tot': 'N',
                  'P_tot': 'P',
                  'Si_tot': 'Si',
                  'POC': 'C',
                  'PON': 'N',
                  'POP': 'P'
                  }
             },
    'DIC': {'class': dim.DIM,
            'coupling': {
                'coupled_resp_sources': ['Phy'] + heterotrophs_list,
                'coupled_uptake_sinks': 'Phy',
                'coupled_sloppy_feeding_sources': heterotrophs_list
            },
            'initialization': {'concentration': 1000},
            'aggregate':
                {'C_tot': 'concentration'}
            },
    'NH4': {'class': dim.DIM,
            'coupling': {
                'coupled_uptake_sinks': 'Phy',
                'coupled_NO3': 'NO3',
                'coupled_sloppy_feeding_sources': heterotrophs_list,
            },
            'initialization': {'concentration': 0.},
            'aggregate':
                {'N_tot': 'concentration'}
            },
    'NO3': {'class': dim.DIM,
            'coupling': {
                'coupled_uptake_sinks': 'Phy',
                'coupled_NH4': 'NH4',
                'coupled_sloppy_feeding_sources': heterotrophs_list
            },
            'initialization': {'concentration': 58.},
            'aggregate':
                {'N_tot': 'concentration'}
            },
    'DIP': {'class': dim.DIM,
            'coupling': {
                'coupled_uptake_sinks': 'Phy',
                'coupled_sloppy_feeding_sources': heterotrophs_list,
            },
            'initialization': {'concentration': 1.6},
            'aggregate':
                {'P_tot': 'concentration'}
            },
    'DSi': {'class': dim.DIM,
            'coupling': {
                'coupled_uptake_sinks': 'Phy',
                'coupled_exud_sources_phyto': 'Phy',
                'coupled_sloppy_feeding_sources': heterotrophs_list
            },
            'initialization': {'concentration': 22.},
            'aggregate':
                {'Si_tot': 'concentration'}
            },

}

# ===============================================================================
# CREATE OPTIMIZATION-OPTIMIZED CONFIGURATIONS WITH MINIMAL DIAGNOSTICS
# ===============================================================================

# Main configuration for classical simulation with full diagnostics output
Onur_full = fns.deep_update(Onur, Phy_diagnostics_full)
