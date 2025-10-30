import numpy as np

from ..core.base import BaseStateVar
from ..utils import functions as fns


class DIM(BaseStateVar):
    def __init__(self,
                 name,
                 r_nit=0.2,  # [d-1] Specific nitrification rate (Onur22)
                 A_E=0.65,  # [-] Activation Energy for temperature scaling (denitrification) (Onur22)
                 T_ref=283.15,  # [K] Reference temperature (denitrification) (Onur22)
                 k_remin=0.0, # [d-1] "external" remineralization rate (i.e., not related to simulated processes and state variables)
                 dt2=False,
                 dtype=np.float64,
                 bound_temp_to_1=True,  # Whether to bound temperature limitation to [0,1]
                 ):

        super().__init__(dtype=dtype)

        self.source_exudation = None
        self.source_sloppy_feeding = None
        self.formulation = None
        self.sinks = None
        self.sources = None
        self.coupled_sloppy_feeding_sources = None
        self.coupled_exud_sources_phyto = None
        self.coupled_NO3 = None
        self.coupled_NH4 = None
        self.coupled_uptake_sinks = None
        self.coupled_resp_sources = None
        self.coupled_remin_sources = None
        self.ICs = None
        self.concentration = None
        self.classname = 'DIM'  # Name used as prefix for variables (used in Model.finalizeres vs varinfos)

        self.r_nit = r_nit
        self.A_E = A_E
        self.T_ref = T_ref

        self.k_remin = k_remin
        self.remineralization_rate = 0

        self.name = name
        self.dt2 = dt2
        self.bound_temp_to_1 = bound_temp_to_1

        # Source and sink terms
        self.source_remineralization = None
        self.source_respiration = None
        self.source_airseaexchange = None
        self.source_redox = None
        self.source_riverine_loads = None
        self.sink_uptake = None
        self.sink_redox = None

    def set_ICs(self,
                concentration,  # Initial concentration
                ):
        self.concentration = concentration
        self.ICs = [self.concentration]

    def set_coupling(self,
                     coupled_remin_sources=None,
                     coupled_resp_sources=None,  # can be a list to sum over multiple respiration source terms
                     coupled_uptake_sinks=None,
                     coupled_NH4=None, coupled_NO3=None,
                     coupled_exud_sources_phyto=None,
                     coupled_sloppy_feeding_sources=None
                     ):
        # Coupling links
        self.coupled_remin_sources = coupled_remin_sources
        self.coupled_resp_sources = coupled_resp_sources
        self.coupled_uptake_sinks = coupled_uptake_sinks
        self.coupled_NH4 = coupled_NH4
        self.coupled_NO3 = coupled_NO3
        self.coupled_exud_sources_phyto = coupled_exud_sources_phyto
        self.coupled_sloppy_feeding_sources = coupled_sloppy_feeding_sources

        # Optimization: Pre-compute temperature limitation array for entire simulation
        if self.setup is not None:
            self._precompute_temp_limitation(
                A_E=self.A_E, T_ref=self.T_ref, boltz=True,
                bound_temp_to_1=self.bound_temp_to_1, suffix=''
            )

    def update_val(self, concentration,
                   t=None,
                   t_idx=None,
                   debugverbose=False):
        self.concentration = concentration
        # Calculate remineralization rate (available for all subsequent steps)
        if self.k_remin > 0 and t is not None:
            # Optimization: Use pre-computed temperature limitation
            self.remineralization_rate = self.k_remin * self.limT_array[t_idx]
        else:
            self.remineralization_rate = 0.

    def get_sources(self, t=None, t_idx=None):
        # SOURCES
        self.get_source_remineralization(t)
        self.get_source_respiration()
        self.get_source_airseaexchange()
        self.get_source_redox(t, t_idx=t_idx)
        self.get_source_sloppy_feeding()
        self.get_source_exudation()
        self.get_source_riverine_loads(t)

        # SOURCE terms of the state equation
        self.sources = (self.source_respiration +
                        self.source_remineralization +
                        self.source_airseaexchange +
                        self.source_redox +
                        self.source_sloppy_feeding +
                        self.source_exudation +
                        self.source_riverine_loads)

        return np.array(self.sources, dtype=self.dtype)

    def get_sinks(self, t=None, t_idx=None):
        # SINKS
        self.get_sink_uptake()
        self.get_sink_redox()

        # SINK terms of the state equation
        self.sinks = (self.sink_uptake +
                      self.sink_redox)

        return np.array(self.sinks, dtype=self.dtype)

    def get_source_remineralization(self, t=None):
        """Calculate remineralization sources for Onur22 formulation."""
        if self.coupled_remin_sources is not None:
            nutrient_map = {'NH4': 'N', 'DIP': 'P', 'DSi': 'Si'}
            self.source_remineralization = fns.get_all_contributors(
                self.coupled_remin_sources, 'sink_remineralization', nutrient_map[self.name])
        else:
            self.source_remineralization = 0.

    def get_source_respiration(self):
        if self.coupled_resp_sources is not None:
            if self.name == "DIC":
                self.source_respiration = fns.get_all_contributors(self.coupled_resp_sources,
                                                                   'sink_respiration', 'C')
        else:
            self.source_respiration = 0.

    def get_source_airseaexchange(self):
        # Sch07
        if self.name == 'DIC':
            self.source_airseaexchange = 0.
        else:
            self.source_airseaexchange = 0.

    def get_source_redox(self, t, t_idx=None):
        # Onur22
        if self.name == 'NO3':
            # Optimization: Use pre-computed temperature limitation
            self.source_redox = self.r_nit * self.limT_array[t_idx] * self.coupled_NH4.concentration
        else:
            self.source_redox = 0.

    def get_source_sloppy_feeding(self):
        if self.name == 'DIC':
            self.source_sloppy_feeding = np.sum(
                [sf.source_ing_C_unassimilated_to_dim for sf in self.coupled_sloppy_feeding_sources])
        elif self.name == 'NH4':
            self.source_sloppy_feeding = np.sum(
                [sf.source_ing_N_unassimilated_to_dim for sf in self.coupled_sloppy_feeding_sources])
        elif self.name == 'NO3':
            self.source_sloppy_feeding = 0.
        elif self.name == 'DIP':
            self.source_sloppy_feeding = np.sum(
                [sf.source_ing_P_unassimilated_to_dim for sf in self.coupled_sloppy_feeding_sources])
        elif self.name == 'DSi':
            self.source_sloppy_feeding = np.sum([sum(sf.source_ingestion.Si.values()) * (1 - sf.f_unass_Si)
                                                 for sf in self.coupled_sloppy_feeding_sources])


    def get_source_exudation(self):
        if self.name == 'DSi':
            self.source_exudation = self.coupled_exud_sources_phyto.sink_exudation.Si
            self.source_exudation = self.source_exudation + self.coupled_exud_sources_phyto.sink_lysis.Si
        else:
            self.source_exudation = 0.

    def get_source_riverine_loads(self, t=None):
        """Get riverine nutrient loads for current timestep."""
        if not self.setup.riverine_loads:
            self.source_riverine_loads = 0.
            return

        if self.name in self.setup.loads.columns:
            self.source_riverine_loads = self.setup.loads.loc[t][self.name]
        else:
            self.source_riverine_loads = 0.


    def get_sink_uptake(self):
        if self.name == 'DIC':  # Onur22
            self.sink_uptake = fns.get_all_contributors(self.coupled_uptake_sinks, 'source_PP', 'C')

        # 202403 TODO priority2 - make this more generic (e.g. dict)
        elif self.name == 'NH4':
            self.sink_uptake = fns.get_all_contributors(self.coupled_uptake_sinks, 'source_uptake', 'NH4')
        elif self.name == 'NO3':
            self.sink_uptake = fns.get_all_contributors(self.coupled_uptake_sinks, 'source_uptake', 'NO3')
        elif self.name == 'DIP':
            self.sink_uptake = fns.get_all_contributors(self.coupled_uptake_sinks, 'source_uptake', 'P')
        elif self.name == 'DSi':
            self.sink_uptake = fns.get_all_contributors(self.coupled_uptake_sinks, 'source_uptake', 'Si')

    def get_sink_redox(self):
        if self.name == 'NH4':
            self.sink_redox = self.coupled_NO3.source_redox
        else:
            self.sink_redox = 0.
