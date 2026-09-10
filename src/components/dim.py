import numpy as np

from ..core.base import BaseStateVar
from ..utils import functions as fns


# Currency carried by each dissolved-inorganic pool: which field of the emitting
# component's Elms it reads. Two tables, because the processes do not split N the same way:
# remineralisation produces bulk N (routed to NH4), while uptake is already split between
# NH4 and NO3 by the phytoplankton. A pool absent from a table has no term for that process.
REMIN_CURRENCY = {'NH4': 'N', 'DIP': 'P', 'DSi': 'Si'}
UPTAKE_CURRENCY = {'NH4': 'NH4', 'NO3': 'NO3', 'DIP': 'P', 'DSi': 'Si'}   # 'DIC': 'C'


class DIM(BaseStateVar):
    def __init__(self,
                 name,
                 r_nit=0.2,  # [d-1] Specific nitrification rate (Kerimoglu22)
                 A_E=0.65,  # [-] Activation energy for temperature scaling (Kerimoglu22)
                 T_ref=283.15,  # [K] Reference temperature (Kerimoglu22)
                 k_remin=0.0, # [d-1] "external" remineralization rate (i.e., not related to simulated processes and state variables)
                 dt2=False,
                 dtype=np.float64,
                 bound_temp_to_1=True,  # Whether to bound temperature limitation to [0,1]
                 ):

        super().__init__(dtype=dtype)

        self.source_exudation = None
        self.source_sloppy_feeding = None
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
        self.source_airseaexchange = 0.   # non resolu, cf. get_source_airseaexchange
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
        # self.get_source_airseaexchange()
        self.get_source_redox(t, t_idx=t_idx)
        self.get_source_sloppy_feeding()
        self.get_source_exudation()
        self.get_source_riverine_loads(t, t_idx=t_idx)

        # SOURCE terms of the state equation
        self.sources = (self.source_respiration +
                        self.source_remineralization +
                        # self.source_airseaexchange +
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
        """Calculate remineralization sources for Kerimoglu22 formulation."""
        if self.coupled_remin_sources is not None:
            self.source_remineralization = fns.get_all_contributors(
                self.coupled_remin_sources, 'sink_remineralization', REMIN_CURRENCY[self.name])
        else:
            self.source_remineralization = 0.

    def get_source_respiration(self):
        # Carbon only: respiration releases C, the other currencies go through
        # remineralisation. Written unconditionally -- the previous form left the
        # attribute untouched (hence stale, or None on the first call) for a non-DIC pool
        # that had coupled_resp_sources set.
        if self.coupled_resp_sources is not None and self.name == "DIC":
            self.source_respiration = fns.get_all_contributors(self.coupled_resp_sources,
                                                               'sink_respiration', 'C')
        else:
            self.source_respiration = 0.

    # Requires a full description of the DIC source and sink dynamics.
    # def get_source_airseaexchange(self):
    #     self.source_airseaexchange = ...

    def get_source_redox(self, t, t_idx=None):
        # Kerimoglu22
        if self.name == 'NO3':
            # Optimization: Use pre-computed temperature limitation
            self.source_redox = self.r_nit * self.limT_array[t_idx] * self.coupled_NH4.concentration
        else:
            self.source_redox = 0.

    def get_source_sloppy_feeding(self):
        if self.name == 'DIC':
            self.source_sloppy_feeding = 0.
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

    def get_source_riverine_loads(self, t=None, t_idx=None):
        """Get riverine nutrient loads for current timestep.

        Optimization: Uses pre-computed arrays indexed by t_idx instead of DataFrame .loc lookups.
        """
        if not self.setup.riverine_loads:
            self.source_riverine_loads = 0.
            return

        loads_array = getattr(self.setup, f'loads_{self.name}_array', None)
        self.source_riverine_loads = loads_array[t_idx] if loads_array is not None else 0.


    def get_sink_uptake(self):
        """Uptake of this nutrient by the coupled autotrophs (Kerimoglu22)."""
        currency = UPTAKE_CURRENCY.get(self.name)
        self.sink_uptake = 0. if currency is None else fns.get_all_contributors(
            self.coupled_uptake_sinks, 'source_uptake', currency)
        # DIC (hence absent from UPTAKE_CURRENCY). Requires a full description of the
        # DIC source and sink dynamics.
        #     self.sink_uptake = fns.get_all_contributors(
        #         self.coupled_uptake_sinks, 'source_PP', 'C')

    def get_sink_redox(self):
        if self.name == 'NH4':
            self.sink_redox = self.coupled_NO3.source_redox
        else:
            self.sink_redox = 0.
