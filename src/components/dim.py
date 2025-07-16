import numpy as np

from ..core.base import BaseStateVar
from ..utils import functions as fns


class DIM(BaseStateVar):
    def __init__(self,
                 name,
                 r_nit=0.2,  # [d-1] Specific nitrification rate (Onur22)
                 A_E=0.65,  # [-] Activation Energy for temperature scaling (denitrification) (Onur22)
                 T_ref=283.,  # [K] Reference temperature (denitrification) (Onur22)
                 dt2=False,
                 ):

        super().__init__()

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

        self.name = name
        self.dt2 = dt2

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

    def update_val(self, concentration,
                   debugverbose=False):
        self.concentration = concentration

    def get_sources(self, t=None):
        # SOURCES
        self.get_source_remineralization()
        self.get_source_respiration()
        self.get_source_airseaexchange()
        self.get_source_redox(t)
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

        return np.array(self.sources)

    def get_sinks(self, t=None):
        # SINKS
        self.get_sink_uptake()
        self.get_sink_redox()

        # SINK terms of the state equation
        self.sinks = (self.sink_uptake +
                      self.sink_redox)

        return np.array(self.sinks)

    def get_source_remineralization(self):
        # Sch07
        if self.coupled_remin_sources is not None:
            if self.name == 'DIC' or self.name == 'DIN':
                self.source_remineralization = fns.get_all_contributors(self.coupled_remin_sources,
                                                                        'sink_remineralization')
            elif self.name == 'TA':
                # /!\ Then remineralization is actually a sink term for Total Alkalinity!
                self.source_remineralization = - (1 + 1 / 16.) * fns.get_all_contributors(self.coupled_remin_sources,
                                                                                          'sink_remineralization')
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

    def get_source_redox(self, t):
        # Onur22
        if self.name == 'NO3':
            self.source_redox = self.r_nit * fns.getlimT(self.setup.T.loc[t]['T'], A_E=self.A_E, T_ref=self.T_ref,
                                                         boltz=True) * self.coupled_NH4.concentration
        else:
            self.source_redox = 0.

    def get_source_sloppy_feeding(self):
        if self.formulation == "Onur22":
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

        else:
            self.source_sloppy_feeding = 0.

    def get_source_exudation(self):
        if self.formulation == 'Onur22':
            if self.name == 'DSi':
                self.source_exudation = self.coupled_exud_sources_phyto.sink_exudation.Si
                self.source_exudation = self.source_exudation + self.coupled_exud_sources_phyto.sink_lysis.Si
            else:
                self.source_exudation = 0.
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
        if self.name == 'DIC':  # Onur22 and Sch07
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


        elif self.name == 'DIN':  # Sch07
            self.sink_uptake = fns.get_all_contributors(self.coupled_uptake_sinks, 'source_uptake', 'N')
        elif self.name == 'TA':  # Sch07
            # /!\ Then uptake is actually a source term for Total Alkalinity!
            self.sink_uptake = - (1 + 1 / 16.) * fns.get_all_contributors(self.coupled_uptake_sinks, 'source_uptake',
                                                                          'N')

    def get_sink_redox(self):
        if self.name == 'NH4':
            self.sink_redox = self.coupled_NO3.source_redox
        else:
            self.sink_redox = 0.
