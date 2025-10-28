import numpy as np

from ..core.base import BaseOrg, Elms
from ..utils import functions as fns
from src.config_model import varinfos


class Detritus(BaseOrg):
    def __init__(self,
                 name,
                 eps_kd=2e-4 * varinfos.molmass_C,  # [m2 mgC-1] Diffuse attenuation cross section
                 kleak=0., #0.25,  # [d-1] Specific leakage rate (out of the system) (Onur22)
                 beta_max=0.033,  # [m3 mmolC-1 d-1] Max collision kernel for A2 (Onur22)
                 KA2=57.48,  # [mmolC m-3] Half saturation cst for TEP dependence of A2 (Onur22)
                 dt2=False,
                 bound_temp_to_1=True,  # Whether to bound temperature limitation to [0,1]
                 ):

        super().__init__()

        self.formulation = None
        self.aggTEP_C = None
        self.aggDetS_Si = None
        self.aggDetS_P = None
        self.aggDetS_N = None
        self.aggDetS_C = None
        self.aggPhy_Chl = None
        self.aggPhy_Si = None
        self.aggPhy_P = None
        self.aggPhy_N = None
        self.sumc = None
        self.aggPhy_C = None
        self.aggPD_tot = None
        self.classname = 'Det'  # Name used as prefix for variables (used in Model.finalizeres vs varinfos)

        self.name = name
        self.eps_kd = eps_kd
        self.kleak = kleak
        self.beta_max = beta_max
        self.KA2 = KA2
        self.dt2 = dt2
        self.bound_temp_to_1 = bound_temp_to_1

        # Source and sink terms
        self.source_aggregation = Elms()
        self.source_mortality = Elms()
        self.source_sloppy_feeding = Elms()
        self.sink_breakdown = Elms()
        self.sink_ingestion = Elms()
        self.sink_aggregation = Elms()
        self.sink_leakage_out = Elms()
        self.sink_remineralization = Elms()
        self.sink_vertical_loss = Elms()  # Vertical loss coupled to mineral flocs


        # Coupling links
        self.coupled_aggreg_sources = None
        self.coupled_mortality_sources = None
        self.coupled_sloppy_feeding_sources = None
        self.coupled_Phy = None
        self.coupled_TEPC = None
        self.coupled_larger_Det = None
        self.coupled_smaller_Det = None
        self.coupled_consumers = None
        self.coupled_aggregate = None
        self.coupled_remin_products = None

    def set_coupling(self,
                     coupled_aggreg_sources=None,
                     coupled_mortality_sources=None,
                     coupled_sloppy_feeding_sources=None,
                     coupled_Phy=None, coupled_TEPC=None, coupled_larger_Det=None, coupled_smaller_Det=None,
                     coupled_consumers=None,
                     coupled_aggregate=None,
                     coupled_remin_products=None
                     ):
        # Coupling links
        self.coupled_aggreg_sources = coupled_aggreg_sources
        self.coupled_mortality_sources = coupled_mortality_sources
        self.coupled_sloppy_feeding_sources = coupled_sloppy_feeding_sources
        self.coupled_Phy = coupled_Phy
        self.coupled_TEPC = coupled_TEPC
        self.coupled_larger_Det = coupled_larger_Det
        self.coupled_smaller_Det = coupled_smaller_Det
        self.coupled_consumers = coupled_consumers
        self.coupled_aggregate = coupled_aggregate
        self.coupled_remin_products = coupled_remin_products

    def get_coupled_processes_indepent_sinks_sources(self, t=None):
        """
        Get those sinks and sources that do not depend on sinks/sources from couples state variables.
        --> they can be computed first and then be available elsewhere for coupling.
        This is a first step to solve the issue of the order in which state variables are given in Config_model.
        :return:
        """
        # SOURCES

        # SINKS
        self.get_sink_remineralization()

    def get_sources(self, t=None):
        # SOURCES
        self.get_source_aggregation()
        self.get_source_mortality()
        self.get_source_sloppy_feeding()

        # SOURCE terms of the state equation
        self.C_sources = (self.source_aggregation.C +
                          self.source_mortality.C +
                          self.source_sloppy_feeding.C)

        self.N_sources = (self.source_aggregation.N +
                          self.source_mortality.N +
                          self.source_sloppy_feeding.N)
        if self.P is not None:
            self.P_sources = (self.source_aggregation.P +
                              self.source_mortality.P +
                              self.source_sloppy_feeding.P)
        if self.Si is not None:
            self.Si_sources = (self.source_aggregation.Si +
                               self.source_mortality.Si +
                               self.source_sloppy_feeding.Si)
        return np.array(
            [sources for sources in (self.C_sources, self.N_sources, self.P_sources, self.Si_sources) if
             sources is not None])

    def get_sinks(self, t=None):
        # SINKS
        self.get_sink_breakdown()
        self.get_sink_aggregation()
        self.get_sink_leakage_out()
        self.get_sink_ingestion()
        self.get_sink_vertical_loss()
        # self.get_sink_remineralization()

        # SINK terms of the state equation
        self.C_sinks = (self.sink_breakdown.C +
                        self.sink_aggregation.C +
                        self.sink_leakage_out.C +
                        self.sink_ingestion.C +
                        self.sink_remineralization.C +
                        self.sink_vertical_loss.C)

        self.N_sinks = (self.sink_breakdown.N +
                        self.sink_aggregation.N +
                        self.sink_leakage_out.N +
                        self.sink_ingestion.N +
                        self.sink_remineralization.N +
                        self.sink_vertical_loss.N)

        if self.P is not None:
            self.P_sinks = (self.sink_breakdown.P +
                            self.sink_aggregation.P +
                            self.sink_leakage_out.P +
                            self.sink_ingestion.P +
                            self.sink_remineralization.P +
                            self.sink_vertical_loss.P)

        if self.Si is not None:
            self.Si_sinks = (self.sink_breakdown.Si +
                             self.sink_aggregation.Si +
                             self.sink_leakage_out.Si +
                             self.sink_ingestion.Si +
                             self.sink_remineralization.Si +
                             self.sink_vertical_loss.Si)

        return np.array(
            [sinks for sinks in (self.C_sinks, self.N_sinks, self.P_sinks, self.Si_sinks) if sinks is not None])

    def get_source_aggregation(self):
        """Calculate aggregation sources for Onur22 formulation."""
        if self.name == 'DetL':
            betaA2 = self.beta_max * self.coupled_TEPC.C / (self.coupled_TEPC.C + self.KA2)
            self.aggPD_tot = betaA2 * (
                    self.coupled_Phy.C + self.coupled_TEPC.C + self.coupled_smaller_Det.C) * self.C

            self.sumc = self.coupled_Phy.C + self.coupled_TEPC.C + self.coupled_smaller_Det.C + self.C

            # PD is Phyto+Detritus
            self.aggPhy_C = self.aggPD_tot * self.coupled_Phy.C / self.sumc
            self.aggPhy_N = self.aggPhy_C * self.coupled_Phy.QN
            self.aggPhy_P = self.aggPhy_C * self.coupled_Phy.QP
            self.aggPhy_Si = self.aggPhy_C * self.coupled_Phy.QSi
            self.aggPhy_Chl = self.aggPhy_C * self.coupled_Phy.thetaC

            self.aggDetS_C = self.aggPD_tot * self.coupled_smaller_Det.C / self.sumc
            self.aggDetS_N = self.aggDetS_C * self.coupled_smaller_Det.QN
            self.aggDetS_P = self.aggDetS_C * self.coupled_smaller_Det.QP
            self.aggDetS_Si = self.aggDetS_C * self.coupled_smaller_Det.QSi

            self.aggTEP_C = self.aggPD_tot * self.coupled_TEPC.C / self.sumc

            self.source_aggregation.C = (self.aggPhy_C +
                                         self.aggDetS_C +
                                         self.aggTEP_C)
            self.source_aggregation.N = (self.aggPhy_N +
                                         self.aggDetS_N)
            self.source_aggregation.P = (self.aggPhy_P +
                                         self.aggDetS_P)
            self.source_aggregation.Si = (self.aggPhy_Si +
                                          self.aggDetS_Si)

        elif self.name == 'DetS':
            self.source_aggregation.C = 0.
            self.source_aggregation.N = 0.
            self.source_aggregation.P = 0.
            self.source_aggregation.Si = 0.

    def get_source_mortality(self):
        if self.name == "DetS":
            self.source_mortality.C = fns.get_all_contributors(self.coupled_mortality_sources, 'sink_mortality', 'C')
            coupled_mortsources_N = [c for c in self.coupled_mortality_sources if c.sink_mortality.N is not None]
            self.source_mortality.N = fns.get_all_contributors(coupled_mortsources_N, 'sink_mortality', 'N')
            coupled_mortsources_P = [c for c in self.coupled_mortality_sources if c.sink_mortality.P is not None]
            self.source_mortality.P = fns.get_all_contributors(coupled_mortsources_P, 'sink_mortality', 'P')
            coupled_mortsources_Si = [c for c in self.coupled_mortality_sources if c.sink_mortality.Si is not None]
            self.source_mortality.Si = fns.get_all_contributors(coupled_mortsources_Si, 'sink_mortality', 'Si')
        else:
            self.source_mortality.C = 0.
            self.source_mortality.N = 0.
            self.source_mortality.P = 0.
            self.source_mortality.Si = 0.

    def get_source_sloppy_feeding(self):
        self.source_sloppy_feeding.C = 0.
        self.source_sloppy_feeding.N = 0.
        self.source_sloppy_feeding.P = 0.
        if self.coupled_sloppy_feeding_sources is not None:
            self.source_sloppy_feeding.Si = np.sum([sum(sf.source_ingestion.Si.values()) * sf.f_unass_Si
                                                    for sf in self.coupled_sloppy_feeding_sources])
        else:
            self.source_sloppy_feeding.Si = 0.

    def get_sink_breakdown(self):
        """Calculate breakdown sinks for Onur22 formulation."""
        self.sink_breakdown.C = 0.
        self.sink_breakdown.N = 0.
        self.sink_breakdown.P = 0.
        self.sink_breakdown.Si = 0.

    def get_sink_aggregation(self):
        if self.name == 'DetL':
            self.sink_aggregation.C = 0.
            self.sink_aggregation.N = 0.
            self.sink_aggregation.P = 0.
            self.sink_aggregation.Si = 0.
        elif self.name == 'DetS':
            self.sink_aggregation.C = self.coupled_larger_Det.aggDetS_C
            self.sink_aggregation.N = self.coupled_larger_Det.aggDetS_N
            self.sink_aggregation.P = self.coupled_larger_Det.aggDetS_P
            self.sink_aggregation.Si = self.coupled_larger_Det.aggDetS_Si

    def get_sink_leakage_out(self):
        # Independant leakage (initially: for Mesocosm loss in Kerimoglu et al., 2022)
        self.sink_leakage_out.C = self.kleak * self.C
        self.sink_leakage_out.N = self.kleak * self.N if self.N else 0.
        self.sink_leakage_out.P = self.kleak * self.P if self.P else 0.
        self.sink_leakage_out.Si = self.kleak * self.Si if self.Si else 0.

    def get_sink_ingestion(self):
        """Calculate ingestion sinks for Onur22 formulation."""
        self.sink_ingestion.C = np.sum([c.source_ingestion.C[self.name] for c in self.coupled_consumers])
        self.sink_ingestion.N = np.sum([c.source_ingestion.N[self.name] for c in self.coupled_consumers])
        self.sink_ingestion.P = np.sum([c.source_ingestion.P[self.name] for c in self.coupled_consumers])
        self.sink_ingestion.Si = np.sum([c.source_ingestion.Si[self.name] for c in self.coupled_consumers])

    def get_sink_remineralization(self, t=None):
        self.sink_remineralization.C = 0.
        self.sink_remineralization.N = 0.
        self.sink_remineralization.P = 0.
        self.sink_remineralization.Si = 0.

        if self.coupled_remin_products is not None:
            for product in self.coupled_remin_products:
                if hasattr(product, 'remineralization_rate') and product.remineralization_rate > 0:
                    if product.name == 'NH4' and hasattr(self, 'N'):
                        self.sink_remineralization.N = product.remineralization_rate * self.N
                    elif product.name == 'DIP' and hasattr(self, 'P'):
                        self.sink_remineralization.P = product.remineralization_rate * self.P
                    elif product.name == 'DSi' and hasattr(self, 'Si'):
                        self.sink_remineralization.Si = product.remineralization_rate * self.Si


    def get_sink_vertical_loss(self):
        """Vertical loss coupled to mineral floc dynamics (sedimentation - resuspension)"""
        if self.coupled_aggregate is not None:
            # Convert s-1 to d-1 when coupling to biogeochemistry
            rate = (self.coupled_aggregate.net_vertical_loss_rate *
                    self.coupled_aggregate.time_conversion_factor)
            self.sink_vertical_loss.C = rate * self.C
            self.sink_vertical_loss.N = rate * self.N
            if self.P is not None:
                self.sink_vertical_loss.P = rate * self.P
            else:
                self.sink_vertical_loss.P = 0.
            if self.Si is not None:
                self.sink_vertical_loss.Si = rate * self.Si
            else:
                self.sink_vertical_loss.Si = 0.
        else:
            self.sink_vertical_loss.C = 0.
            self.sink_vertical_loss.N = 0.
            self.sink_vertical_loss.P = 0.
            self.sink_vertical_loss.Si = 0.
