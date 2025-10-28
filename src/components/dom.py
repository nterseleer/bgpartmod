import numpy as np

from ..core.base import BaseOrg, Elms
from ..utils import functions as fns
from src.config_model import varinfos

class DOM(BaseOrg):
    def __init__(self,
                 name,
                 rho_TEP=0.01,
                 alpha=0.001,  # [-] Self attachment probability (Onur22)
                 alpha_TEPC=0.4,  # [-] Particle stickiness DOC-TEPC (Onur22)
                 beta=0.86,  # [m3 mmolC-1 d-1] Self collision kernel (Onur22)
                 beta_TEPC=0.064,  # [m3 mmolC-1 d-1] Collision kernel DOC-TEPC (Onur22)
                 kleak=0.,           # [d-1] Specific leakage rate (out of the system)
                 dt2=False,
                 bound_temp_to_1=True,  # Whether to bound temperature limitation to [0,1]
                 ):

        super().__init__()

        self.formulation = None
        self.classname = 'DOM'  # Name used as prefix for variables (used in Model.finalizeres vs varinfos)
        self.name = name
        self.rho_TEP = rho_TEP
        self.alpha = alpha
        self.alpha_TEPC = alpha_TEPC
        self.beta = beta
        self.beta_TEPC = beta_TEPC
        self.kleak = kleak
        self.dt2 = dt2
        self.bound_temp_to_1 = bound_temp_to_1

        # Source and sink terms
        self.source_exudation = Elms()
        self.source_breakdown = Elms()
        self.source_aggregation = Elms()
        self.source_sloppy_feeding = Elms()
        self.sink_ingestion = Elms()
        self.sink_remineralization = Elms()
        self.sink_breakdown = Elms()
        self.sink_aggregation = Elms()
        self.sink_vertical_loss = Elms()  # Vertical loss coupled to mineral flocs
        self.sink_leakage_out = Elms()

        # Coupling links
        self.coupled_aggregate = None  # Link to Macroflocs for vertical dynamics
        self.coupled_exud_sources_phyto = None
        self.coupled_exud_sources_zoo = None
        self.coupled_breakdown_sources = None
        self.coupled_lysis_sources = None
        self.coupled_TEPC = None
        self.coupled_aggreg_sources = None
        self.coupled_aggreg_target = None
        self.coupled_sloppy_feeding_sources = None
        self.coupled_consumers = None
        self.coupled_remin_products = None

    def set_coupling(self,
                     coupled_exud_sources_phyto=None, coupled_exud_sources_zoo=None,
                     coupled_breakdown_sources=None, coupled_lysis_sources=None,
                     coupled_aggreg_sources=None,
                     coupled_aggreg_target=None,
                     coupled_TEPC=None,
                     coupled_sloppy_feeding_sources=None,
                     coupled_consumers=None,
                     coupled_remin_products=None,
                     coupled_aggregate=None
                     ):
        self.coupled_exud_sources_phyto = coupled_exud_sources_phyto
        self.coupled_exud_sources_zoo = coupled_exud_sources_zoo
        self.coupled_breakdown_sources = coupled_breakdown_sources
        self.coupled_lysis_sources = coupled_lysis_sources
        self.coupled_TEPC = coupled_TEPC
        self.coupled_aggreg_sources = coupled_aggreg_sources
        self.coupled_aggreg_target = coupled_aggreg_target
        self.coupled_sloppy_feeding_sources = coupled_sloppy_feeding_sources
        try:
            iterator = iter(coupled_consumers)
        except TypeError:
            self.coupled_consumers = [coupled_consumers]
        else:
            self.coupled_consumers = coupled_consumers
        self.coupled_remin_products = coupled_remin_products
        self.coupled_aggregate = coupled_aggregate


    def get_coupled_processes_indepent_sinks_sources(self, t=None):
        """
        Get those sinks and sources that do not depend on sinks/sources from couples state variables.
        --> they can be computed first and then be available elsewhere for coupling.
        This is a first step to solve the issue of the order in which state variables are given in Config_model.
        :return:
        """
        # SOURCES

        # SINKS
        self.get_sink_breakdown()
        self.get_sink_remineralization()
        self.get_sink_leakage_out()


    def get_sources(self, t=None):
        # SOURCES
        self.get_source_exudation()
        self.get_source_breakdown()
        self.get_source_aggregation()
        self.get_source_sloppy_feeding()

        # SINKS
        self.get_sink_aggregation()

        # SOURCE terms of the state equation
        self.C_sources = (self.source_exudation.C +
                          self.source_breakdown.C +
                          self.source_aggregation.C +
                          self.source_sloppy_feeding.C)
        if self.N is not None:
            self.N_sources = (self.source_exudation.N +
                              self.source_breakdown.N +
                              self.source_aggregation.N +
                              self.source_sloppy_feeding.N)
        if self.P is not None:
            self.P_sources = (self.source_exudation.P +
                              self.source_breakdown.P +
                              self.source_aggregation.P +
                              self.source_sloppy_feeding.P)

        return np.array(
            [sources for sources in (self.C_sources, self.N_sources, self.P_sources) if sources is not None])

    def get_sinks(self, t=None):
        # SINKS
        self.get_sink_ingestion()
        self.get_sink_vertical_loss()
        # self.get_sink_remineralization()
        # self.get_sink_breakdown()
        # self.get_sink_aggregation()
        # self.get_sink_leakage_out()


        # SINK terms of the state equation
        self.C_sinks = (self.sink_ingestion.C +
                        self.sink_remineralization.C +
                        self.sink_breakdown.C +
                        self.sink_aggregation.C +
                        self.sink_vertical_loss.C +
                        self.sink_leakage_out.C)
        if self.N is not None:
            self.N_sinks = (self.sink_ingestion.N +
                            self.sink_remineralization.N +
                            self.sink_breakdown.N +
                            self.sink_aggregation.N +
                            self.sink_vertical_loss.N +
                        self.sink_leakage_out.N)
        if self.P is not None:
            self.P_sinks = (self.sink_ingestion.P +
                            self.sink_remineralization.P +
                            self.sink_breakdown.P +
                            self.sink_aggregation.P +
                            self.sink_vertical_loss.P +
                        self.sink_leakage_out.N)

        return np.array(
            [sinks for sinks in (self.C_sinks, self.N_sinks, self.P_sinks) if sinks is not None])

    def get_source_exudation(self):
        """Calculate exudation sources for Onur22 formulation."""
        if self.name == "DOCS":
            self.source_exudation.C = (self.coupled_exud_sources_phyto.sink_exudation.C *
                                       self.coupled_exud_sources_phyto.frac_exud_small)
            self.source_exudation.N = self.coupled_exud_sources_phyto.sink_exudation.N
            self.source_exudation.P = self.coupled_exud_sources_phyto.sink_exudation.P
        elif self.name == "DOCL":
            self.source_exudation.C = (self.coupled_exud_sources_phyto.sink_exudation.C *
                                       (1. - self.coupled_exud_sources_phyto.frac_exud_small))
        else:
            self.source_exudation.C = 0.

    def get_source_breakdown(self):
        """Calculate breakdown sources for Onur22 formulation."""
        if self.name == "DOCS":
            # No C source, but lysed material from phytoplankton and Heterotrophs go to DON and DOP
            self.source_breakdown.C = 0.  # All C lysis/breakdown goes to DOCL
            self.source_breakdown.N = fns.get_all_contributors(self.coupled_lysis_sources, 'sink_lysis', 'N')
            self.source_breakdown.P = fns.get_all_contributors(self.coupled_lysis_sources, 'sink_lysis', 'P')
        elif self.name == "DOCL":
            self.source_breakdown.C = (fns.get_all_contributors(self.coupled_lysis_sources, 'sink_lysis', 'C') +
                                       fns.get_all_contributors(self.coupled_breakdown_sources, 'sink_breakdown', 'C'))
        else:
            self.source_breakdown.C = 0.
            self.source_breakdown.N = 0.
            self.source_breakdown.P = 0.

    def get_source_aggregation(self):
        """Calculate aggregation sources for Onur22 formulation."""
        if self.name == "TEPC":
            self.source_aggregation.C = fns.get_all_contributors(self.coupled_aggreg_sources, 'sink_aggregation', 'C')
        else:
            self.source_aggregation.C = 0.
            self.source_aggregation.N = 0.
            self.source_aggregation.P = 0.

    def get_source_sloppy_feeding(self):
        """Calculate sloppy feeding sources (vectorized)."""
        if self.name == 'DOCS':
            # Vectorized: extract all unassimilated contributions at once
            C_unassim = np.array([sf.source_ing_C_unassimilated_to_dom for sf in self.coupled_sloppy_feeding_sources])
            N_unassim = np.array([sf.source_ing_N_unassimilated_to_dom for sf in self.coupled_sloppy_feeding_sources])
            P_unassim = np.array([sf.source_ing_P_unassimilated_to_dom for sf in self.coupled_sloppy_feeding_sources])
            self.source_sloppy_feeding.C = np.sum(C_unassim)
            self.source_sloppy_feeding.N = np.sum(N_unassim)
            self.source_sloppy_feeding.P = np.sum(P_unassim)
        else:
            self.source_sloppy_feeding.C = 0

    def get_sink_ingestion(self):
        """Calculate ingestion sinks for Onur22 formulation (vectorized)."""
        # Vectorized: extract all consumer contributions
        C_ing = np.array([c.source_ingestion.C[self.name] for c in self.coupled_consumers])
        self.sink_ingestion.C = np.sum(C_ing)

        if self.name == "DOCS":
            N_ing = np.array([c.source_ingestion.N[self.name] for c in self.coupled_consumers])
            P_ing = np.array([c.source_ingestion.P[self.name] for c in self.coupled_consumers])
            self.sink_ingestion.N = np.sum(N_ing)
            self.sink_ingestion.P = np.sum(P_ing)
        else:
            self.sink_ingestion.N = 0.
            self.sink_ingestion.P = 0.

    def get_sink_remineralization(self, t=None):
        """Calculate remineralization sinks for Onur22 formulation."""
        self.sink_remineralization.C = 0.
        self.sink_remineralization.N = 0.
        self.sink_remineralization.P = 0.

        if self.coupled_remin_products is not None:
            for product in self.coupled_remin_products:
                if product.name == 'NH4' and hasattr(self, 'N'):
                    self.sink_remineralization.N = product.remineralization_rate * self.N
                elif product.name == 'DIP' and hasattr(self, 'P'):
                    self.sink_remineralization.P = product.remineralization_rate * self.P

    def get_sink_breakdown(self, t=None):
        """Calculate breakdown sinks for Onur22 formulation."""
        if self.name == 'TEPC':
            self.sink_breakdown.C = self.rho_TEP * self.C
        else:
            self.sink_breakdown.C = 0.
        self.sink_breakdown.N = 0.
        self.sink_breakdown.P = 0.

    def get_sink_aggregation(self):
        """Calculate aggregation sinks for Onur22 formulation."""
        if self.name == 'DOCL':
            self.sink_aggregation.C = (self.alpha * self.beta * self.C ** 2 +
                                       self.alpha_TEPC * self.beta_TEPC * self.C * self.coupled_TEPC.C)
        elif self.name == "DOCS":
            self.sink_aggregation.C = self.alpha_TEPC * self.beta_TEPC * self.C * self.coupled_TEPC.C
            # DON and DOP are considered not to aggregate
            self.sink_aggregation.N = 0.
            self.sink_aggregation.P = 0.
        elif self.name == "TEPC":
            self.sink_aggregation.C = self.coupled_aggreg_target.aggTEP_C

    def get_sink_vertical_loss(self):
        """Vertical loss coupled to mineral floc dynamics (sedimentation - resuspension)"""
        if self.coupled_aggregate is not None:
            # Convert s-1 to d-1 when coupling to biogeochemistry
            rate = (self.coupled_aggregate.net_vertical_loss_rate *
                    self.coupled_aggregate.time_conversion_factor)
            self.sink_vertical_loss.C = rate * self.C
            # DOM components (TEPC, DOCS, DOCL) are C-only
            if self.N is not None:
                self.sink_vertical_loss.N = rate * self.N
            else:
                self.sink_vertical_loss.N = 0.
            if self.P is not None:
                self.sink_vertical_loss.P = rate * self.P
            else:
                self.sink_vertical_loss.P = 0.
        else:
            self.sink_vertical_loss.C = 0.
            self.sink_vertical_loss.N = 0.
            self.sink_vertical_loss.P = 0.

    def get_sink_leakage_out(self):
        # Independent leakage = pure loss term (open system!) ~ loss of TEP to TEP_non_reactive
        self.sink_leakage_out.C = self.kleak * self.C
        self.sink_leakage_out.N = self.kleak * self.N if self.N else 0.
        self.sink_leakage_out.P = self.kleak * self.P if self.P else 0.
        self.sink_leakage_out.Si = self.kleak * self.Si if self.Si else 0.
