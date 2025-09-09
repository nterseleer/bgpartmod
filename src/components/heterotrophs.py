import numpy as np

from ..core.base import BaseOrg, Elms
from ..utils import functions as fns


class Heterotrophs(BaseOrg):
    def __init__(self,
                 name,
                 g_max=0.230,  # [-]?? or [d-1] Maximum ingestion rate
                 K_i=10.,  # [-]?? or [mmolC m-3] Half-saturation cst for ingestion
                 eff_C=0.6,  # [-] C assimilation efficiency
                 eff_N=1.,  # [-] N assimilation efficiency
                 eff_P=1.,  # [-] N assimilation efficiency
                 zeta_resp=0.05,  # [d-1] Basal respiration rate
                 lysrate_lin=0.1,  # [d-1] Linear lysis rate
                 lysrate_quad=0.,  # [d-1] Quadratic lysis rate
                 mortrate_lin=0.,  # [d-1] Linear mortality rate
                 mortrate_quad=0.,  # [d-1] Quadratic mortality rate
                 f_unass_excr=0.8,  # [-] Organic fraction of unassimilated excretes
                 f_unass_Si=0.9,  # [-] Detrital fraction of unassimilated silicate
                 A_E=0.65,  # [-] Activation Energy for T° scaling
                 T_ref=283,  # [K] Reference T°
                 eps_kd=0.012,  # [m2 mmolC-1] Specific attenuation coefficient

                 dt2=False,
                 **kwargs,
                 ):

        super().__init__()

        self.source_ing_P_unassimilated_to_dim = None
        self.source_ing_N_unassimilated_to_dim = None
        self.source_ing_C_unassimilated_to_dim = None
        self.source_ing_P_unassimilated_to_dom = None
        self.source_ing_N_unassimilated_to_dom = None
        self.source_ing_P_unassimilated = None
        self.source_ing_C_unassimilated_to_dom = None
        self.source_ing_N_unassimilated = None
        self.source_ing_C_unassimilated = None
        self.source_ing_N_assimilated = None
        self.source_ing_P_assimilated = None
        self.source_ing_C_assimilated = None
        self.ing_P = None
        self.ing_N = None
        self.ing_C = None
        self.lim_T = None
        self.coupled_consumers = None
        self.coupled_targets = None
        self.classname = 'Het'  # Name used as prefix for variables (used in Model.finalizeres vs varinfos)
        self.name = name
        self.g_max = g_max
        self.K_i = K_i
        self.eff_C = eff_C
        self.eff_N = eff_N
        self.eff_P = eff_P
        self.zeta_resp = zeta_resp
        self.lysrate_lin = lysrate_lin
        self.lysrate_quad = lysrate_quad
        self.mortrate_lin = mortrate_lin
        self.mortrate_quad = mortrate_quad
        self.f_unass_excr = f_unass_excr
        self.f_unass_Si = f_unass_Si
        self.A_E = A_E
        self.T_ref = T_ref
        self.eps_kd = eps_kd

        self.pref = {}
        for k, v in kwargs.items():
            if k.startswith('pref'):
                self.pref[k.split("_", maxsplit=1)[1]] = v
            # else:
            #     setattr(self, k, v)

        self.dt2 = dt2

        # Source and sink terms
        self.source_ingestion = Elms(dict=True)
        self.sink_ingestion = Elms()
        self.sink_respiration = Elms()
        self.sink_lysis = Elms()
        self.sink_mortality = Elms()

        # Initialize SMS terms
        self.C_SMS = 0.
        self.N_SMS = 0.
        self.P_SMS = 0.

    def set_coupling(self, coupled_targets,
                     coupled_consumers=None):
        # Coupling links
        self.coupled_targets = coupled_targets  # preys (for zooplankton), but can also be OM
        self.coupled_consumers = coupled_consumers
        for t in self.coupled_targets:
            if t.name not in self.pref.keys():
                self.pref[t.name] = 0.
        if sum(self.pref.values()) != 1:
            print('WARNING PREF TROUBLE for {} with sum preferences = {}'.format(self.name, sum(self.pref.values())))
            input()

    def get_sources(self, t=None):
        # Limitation functions
        self.lim_T = fns.getlimT(self.setup.T.loc[t]['T'], A_E=self.A_E, T_ref=self.T_ref, boltz=True)

        # SOURCES
        self.get_source_ingestion()

        # SINKS
        self.get_sink_mortality()
        self.get_sink_lysis()
        self.get_sink_respiration()

        # Correct the eff_X values in order to keep a constant QN and QP (given N, P are diagnostic variables)
        #       as in Kerimoglu et al. (2020)
        eff_C = self.eff_C
        eff_N = self.eff_N
        eff_P = self.eff_P

        # Respiration.C is removed here as this is the only other term that has no N and P counterpart that would
        #   be computed from QN and QP. Consequently, it needs to be removed from C_ingestion here to take this
        #   C loss into account when updating eff_C, eff_N and eff_P to conserve QN and QP
        self.ing_C = sum(self.source_ingestion.C.values()) - self.sink_respiration.C
        self.ing_N = sum(self.source_ingestion.N.values())
        self.ing_P = sum(self.source_ingestion.P.values())

        if eff_P * self.ing_P > eff_C * self.ing_C * self.QP:
            eff_P = eff_C * self.ing_C * self.QP / self.ing_P
        elif eff_P * self.ing_P < eff_C * self.ing_C * self.QP:
            eff_C = eff_P * self.ing_P / (self.ing_C * self.QP)

        if eff_N * self.ing_N > eff_C * self.ing_C * self.QN:
            eff_N = eff_C * self.ing_C * self.QN / self.ing_N
        elif eff_N * self.ing_N < eff_C * self.ing_C * self.QN:
            eff_C = eff_N * self.ing_N / (self.ing_C * self.QN)

        if eff_P * self.ing_P > eff_C * self.ing_C * self.QP:
            eff_P = eff_C * self.ing_C * self.QP / self.ing_P

        # SOURCE terms of the state equation
        self.source_ing_C_assimilated = self.ing_C * eff_C
        self.source_ing_N_assimilated = self.ing_N * eff_N
        self.source_ing_P_assimilated = self.ing_P * eff_P

        # Remaining part to be distributed elsewhere
        self.source_ing_C_unassimilated = self.ing_C - self.source_ing_C_assimilated
        self.source_ing_N_unassimilated = self.ing_N - self.source_ing_N_assimilated
        self.source_ing_P_unassimilated = self.ing_P - self.source_ing_P_assimilated
        self.source_ing_C_unassimilated_to_dom = self.source_ing_C_unassimilated * self.f_unass_excr
        self.source_ing_N_unassimilated_to_dom = self.source_ing_N_unassimilated * self.f_unass_excr
        self.source_ing_P_unassimilated_to_dom = self.source_ing_P_unassimilated * self.f_unass_excr
        self.source_ing_C_unassimilated_to_dim = self.source_ing_C_unassimilated - self.source_ing_C_unassimilated_to_dom
        self.source_ing_N_unassimilated_to_dim = self.source_ing_N_unassimilated - self.source_ing_N_unassimilated_to_dom
        self.source_ing_P_unassimilated_to_dim = self.source_ing_P_unassimilated - self.source_ing_P_unassimilated_to_dom

        self.C_sources = self.source_ing_C_assimilated
        self.N_sources = self.source_ing_N_assimilated
        self.P_sources = self.source_ing_P_assimilated

        return np.array((self.C_sources, self.N_sources, self.P_sources))

    def get_sinks(self, t=None):
        # SINKS
        self.get_sink_ingestion()
        # self.get_sink_respiration()
        # self.get_sink_lysis()
        # self.get_sink_mortality()

        # SINK terms of the state equation
        self.C_sinks = (self.sink_ingestion.C +
                        # self.sink_respiration.C + # commented out because it is already removed from ingestion (conserved QN and QP)
                        self.sink_lysis.C +
                        self.sink_mortality.C)

        self.N_sinks = (self.sink_ingestion.N +
                        self.sink_respiration.N +
                        self.sink_lysis.N +
                        self.sink_mortality.N)

        self.P_sinks = (self.sink_ingestion.P +
                        self.sink_respiration.P +
                        self.sink_lysis.P +
                        self.sink_mortality.P)

        return np.array((self.C_sinks, self.N_sinks, self.P_sinks))

    def get_source_ingestion(self):
        realpref = {}
        sumpc = np.sum([self.pref[t.name] * t.C for t in self.coupled_targets])
        for t in self.coupled_targets:
            realpref[t.name] = (self.pref[t.name] * t.C) / sumpc

        sumrpc = np.sum([realpref[t.name] * t.C for t in self.coupled_targets])
        for t in self.coupled_targets:
            self.source_ingestion.C[t.name] = self.C * self.g_max * self.lim_T * realpref[t.name] * t.C / (self.K_i + sumrpc)
            if t.N is not None:
                self.source_ingestion.N[t.name] = self.source_ingestion.C[t.name] * t.N / t.C
            else:
                self.source_ingestion.N[t.name] = 0.
            if t.P is not None:
                self.source_ingestion.P[t.name] = self.source_ingestion.C[t.name] * t.P / t.C
            else:
                self.source_ingestion.P[t.name] = 0.
            if t.Si is not None:
                self.source_ingestion.Si[t.name] = self.source_ingestion.C[t.name] * t.Si / t.C
            else:
                self.source_ingestion.Si[t.name] = 0.


    def get_sink_ingestion(self):
        self.sink_ingestion.C = fns.get_all_contributors(self.coupled_consumers, 'source_ingestion', 'C', dkey=self.name)
        self.sink_ingestion.N = fns.get_all_contributors(self.coupled_consumers, 'source_ingestion', 'N', dkey=self.name)
        self.sink_ingestion.P = fns.get_all_contributors(self.coupled_consumers, 'source_ingestion', 'P', dkey=self.name)

    def get_sink_respiration(self):
        self.sink_respiration.C = self.C * self.zeta_resp * self.lim_T
        self.sink_respiration.N = 0.
        self.sink_respiration.P = 0.

    def get_sink_lysis(self):
        self.sink_lysis.C = self.C * (self.lysrate_lin + self.C * self.lysrate_quad) * self.lim_T
        self.sink_lysis.N = self.sink_lysis.C * self.QN
        self.sink_lysis.P = self.sink_lysis.C * self.QP

    def get_sink_mortality(self):
        self.sink_mortality.C = self.C * (self.mortrate_lin + self.C * self.mortrate_quad) * self.lim_T
        self.sink_mortality.N = self.sink_mortality.C * self.QN
        self.sink_mortality.P = self.sink_mortality.C * self.QP
