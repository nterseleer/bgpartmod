import numpy as np

from ..core.base import BaseOrg, Elms
from src.config_model import varinfos
from ..utils import functions as fns


class Phyto(BaseOrg):
    def __init__(self,
                 mu_max=3.,  # [d-1]
                 v_max=0.5,  # [gN gC-1 d-1]
                 v_max_N=0.78,  # [molN mmolC-1 d-1] # Onur22
                 v_max_P=0.075,  # [molP mmolC-1 d-1] # Onur22
                 v_max_Si=1.2,  # [molSi mmolC-1 d-1] # Onur22
                 alpha=0.46 * 1e-5,  # [gC m-2 (g Chla µmol quanta)-1
                 QN_min=0.04,  # [gN gC-1]
                 QN_max=0.167,  # [gN gC-1]
                 QP_min=0.003,
                 QP_max=0.012,
                 QSi_min=0.06,
                 QSi_max=0.18,
                 boundQN=True,  # [Boolean]
                 n=0.01,  # [-]
                 KNH4=3 * varinfos.molmass_N / 1e3,  # [g m-3] converted from 3 µM
                 KNO3=2.,
                 KDIP=0.05,
                 KDSi=0.43,
                 thetaN_max=0.3,  # [gChla gN-1]

                 theta_max=0.07,  # [gChl gC -1] (Onur22)

                 T_ref=18. + varinfos.degCtoK,  # [K]
                 A_E=0.32,  # [-] Activation energy for T° scaling (Onur22)
                 r_ref=0.025,  # [d-1]

                 sigmaN_C=1000 * varinfos.molmass_N ** 2 / varinfos.molmass_C ** 2,
                 # [mmolN2 mmolC-2] Slope parameter for DIN-uptake regulation (Sch07)
                 gamma_C=0.,  # 0.250,  # [d-1] Phytoplankton linear C loss rate (Sch07)
                 gamma_N=0.,  # 0.180,  # [d-1] Phytoplankton linear N loss rate (Sch07)
                 gamma_chl=0.,  # 0.001,  # [d-1] Chla degradation rate (Sch07)

                 gamma_C_exud_base=.02,  # [d-1] Basal exudation rate (Onur22)
                 gamma_C_exud_prod=.11,  # [-] Production sp. exudation rate (Onur22)

                 phi_PP=0.02 / (varinfos.molmass_N ** 2) * 1e6,  # [m6 mmolN-2 d-1]  (Sch07)
                 phi_PD=0.,  # [?] Deviation from a priori guess on Phy-Det aggregation
                 zeta=2.,  # [gC gN-1]

                 zeta_resp_base=.01,  # [d-1] Basal respiration rate (Onur22)
                 zeta_resp_prod=.01,  # [-] Production sp. respiration rate (Onur22)

                 beta_fact=0.033,  # see Sch07 Equ B13
                 k_beta=53.125,  # half-sat cst for aggregation vs TEPC cf Ruiz et al 2002

                 lysrate=0.1,  # [d-1] Onur22
                 mortrate=0.05,  # [d-1] Onur22

                 fixedstoichiometry=False,  # [Boolean]
                 boundrhoChl=False,  # [Boolean]
                 r_ref_tweakGMK=False,

                 kdvar=False,  # Whether to apply an extinction coefficient in the water column to incident light
                 eps_kd=8e-4 * varinfos.molmass_C,  # Diffuse attenuation cross section of phytoplankton [m2 mgC-1]

                 dt2=False,
                 name='DefaultPhyto',
                 # formulation='Sch07'
                 ):

        super().__init__()

        self.frac_exud_small = None
        self.rho_Chl = None
        self.PC = None
        self.PC_max = None
        self.kd = None
        self.t = None
        self.formulation = None
        self.diagnostics = None
        self.classname = 'Phy'  # Name used as prefix for variables (used in Model.finalizeres vs varinfos)
        self.mu_max = mu_max
        self.v_max = v_max
        self.v_max_N = v_max_N
        self.v_max_P = v_max_P
        self.v_max_Si = v_max_Si
        self.alpha = alpha
        self.QN_min = QN_min
        self.QN_max = QN_max
        self.QP_min = QP_min
        self.QP_max = QP_max
        self.QSi_min = QSi_min
        self.QSi_max = QSi_max
        self.boundQN = boundQN
        self.n = n
        self.KNH4 = KNH4
        self.KNO3 = KNO3
        self.KDIP = KDIP
        self.KDSi = KDSi
        self.thetaN_max = thetaN_max
        self.theta_max = theta_max
        self.T_ref = T_ref
        self.A_E = A_E
        self.r_ref = r_ref
        self.sigmaN_C = sigmaN_C
        self.gamma_C = gamma_C
        self.gamma_N = gamma_N
        self.gamma_chl = gamma_chl
        self.gamma_C_exud_base = gamma_C_exud_base
        self.gamma_C_exud_prod = gamma_C_exud_prod
        self.phi_PP = phi_PP
        self.phi_PD = phi_PD
        self.zeta = zeta
        self.zeta_resp_base = zeta_resp_base
        self.zeta_resp_prod = zeta_resp_prod
        self.beta_fact = beta_fact
        self.lysrate = lysrate
        self.mortrate = mortrate
        self.k_beta = k_beta
        self.fixedstoichiometry = fixedstoichiometry
        self.boundrhoChl = boundrhoChl
        self.r_ref_tweakGMK = r_ref_tweakGMK
        self.kdvar = kdvar
        self.eps_kd = eps_kd
        self.dt2 = dt2
        self.name = name

        self.lim_N = None
        self.lim_P = None
        self.lim_Si = None
        self.lim_QN = None
        self.lim_QN2 = None
        self.limNUT = None
        self.limT = None
        self.fnut = None
        self.mmNH4 = None
        self.mmNO3 = None
        self.mmDIP = None
        self.mmDSi = None
        self.limI = None


        self.limQUOTA = Elms()
        self.limQUOTAmin = Elms()

        # Source and sink terms
        self.source_PP = Elms()
        self.source_uptake = Elms()
        self.source_Chlprod = Elms()
        self.sink_lysis = Elms()
        self.sink_mortality = Elms()
        self.sink_exudation = Elms()
        self.sink_respiration = Elms()
        self.sink_ingestion = Elms()
        self.sink_aggregation = Elms()

        # Coupling links
        self.coupled_consumer = None
        self.coupled_NH4 = None
        self.coupled_NO3 = None
        self.coupled_DIP = None
        self.coupled_DSi = None
        self.coupled_TEPC = None
        self.coupled_Det = None
        self.coupled_SPM = None
        self.coupled_aggreg_target = None
        self.coupled_light_attenuators = None

    def set_coupling(self, coupled_consumer=None,
                     coupled_NH4=None, coupled_NO3=None, coupled_DIP=None, coupled_DSi=None,
                     coupled_TEPC=None, coupled_Det=None, coupled_SPM=None,
                     coupled_aggreg_target=None,
                     coupled_light_attenuators=None):
        # Coupling links
        self.coupled_consumer = coupled_consumer
        self.coupled_NH4 = coupled_NH4
        self.coupled_NO3 = coupled_NO3
        self.coupled_DIP = coupled_DIP
        self.coupled_DSi = coupled_DSi
        self.coupled_TEPC = coupled_TEPC
        self.coupled_Det = coupled_Det
        self.coupled_SPM = coupled_SPM
        self.coupled_aggreg_target = coupled_aggreg_target
        self.coupled_light_attenuators = coupled_light_attenuators

    def get_sources(self, t):
        # Limitation functions
        self.get_limNUT()
        if self.formulation == 'Onur22':
            self.limT = fns.getlimT(self.setup.T, A_E=self.A_E, T_ref=self.T_ref, boltz=True)
        else:
            self.limT = fns.getlimT(self.setup.T)
        if self.P is not None and self.Si is not None:
            self.fnut = min(self.QN, self.QP, self.QSi)

        # SOURCES
        self.get_source_PP(t)
        self.get_source_uptake()
        self.get_source_Chlprod()

        # SINKS
        self.get_sink_lysis()
        self.get_sink_mortality()
        self.get_sink_exudation()
        self.get_sink_respiration()

        # SOURCE terms of the state equation
        self.C_sources = self.source_PP.C

        self.N_sources = self.source_uptake.N
        self.Chl_sources = self.source_Chlprod.Chl
        if self.P is not None:
            self.P_sources = self.source_uptake.P
        if self.Si is not None:
            self.Si_sources = self.source_uptake.Si

        return np.array(
            [sources for sources in (self.C_sources, self.N_sources, self.Chl_sources, self.P_sources, self.Si_sources)
             if sources is not None])

    # TODO priority1 this is dangerous: if for some reason self.P is not None but self.P_sources is None then
    #       it will not be passed. Normally this should cause a fatal issue but still it can continue silently!

    def get_sinks(self, t):
        # SINKS
        # self.get_sink_lysis()
        # self.get_sink_mortality()
        # self.get_sink_exudation()
        # self.get_sink_respiration()
        self.get_sink_ingestion()
        self.get_sink_aggregation()

        # SINK terms of the state equation
        self.C_sinks = (self.sink_respiration.C +
                        self.sink_exudation.C +
                        self.sink_aggregation.C +
                        self.sink_ingestion.C +
                        self.sink_lysis.C +
                        self.sink_mortality.C)

        self.N_sinks = (self.sink_respiration.N +
                        self.sink_exudation.N +
                        self.sink_aggregation.N +
                        self.sink_ingestion.N +
                        self.sink_lysis.N +
                        self.sink_mortality.N)

        self.Chl_sinks = (self.sink_respiration.Chl +
                          self.sink_exudation.Chl +
                          self.sink_aggregation.Chl +
                          self.sink_ingestion.Chl +
                          self.sink_lysis.Chl +
                          self.sink_mortality.Chl)

        if self.P is not None:
            self.P_sinks = (self.sink_respiration.P +
                            self.sink_exudation.P +
                            self.sink_aggregation.P +
                            self.sink_ingestion.P +
                            self.sink_lysis.P +
                            self.sink_mortality.P)

        if self.Si is not None:
            self.Si_sinks = (self.sink_respiration.Si +
                             self.sink_exudation.Si +
                             self.sink_aggregation.Si +
                             self.sink_ingestion.Si +
                             self.sink_lysis.Si +
                             self.sink_mortality.Si)

        return np.array(
            [sinks for sinks in (self.C_sinks, self.N_sinks, self.Chl_sinks, self.P_sinks, self.Si_sinks) if
             sinks is not None])

    def get_diagnostic_variables(self):
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics])


    def get_limNUT(self, verbose=True):
        if self.formulation == 'Onur22':
            self.lim_N = 1 - self.QN_min / self.QN
            self.lim_P = 1 - self.QP_min / self.QP
            self.lim_Si = 1 - self.QSi_min / self.QSi
            self.limNUT = min(self.lim_N, self.lim_P, self.lim_Si)
            self.limQUOTAmin.N = (self.QN - self.QN_min) / (self.QN_max - self.QN_min)
            self.limQUOTA.N = 1 - self.limQUOTAmin.N
            self.limQUOTAmin.P = (self.QP - self.QP_min) / (self.QP_max - self.QP_min)
            self.limQUOTA.P = 1 - self.limQUOTAmin.P
            self.limQUOTAmin.Si = (self.QSi - self.QSi_min) / (self.QSi_max - self.QSi_min)
            self.limQUOTA.Si = 1 - self.limQUOTAmin.Si


        else:
            self.lim_N = self.coupled_NH4.concentration / (self.coupled_NH4.concentration + self.KNH4)
            if self.boundQN:
                self.lim_QN = max((self.QN - self.QN_min), 0) / (self.QN_max - self.QN_min)
                self.lim_QN2 = (max((self.QN_max - self.QN), 0) / (self.QN_max - self.QN_min)) ** self.n
            else:
                self.lim_QN = (self.QN - self.QN_min) / (self.QN_max - self.QN_min)
                self.lim_QN2 = ((self.QN_max - self.QN) / (self.QN_max - self.QN_min)) ** self.n
            if verbose and self.formulation == 'GMK98' and (self.QN_max - self.QN) < 0:
                # Then lim_QN2 is np.nan and so is self.VC_N etc...
                # Note: lim_QN2 not used in Schartau
                print('(self.QN_max - self.QN)<0 !!! at t=', self.t)
            if verbose and (self.QN - self.QN_min) < 0:
                print('(self.QN - self.QN_min)<0 !!! at t=', self.t)
            self.limNUT = self.lim_QN

    def get_kd(self):
        if self.kdvar:
            if self.formulation == 'GMK98':
                self.kd = self.setup.k_att * self.setup.Chl_tot  # TBC !!!

            elif self.formulation == "Onur22":
                self.kd = self.eps_kd * self.C + np.sum([x.eps_kd * x.C for x in self.coupled_light_attenuators])
                if self.coupled_SPM:
                    # print(["{} with eps {}".format(x.name, x.eps_kd) for x in self.coupled_SPM])
                    # print(self.kd)
                    self.kd += np.sum([x.eps_kd * x.massconcentration for x in self.coupled_SPM])
                    # print(self.kd)

            else:
                self.kd = (self.setup.kb +
                           self.coupled_SPM.eps_kd * np.sqrt(self.coupled_SPM.concentration) *1e9 +
                           self.coupled_Det.eps_kd * self.coupled_Det.C +
                           self.eps_kd * self.setup.Cphy_tot)

    def get_source_PP(self, t):
        I_t = self.setup.PAR.loc[t]['PAR']
        if self.kdvar:
            self.get_kd()
            I_t = I_t * np.exp(-self.kd * self.setup.z)

        if not self.fixedstoichiometry:
            self.PC_max = self.mu_max * self.limNUT * self.limT  # OK
            self.limI = 1 - np.exp((-self.alpha * self.thetaC * I_t) / self.PC_max) if self.PC_max > 0 else 0  # Eq. 4
            self.PC = self.PC_max * self.limI
            self.source_PP.C = self.PC * self.C


            if self.boundrhoChl:
                self.rho_Chl = 0. if not light else self.thetaN_max * self.PC / (self.alpha * self.thetaC * I_t_LD)
            else:
                # TODO check this
                if self.formulation == "Onur22":
                    self.rho_Chl = 0. if I_t<0.01 else self.theta_max / self.QN_max * self.PC / (
                            self.alpha * self.thetaC * I_t) * 12.

                else:
                    self.rho_Chl = self.thetaN_max * self.PC / (
                            self.alpha * self.thetaC * I_t)  # TBC HERE I_t_LD or I ?? # TROUBLE
        else:
            self.limI = 1 - np.exp((-self.alpha * self.thetaC * I_t_LD) / self.mu_max)
            self.source_PP.C = self.mu_max * self.lim_N * self.limI * self.C
            self.rho_Chl = None

    def get_source_uptake(self):
        if not self.fixedstoichiometry:
            # GMK1998
            if self.formulation == 'GMK98':
                VC_max = self.v_max * self.lim_QN2 * self.limT
                self.VC_N = VC_max * self.lim_N
                self.source_uptake.N = self.VC_N / self.QN * self.N
            elif self.formulation == 'Sch07':
                RN_C = 1 - np.exp(-self.sigmaN_C * (np.abs(self.QN - self.QN_max) - (self.QN - self.QN_max)) ** 2)
                self.VC_N = self.PC_max * self.QN_max * RN_C * self.lim_N
                self.source_uptake.N = self.VC_N / self.QN * self.N
            elif self.formulation == 'Onur22':
                wNH4 = self.coupled_NH4.concentration / self.KNH4
                wNO3 = self.coupled_NO3.concentration / self.KNO3
                self.mmNH4 = wNH4 / (1 + wNH4 + wNO3)
                self.mmNO3 = wNO3 / (1 + wNH4 + wNO3)
                self.mmDIP = self.coupled_DIP.concentration / (self.coupled_DIP.concentration + self.KDIP)
                self.mmDSi = self.coupled_DSi.concentration / (self.coupled_DSi.concentration + self.KDSi)
                self.source_uptake.NH4 = self.v_max_N * self.limQUOTA.N * self.mmNH4 * self.limT * self.C
                self.source_uptake.NO3 = self.v_max_N * self.limQUOTA.N * self.mmNO3 * self.limT * self.C
                self.source_uptake.N = self.source_uptake.NH4 + self.source_uptake.NO3
                self.source_uptake.P = self.v_max_P * self.limQUOTA.P * self.mmDIP * self.limT * self.C
                self.source_uptake.Si = self.v_max_Si * self.limQUOTA.Si * self.mmDSi * self.limT * self.C
        else:
            self.source_uptake.N = self.source_PP.C * self.QN

    def get_source_Chlprod(self):
        if not self.fixedstoichiometry:
            if self.formulation == 'Onur22':
                """
                20240702 - dPhyChl = VN*rho_chl - Without self.C
                Indeed, rho_Chl is expressed as gChl/molN!
                """
                self.source_Chlprod.Chl = self.rho_Chl * self.source_uptake.N
            else:
                # GMK1998, Sch07
                self.source_Chlprod.Chl = self.C * self.rho_Chl * self.VC_N  # equivalent to self.VC_N / self.thetaC * self.rho_Chl * self.Chl
        else:
            self.source_Chlprod.Chl = self.source_PP.C * self.thetaC

    def get_sink_lysis(self):
        if self.formulation == 'Onur22':
            flysis = 0.05 + 0.95 * (1. + np.exp(-(self.fnut - 0.2) * 30)) ** (-1)
            self.sink_lysis.C = self.C * flysis * self.lysrate
            self.sink_lysis.N = self.sink_lysis.C * self.QN
            self.sink_lysis.Chl = self.sink_lysis.C * self.thetaC
            self.sink_lysis.P = self.sink_lysis.C * self.QP
            self.sink_lysis.Si = self.sink_lysis.C * self.QSi
        else:
            self.sink_lysis.C = 0
            self.sink_lysis.N = 0
            self.sink_lysis.Chl = 0
            self.sink_lysis.P = 0
            self.sink_lysis.Si = 0

    def get_sink_mortality(self):
        if self.formulation == 'Onur22':
            self.sink_mortality.C = self.C * self.mortrate
            self.sink_mortality.N = self.sink_mortality.C * self.QN
            self.sink_mortality.Chl = self.sink_mortality.C * self.thetaC
            self.sink_mortality.P = self.sink_mortality.C * self.QP
            self.sink_mortality.Si = self.sink_mortality.C * self.QSi
        else:
            self.sink_mortality.C = 0
            self.sink_mortality.N = 0
            self.sink_mortality.Chl = 0
            self.sink_mortality.P = 0
            self.sink_mortality.Si = 0

    def get_sink_exudation(self):
        if self.formulation == 'Onur22':
            self.sink_exudation.C = self.gamma_C_exud_base * self.C + self.gamma_C_exud_prod * self.PC
            self.sink_exudation.Chl = 0.
            xquota = 0.99
            self.sink_exudation.N = self.sink_exudation.C * self.QN if self.limQUOTAmin.N > xquota else 0.
            self.sink_exudation.P = self.sink_exudation.C * self.QP if self.limQUOTAmin.P > xquota else 0.
            self.sink_exudation.Si = self.sink_exudation.C * self.QSi if self.limQUOTAmin.Si > xquota else 0.
            self.frac_exud_small = 1 / (1 + np.exp(-(self.fnut - 0.2) * 30))
        else:
            self.sink_exudation.C = self.gamma_C * self.C
            self.sink_exudation.N = self.gamma_N * self.N
            self.sink_exudation.Chl = self.gamma_chl * self.Chl

    def get_sink_respiration(self):
        if self.formulation == 'Onur22':
            self.sink_respiration.C = self.zeta_resp_base * self.C + self.zeta_resp_prod * self.PC
            self.sink_respiration.Chl = self.sink_respiration.C * self.thetaC  # "Chl degradation"
            self.sink_respiration.N = 0.
            self.sink_respiration.P = 0.
            self.sink_respiration.Si = 0.
        else:
            # Respiration (maintenance metabolic cost) + biosynthetic costs
            if self.r_ref_tweakGMK:
                self.r_ref = self.r_ref * self.source_PP.C
            self.sink_respiration.C = (self.r_ref + self.zeta * self.VC_N) * self.C

            if self.formulation == 'GMK98':
                self.sink_respiration.N = self.r_ref * self.N
                self.sink_respiration.Chl = self.r_ref * self.Chl
            else:
                self.sink_respiration.N = 0.
                self.sink_respiration.Chl = 0.

    def get_sink_ingestion(self):
        if self.coupled_consumer is not None:
            self.sink_ingestion.C = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion', 'C')
            self.sink_ingestion.N = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion', 'N')
            self.sink_ingestion.Chl = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion', 'C') * self.thetaC

        else:
            self.sink_ingestion.C = 0.
            self.sink_ingestion.N = 0.
            self.sink_ingestion.Chl = 0.
            self.sink_ingestion.P = 0.
            self.sink_ingestion.Si = 0.

    def get_sink_aggregation(self, ):
        if self.formulation == 'Onur22':
            self.sink_aggregation.C = self.coupled_aggreg_target.aggPhy_C
            self.sink_aggregation.N = self.coupled_aggreg_target.aggPhy_N
            self.sink_aggregation.P = self.coupled_aggreg_target.aggPhy_P
            self.sink_aggregation.Si = self.coupled_aggreg_target.aggPhy_Si
            self.sink_aggregation.Chl = self.coupled_aggreg_target.aggPhy_Chl

        else:
            if self.coupled_TEPC is not None:
                # TBC - 20231009 - the 1/self.QN is TBC (partially present in Sch07 Equ B13)
                beta = self.beta_fact / self.QN * self.coupled_TEPC.concentration / (
                        self.k_beta + self.coupled_TEPC.concentration)
                self.sink_aggregation.N = self.phi_PP * self.N ** 2 + self.phi_PD * beta * self.N * self.coupled_Det.N
                self.sink_aggregation.C = self.sink_aggregation.N / self.QN
                self.sink_aggregation.Chl = self.sink_aggregation.C * self.thetaC
            else:
                self.sink_aggregation.N = 0.
                self.sink_aggregation.C = 0.
                self.sink_aggregation.Chl = 0.
