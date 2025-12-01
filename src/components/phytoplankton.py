import numpy as np

from ..core.base import BaseOrg, Elms
from src.config_model import varinfos
from ..utils import functions as fns


class Phyto(BaseOrg):
    def __init__(self,
                 mu_max=3.,  # [d-1]
                 v_max_N=0.78,  # [molN molC-1 d-1] # Onur22
                 v_max_P=0.075,  # [molP molC-1 d-1] # Onur22
                 v_max_Si=1.2,  # [molSi molC-1 d-1] # Onur22
                 alpha=0.46 * 1e-5,  # [gC m-2 (g Chla µmol quanta)-1
                 QN_min=0.04,  # [gN gC-1]
                 QN_max=0.167,  # [gN gC-1]
                 QP_min=0.003,
                 QP_max=0.012,
                 QSi_min=0.06,
                 QSi_max=0.18,
                 KNH4=3 * varinfos.molmass_N / 1e3,  # [g m-3] converted from 3 µM
                 KNO3=2.,
                 KDIP=0.05,
                 KDSi=0.43,
                 thetaN_max=0.3,  # [gChla gN-1]
                 theta_max=0.07 * varinfos.molmass_C,  # [gChl molC-1] from [gChl gC-1] (Onur22)
                 T_ref=18. + varinfos.degCtoK,  # [K]
                 A_E=0.32,  # [-] Activation energy for T° scaling (Onur22)
                 QratioIC=None,  # Ratio to correct for the ICs internal pools vs Q_max (in BaseOrg.set_ICs)
                 checkQmax=True,  # Boolean: whether to check that IC are OK wrt Q_max
                 gamma_C_exud_base=.02,  # [d-1] Basal exudation rate (Onur22)
                 gamma_C_exud_prod=.11,  # [-] Production sp. exudation rate (Onur22)
                 zeta_resp_base=.01,  # [d-1] Basal respiration rate (Onur22)
                 zeta_resp_prod=.01,  # [-] Production sp. respiration rate (Onur22)
                 lysrate=0.1,  # [d-1] Onur22
                 mortrate=0.05,  # [d-1] Onur22
                 # Implicit copepod grazing closure term
                 A_E_grazing=0.65,  # [-] Activation energy for temperature scaling of grazing
                 T_ref_grazing=283.15,  # [K] Reference temperature for grazing (10°C)
                 K_grazing=20.0,  # [mmol C m-3] Half-saturation constant for grazing (Michaelis-Menten)
                 grazing_loss_max=0.0,  # [mmol C m-3 d-1] Maximum grazing loss rate (default: disabled)
                 kdvar=False,  # Whether to apply an extinction coefficient in the water column to incident light
                 eps_kd=8e-4 * varinfos.molmass_C,  # Diffuse attenuation cross section of phytoplankton [m2 mgC-1]
                 divide_water_depth_ratio=1.,  # ratio to divide water depth to restrain light attenuation
                 dt2=False,
                 dtype=np.float64,
                 name='DefaultPhyto',
                 bound_temp_to_1=True,  # Whether to bound temperature limitation to [0,1]
                 ):

        super().__init__(dtype=dtype)

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
        self.KNH4 = KNH4
        self.KNO3 = KNO3
        self.KDIP = KDIP
        self.KDSi = KDSi
        self.thetaN_max = thetaN_max
        self.theta_max = theta_max
        self.T_ref = T_ref
        self.A_E = A_E
        self.QratioIC = QratioIC
        self.checkQmax = checkQmax
        self.gamma_C_exud_base = gamma_C_exud_base
        self.gamma_C_exud_prod = gamma_C_exud_prod
        self.zeta_resp_base = zeta_resp_base
        self.zeta_resp_prod = zeta_resp_prod
        self.lysrate = lysrate
        self.mortrate = mortrate
        self.A_E_grazing = A_E_grazing
        self.T_ref_grazing = T_ref_grazing
        self.K_grazing = K_grazing
        self.grazing_loss_max = grazing_loss_max
        self.kdvar = kdvar
        self.eps_kd = eps_kd
        self.divide_water_depth_ratio = divide_water_depth_ratio
        self.dt2 = dt2
        self.name = name
        self.bound_temp_to_1 = bound_temp_to_1

        self.lim_N = None
        self.lim_P = None
        self.lim_Si = None
        self.limNUT = None
        self.limT = None
        self.fnut = None
        self.mmNH4 = None
        self.mmNO3 = None
        self.mmDIP = None
        self.mmDSi = None
        self.limI = None
        self.PAR_t = None
        self.PAR_t_water_column = None
        self.PAR_t_water_column_theoretical = None

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

        # Optimization: Lazy initialization flag for eps_kd cache
        self._eps_kds_cached = False

        # Optimization: Pre-compute temperature limitation arrays for entire simulation
        if self.setup is not None:
            self._precompute_temp_limitation(
                A_E=self.A_E, T_ref=self.T_ref, boltz=True,
                bound_temp_to_1=self.bound_temp_to_1, suffix='_growth'
            )
            self._precompute_temp_limitation(
                A_E=self.A_E_grazing, T_ref=self.T_ref_grazing, boltz=True,
                bound_temp_to_1=self.bound_temp_to_1, suffix='_grazing'
            )

    def get_sources(self, t, t_idx):
        # Optimization: Use pre-computed arrays for faster access

        # Limitation functions
        self.get_limNUT()
        self.limT = self.limT_growth_array[t_idx]
        if self.P is not None and self.Si is not None:
            self.fnut = min(self.QN, self.QP, self.QSi)

        # SOURCES
        self.get_source_PP(t, t_idx=t_idx)
        self.get_source_uptake()
        self.get_source_Chlprod()

        # SINKS
        self.get_sink_lysis()
        self.get_sink_mortality(t, t_idx=t_idx)
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
             if sources is not None], dtype=self.dtype)

    # TODO priority1 this is dangerous: if for some reason self.P is not None but self.P_sources is None then
    #       it will not be passed. Normally this should cause a fatal issue but still it can continue silently!

    def get_sinks(self, t, t_idx=None):
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
             sinks is not None], dtype=self.dtype)

    def get_diagnostic_variables(self):
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics], dtype=self.dtype)

    def get_limNUT(self):
        """Calculate nutrient limitations for Onur22 formulation."""
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

    def get_kd(self):
        """Calculate light attenuation coefficient using cached eps_kd arrays (lazy init)."""
        # Lazy initialization: compute cache on first call (after all set_coupling done)
        if not self._eps_kds_cached:
            self._cached_eps_kds = np.array([x.eps_kd for x in self.coupled_light_attenuators])
            self._cached_spm_eps_kds = np.array([x.eps_kd for x in self.coupled_SPM]) if self.coupled_SPM else None
            self._eps_kds_cached = True

        # Compute kd using cached eps_kd values
        Cs = np.array([x.C for x in self.coupled_light_attenuators])
        self.kd = self.eps_kd * self.C + np.sum(self._cached_eps_kds * Cs)

        # SPM contribution (if applicable)
        if self.coupled_SPM:
            spm_masses = np.array([x.massconcentration for x in self.coupled_SPM])
            self.kd += np.sum(self._cached_spm_eps_kds * spm_masses)

    def get_source_PP(self, t, t_idx=None):
        """Calculate primary production for Onur22 formulation."""
        # Optimization: Use pre-computed arrays for faster access
        self.PAR_t = self.setup.PAR_array[t_idx]
        if self.kdvar:
            self.get_kd()
            self.PAR_t_water_column_theoretical = self.PAR_t * np.exp(-self.kd * self.setup.water_depth_array[t_idx])
            self.PAR_t_water_column = self.PAR_t * np.exp(-self.kd * self.setup.water_depth_array[t_idx] /
                                                          self.divide_water_depth_ratio)

        self.PC_max = self.mu_max * self.limNUT * self.limT
        self.limI = 1 - np.exp((-self.alpha * self.thetaC * self.PAR_t_water_column) / (
                self.PC_max * varinfos.molmass_C)) if self.PC_max > 0 else 0
        self.PC = self.PC_max * self.limI
        self.source_PP.C = self.PC * self.C

        # Calculate rho_Chl for chlorophyll production
        # rho_chl must be in [mgChl mmolN-1] (see source_Chlprod.Chl)
        # QN_max is in [molN:molC], theta_max is in [mgChl molC-1]
        self.rho_Chl = 0. if self.PAR_t_water_column < 0.01 else self.theta_max / self.QN_max * (
                self.PC * varinfos.molmass_C / (self.alpha * self.thetaC * self.PAR_t_water_column))

    def get_source_uptake(self):
        """Calculate nutrient uptake for Onur22 formulation."""
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

    def get_source_Chlprod(self):
        """
        Calculate chlorophyll production for Onur22 formulation.
        dPhyChl = VN*rho_chl (without self.C)
        rho_Chl is in [mgChl mmolN-1], source_uptake.N is in [mmolN m-3 d-1]
        Hence [mgChl m-3 d-1] = [mgChl mmolN-1] * [mmolN m-3 d-1]
        """
        self.source_Chlprod.Chl = self.rho_Chl * self.source_uptake.N

    def get_sink_lysis(self):
        flysis = 0.05 + 0.95 * (1. + np.exp(-(self.fnut - 0.2) * 30)) ** (-1)
        self.sink_lysis.C = self.C * flysis * self.lysrate
        self.sink_lysis.N = self.sink_lysis.C * self.QN
        self.sink_lysis.Chl = self.sink_lysis.C * self.thetaC
        self.sink_lysis.P = self.sink_lysis.C * self.QP
        self.sink_lysis.Si = self.sink_lysis.C * self.QSi

    def get_sink_mortality(self, t, t_idx=None):
        # Optimization: Use pre-computed temperature limitation for grazing
        grazing_loss = (self.grazing_loss_max *
                        self.limT_grazing_array[t_idx] *
                        self.C / (self.K_grazing + self.C))
        self.sink_mortality.C = self.C * self.mortrate + grazing_loss
        self.sink_mortality.N = self.sink_mortality.C * self.QN
        self.sink_mortality.Chl = self.sink_mortality.C * self.thetaC
        self.sink_mortality.P = self.sink_mortality.C * self.QP
        self.sink_mortality.Si = self.sink_mortality.C * self.QSi

    def get_sink_exudation(self):
        """Calculate exudation for Onur22 formulation."""
        self.sink_exudation.C = self.gamma_C_exud_base * self.C + self.gamma_C_exud_prod * self.PC
        self.sink_exudation.Chl = 0.
        xquota = 0.99
        self.sink_exudation.N = self.sink_exudation.C * self.QN if self.limQUOTAmin.N > xquota else 0.
        self.sink_exudation.P = self.sink_exudation.C * self.QP if self.limQUOTAmin.P > xquota else 0.
        self.sink_exudation.Si = self.sink_exudation.C * self.QSi if self.limQUOTAmin.Si > xquota else 0.
        self.frac_exud_small = 1 / (1 + np.exp(-(self.fnut - 0.2) * 30))

    def get_sink_respiration(self):
        """Calculate respiration for Onur22 formulation."""
        self.sink_respiration.C = self.zeta_resp_base * self.C + self.zeta_resp_prod * self.PC
        self.sink_respiration.Chl = self.sink_respiration.C * self.thetaC  # Chl degradation
        self.sink_respiration.N = 0.
        self.sink_respiration.P = 0.
        self.sink_respiration.Si = 0.

    def get_sink_ingestion(self):
        if self.coupled_consumer is not None:
            self.sink_ingestion.C = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion', 'C')
            self.sink_ingestion.N = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion', 'N')
            self.sink_ingestion.Chl = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion',
                                                               'C') * self.thetaC

        else:
            self.sink_ingestion.C = 0.
            self.sink_ingestion.N = 0.
            self.sink_ingestion.Chl = 0.
            self.sink_ingestion.P = 0.
            self.sink_ingestion.Si = 0.

    def get_sink_aggregation(self):
        """Calculate aggregation for Onur22 formulation."""
        self.sink_aggregation.C = self.coupled_aggreg_target.aggPhy_C
        self.sink_aggregation.N = self.coupled_aggreg_target.aggPhy_N
        self.sink_aggregation.P = self.coupled_aggreg_target.aggPhy_P
        self.sink_aggregation.Si = self.coupled_aggreg_target.aggPhy_Si
        self.sink_aggregation.Chl = self.coupled_aggreg_target.aggPhy_Chl
