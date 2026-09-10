import numpy as np
from scipy.special import exp1, gamma, gammaincc

from ..core.base import BaseOrg, Elms
from src.config_model import varinfos
from ..utils import functions as fns


def _floor_clip(x):
    """Clip a limitation factor to [1e-6, 1] (numerical-protection mode)."""
    return np.clip(x, 1e-6, 1.)


def _floor_min(x):
    """Floor a limitation factor at 1e-6, no upper bound (default mode)."""
    return max(1e-6, x)

_M_EPS = 1e-6  # slope below which the exponential biomass profile is numerically
               # indistinguishable from homogeneous: use the exact E1 form there
               # (gamma(s)*gammaincc(s,x) loses precision as s = m/kd -> 0).

class Phyto(BaseOrg):
    def __init__(self,
                 mu_max=5.2,  # [d-1] (Kerimoglu22)
                 v_max_N=0.78,  # [molN molC-1 d-1] (Kerimoglu22)
                 v_max_P=0.075,  # [molP molC-1 d-1] (Kerimoglu22)
                 v_max_Si=1.2,  # [molSi molC-1 d-1] (Kerimoglu22)
                 alpha=7.e-6,  # [mgC mgChl-1 (µmol quanta m-2)-1] (Kerimoglu22)
                 QN_min=0.05,  # [molN molC-1] (Kerimoglu22)
                 QN_max=0.15,  # [molN molC-1] (Kerimoglu22)
                 QP_min=0.003,  # [molP molC-1] (Kerimoglu22)
                 QP_max=0.012,  # [molP molC-1] (Kerimoglu22)
                 QSi_min=0.06,  # [molSi molC-1] (Kerimoglu22)
                 QSi_max=0.18,  # [molSi molC-1] (Kerimoglu22)
                 KNH4=0.5,  # [mmolN m-3] (Kerimoglu22)
                 KNO3=2.,  # [mmolN m-3] (Kerimoglu22)
                 KDIP=0.05,  # [mmolP m-3] (Kerimoglu22)
                 KDSi=0.43,  # [mmolSi m-3] (Kerimoglu22)
                 theta_max=0.07 * varinfos.molmass_C,  # [gChl molC-1] from [gChl gC-1] (Kerimoglu22)
                 T_ref=283.15,  # [K] (Kerimoglu22)
                 A_E=0.32,  # [-] Activation energy for T° scaling (Kerimoglu22)
                 QratioIC=None,  # Ratio to correct for the ICs internal pools vs Q_max (in BaseOrg.set_ICs)
                 checkQmax=True,  # Boolean: whether to check that IC are OK wrt Q_max
                 gamma_C_exud_base=.02,  # [d-1] Basal exudation rate (Kerimoglu22)
                 gamma_C_exud_prod=.11,  # [-] Production sp. exudation rate (Kerimoglu22)
                 zeta_resp_base=.01,  # [d-1] Basal respiration rate (Kerimoglu22)
                 zeta_resp_prod=.01,  # [-] Production sp. respiration rate (Kerimoglu22)
                 lysrate=0.1,  # [d-1] (Kerimoglu22)
                 mortrate=0.05,  # [d-1] (Kerimoglu22)
                 # Implicit copepod grazing closure term
                 A_E_grazing=0.65,  # [-] Activation energy for temperature scaling of grazing
                 T_ref_grazing=283.15,  # [K] Reference temperature for grazing (10°C)
                 K_grazing=20.0,  # [mmol C m-3] Half-saturation constant for grazing (Michaelis-Menten)
                 grazing_loss_max=0.0,  # [mmol C m-3 d-1] Maximum grazing loss rate (default: disabled)
                 grazing_holling_exponent=1.0,  # [-] Holling exponent: 1.0=type II (MM), 2.0=type III (sigmoidal)
                 eps_kd=0.024,  # [m2 mmolC-1] Diffuse attenuation cross section of phytoplankton (Kerimoglu22)
                 prescribe_SPM_from_setup=False,  # Whether to use prescribed SPM from Setup (for BGC-only runs)
                 prescribed_SPM_eps_kd=0.066e3,  # [m² kg⁻¹] Light attenuation cross section for prescribed SPM
                 dt2=False,
                 dtype=np.float64,
                 name='DefaultPhyto',
                 bound_temp_to_1=True,  # Whether to bound temperature limitation to [0,1]
                 apply_numerical_protections=False,  # Whether to apply numerical safeguards
                 I_min = 1, # [µmol photon m-2 d-1] Non-zero irradiance floor bounding the
                            # light-integration depth z_upper_layer (numerical stability,
                            # NOT the photic-zone limit). See get_source_PP.
                 biomass_profile_slope_base=0.0,  # [m-1] m0, baseline slope of the exponential biomass
                                             # profile B(z)=B0*exp(-m*z); 0 = homogeneous (legacy m=0)
                 biomass_profile_slope_seasonal_amp=0.0,  # [m-1] seasonal amplitude of the slope; 0 = none
                 biomass_profile_slope_peak_doy=165,  # day-of-year of the seasonal max (~mid-June,
                                                      # Silori 2025 Fig 4b)
                 compound_production_factor=1.0,  # [-] lumped multiplier on realised production for
                                                  # 0D-unresolved processes; 1.0 = off
                 ):

        super().__init__(dtype=dtype)

        self.frac_exud_small = None
        self.rho_Chl = None
        self.PC = None
        self.PC_max = None
        self.kd = None
        self.t = None
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
        self.grazing_holling_exponent = grazing_holling_exponent
        self.eps_kd = eps_kd
        self.prescribe_SPM_from_setup = prescribe_SPM_from_setup
        self.prescribed_SPM_eps_kd = prescribed_SPM_eps_kd
        self.dt2 = dt2
        self.name = name
        self.bound_temp_to_1 = bound_temp_to_1
        self.apply_numerical_protections = apply_numerical_protections
        # Floor policy fixed once at construction (the flag never changes at runtime):
        # clip to [1e-6, 1] under numerical protection, else just floor at 1e-6.
        self._floor = _floor_clip if apply_numerical_protections else _floor_min
        self._C_MIN = 1e-12  # [mmol m-3] Minimum pool for safe calculations
        self.I_min = I_min
        self.biomass_profile_slope_base = biomass_profile_slope_base
        self.biomass_profile_slope_seasonal_amp = biomass_profile_slope_seasonal_amp
        self.biomass_profile_slope_peak_doy = biomass_profile_slope_peak_doy
        self.compound_production_factor = compound_production_factor
        # Effective biomass-profile slope per timestep (precomputed in set_coupling) and its
        # current value (diagnostic, set in get_source_PP).
        self.biomass_profile_slope_array = None
        self.biomass_profile_slope_t = None

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
        self.limI_upper_layer = None
        self.limI_theoretical = None
        self.z_upper_layer = None
        self.PAR_t = None
        self.I_H = None  # Irradiance reaching the seabed [µE m-2 d-1] (diagnostic)

        self.limQUOTA = Elms()
        self.limQUOTAmin = Elms()

        # Source and sink terms
        self.source_PP = Elms()
        self.source_uptake = Elms()
        self.source_Chlprod = Elms()
        self.sink_lysis = Elms()
        self.sink_mortality = Elms()
        self.sink_grazing = Elms()  # DIAGNOSTIC-ONLY mirror of the grazing-closure flux (see get_sink_mortality)
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

        self.coupled_Micro_in_Macro_massconcentration = None
        self.coupled_Microflocs_massconcentration = None

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

        # Mineral pools attenuating light. Empty in a configuration without the
        # flocculation module (the mineral term of kd is then simply 0), which is the
        # case of the reference biogeochemical configuration base_config.Kerimoglu2022.
        if self.prescribe_SPM_from_setup:
            from ..components.flocs import PrescribedFlocs
            self.coupled_SPM = [
                PrescribedFlocs(name="Microflocs", eps_kd=self.prescribed_SPM_eps_kd),
                PrescribedFlocs(name="Micro_in_Macro", eps_kd=self.prescribed_SPM_eps_kd)
            ]
        else:
            self.coupled_SPM = coupled_SPM or []
        self._microflocs_idx = self._spm_index('Microflocs')
        self._micro_in_macro_idx = self._spm_index('Micro_in_Macro')

        self.coupled_aggreg_target = coupled_aggreg_target
        self.coupled_light_attenuators = coupled_light_attenuators or []

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

            # Pre-compute the effective biomass-profile slope m(t) for the biomass-weighted
            # light-limitation path. doy is taken from setup.dates (same indexing as PAR_array
            # and the limT arrays). With seasonal_amp=0 this is a constant array (= baseline
            # slope), so biomass_profile_slope_base=0 reproduces the homogeneous (legacy m=0) case.
            # np.maximum(0, .) enforces m >= 0.
            # A fractional (intra-day) component is added to dayofyear so that m(t) varies
            # smoothly at the setup timestep instead of stepping once per calendar day.
            intraday_frac = (
                self.setup.dates.hour * 3600
                + self.setup.dates.minute * 60
                + self.setup.dates.second
            ).to_numpy() / 86400.0
            doy = self.setup.dates.dayofyear.to_numpy().astype(self.dtype) + intraday_frac.astype(self.dtype)

            seasonal_factor = (1.0 + np.cos(2 * np.pi * (doy - self.biomass_profile_slope_peak_doy) / 365.25)) / 2.0
            self.biomass_profile_slope_array = np.maximum(
                0.0,
                self.biomass_profile_slope_base + self.biomass_profile_slope_seasonal_amp * seasonal_factor
            )

    def _spm_index(self, name):
        """Index of a mineral pool in coupled_SPM, or None when it is not coupled."""
        return next((i for i, spm in enumerate(self.coupled_SPM) if spm.name == name), None)

    def get_sources(self, t, t_idx):
        # Optimization: Use pre-computed arrays for faster access

        # Update prescribed SPM if applicable
        if self.prescribe_SPM_from_setup and self.coupled_SPM:
            self.coupled_SPM[0].massconcentration = self.setup.Microflocs_massconc_array[t_idx]
            self.coupled_SPM[1].massconcentration = self.setup.Micro_in_Macro_massconc_array[t_idx]
        self.coupled_Microflocs_massconcentration = (
            self.coupled_SPM[self._microflocs_idx].massconcentration
            if self._microflocs_idx is not None else 0.)
        self.coupled_Micro_in_Macro_massconcentration = (
            self.coupled_SPM[self._micro_in_macro_idx].massconcentration
            if self._micro_in_macro_idx is not None else 0.)

        # Limitation functions
        self.get_limNUT()
        self.limT = self.limT_growth_array[t_idx]

        # Calculate fnut from active currencies only
        active_quotas = []
        if self.N is not None:
            active_quotas.append(self.QN)
        if self.P is not None:
            active_quotas.append(self.QP)
        if self.Si is not None:
            active_quotas.append(self.QSi)
        self.fnut = min(active_quotas) if active_quotas else 1.0

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

    def get_sinks(self, t, t_idx=None):
        # SINKS
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

        if self.debug_mode and self.C < 1e-3:
            print(f'DEBUG C={self.C:.2e} | source_PP={self.source_PP.C:.2e} | sinks={self.C_sinks:.2e}')
            print(f'  resp={self.sink_respiration.C:.2e}, exud={self.sink_exudation.C:.2e}, '
                  f'agg={self.sink_aggregation.C:.2e}, ing={self.sink_ingestion.C:.2e}, '
                  f'lys={self.sink_lysis.C:.2e}, mort={self.sink_mortality.C:.2e}')

        return np.array(
            [sinks for sinks in (self.C_sinks, self.N_sinks, self.Chl_sinks, self.P_sinks, self.Si_sinks) if
             sinks is not None], dtype=self.dtype)

    def get_diagnostic_variables(self):
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics], dtype=self.dtype)

    def get_limNUT(self):
        """Calculate nutrient limitations for Kerimoglu22 formulation."""
        active_lims = []
        _clip = np.clip if self.apply_numerical_protections else lambda x, a, b: x

        if self.N is not None:
            QN_safe = max(self.QN, self.QN_min * 1.001) if self.apply_numerical_protections else self.QN
            self.lim_N = _clip(1 - self.QN_min / QN_safe, 0., 1.)
            active_lims.append(self.lim_N)
            self.limQUOTAmin.N = _clip((self.QN - self.QN_min) / (self.QN_max - self.QN_min), 0., 1.)
            self.limQUOTA.N = 1 - self.limQUOTAmin.N
        else:
            self.lim_N = 1.0

        if self.P is not None:
            QP_safe = max(self.QP, self.QP_min * 1.001) if self.apply_numerical_protections else self.QP
            self.lim_P = _clip(1 - self.QP_min / QP_safe, 0., 1.)
            active_lims.append(self.lim_P)
            self.limQUOTAmin.P = _clip((self.QP - self.QP_min) / (self.QP_max - self.QP_min), 0., 1.)
            self.limQUOTA.P = 1 - self.limQUOTAmin.P
        else:
            self.lim_P = 1.0

        if self.Si is not None:
            QSi_safe = max(self.QSi, self.QSi_min * 1.001) if self.apply_numerical_protections else self.QSi
            self.lim_Si = _clip(1 - self.QSi_min / QSi_safe, 0., 1.)
            active_lims.append(self.lim_Si)
            self.limQUOTAmin.Si = _clip((self.QSi - self.QSi_min) / (self.QSi_max - self.QSi_min), 0., 1.)
            self.limQUOTA.Si = 1 - self.limQUOTAmin.Si
        else:
            self.lim_Si = 1.0

        # limNUT = minimum of active currency limitations only
        self.limNUT = max(0., min(active_lims)) if active_lims else 1.0

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
        """Calculate primary production for Kerimoglu22 formulation.

        Light limitation is integrated vertically over the sunlit upper layer
        [0, z_upper_layer]; rho_Chl (Chl synthesis regulation) is derived from the
        same integral so the two stay mutually consistent. The depth-mean P-I limitation
        is weighted by an exponential biomass profile B(z)=B0*exp(-m*z) and diluted over
        the column by the biomass fraction; m=0 reduces exactly to the homogeneous case.
        Realised production may additionally carry compound_production_factor
        (production only).
        When there is no usable light (night, or null PC_max/kd) production and every
        light-related term collapse to exactly 0.
        """
        # Optimization: Use pre-computed arrays for faster access
        self.PAR_t = self.setup.PAR_array[t_idx]
        self.get_kd()

        self.PC_max = self.mu_max * self.limNUT * self.limT
        H = self.setup.water_depth_array[t_idx]       # total water depth
        self.I_H = self.PAR_t * np.exp(-self.kd * H)  # irradiance reaching the seabed (diagnostic)

        # Diagnostic: effective biomass-profile slope at this step.
        self.biomass_profile_slope_t = (self.biomass_profile_slope_array[t_idx]
                                        if self.biomass_profile_slope_array is not None else 0.0)

        if self.PC_max > 0 and self.PAR_t > self.I_min and self.kd > 0:
            # The depth-mean P-I limitation is weighted by an exponential biomass
            # profile B(z)=B0*exp(-m*z) (m>=0; m=0 = homogeneous). The weighted mean
            # of 1-exp(-a*I(z)) over [0, z_upper_layer], with I(z)=I0*exp(-kd*z),
            # closes in the upper incomplete gamma function Γ(s,x), s=m/kd
            # (s=0 -> Γ(0,x)=E1(x), the homogeneous case).
            a = self.alpha * self.thetaC / (self.PC_max * varinfos.molmass_C)
            I_0 = self.PAR_t                          # incident irradiance at the surface
            # I_base: irradiance at the base of the upper layer. I_min is a small
            # non-zero irradiance floor, NOT the photic-zone limit (far shallower given
            # the turbidity): it caps the integration depth so the 1/(kd*z) and exp1()
            # terms stay numerically stable in turbid/deep water.
            I_base = max(self.I_H, self.I_min)
            z_up = min(np.log(I_0 / I_base) / self.kd, H)
            m = self.biomass_profile_slope_array[t_idx]

            # Biomass-weighted mean light limitation over [0, z_up] (E1 fallback for
            # m <= _M_EPS (numerically homogeneous),
            # exact and numerically cleanest). The full-column dilution uses the biomass
            # fraction in the upper layer (NOT z_up/H, which under-estimates limI once m>0).
            if m <= _M_EPS:
                limI_up = 1.0 - (exp1(a * I_base) - exp1(a * I_0)) / (self.kd * z_up)
                frac = z_up / H
            else:
                s = m / self.kd
                # gammaincc is the *regularised* upper incomplete gamma (needs s>0), so
                # Γ(s,x) = gamma(s)*gammaincc(s,x); G = Γ(s,a*I_base) - Γ(s,a*I_0).
                G = gamma(s) * (gammaincc(s, a * I_base) - gammaincc(s, a * I_0))
                limI_up = 1.0 - s * (a * I_0) ** (-s) * G / (1.0 - np.exp(-m * z_up))
                frac = (1.0 - np.exp(-m * z_up)) / (1.0 - np.exp(-m * H))

            self.z_upper_layer = z_up
            self.limI_upper_layer = self._floor(limI_up)
            self.limI = self._floor(limI_up * frac)   # geometric, biomass-weighted

            # Homogeneous reference (forced m=0) for the limI/limI_theoretical ratio.
            limI_up_m0 = 1.0 - (exp1(a * I_base) - exp1(a * I_0)) / (self.kd * z_up)
            self.limI_theoretical = self._floor(limI_up_m0 * (z_up / H))

            # Realised production carries the compound factor (production ONLY), capped at
            # PC_max. PC_geom is the geometric rate feeding the Chl:C acclimation (rho_Chl),
            # so the compound factor does not contaminate Chl:C.
            PC_geom = self.PC_max * self.limI
            self.PC = self.PC_max * min(self.limI * self.compound_production_factor, 1.0)

            # rho_Chl (Chl synthesis regulation) coherent with the vertical integral, tied
            # to the GEOMETRIC field (PC_geom, NOT self.PC). Invert the limitation to the
            # effective irradiance Ī: 1 - exp(-a·Ī) = limI.
            # rho_Chl in [mgChl mmolN-1]; QN_max in [molN molC-1]; theta_max in
            # [mgChl mmolC-1], so theta_max/QN_max is in [mgChl mmolN-1].
            if self.limI > 1e-6:
                I_eff = -np.log(1.0 - self.limI) / a if self.limI < 1.0 else I_0
                thetaC_safe = max(self.thetaC, 1e-9) if self.apply_numerical_protections else self.thetaC
                denom = self.alpha * thetaC_safe * I_eff
                rho_Chl_raw = self.theta_max / self.QN_max * (PC_geom * varinfos.molmass_C / denom)
                self.rho_Chl = np.clip(rho_Chl_raw, 0., self.theta_max * 10) \
                    if self.apply_numerical_protections else rho_Chl_raw
            else:
                self.rho_Chl = 0.

        else:
            # No usable light: production and all light-related terms vanish. The
            # diagnostic ratios limI/limI_* become NaN here (0/0); plotting.py
            # sanitizes those to gaps.
            self.z_upper_layer = 0.0
            self.limI = 0.0
            self.limI_upper_layer = 0.0
            self.limI_theoretical = 0.0
            self.PC = 0.0
            self.rho_Chl = 0.0

        C_safe = max(self.C, self._C_MIN) if self.apply_numerical_protections else self.C
        self.source_PP.C = self.PC * C_safe

        if self.debug_mode and any(p is not None and p < 0 for p in (self.C, self.N, self.P, self.Si, self.Chl)):
            print(f'DEBUG: Negative pools - C={self.C}, N={self.N}, P={self.P}, Si={self.Si}, Chl={self.Chl}')
            print(f'  limI={self.limI}, PC_max={self.PC_max}, PC={self.PC}')



    def get_source_uptake(self):
        """Calculate nutrient uptake for Kerimoglu22 formulation."""
        if self.N is not None:
            wNH4 = self.coupled_NH4.concentration / self.KNH4
            wNO3 = self.coupled_NO3.concentration / self.KNO3
            self.mmNH4 = wNH4 / (1 + wNH4 + wNO3)
            self.mmNO3 = wNO3 / (1 + wNH4 + wNO3)
            self.source_uptake.NH4 = self.v_max_N * self.limQUOTA.N * self.mmNH4 * self.limT * self.C
            self.source_uptake.NO3 = self.v_max_N * self.limQUOTA.N * self.mmNO3 * self.limT * self.C
            self.source_uptake.N = self.source_uptake.NH4 + self.source_uptake.NO3
        if self.P is not None:
            self.mmDIP = self.coupled_DIP.concentration / (self.coupled_DIP.concentration + self.KDIP)
            self.source_uptake.P = self.v_max_P * self.limQUOTA.P * self.mmDIP * self.limT * self.C
        if self.Si is not None:
            self.mmDSi = self.coupled_DSi.concentration / (self.coupled_DSi.concentration + self.KDSi)
            self.source_uptake.Si = self.v_max_Si * self.limQUOTA.Si * self.mmDSi * self.limT * self.C

    def get_source_Chlprod(self):
        """
        Calculate chlorophyll production for Kerimoglu22 formulation.
        dPhyChl = VN*rho_chl (without self.C)
        rho_Chl is in [mgChl mmolN-1], source_uptake.N is in [mmolN m-3 d-1]
        Hence [mgChl m-3 d-1] = [mgChl mmolN-1] * [mmolN m-3 d-1]
        """
        if self.N is not None:
            self.source_Chlprod.Chl = self.rho_Chl * self.source_uptake.N
        else:
            self.source_Chlprod.Chl = 0.

    def get_sink_lysis(self):
        flysis = 0.05 + 0.95 * (1. + np.exp(-(self.fnut - 0.2) * 30)) ** (-1)
        self.sink_lysis.C = self.C * flysis * self.lysrate
        if self.N is not None:
            self.sink_lysis.N = self.sink_lysis.C * self.QN
        self.sink_lysis.Chl = self.sink_lysis.C * self.thetaC
        if self.P is not None:
            self.sink_lysis.P = self.sink_lysis.C * self.QP
        if self.Si is not None:
            self.sink_lysis.Si = self.sink_lysis.C * self.QSi

    def get_sink_mortality(self, t, t_idx=None):
        # Optimization: Use pre-computed temperature limitation for grazing
        n = self.grazing_holling_exponent
        grazing_loss = (self.grazing_loss_max *
                        self.limT_grazing_array[t_idx] *
                        self.C**n / (self.K_grazing**n + self.C**n))
        self.sink_mortality.C = self.C * self.mortrate + grazing_loss
        if self.N is not None:
            self.sink_mortality.N = self.sink_mortality.C * self.QN
        self.sink_mortality.Chl = self.sink_mortality.C * self.thetaC
        if self.P is not None:
            self.sink_mortality.P = self.sink_mortality.C * self.QP
        if self.Si is not None:
            self.sink_mortality.Si = self.sink_mortality.C * self.QSi

        # DIAGNOSTIC-ONLY split of the implicit grazing-closure flux.
        # grazing_loss is *kept inside* sink_mortality.C above, hence still routed to DetS
        # via DetS.coupled_mortality_sources (routing and mass balance unchanged). sink_grazing
        # is a passive mirror so the grazing term can be isolated in the analysis
        # (Phy_sink_grazing.*). It is therefore deliberately NOT added to the *_sinks sums in
        # get_sinks() — doing so would double-count the loss.
        #
        # TODO (next optimization round, when backward-compat with running/previous sims is no
        #   longer required): promote this to a first-class flux — move grazing_loss OUT of
        #   get_sink_mortality into a dedicated get_sink_grazing(t, t_idx), add sink_grazing to
        #   the *_sinks sums in get_sinks(), and add the symmetric get_source_grazing() in
        #   detritus.Detritus fed by a new DetS 'coupled_grazing_sources': ['Phy'] in
        #   base_config (plus matching config_diagnostics entries). Mass-neutral, but it changes
        #   the model/config structure, hence deferred.
        self.sink_grazing.C = grazing_loss
        if self.N is not None:
            self.sink_grazing.N = grazing_loss * self.QN
        if self.P is not None:
            self.sink_grazing.P = grazing_loss * self.QP
        if self.Si is not None:
            self.sink_grazing.Si = grazing_loss * self.QSi

    def get_sink_exudation(self):
        """Calculate exudation for Kerimoglu22 formulation."""
        C_safe = max(self.C, self._C_MIN) if self.apply_numerical_protections else self.C
        self.sink_exudation.C = self.gamma_C_exud_base * C_safe + self.gamma_C_exud_prod * self.source_PP.C
        self.sink_exudation.Chl = 0.
        xquota = 0.99
        if self.N is not None:
            self.sink_exudation.N = self.sink_exudation.C * self.QN if self.limQUOTAmin.N > xquota else 0.
        if self.P is not None:
            self.sink_exudation.P = self.sink_exudation.C * self.QP if self.limQUOTAmin.P > xquota else 0.
        if self.Si is not None:
            self.sink_exudation.Si = self.sink_exudation.C * self.QSi if self.limQUOTAmin.Si > xquota else 0.
        self.frac_exud_small = 1 / (1 + np.exp(-(self.fnut - 0.2) * 30))

    def get_sink_respiration(self):
        """Calculate respiration for Kerimoglu22 formulation."""
        C_safe = max(self.C, self._C_MIN) if self.apply_numerical_protections else self.C
        self.sink_respiration.C = self.zeta_resp_base * C_safe + self.zeta_resp_prod * self.source_PP.C
        self.sink_respiration.Chl = self.sink_respiration.C * self.thetaC  # Chl degradation
        self.sink_respiration.N = 0.
        self.sink_respiration.P = 0.
        self.sink_respiration.Si = 0.

    def get_sink_ingestion(self):
        if self.coupled_consumer is not None:
            self.sink_ingestion.C = fns.get_all_contributors(self.coupled_consumer, 'source_ingestion', 'C')
            if self.N is not None:
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
        """Calculate aggregation for Kerimoglu22 formulation."""
        self.sink_aggregation.C = self.coupled_aggreg_target.aggPhy_C
        if self.N is not None:
            self.sink_aggregation.N = self.coupled_aggreg_target.aggPhy_N
        if self.P is not None:
            self.sink_aggregation.P = self.coupled_aggreg_target.aggPhy_P
        if self.Si is not None:
            self.sink_aggregation.Si = self.coupled_aggreg_target.aggPhy_Si
        self.sink_aggregation.Chl = self.coupled_aggreg_target.aggPhy_Chl
