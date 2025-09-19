import numpy as np

from ..core.base import BaseStateVar
from src.config_model import varinfos
from ..utils import functions as fns


class Flocs(BaseStateVar):
    def __init__(self,
                 name,
                 alpha_PP=0.2,
                 alpha_PF=0.2,
                 alpha_FF=0.3,
                 p_exp=0.4,
                 q_exp=0.1,
                 f_frac_floc_break=0.1,
                 efficiency_break=2e-4,
                 fyflocstrength=1e-10,
                 mu_viscosity=1e-6,

                 d_p_microflocdiam=18e-6,
                 nf_fractal_dim=2,

                 density=2500,  # [kg/m3]

                 # Base values for additive TEP formulation - values in the absence of organic TEP (= purely mineral)
                 alpha_FF_base = 0.02,     # [-] Base FF collision efficiency (mineral only)
                 alpha_PP_base = 0.02,     # [-] Base PP collision efficiency (mineral only)
                 alpha_PF_base = 0.02,     # [-] Base PF collision efficiency (mineral only)
                 fyflocstrength_base = 1e-10,  # [N] Base floc strength (mineral only)
                 tau_cr_base = 0.1,        # [Pa] Base critical shear stress (mineral only)

                 # Delta values for TEP effect (additive increments)
                 delta_alpha_FF = 0.03,    # [-] TEP increment for FF collision efficiency
                 delta_alpha_PP = 0.03,    # [-] TEP increment for PP collision efficiency
                 delta_alpha_PF = 0.03,    # [-] TEP increment for PF collision efficiency
                 deltaFymax = 1e-9,        # [N] TEP increment for floc strength
                 delta_tau_cr = 0.2,       # [Pa] TEP increment for critical shear stress

                 # TEP coupling parameters
                 K_glue = None,            # [mmol m-3] Half-saturation for TEP effect

                 # Legacy parameters (for backward compatibility, will be deprecated)
                 alpha_FF_ref = None,
                 alpha_PP_ref = None,
                 alpha_PF_ref = None,
                 delta_tau = 2.,       # [-] Legacy parameter

                 #
                 resuspension_rate = 0.,       # [kg/m²/s/Pa]
                 settling_vel_min_fraction = 0.1,  # [-] Minimum fraction of settling velocity retained at max shear
                 settling_vel_max_fraction = 1.0,  # [-] Maximum fraction of settling velocity retained at min shear
                 apply_settling = True, # Boolean. Whether to apply sediment settling
                 counter_settling_by_turbulence = False, # Boolean, whether to account for settling by turbulent water

                 # eps_kd=2e-5 * varinfos.molmass_C,  # [m2 mgC-1] Diffuse attenuation cross section
                 # # value from
                 # eps_kd=2e-5,  # [m2 mg-1] Diffuse attenuation cross section of SPM with kd = ... + eps_kd * sqrt(SPM)
                 # # value from Tian et al., 2009
                 eps_kd = 0.066* 1e3, # from 0.066 [m-1 (mg l-1)-1] to [m-1 (g l-1)-1] = [m-1 (kg m-3)-1] (model units) Light attenuation due to SPM
                 # model is in kg m-3

                 sinking_leak=0,

                 dt2=True,
                 time_conversion_factor = 86400,   # Flocs run in s-1 while the rest of the model is in d-1
                 spinup_days = 0,   # Days of spin-up before interactions with biological components
                 ):

        super().__init__()

        # Additive formulation parameters
        self.alpha_FF_base = alpha_FF_base
        self.alpha_PP_base = alpha_PP_base
        self.alpha_PF_base = alpha_PF_base
        self.fyflocstrength_base = fyflocstrength_base
        self.tau_cr_base = tau_cr_base

        self.delta_alpha_FF = delta_alpha_FF
        self.delta_alpha_PP = delta_alpha_PP
        self.delta_alpha_PF = delta_alpha_PF
        self.deltaFymax = deltaFymax
        self.delta_tau_cr = delta_tau_cr

        self.K_glue = K_glue

        # Legacy parameters (for backward compatibility)
        self.alpha_FF_ref = alpha_FF_ref
        self.alpha_PP_ref = alpha_PP_ref
        self.alpha_PF_ref = alpha_PF_ref
        self.delta_tau = delta_tau

        self.SMS = None
        self.settling_vel = None
        self.settling_vel_base = None
        self.breaksource = None
        self.ffloss = None
        self.ppsource = None
        self.settling_loss = None
        self.g_shear_rate_at_t = None
        self.bed_shear_stress_at_t = None
        self.erosion_factor = None
        self.Ncnum = None
        self.coupled_glue = None
        self.coupled_Nt = None
        self.coupled_Nf = None
        self.coupled_Np = None
        self.sink_breakdown = None
        self.sink_aggregation = None
        self.source_breakdown = None
        self.ICs = None
        self.source_aggregation = None
        self.numconc = None
        self.massconcentration = None
        self.volconcentration = None
        self.sedimentation = None
        self.resuspension = None
        self.diagnostics = None
        self.classname = 'Floc'
        self.name = name
        # Initialize alphas with base values (will be updated by TEP coupling if active)
        self.alpha_PP = alpha_PP_base
        self.alpha_PF = alpha_PF_base
        self.alpha_FF = alpha_FF_base
        self.p_exp = p_exp
        self.q_exp = q_exp
        self.f_frac_floc_break = f_frac_floc_break
        self.efficiency_break = efficiency_break
        # Initialize with base value (will be updated by TEP coupling if active)
        self.fyflocstrength = fyflocstrength_base
        self.mu_viscosity = mu_viscosity
        self.sinking_leak = sinking_leak

        self.resuspension_rate = resuspension_rate
        self.settling_vel_min_fraction = settling_vel_min_fraction
        self.settling_vel_max_fraction = settling_vel_max_fraction
        self.counter_settling_by_turbulence = counter_settling_by_turbulence
        self.apply_settling = apply_settling

        # Initialize tau_cr with base value (will be updated by TEP coupling if active)
        self.tau_cr = tau_cr_base

        self.d_p_microflocdiam = d_p_microflocdiam
        self.diam = d_p_microflocdiam
        self.nf_fractal_dim = nf_fractal_dim
        self.density = density
        self.eps_kd = eps_kd
        self.dt2 = dt2
        self.time_conversion_factor = time_conversion_factor
        self.spinup_days = spinup_days

        self.setup = None


    def _calculate_macrofloc_diameter(self):
        """
        Calculate macrofloc diameter using Lee et al. (2011) fractal theory.
        DF = NC^(1/nf) * DP where NC is microflocs per macrofloc

        Returns:
            float: Calculated macrofloc diameter
        """
        if self.name != 'Macroflocs':
            return self.diam

        # Calculate NC (microflocs per macrofloc)
        if hasattr(self, 'coupled_Nt') and hasattr(self, 'coupled_Nf') and \
                self.coupled_Nt is not None and self.coupled_Nf is not None:

            # Guard against division by zero
            if self.coupled_Nf.numconc > 0:
                nc = self.coupled_Nt.numconc / self.coupled_Nf.numconc
            else:
                nc = 1.0  # Default to 1 if no macroflocs present

            # Calculate diameter based on fractal theory
            macrofloc_diam = nc ** (1 / self.nf_fractal_dim) * self.coupled_Np.diam

            # Ensure minimum size is at least double the microfloc diameter
            min_diam = 2 * self.coupled_Np.diam
            return max(macrofloc_diam, min_diam)

        # If couplings aren't established yet, return the default diameter
        return self.diam

    def set_ICs(self,
                numconc
                ):
        self.numconc = numconc
        self.ICs = [self.numconc]

        if self.name == "Microflocs" or self.name == "Micro_in_Macro":
            self.massconcentration = numconc * np.pi / 6. * self.diam*self.diam*self.diam * self.density

        # Source and sink terms
        self.source_aggregation = None
        self.source_breakdown = None
        self.sink_aggregation = None
        self.sink_breakdown = None

    def set_coupling(self,
                     coupled_Np=None,
                     coupled_Nf=None,
                     coupled_Nt=None,
                     coupled_glue=None,
                     ):
        self.coupled_Np = coupled_Np if self.name != 'Microflocs' else self
        self.coupled_Nf = coupled_Nf if self.name != 'Macroflocs' else self
        self.coupled_Nt = coupled_Nt if self.name != 'Micro_in_Macro' else self
        self.coupled_glue = coupled_glue

        if self.name != 'Microflocs':
            self.nf_fractal_dim = self.coupled_Np.nf_fractal_dim
            # Initialize with base values (will be updated by TEP coupling if active)
            self.alpha_PP = self.coupled_Np.alpha_PP_base
            self.alpha_PF = self.coupled_Np.alpha_PF_base
            self.alpha_FF = self.coupled_Np.alpha_FF_base
            self.f_frac_floc_break = self.coupled_Np.f_frac_floc_break
            self.efficiency_break = self.coupled_Np.efficiency_break
            # Initialize with base value (will be updated by TEP coupling if active)
            self.fyflocstrength = self.coupled_Np.fyflocstrength_base

            self.eps_kd = self.coupled_Np.eps_kd

            # Copy new additive formulation parameters
            self.K_glue = self.coupled_Np.K_glue
            self.alpha_FF_base = self.coupled_Np.alpha_FF_base
            self.alpha_PP_base = self.coupled_Np.alpha_PP_base
            self.alpha_PF_base = self.coupled_Np.alpha_PF_base
            self.delta_alpha_FF = self.coupled_Np.delta_alpha_FF
            self.delta_alpha_PP = self.coupled_Np.delta_alpha_PP
            self.delta_alpha_PF = self.coupled_Np.delta_alpha_PF
            self.fyflocstrength_base = self.coupled_Np.fyflocstrength_base
            self.deltaFymax = self.coupled_Np.deltaFymax
            self.tau_cr_base = self.coupled_Np.tau_cr_base
            self.delta_tau_cr = self.coupled_Np.delta_tau_cr

            # Copy legacy parameters for backward compatibility
            self.alpha_FF_ref = self.coupled_Np.alpha_FF_ref
            self.alpha_PP_ref = self.coupled_Np.alpha_PP_ref
            self.alpha_PF_ref = self.coupled_Np.alpha_PF_ref
            self.delta_tau = self.coupled_Np.delta_tau

            self.spinup_days = self.coupled_Np.spinup_days

        if self.name == 'Micro_in_Macro':
            # 20250910 - This should be obsolete and deleted (only computed for Macroflocs, i.e. coupled_Nf)
            self.sinking_leak = self.coupled_Nf.sinking_leak
            self.resuspension_rate = self.coupled_Nf.resuspension_rate
            self.tau_cr_base = self.coupled_Nf.tau_cr_base
            self.delta_tau_cr = self.coupled_Nf.delta_tau_cr
            # Legacy parameter
            self.delta_tau = self.coupled_Nf.delta_tau
            self.apply_settling = self.coupled_Nf.apply_settling

        # Initialize macrofloc diameter if couplings are now established
        if self.name == 'Macroflocs':
            self.diam = self._calculate_macrofloc_diameter()


    def update_val(self, numconc,
                   t=None,
                   debugverbose=False):
        self.numconc = numconc

        # For macroflocs, recalculate diameter based on fractal theory
        if self.name == 'Macroflocs':
            self.diam = self._calculate_macrofloc_diameter()

        self.volconcentration = numconc * np.pi / 6. * self.diam * self.diam * self.diam
        
        if self.name == "Microflocs" or self.name == "Micro_in_Macro":
            self.massconcentration = numconc * np.pi / 6. * self.diam*self.diam*self.diam * self.density


    def get_diagnostic_variables(self):
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics])


    def get_sources(self, t=None):
        """
        Note: the "FABM" and "get_sources/get_sinks" approach is not adopted here - for simplicity and also
        with some mathematical simplifications within the SMS equations to enhance computational efficiency
        (since the floc module is run at smaller time step)
        """

        self.Ncnum = self.coupled_Nt.numconc / self.coupled_Nf.numconc
        self.g_shear_rate_at_t = self.setup.g_shear_rate.loc[t]['ShearRate']
        self.bed_shear_stress_at_t = self.setup.bed_shear_stress.loc[t]['BedShearStress']
        self.water_depth_at_t = self.setup.water_depth.loc[t]['WaterDepth']

        # Check if we're in spin-up phase
        is_spinup = hasattr(self.setup, 'in_spinup_phase') and self.setup.in_spinup_phase

        if is_spinup:
            # Spin-up mode: disable TEP coupling and settling/resuspension (only mineral flocculation dynamics)
            use_tep_coupling = False

            # Store original values and disable settling
            original_apply_settling = self.apply_settling
            original_resuspension_rate = self.resuspension_rate
            self.apply_settling = False
            self.resuspension_rate = 0.0
        else:
            # Normal mode: use TEP coupling if available
            use_tep_coupling = self.coupled_glue and self.K_glue
            # Initialize variables to avoid reference errors
            original_apply_settling = None
            original_resuspension_rate = None

        if use_tep_coupling:
            mm_TEP = self.coupled_glue.C / (self.K_glue + self.coupled_glue.C)

            # New unified additive formulation: variable = base + delta * mm_TEP
            self.alpha_FF = self.alpha_FF_base + self.delta_alpha_FF * mm_TEP
            self.alpha_PP = self.alpha_PP_base + self.delta_alpha_PP * mm_TEP
            self.alpha_PF = self.alpha_PF_base + self.delta_alpha_PF * mm_TEP
            self.fyflocstrength = self.fyflocstrength_base + self.deltaFymax * mm_TEP

            # Apply tau_cr modification only for settling flocs
            if self.name == 'Macroflocs' or self.name == 'Micro_in_Macro':
                self.tau_cr = self.tau_cr_base + self.delta_tau_cr * mm_TEP
        else:
            # No TEP coupling: use base values
            self.alpha_FF = self.alpha_FF_base
            self.alpha_PP = self.alpha_PP_base
            self.alpha_PF = self.alpha_PF_base
            self.fyflocstrength = self.fyflocstrength_base
            if self.name == 'Macroflocs' or self.name == 'Micro_in_Macro':
                self.tau_cr = self.tau_cr_base



        if self.name == 'Microflocs':
            self.SMS = (
                        # Loss from microfloc-microfloc collisions forming macroflocs
                        -2. / 3. * self.alpha_PP * self.coupled_Np.diam * self.coupled_Np.diam * self.coupled_Np.diam *
                        self.g_shear_rate_at_t * self.coupled_Np.numconc * self.coupled_Np.numconc * (
                                    self.Ncnum / (self.Ncnum - 1.)) -

                        # Loss from microfloc-macrofloc collisions
                        1. / 6. * self.alpha_PF * (
                                    self.Ncnum ** (1. / self.nf_fractal_dim) + 1.) ** 3. *
                        self.coupled_Np.diam * self.coupled_Np.diam * self.coupled_Np.diam *
                        self.g_shear_rate_at_t * self.coupled_Np.numconc * self.coupled_Nf.numconc +

                        # Gain from macrofloc breakage
                        self.f_frac_floc_break * self.Ncnum * self.efficiency_break * self.g_shear_rate_at_t *
                        (self.Ncnum ** (1. / self.nf_fractal_dim) - 1.) ** self.p_exp *
                        (self.mu_viscosity * self.g_shear_rate_at_t * self.Ncnum ** (2. / self.nf_fractal_dim) *
                         self.coupled_Np.diam * self.coupled_Np.diam / self.fyflocstrength) ** self.q_exp * self.coupled_Nf.numconc)

        elif self.name == 'Macroflocs':
            # Winterwerp formula: W_s,F = (1/18) × (ρ_s - ρ_w) / μ × g × D_p^(3-nf) × D_F^(nf-1)
            self.settling_vel_base = (1/18) * (self.density - self.setup.rho_water) / self.setup.mu_water * 9.81 * (
                    self.d_p_microflocdiam ** (3 - self.nf_fractal_dim)) * (self.diam ** (self.nf_fractal_dim - 1))

            self.settling_vel_base = self.settling_vel_base * self.apply_settling

            if self.counter_settling_by_turbulence:
                # # Apply shear-dependent modulation
                normalized_shear = (self.g_shear_rate_at_t - self.setup.g_shear_rate_min) / (
                    self.setup.delta_g_shear_rate)
                shear_factor = self.settling_vel_min_fraction + (self.settling_vel_max_fraction - self.settling_vel_min_fraction) * 0.5 * (
                            1 + np.cos(normalized_shear * np.pi))
                self.settling_vel = self.settling_vel_base * shear_factor
            else:
                self.settling_vel = self.settling_vel_base
            # print('BLA base vs countered', self.diam, self.settling_vel_base, self.settling_vel, self.settling_vel/self.settling_vel_base)


            if self.resuspension_rate > 0:
                # Physical settling and resuspension
                self.sedimentation = self.settling_vel * self.numconc / self.water_depth_at_t

                self.erosion_factor = max(0, (self.bed_shear_stress_at_t / self.tau_cr - 1))
                self.resuspension = self.resuspension_rate * self.erosion_factor / self.water_depth_at_t

                self.settling_loss = self.sedimentation - self.resuspension
            else:
                # Fallback to original formulation
                self.settling_loss = self.sinking_leak * self.settling_vel * self.numconc


            # Formation from microfloc-microfloc collisions
            self.ppsource = (2. / 3. * self.alpha_PP * self.coupled_Np.diam * self.coupled_Np.diam * self.coupled_Np.diam *
                        self.g_shear_rate_at_t * self.coupled_Np.numconc * self.coupled_Np.numconc * (
                                    1. / (self.Ncnum - 1.)))
            # Loss from macrofloc-macrofloc collisions
            self.ffloss = (2. / 3. * self.alpha_FF * self.Ncnum ** (3. / self.nf_fractal_dim) *
                        self.coupled_Np.diam * self.coupled_Np.diam * self.coupled_Np.diam *
                        self.g_shear_rate_at_t * self.coupled_Nf.numconc * self.coupled_Nf.numconc)
            # Source from breakage (increases macroflocs numbers)
            self.breaksource = (self.efficiency_break * self.g_shear_rate_at_t * (
                                self.Ncnum ** (1. / self.nf_fractal_dim) - 1.) ** self.p_exp *
                        (self.mu_viscosity * self.g_shear_rate_at_t * self.Ncnum ** (2. / self.nf_fractal_dim) *
                         self.coupled_Np.diam * self.coupled_Np.diam / self.fyflocstrength) ** self.q_exp * self.coupled_Nf.numconc)

            self.SMS = (self.ppsource -
                        self.ffloss +
                        self.breaksource -
                        self.settling_loss
                        )

        else:

            self.settling_loss = self.coupled_Nf.settling_loss * self.Ncnum

            self.SMS = (
                        # Gain from microfloc-microfloc collisions
                        2. / 3. * self.alpha_PP * self.coupled_Np.diam * self.coupled_Np.diam * self.coupled_Np.diam *
                        self.g_shear_rate_at_t * self.coupled_Np.numconc * self.coupled_Np.numconc * (
                                    self.Ncnum / (self.Ncnum - 1.)) +

                        # Gain from microfloc-macrofloc collisions
                        1. / 6. * self.alpha_PF * (
                                    self.Ncnum ** (1. / self.nf_fractal_dim) + 1.) ** 3 *
                        self.coupled_Np.diam * self.coupled_Np.diam * self.coupled_Np.diam *
                        self.g_shear_rate_at_t * self.coupled_Np.numconc * self.coupled_Nf.numconc -

                        # Loss from macrofloc breakage
                        self.f_frac_floc_break * self.Ncnum * self.efficiency_break * self.g_shear_rate_at_t *
                        (self.Ncnum ** (1. / self.nf_fractal_dim) - 1.) ** self.p_exp *
                        (self.mu_viscosity * self.g_shear_rate_at_t * self.Ncnum ** (2. / self.nf_fractal_dim) *
                         self.coupled_Np.diam * self.coupled_Np.diam / self.fyflocstrength) ** self.q_exp * self.coupled_Nf.numconc -

                        # Loss from sinking
                        self.settling_loss

                        )

        # Restore original values if we were in spin-up mode
        if hasattr(self.setup, 'in_spinup_phase') and self.setup.in_spinup_phase:
            self.apply_settling = original_apply_settling
            self.resuspension_rate = original_resuspension_rate

        return np.array(self.SMS)


    def get_sinks(self, t=None):
        return np.array([0])

    """
    def get_beta(self, i, j):
        return 1 / 6 * (i.diam + j.diam) ** 3 * self.g_shear_rate_at_t

    def get_a_breakage(self):
        return (self.efficiency_break * self.g_shear_rate_at_t *
                ((self.coupled_Nf.diam - self.coupled_Np.diam) / self.coupled_Np.diam) ** self.p_exp *
                (
                        self.mu_viscosity * self.g_shear_rate_at_t * self.coupled_Nf.numconc ** 2 / self.fyflocstrength) ** self.q_exp)

    def get_source_breakdown(self):
        if self.name == 'Microflocs':
            self.source_breakdown = self.macro_breakage * self.fyflocstrength * self.Ncnum
        elif self.name == 'Macroflocs':
            self.source_breakdown = self.macro_breakage
        elif self.name == 'Micro_in_Macro':
            self.source_breakdown = 0.

    def get_source_aggregation(self):
        if self.name == 'Microflocs':
            self.source_aggregation = 0
        elif self.name == 'Macroflocs':
            self.source_aggregation = self.microflocs_collision * 0.5 * (1 / (self.Ncnum - 1))
        elif self.name == 'Micro_in_Macro':
            self.source_aggregation = (self.microflocs_collision * 0.5 * (self.Ncnum / (self.Ncnum - 1)) +
                                       self.micromacro_collision)

    def get_sink_breakdown(self):
        if self.name == 'Microflocs' or self.name == 'Macroflocs':
            self.sink_breakdown = 0.
        elif self.name == 'Micro_in_Macro':
            self.sink_breakdown = self.macro_breakage * self.fyflocstrength * self.Ncnum

    def get_sink_aggregation(self):
        if self.name == 'Microflocs':
            self.sink_aggregation = (self.microflocs_collision * 0.5 * (self.Ncnum / (self.Ncnum - 1)) +
                                     self.micromacro_collision)
        elif self.name == 'Macroflocs':
            self.sink_aggregation = self.macroflocs_collision * 0.5
        elif self.name == 'Micro_in_Macro':
            self.sink_aggregation = 0.
    """

