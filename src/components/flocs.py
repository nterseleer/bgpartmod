import numpy as np

from ..core.base import BaseStateVar
from ..utils import functions as fns


class PrescribedTEP:
    """
    Lightweight class to hold prescribed TEP concentration for Flocs-only optimization.
    Mimics the interface of a full TEP state variable but only stores C concentration.
    The value of C is updated at each timestep from Setup.TEP_array.
    """
    def __init__(self):
        self.C = 0.0  # Will be updated at each timestep


class PrescribedFlocs:
    """
    Lightweight class to hold prescribed Flocs data for BGC-only optimization.
    Mimics the interface of full Flocs instances but loads data from Setup.
    Updated at each timestep from Setup flocs arrays.

    Provides:
    - massconcentration (Microflocs, Macroflocs) for light attenuation
    - sink_sedimentation, source_resuspension, numconc for vertical coupling
    """
    def __init__(self, name, eps_kd=0.066e3, time_conversion_factor=86400,
                 resusp_ewma_alpha=0.0, organomin_coupling_fraction=1.0):
        self.name = name
        self.eps_kd = eps_kd
        self.time_conversion_factor = time_conversion_factor
        self.resusp_ewma_alpha = resusp_ewma_alpha
        self.organomin_coupling_fraction = organomin_coupling_fraction

        # Light attenuation (for all floc types)
        self.massconcentration = 0.0  # [kg m⁻³] - updated at each timestep

        # Vertical coupling (Macroflocs only)
        if name == "Macroflocs":
            self.sink_sedimentation = 0.0  # [# m⁻³ s⁻¹] - updated at each timestep
            self.source_resuspension = 0.0  # [# m⁻³ s⁻¹] - updated at each timestep
            self.numconc = 0.0  # [# m⁻³] - updated at each timestep


class SharedFlocTEPParameters:
    """
    Centralized computation of TEP-dependent parameters for flocculation.
    Shared by all Flocs instances (Microflocs, Macroflocs, Micro_in_Macro) to avoid
    redundant calculations and ensure consistency.
    Includes: alphas (FF, PP, PF), fyflocstrength, tau_cr, nf_fractal_dim
    """

    def __init__(self, microfloc_instance, unified_alphas=True):
        """
        Initialize shared alpha parameters from the first Flocs instance (Microflocs).

        Args:
            microfloc_instance: Microflocs instance to extract configuration from
            unified_alphas: If True, use single alpha for FF=PP=PF (optimized)
        """
        self.unified_alphas = unified_alphas

        # Extract configuration from Microflocs instance
        self.alpha_PP_base = microfloc_instance.alpha_PP_base
        self.delta_alpha_PP = microfloc_instance.delta_alpha_PP
        if self.unified_alphas:
            self.alpha_FF_base = self.alpha_PP_base
            self.alpha_PF_base = self.alpha_PP_base
            self.delta_alpha_FF = self.delta_alpha_PP
            self.delta_alpha_PF = self.delta_alpha_PP
        else:
            self.alpha_FF_base = microfloc_instance.alpha_FF_base
            self.alpha_PF_base = microfloc_instance.alpha_PF_base
            self.delta_alpha_FF = microfloc_instance.delta_alpha_FF
            self.delta_alpha_PF = microfloc_instance.delta_alpha_PF

        self.fyflocstrength_base = microfloc_instance.fyflocstrength_base
        self.deltaFymax = microfloc_instance.deltaFymax

        self.tau_cr_base = microfloc_instance.tau_cr_base
        self.delta_tau_cr = microfloc_instance.delta_tau_cr

        self.nf_fractal_dim_base = microfloc_instance.nf_fractal_dim
        self.delta_nf_fractal_dim = microfloc_instance.delta_nf_fractal_dim

        self.K_glue = getattr(microfloc_instance, 'K_glue', None)

    def get_tep_parameters(self, coupled_glue):
        """
        Calculate all TEP-dependent parameters in one go.
        Returns: (alpha_FF, alpha_PP, alpha_PF), fyflocstrength, tau_cr, nf_fractal_dim, mm_TEP
        """
        # Determine TEP coupling status and calculate mm_TEP
        use_tep_coupling = coupled_glue is not None and self.K_glue is not None

        if use_tep_coupling:
            mm_TEP = coupled_glue.C / (self.K_glue + coupled_glue.C)
        else:
            mm_TEP = 0.0

        # Calculate alphas
        if self.unified_alphas:
            alpha = self.alpha_PP_base + (self.delta_alpha_PP * mm_TEP if use_tep_coupling else 0)
            alphas = (alpha, alpha, alpha)
        else:
            alphas = (
                self.alpha_FF_base + (self.delta_alpha_FF * mm_TEP if use_tep_coupling else 0),
                self.alpha_PP_base + (self.delta_alpha_PP * mm_TEP if use_tep_coupling else 0),
                self.alpha_PF_base + (self.delta_alpha_PF * mm_TEP if use_tep_coupling else 0)
            )

        # Calculate fyflocstrength
        fyflocstrength = self.fyflocstrength_base + (self.deltaFymax * mm_TEP if use_tep_coupling else 0)

        # Calculate tau_cr
        tau_cr = self.tau_cr_base + (self.delta_tau_cr * mm_TEP if use_tep_coupling else 0)

        # Calculate nf_fractal_dim
        nf_fractal_dim = self.nf_fractal_dim_base + (self.delta_nf_fractal_dim * mm_TEP if use_tep_coupling else 0)

        # Store calculated values for reuse by other instances
        self.current_alphas = alphas
        self.current_fyflocstrength = fyflocstrength
        self.current_tau_cr = tau_cr
        self.current_nf_fractal_dim = nf_fractal_dim
        self.current_mm_TEP = mm_TEP

        return alphas, fyflocstrength, tau_cr, nf_fractal_dim, mm_TEP


class Flocs(BaseStateVar):
    # Class-level shared alpha instance for performance optimization
    # Recreated for each new simulation to ensure fresh parameters
    _shared_alphas_instance = None

    def __init__(self,
                 name,
                 p_exp=1.0,  # [-] Exponent of the breakage kinetics (Lee11)
                 q_exp=1.0,  # [-] Exponent of the breakage kinetics; = 3 - nf at nf = 2.0 (Lee11)
                 dynamic_q_exp=True,  # [-] If True, q = 3 - nf(t) as in Lee11/Lee14, instead of the fixed q_exp
                 d_crit_growth=None,   # [m] Critical diameter above which breakage ramps up (Lee11: 450e-6). None disables it
                 d_crit_exponent=10.,  # [-] Steepness of that ramp
                 d_crit_max_factor=100.,  # [-] Cap on the ramp, so an explicit Euler step cannot overshoot
                 f_frac_floc_break=0.1,  # [-] Fraction of microflocs released by breakage (Lee11)
                 efficiency_break=1.0e-4,  # [s^0.5 m-1] Efficiency factor for breakage (Lee11)

                 d_p_microflocdiam=18e-6,  # [m] Diameter of the flocculi (Lee11)
                 nf_fractal_dim=2.0,  # [-] Fractal dimension of the macroflocs (Lee11)

                 density=1600,  # [kg m-3] Density of the flocculi (Lee11)

                 # Base values for additive TEP formulation - values in the absence of organic TEP (= purely mineral)
                 alpha_FF_base = 0.10,     # [-] Base FF collision efficiency, mineral only (Lee11, TCPBE)
                 alpha_PP_base = 0.10,     # [-] Base PP collision efficiency, mineral only (Lee11, TCPBE)
                 alpha_PF_base = 0.10,     # [-] Base PF collision efficiency, mineral only (Lee11, TCPBE)
                 fyflocstrength_base = 1e-10,  # [N] Base floc yield strength, mineral only (Lee11; their
                                               # Table 3 reads [Pa], but the breakage kernel needs a force
                                               # for mu*G*D_F^2/F_y to be dimensionless)
                 tau_cr_base = 0.5,        # [Pa] Base critical shear stress for erosion, mineral only.
                                           # Not in Lee11: resuspension is specific to this model.

                 # Delta values for TEP effect (additive increments). None of these are in Lee11:
                 # the TEP coupling is specific to this model.
                 delta_alpha_FF = 0.03,    # [-] TEP increment for FF collision efficiency
                 delta_alpha_PP = 0.03,    # [-] TEP increment for PP collision efficiency
                 delta_alpha_PF = 0.03,    # [-] TEP increment for PF collision efficiency
                 deltaFymax = 1e-9,        # [N] TEP increment for floc strength

                 # Optimization parameter for alpha sharing
                 unified_alphas = True,    # [-] Use single alpha for FF=PP=PF (performance optimization)
                 include_PP_collision = True,  # [-] Include PP collision term (microfloc-microfloc)
                 delta_tau_cr = 0.2,       # [Pa] TEP increment for critical shear stress
                 delta_nf_fractal_dim = 0.0,  # [-] TEP increment for fractal dimension

                 # TEP coupling parameters
                 K_glue = None,            # [mmol m-3] Half-saturation for TEP effect
                 prescribe_tep_from_setup = False,  # [-] Use prescribed TEP from Setup instead of coupled_glue
                 #
                 resuspension_rate = 0.,   # [# m-2 s-1 Pa-1] Erosion rate constant (0 = no resuspension)
                 apply_settling = True, # Boolean. Whether to apply sediment settling
                 settling_velocity_factor = None,  # [-] Constant fraction of base settling velocity (None = w_s natif)

                 # Vertical coupling parameters for BGC components
                 resusp_ewma_alpha = 0.0,  # [-] α=0: fixed ratio (Option A), α>0: adaptive smoothing (Option B1)
                 organomin_coupling_fraction = 1.0,  # [-] Fraction forming organo-mineral aggregates

                 # Light attenuation by SPM: 0.066 [m-1 (mg l-1)-1] converted to model units
                 # [m-1 (kg m-3)-1]. Tian et al. (2009) give an alternative in sqrt(SPM).
                 eps_kd = 0.066 * 1e3,  # [m2 kg-1]

                 sinking_leak=0,

                 dt2=True,
                 dtype=np.float64,
                 time_conversion_factor = 86400,   # Flocs run in s-1 while the rest of the model is in d-1
                 spinup_days = 0,   # Days of spin-up before interactions with biological components
                 ):

        super().__init__(dtype=dtype)

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
        self.delta_nf_fractal_dim = delta_nf_fractal_dim

        self.K_glue = K_glue
        self.include_PP_collision = include_PP_collision
        self.prescribe_tep_from_setup = prescribe_tep_from_setup

        # SharedFlocTEPParameters instance (injected in set_coupling)
        self.shared_alphas = None
        self.unified_alphas = unified_alphas

        self.SMS = None
        self.settling_vel = None
        self.settling_vel_base = None
        self.breaksource = None
        self.ffloss = None
        self.ppsource = None
        self.settling_loss = None
        self.net_vertical_loss_rate = None  # [s-1] Net vertical loss rate for organic coupling
        self.g_shear_rate_at_t = None
        self.bed_shear_stress_at_t = None
        self.mu_water_at_t = None
        self.erosion_factor = None
        self.water_depth_at_t = None
        self.Ncnum = None
        self.coupled_glue = None
        self.coupled_glue_C = None
        self.coupled_Nt = None
        self.coupled_Nf = None
        self.coupled_Np = None
        self.sink_breakdown = None
        self.sink_aggregation = None
        self.source_breakdown = None
        self.ICs = None
        self.source_aggregation = None
        self.sink_breakage = None
        self.source_PF_collision = None
        self.source_PP_collision = None
        self.source_breakage = None
        self.sink_FF_collision = None
        self.sink_PF_collision = None
        self.sink_PP_collision = None
        self.numconc = None
        self.massconcentration = None
        self.volconcentration = None
        self.sink_sedimentation = None
        self.source_resuspension = None
        self.diagnostics = None
        self.classname = 'Floc'
        self.name = name
        self._seed_tep_diagnostics()
        self.p_exp = p_exp
        self.q_exp = q_exp
        self.dynamic_q_exp = dynamic_q_exp
        self.q_exp_at_t = q_exp  # [-] Exponent actually used at this step (diagnostic)
        self.d_crit_growth = d_crit_growth
        self.d_crit_exponent = d_crit_exponent
        self.d_crit_max_factor = d_crit_max_factor
        self.growth_limiter = 1.0  # [-] Breakage enhancement applied at this step (diagnostic)
        self.f_frac_floc_break = f_frac_floc_break
        self.efficiency_break = efficiency_break
        self.sinking_leak = sinking_leak

        self.resuspension_rate = resuspension_rate
        self.settling_velocity_factor = settling_velocity_factor
        self.apply_settling = apply_settling
        self.resusp_ewma_alpha = resusp_ewma_alpha
        self.organomin_coupling_fraction = organomin_coupling_fraction

        self.d_p_microflocdiam = d_p_microflocdiam
        self.diam = d_p_microflocdiam
        self.nf_fractal_dim = nf_fractal_dim
        self.density = density
        self.eps_kd = eps_kd
        self.dt2 = dt2
        self.time_conversion_factor = time_conversion_factor
        self.spinup_days = spinup_days

        self.setup = None


    def _seed_tep_diagnostics(self):
        """Seed the TEP-dependent parameters with their base (mineral-only) values.

        Only the row-0 diagnostics depend on this: get_tep_parameters() overwrites them at
        every derivative evaluation, so the trajectory is unaffected. Called again from
        set_coupling() on Macroflocs and Micro_in_Macro, once they have inherited the bases
        from the Microflocs "master" -- otherwise row 0 would report the class defaults.
        """
        if self.unified_alphas:
            self.alpha_PP = self.alpha_PF = self.alpha_FF = self.alpha_PP_base
        else:
            self.alpha_PP = self.alpha_PP_base
            self.alpha_PF = self.alpha_PF_base
            self.alpha_FF = self.alpha_FF_base
        self.fyflocstrength = self.fyflocstrength_base
        self.tau_cr = self.tau_cr_base

    def _calculate_macrofloc_diameter(self):
        """
        Calculate macrofloc diameter using Lee et al. (2011) fractal theory.
        DF = NC^(1/nf) * DP where NC is microflocs per macrofloc

        Returns:
            float: Calculated macrofloc diameter
        """
        return (self.coupled_Nt.numconc / self.coupled_Nf.numconc) ** (1 / self.nf_fractal_dim) * self.coupled_Np.diam


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
                     coupled_Np=None, # Microflocs
                     coupled_Nf=None, # Macroflocs
                     coupled_Nt=None, # Micro_in_Macro
                     coupled_glue=None,
                     ):
        self.coupled_Np = coupled_Np if self.name != 'Microflocs' else self
        self.coupled_Nf = coupled_Nf if self.name != 'Macroflocs' else self
        self.coupled_Nt = coupled_Nt if self.name != 'Micro_in_Macro' else self

        if self.name != 'Microflocs':
            self.nf_fractal_dim = self.coupled_Np.nf_fractal_dim
            self.f_frac_floc_break = self.coupled_Np.f_frac_floc_break

            self.p_exp = self.coupled_Np.p_exp
            self.q_exp = self.coupled_Np.q_exp
            self.dynamic_q_exp = self.coupled_Np.dynamic_q_exp
            self.d_crit_growth = self.coupled_Np.d_crit_growth
            self.d_crit_exponent = self.coupled_Np.d_crit_exponent
            self.d_crit_max_factor = self.coupled_Np.d_crit_max_factor
            self.d_p_microflocdiam = self.coupled_Np.d_p_microflocdiam
            self.sinking_leak = self.coupled_Np.sinking_leak
            self.diam = self.coupled_Np.diam

            self.efficiency_break = self.coupled_Np.efficiency_break
            # density must follow d_p_microflocdiam: both define the flocculus building block
            # and are only ever configured on the Microflocs "master". Without this line,
            # Micro_in_Macro computes its mass concentration and Macroflocs its settling
            # constant with the class default (2500) whatever the master carries.
            self.density = self.coupled_Np.density

            self.eps_kd = self.coupled_Np.eps_kd

            self.fyflocstrength_base = self.coupled_Np.fyflocstrength_base
            self.deltaFymax = self.coupled_Np.deltaFymax

            # Bases of the TEP-dependent parameters: configured on the master only.
            self.unified_alphas = self.coupled_Np.unified_alphas
            self.alpha_FF_base = self.coupled_Np.alpha_FF_base
            self.alpha_PP_base = self.coupled_Np.alpha_PP_base
            self.alpha_PF_base = self.coupled_Np.alpha_PF_base
            self.tau_cr_base = self.coupled_Np.tau_cr_base
            self._seed_tep_diagnostics()

            self.spinup_days = self.coupled_Np.spinup_days

            self.prescribe_tep_from_setup = self.coupled_Np.prescribe_tep_from_setup

        # set_ICs() runs before the couplings are wired, so Micro_in_Macro seeded its
        # mass concentration with its own class defaults for diam and density. Refresh it
        # now that both are inherited, otherwise row 0 of the SPMC diagnostic is wrong
        # (the trajectory itself is fine: update_val recomputes it at every step).
        if self.name == 'Micro_in_Macro' and self.numconc is not None:
            self.massconcentration = (self.numconc * np.pi / 6.
                                      * self.diam * self.diam * self.diam * self.density)

        # Initialize macrofloc diameter if couplings are now established
        if self.name == 'Macroflocs':
            self.diam = self._calculate_macrofloc_diameter()

        # Handle TEP coupling: either dynamic (coupled_glue) or prescribed (from Setup)
        if self.prescribe_tep_from_setup:
            # Create a PrescribedTEP instance that will be updated at each timestep
            self.coupled_glue = PrescribedTEP()
        else:
            # Standard dynamic coupling
            self.coupled_glue = coupled_glue

        # Setup shared alpha computation (performance optimization)
        self._setup_shared_alphas()

        # Optimization: Pre-compute settling velocity constants (Macroflocs only)
        if self.name == 'Macroflocs' and self.setup is not None:
            self._precompute_settling_constants()

    def _setup_shared_alphas(self):
        """
        Setup shared alpha computation for Flocs instances (performance optimization).

        SharedFlocTEPParameters is shared across Microflocs, Macroflocs, and Micro_in_Macro
        within a single simulation to avoid redundant calculations.

        This method is called during set_coupling(), which occurs once per Flocs instance
        during model initialization. Microflocs (first to be created) initializes the shared
        instance, and Macroflocs/Micro_in_Macro reuse it.
        """
        # Only Microflocs creates the shared instance
        # Always recreate to ensure each new simulation has fresh parameters
        if self.name == 'Microflocs':
            Flocs._shared_alphas_instance = SharedFlocTEPParameters(self, self.unified_alphas)

        # All Flocs instances (including Microflocs) get a reference to the shared instance
        if Flocs._shared_alphas_instance is not None:
            self.shared_alphas = Flocs._shared_alphas_instance

    def _precompute_settling_constants(self):
        """
        Pre-compute settling velocity constants for Macroflocs (performance optimization).
        Called once during setup to avoid recalculating constants at each timestep.
        Note: mu_water excluded as it may vary with temperature (Sharqawy et al. 2010).
        """
        if self.name != 'Macroflocs' or self.setup is None:
            return

        # Pre-compute mu_water-independent part of Winterwerp formula
        # W_s,F = (1/18) × (ρ_s - ρ_w) / μ × g × D_p^(3-nf) × D_F^(nf-1)
        # _settling_constant_base excludes μ to allow T-dependent viscosity
        self._settling_constant_base = (1.0/18.0) * 9.81 * (self.density - self.setup.rho_water)

    def update_val(self, numconc,
                   t=None,
                   t_idx=None,
                   debugverbose=False):
        self.numconc = numconc

        # For macroflocs, recalculate diameter based on fractal theory
        if self.name == 'Macroflocs':
            self.diam = self._calculate_macrofloc_diameter()

        self.volconcentration = numconc * np.pi / 6. * self.diam * self.diam * self.diam


        if self.name == "Microflocs" or self.name == "Micro_in_Macro":
            self.massconcentration = numconc * np.pi / 6. * self.diam*self.diam*self.diam * self.density


    def get_diagnostic_variables(self):
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics], dtype=self.dtype)

    def _compute_net_vertical_flux(self):
        """
        Compute net vertical flux for Macroflocs (sedimentation and resuspension).
        Called once per timestep, results stored in self.settling_loss.
        Reused by Micro_in_Macro via coupled_Nf.settling_loss.
        """
        # Settling velocity (Winterwerp formula). _settling_constant_base is the
        # mu-independent part, computed once in _precompute_settling_constants; only
        # mu_water_at_t varies with T (Sharqawy et al. 2010).
        self.settling_vel_base = (self._settling_constant_base / self.mu_water_at_t *
                                  self.d_p_microflocdiam ** (3.0 - self.nf_fractal_dim) *
                                  (self.diam ** (self.nf_fractal_dim - 1.0)) * self.apply_settling)

        if self.settling_velocity_factor is not None:
            # Simple approach: constant fraction of base settling velocity
            self.settling_vel = self.settling_vel_base * self.settling_velocity_factor
        else:
            self.settling_vel = self.settling_vel_base

        if self.resuspension_rate > 0:
            # Physical settling and resuspension
            self.sink_sedimentation = self.settling_vel * self.numconc / self.water_depth_at_t
            self.erosion_factor = max(0, (self.bed_shear_stress_at_t / self.tau_cr - 1))
            self.source_resuspension = self.resuspension_rate * self.erosion_factor / self.water_depth_at_t
            self.settling_loss = self.sink_sedimentation - self.source_resuspension

        else:
            # Fallback to original formulation
            self.sink_sedimentation = self.sinking_leak * self.settling_vel * self.numconc
            self.source_resuspension = 0
            self.settling_loss = self.sinking_leak * self.settling_vel * self.numconc

        # Compute net vertical loss rate for organic component coupling
        if self.numconc > 0:
            self.net_vertical_loss_rate = self.settling_loss / self.numconc  # [s-1]
        else:
            self.net_vertical_loss_rate = 0.0

    def get_sources(self, t=None, t_idx=None):
        """
        Note: the "FABM" and "get_sources/get_sinks" approach is not adopted here - for simplicity and also
        with some mathematical simplifications within the SMS equations to enhance computational efficiency
        (since the floc module is run at smaller time step)
        """

        # Update prescribed TEP value at this timestep (if using prescribed mode)
        if self.prescribe_tep_from_setup and isinstance(self.coupled_glue, PrescribedTEP):
            self.coupled_glue.C = self.setup.TEP_array[t_idx]
        self.coupled_glue_C = self.coupled_glue.C if self.coupled_glue else None

        # Check if we're in spin-up phase
        is_spinup = self.setup.in_spinup_phase

        if is_spinup:
            # Spin-up mode: disable TEP coupling and settling/resuspension (only mineral flocculation dynamics)
            use_tep_coupling = False
            # Store original values and disable settling
            original_apply_settling = self.apply_settling
            original_resuspension_rate = self.resuspension_rate
            self.apply_settling = False
            self.resuspension_rate = 0.0
        else:
            # Normal mode: use TEP coupling if available (including prescribed TEP)
            use_tep_coupling = self.coupled_glue and self.K_glue
            # Initialize variables to avoid reference errors
            original_apply_settling = None
            original_resuspension_rate = None


        if self.name != 'Micro_in_Macro':
            # Micro_in_Macro reads none of these: it mirrors them from Macroflocs below.
            self.g_shear_rate_at_t = self.setup.g_shear_rate_array[t_idx]
            self.bed_shear_stress_at_t = self.setup.bed_shear_stress_array[t_idx]
            self.water_depth_at_t = self.setup.water_depth_array[t_idx]
            self.mu_water_at_t = self.setup.mu_water_array[t_idx]

                # Optimization: Calculate once (Microflocs calculates, others reuse)
        # NOTE: Microflocs instance MUST be defined first in model configuration
        if self.name == 'Microflocs':
            # Calculate shared Ncnum once per timestep
            self._shared_ncnum = self.coupled_Nt.numconc / self.coupled_Nf.numconc
            # Compute Ncnum factors (used multiple times)
            self._factor_Ncnum_ratio = self._shared_ncnum / (self._shared_ncnum - 1)
            self._factor_inverse = 1 / (self._shared_ncnum - 1)

            # Get TEP-dependent parameters for Microflocs (used in flux computation)
            alphas, fyflocstrength, tau_cr, nf_fractal_dim, mm_TEP = self.shared_alphas.get_tep_parameters(
                self.coupled_glue if use_tep_coupling else None
            )
            self.alpha_FF, self.alpha_PP, self.alpha_PF = alphas
            self.fyflocstrength = fyflocstrength
            self.nf_fractal_dim = nf_fractal_dim


            # Compute expensive fractal terms once per timestep (optimized: reuse via multiplication)
            frac_inv_nf = 1.0 / self.nf_fractal_dim
            self._ncnum_frac_1_div_nf = self._shared_ncnum ** frac_inv_nf
            self._ncnum_frac_2_div_nf = self._ncnum_frac_1_div_nf * self._ncnum_frac_1_div_nf
            self._ncnum_frac_3_div_nf = self._ncnum_frac_2_div_nf * self._ncnum_frac_1_div_nf

            # Compute derived fractal terms
            self._ncnum_frac_1_div_nf_plus_1_cubed = (self._ncnum_frac_1_div_nf + 1.0) ** 3.0
            self._ncnum_frac_1_div_nf_minus_1 = self._ncnum_frac_1_div_nf - 1.0

            # Geometric terms
            self._np_diam_cubed = self.coupled_Np.diam ** 3.0
            self._np_diam_squared = self.coupled_Np.diam ** 2.0

            # Physical combinations
            self._mu_times_g_shear = self.mu_water_at_t * self.g_shear_rate_at_t

        # All instances: Load essential shared data
        self.Ncnum = self.coupled_Np._shared_ncnum

        # =====================================================
        # FLUX COMPUTATION (by first instance = Microflocs)
        # =====================================================
        if self.name == 'Microflocs':
            # Compute ALL collision flux bases (regardless of who uses them)
            if self.include_PP_collision:
                self._PP_collision_base = (2. / 3. * self.alpha_PP * self._np_diam_cubed *
                                           self.g_shear_rate_at_t * self.coupled_Np.numconc * self.coupled_Np.numconc)
            else:
                self._PP_collision_base = 0.0

            self._PF_collision_base = (1. / 6. * self.alpha_PF * self._ncnum_frac_1_div_nf_plus_1_cubed *
                                       self._np_diam_cubed * self.g_shear_rate_at_t *
                                       self.coupled_Np.numconc * self.coupled_Nf.numconc)

            self._FF_collision_base = (2. / 3. * self.alpha_FF * self._ncnum_frac_3_div_nf *
                                       self._np_diam_cubed * self.g_shear_rate_at_t *
                                       self.coupled_Nf.numconc * self.coupled_Nf.numconc)

            # q = 3 - nf closes the breakage kinetics dimensionally in Lee et al. (2011, 2014).
            # Following nf(t) rather than a fixed q_exp turns it into an extra TEP -> breakage
            # channel, since nf carries the TEP effect (nf_base + delta_nf * mm_TEP).
            self.q_exp_at_t = (3.0 - self.nf_fractal_dim) if self.dynamic_q_exp else self.q_exp

            # Shear-induced breakage alone cannot bound macrofloc growth: less breakage gives a
            # larger D_F, FF aggregation scales as D_F^(3-nf) and so grows too, N_F collapses
            # and the run ends in NaN. Lee et al. (2011) close this with a critical diameter
            # above which "the breakage rate was set sufficiently high to break all flocs".
            # Same idea here as a smooth ramp (an abrupt threshold would ring at dt = 86 s),
            # capped so that one explicit Euler step cannot overshoot the correction.
            # Inactive while D_F < d_crit_growth, hence a no-op for the reference configuration.
            if self.d_crit_growth is not None:
                floc_diam = self._ncnum_frac_1_div_nf * self.coupled_Np.diam
                self.growth_limiter = min(self.d_crit_max_factor,
                                          max(1.0, (floc_diam / self.d_crit_growth)
                                              ** self.d_crit_exponent))
            else:
                self.growth_limiter = 1.0

            self._breakage_base = (self.efficiency_break * self.g_shear_rate_at_t *
                                   self._ncnum_frac_1_div_nf_minus_1 ** self.p_exp *
                                   (self._mu_times_g_shear * self._ncnum_frac_2_div_nf *
                                    self._np_diam_squared / self.fyflocstrength) ** self.q_exp_at_t *
                                   self.coupled_Nf.numconc * self.growth_limiter)

        # =====================================================
        # SMS ASSEMBLY (pool-specific)
        # =====================================================
        if self.name == 'Microflocs':



            # Flux assembly for Microflocs
            self.sink_PP_collision = self.coupled_Np._PP_collision_base * self._factor_Ncnum_ratio
            self.sink_PF_collision = self.coupled_Np._PF_collision_base
            self.source_breakage = self.f_frac_floc_break * self.Ncnum * self.coupled_Np._breakage_base

            self.SMS = (-self.sink_PP_collision -
                        self.sink_PF_collision +
                        self.source_breakage)

        elif self.name == 'Macroflocs':
            # Get TEP-dependent parameters needed for settling
            self.tau_cr = self.shared_alphas.current_tau_cr
            self.nf_fractal_dim = self.shared_alphas.current_nf_fractal_dim
            # Mirror the shared alphas so the per-instance diagnostics report the value
            # actually in use (they are computed once, by Microflocs).
            self.alpha_FF, self.alpha_PP, self.alpha_PF = self.shared_alphas.current_alphas

            # Compute net vertical flux (stored in self.settling_loss)
            self._compute_net_vertical_flux()


            # Flux assembly for Macroflocs
            self.source_PP_collision = self.coupled_Np._PP_collision_base * self.coupled_Np._factor_inverse
            self.sink_FF_collision = self.coupled_Np._FF_collision_base
            self.source_breakage = self.coupled_Np._breakage_base

            self.SMS = (self.source_PP_collision -
                        self.sink_FF_collision +
                        self.source_breakage -
                        self.settling_loss)

        else:  # Micro_in_Macro
            # This pool owns no dynamics of its own: the flocculi it counts are the ones
            # bound inside the macroflocs, so every property below is Macroflocs' (or the
            # shared TEP set computed by Microflocs), mirrored here for the diagnostics.
            self.alpha_FF, self.alpha_PP, self.alpha_PF = self.shared_alphas.current_alphas
            self.fyflocstrength = self.shared_alphas.current_fyflocstrength
            self.tau_cr = self.shared_alphas.current_tau_cr
            self.nf_fractal_dim = self.shared_alphas.current_nf_fractal_dim
            self.settling_vel = self.coupled_Nf.settling_vel
            self.settling_vel_base = self.coupled_Nf.settling_vel_base
            self.erosion_factor = self.coupled_Nf.erosion_factor
            self.g_shear_rate_at_t = self.coupled_Nf.g_shear_rate_at_t
            self.bed_shear_stress_at_t = self.coupled_Nf.bed_shear_stress_at_t
            self.water_depth_at_t = self.coupled_Nf.water_depth_at_t
            self.mu_water_at_t = self.coupled_Nf.mu_water_at_t

            # Flux assembly for Micro_in_Macro. The number fluxes DO scale: this pool
            # counts flocculi, Macroflocs counts flocs, hence the Ncnum factor.
            self.source_PP_collision = self.coupled_Np._PP_collision_base * self.coupled_Np._factor_Ncnum_ratio
            self.source_PF_collision = self.coupled_Np._PF_collision_base
            self.sink_breakage = self.f_frac_floc_break * self.Ncnum * self.coupled_Np._breakage_base
            self.settling_loss = self.coupled_Nf.settling_loss * self.Ncnum
            self.sink_sedimentation = self.coupled_Nf.sink_sedimentation * self.Ncnum
            self.source_resuspension = self.coupled_Nf.source_resuspension * self.Ncnum

            self.SMS = (self.source_PP_collision +
                        self.source_PF_collision -
                        self.sink_breakage -
                        self.settling_loss)

        # Restore original values if we were in spin-up mode
        if self.setup.in_spinup_phase:
            self.apply_settling = original_apply_settling
            self.resuspension_rate = original_resuspension_rate

        return np.array(self.SMS, dtype=self.dtype)


    def get_sinks(self, t=None, t_idx=None):
        return np.array([0], dtype=self.dtype)

