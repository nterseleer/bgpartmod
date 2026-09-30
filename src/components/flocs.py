import numpy as np

from ..core.base import BaseStateVar


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
    """No longer used by the model (replaced by FlocProcesses on 2026-09-30). Kept only so
    that the simulations saved before still load: their Flocs carry an instance of it."""


class FlocProcesses:
    """The processes of the floc population, computed once per timestep for its three pools.

    The three Flocs pools -- Microflocs (free flocculi, N_P), Macroflocs (flocs, N_F) and
    Micro_in_Macro (flocculi bound in flocs, N_T) -- exchange numbers through the same
    processes: a collision or a breakage event is a sink for one pool and a source for
    another. The processes are therefore computed here, once per timestep, from the current
    state of the three pools; each pool then assembles its own sources and sinks from them
    (Flocs.get_sources / get_sinks).

    Shared by the three pools (Flocs._processes). Its parameters are those configured on
    Microflocs, except the vertical exchange ones (resuspension_rate, apply_settling), which
    are configured on Macroflocs.
    """

    # Properties of the floc population that every pool reports as a diagnostic
    COMMON_DIAGNOSTICS = ('alpha_FF', 'alpha_PP', 'alpha_PF', 'fyflocstrength', 'tau_cr',
                          'nf_fractal_dim', 'Ncnum', 'coupled_glue_C', 'q_exp_at_t',
                          'growth_limiter', 'settling_vel_base', 'settling_vel',
                          'erosion_factor', 'g_shear_rate_at_t', 'bed_shear_stress_at_t',
                          'water_depth_at_t', 'mu_water_at_t')

    def __init__(self, microflocs, macroflocs, micro_in_macro):
        self.microflocs = microflocs
        self.macroflocs = macroflocs
        self.micro_in_macro = micro_in_macro
        self.t_idx = None     # timestep of the last computation (None: never computed)
        self.stale = True     # a pool changed since the last computation
        self._bound = False

    def _bind(self):
        """Read the parameters, once, when every coupling is in place."""
        master, macro = self.microflocs, self.macroflocs
        self.setup = master.setup
        self.glue = master.coupled_glue

        self.unified_alphas = master.unified_alphas
        self.alpha_FF_base, self.delta_alpha_FF = master.alpha_FF_base, master.delta_alpha_FF
        self.alpha_PP_base, self.delta_alpha_PP = master.alpha_PP_base, master.delta_alpha_PP
        self.alpha_PF_base, self.delta_alpha_PF = master.alpha_PF_base, master.delta_alpha_PF
        self.fyflocstrength_base, self.deltaFymax = master.fyflocstrength_base, master.deltaFymax
        self.tau_cr_base, self.delta_tau_cr = master.tau_cr_base, master.delta_tau_cr
        self.nf_base, self.delta_nf = master.nf_fractal_dim_base, master.delta_nf_fractal_dim
        self.K_glue = master.K_glue
        self.prescribe_tep = master.prescribe_tep_from_setup

        self.include_PP_collision = master.include_PP_collision
        self.p_exp, self.q_exp, self.dynamic_q_exp = master.p_exp, master.q_exp, master.dynamic_q_exp
        self.d_crit_growth = master.d_crit_growth
        self.d_crit_exponent, self.d_crit_max_factor = master.d_crit_exponent, master.d_crit_max_factor
        self.f_frac_floc_break = master.f_frac_floc_break
        self.efficiency_break = master.efficiency_break
        self.D_p = master.d_p_microflocdiam

        self.resuspension_rate = macro.resuspension_rate
        self.apply_settling = macro.apply_settling
        # Winterwerp settling velocity, the part independent of the viscosity (which
        # varies with T, Sharqawy et al. 2010): (1/18) g (rho_s - rho_w)
        self.settling_constant = (1.0 / 18.0) * 9.81 * (master.density - self.setup.rho_water)
        self._bound = True

    def update(self, t_idx):
        """Compute the processes at timestep t_idx, unless they are already up to date."""
        if not self.stale and t_idx == self.t_idx:
            return
        if not self._bound:
            self._bind()
        self.stale = False
        self.t_idx = t_idx
        setup = self.setup

        N_P = self.microflocs.numconc
        N_F = self.macroflocs.numconc
        N_T = self.micro_in_macro.numconc

        # --- Physical forcing ---
        G = float(setup.g_shear_rate_array[t_idx])          # [s-1] shear rate
        tau_b = float(setup.bed_shear_stress_array[t_idx])  # [Pa] bed shear stress
        h = float(setup.water_depth_array[t_idx])           # [m] water depth
        mu = float(setup.mu_water_array[t_idx])             # [kg m-1 s-1] viscosity

        # --- TEP effect: each parameter = mineral base + delta * TEP saturation ---
        if self.prescribe_tep:
            self.glue.C = float(setup.TEP_array[t_idx])
        self.coupled_glue_C = self.glue.C if self.glue else None
        mm_TEP = (self.coupled_glue_C / (self.K_glue + self.coupled_glue_C)
                  if self.glue and self.K_glue else 0.0)
        if self.unified_alphas:
            self.alpha_FF = self.alpha_PP = self.alpha_PF = \
                self.alpha_PP_base + self.delta_alpha_PP * mm_TEP
        else:
            self.alpha_FF = self.alpha_FF_base + self.delta_alpha_FF * mm_TEP
            self.alpha_PP = self.alpha_PP_base + self.delta_alpha_PP * mm_TEP
            self.alpha_PF = self.alpha_PF_base + self.delta_alpha_PF * mm_TEP
        self.fyflocstrength = F_y = self.fyflocstrength_base + self.deltaFymax * mm_TEP
        self.tau_cr = tau_cr = self.tau_cr_base + self.delta_tau_cr * mm_TEP
        self.nf_fractal_dim = n_f = self.nf_base + self.delta_nf * mm_TEP

        # --- Floc geometry (fractal theory, Lee et al. 2011) ---
        # N_c = flocculi per floc; D_F = N_c^(1/nf) D_P
        self.Ncnum = N_c = N_T / N_F
        x1 = N_c ** (1.0 / n_f)
        x2 = x1 * x1
        x3 = x2 * x1
        D_P = self.D_p
        D_P2 = D_P * D_P
        D_P3 = D_P ** 3.0
        self.D_F = x1 * D_P
        self.volconc_F = N_F * np.pi / 6. * self.D_F * self.D_F * self.D_F

        # --- Collisions [# m-3 s-1] ---
        # flocculus-flocculus (PP), flocculus-floc (PF), floc-floc (FF)
        PP = (2. / 3. * self.alpha_PP * D_P3 * G * N_P * N_P
              if self.include_PP_collision else 0.0)
        PF = 1. / 6. * self.alpha_PF * (x1 + 1.0) ** 3.0 * D_P3 * G * N_P * N_F
        FF = 2. / 3. * self.alpha_FF * x3 * D_P3 * G * N_F * N_F

        # --- Breakage of flocs [# m-3 s-1] ---
        # q = 3 - nf closes the breakage kinetics dimensionally in Lee et al. (2011, 2014).
        # Following nf(t) rather than a fixed q_exp turns it into an extra TEP -> breakage
        # channel, since nf carries the TEP effect.
        self.q_exp_at_t = q = (3.0 - n_f) if self.dynamic_q_exp else self.q_exp
        # Shear-induced breakage alone cannot bound macrofloc growth: less breakage gives a
        # larger D_F, FF aggregation scales as D_F^(3-nf) and so grows too, N_F collapses
        # and the run ends in NaN. Lee et al. (2011) close this with a critical diameter
        # above which "the breakage rate was set sufficiently high to break all flocs".
        # Same idea here as a smooth ramp (an abrupt threshold would ring at dt = 86 s),
        # capped so that one explicit Euler step cannot overshoot the correction.
        # Inactive while D_F < d_crit_growth, hence a no-op for the reference configuration.
        if self.d_crit_growth is not None:
            self.growth_limiter = min(self.d_crit_max_factor,
                                      max(1.0, (self.D_F / self.d_crit_growth) ** self.d_crit_exponent))
        else:
            self.growth_limiter = 1.0
        B = (self.efficiency_break * G * (x1 - 1.0) ** self.p_exp *
             (mu * G * x2 * D_P2 / F_y) ** q * N_F * self.growth_limiter)

        # --- Vertical exchange of flocs [# m-3 s-1] ---
        # Settling velocity (Winterwerp): (1/18) g (rho_s - rho_w) / mu D_P^(3-nf) D_F^(nf-1)
        self.settling_vel_base = self.settling_vel = (
            self.settling_constant / mu * D_P ** (3.0 - n_f) * self.D_F ** (n_f - 1.0)
            * self.apply_settling)
        if self.resuspension_rate > 0:
            self.sedimentation = self.settling_vel * N_F / h
            self.erosion_factor = max(0.0, tau_b / tau_cr - 1.0)
            self.resuspension = self.resuspension_rate * self.erosion_factor / h
        else:
            # No resuspension: no vertical exchange at all
            self.sedimentation = self.resuspension = self.erosion_factor = 0.0
        self.settling_loss = self.sedimentation - self.resuspension
        self.net_vertical_loss_rate = self.settling_loss / N_F if N_F > 0 else 0.0  # [s-1]

        # --- The same processes, counted in flocs or in flocculi ---
        # A PP collision turns N_c/(N_c-1) free flocculi into 1/(N_c-1) new flocs; a broken
        # floc releases a fraction f of its N_c flocculi; a settling floc carries N_c.
        self.PP_flocculi = PP * N_c / (N_c - 1.0)
        self.PP_flocs = PP / (N_c - 1.0)
        self.PF = PF
        self.FF = FF
        self.breakage = B
        self.breakage_flocculi = self.f_frac_floc_break * N_c * B
        self.sedimentation_flocculi = self.sedimentation * N_c
        self.resuspension_flocculi = self.resuspension * N_c
        self.settling_loss_flocculi = self.settling_loss * N_c

        # Row-for-row with the other two pools, the current TEP-dependent properties
        # (the pools only report them as diagnostics: see publish_diagnostics)
        self.mm_TEP = mm_TEP
        self.g_shear_rate_at_t = G
        self.bed_shear_stress_at_t = tau_b
        self.water_depth_at_t = h
        self.mu_water_at_t = mu

    def publish_diagnostics(self, pool):
        """Copy the properties of the floc population onto `pool`, for its diagnostics.

        Done when the diagnostics are read, not at every timestep: most of them are only
        ever recorded on the output steps."""
        if self.t_idx is None:   # not computed yet: the pool keeps its seeded values
            return
        for name in self.COMMON_DIAGNOSTICS:
            setattr(pool, name, getattr(self, name))


class Flocs(BaseStateVar):

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

                 unified_alphas = True,    # [-] Single alpha for FF=PP=PF (alpha_PP_base, delta_alpha_PP)
                 include_PP_collision = True,  # [-] Include PP collision term (microfloc-microfloc)
                 delta_tau_cr = 0.2,       # [Pa] TEP increment for critical shear stress
                 delta_nf_fractal_dim = 0.0,  # [-] TEP increment for fractal dimension

                 # TEP coupling parameters
                 K_glue = None,            # [mmol m-3] Half-saturation for TEP effect
                 prescribe_tep_from_setup = False,  # [-] Use prescribed TEP from Setup instead of coupled_glue
                 #
                 resuspension_rate = 0.,   # [# m-2 s-1 Pa-1] Erosion rate constant (0 = no vertical exchange)
                 apply_settling = True, # Boolean. Whether to apply sediment settling

                 # Vertical coupling parameters for BGC components
                 resusp_ewma_alpha = 0.0,  # [-] α=0: fixed ratio (Option A), α>0: adaptive smoothing (Option B1)
                 organomin_coupling_fraction = 1.0,  # [-] Fraction forming organo-mineral aggregates

                 # Light attenuation by SPM: 0.066 [m-1 (mg l-1)-1] converted to model units
                 # [m-1 (kg m-3)-1]. Tian et al. (2009) give an alternative in sqrt(SPM).
                 eps_kd = 0.066 * 1e3,  # [m2 kg-1]

                 dt2=True,
                 dtype=np.float64,
                 time_conversion_factor = 86400,   # Flocs run in s-1 while the rest of the model is in d-1
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
        self.unified_alphas = unified_alphas

        # Processes shared by the three pools (FlocProcesses, wired in set_coupling).
        # _floc_processes is held by Microflocs, _processes is each pool's reference to it.
        self._floc_processes = None
        self._processes = None

        self.settling_vel = None
        self.settling_vel_base = None
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
        self.ICs = None
        self.sink_breakage = None
        self.source_PF_collision = None
        self.source_PP_collision = None
        self.source_breakage = None
        self.sink_FF_collision = None
        self.sink_PF_collision = None
        self.sink_PP_collision = None
        self.settling_loss = None
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

        self.resuspension_rate = resuspension_rate
        self.apply_settling = apply_settling
        self.resusp_ewma_alpha = resusp_ewma_alpha
        self.organomin_coupling_fraction = organomin_coupling_fraction

        self.d_p_microflocdiam = d_p_microflocdiam
        self.diam = d_p_microflocdiam
        self.nf_fractal_dim_base = nf_fractal_dim   # configured value (mineral base)
        self.nf_fractal_dim = nf_fractal_dim        # value in use at this step (diagnostic)
        self.density = density
        self.eps_kd = eps_kd
        self.dt2 = dt2
        self.time_conversion_factor = time_conversion_factor

        self.setup = None


    def _seed_tep_diagnostics(self):
        """Seed the TEP-dependent parameters with their base (mineral-only) values.

        Only the row-0 diagnostics depend on this: FlocProcesses overwrites them at every
        timestep, so the trajectory is unaffected. Called again from set_coupling() on
        Macroflocs and Micro_in_Macro, once they have inherited the bases from the
        Microflocs "master" -- otherwise row 0 would report the class defaults.
        """
        if self.unified_alphas:
            self.alpha_PP = self.alpha_PF = self.alpha_FF = self.alpha_PP_base
        else:
            self.alpha_PP = self.alpha_PP_base
            self.alpha_PF = self.alpha_PF_base
            self.alpha_FF = self.alpha_FF_base
        self.fyflocstrength = self.fyflocstrength_base
        self.tau_cr = self.tau_cr_base

    def set_ICs(self,
                numconc
                ):
        self.numconc = numconc
        self.ICs = [self.numconc]

        if self.name == "Microflocs" or self.name == "Micro_in_Macro":
            self.massconcentration = numconc * np.pi / 6. * self.diam*self.diam*self.diam * self.density

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
            self.nf_fractal_dim = self.coupled_Np.nf_fractal_dim_base
            self.f_frac_floc_break = self.coupled_Np.f_frac_floc_break

            self.p_exp = self.coupled_Np.p_exp
            self.q_exp = self.coupled_Np.q_exp
            self.dynamic_q_exp = self.coupled_Np.dynamic_q_exp
            self.d_crit_growth = self.coupled_Np.d_crit_growth
            self.d_crit_exponent = self.coupled_Np.d_crit_exponent
            self.d_crit_max_factor = self.coupled_Np.d_crit_max_factor
            self.d_p_microflocdiam = self.coupled_Np.d_p_microflocdiam
            self.diam = self.coupled_Np.diam

            self.efficiency_break = self.coupled_Np.efficiency_break
            # density must follow d_p_microflocdiam: both define the flocculus building block
            # and are only ever configured on the Microflocs "master". Without this line,
            # Micro_in_Macro computes its mass concentration and Macroflocs its settling
            # constant with the class default whatever the master carries.
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

            self.prescribe_tep_from_setup = self.coupled_Np.prescribe_tep_from_setup

        # set_ICs() runs before the couplings are wired, so Micro_in_Macro seeded its
        # mass concentration with its own class defaults for diam and density. Refresh it
        # now that both are inherited, otherwise row 0 of the SPMC diagnostic is wrong
        # (the trajectory itself is fine: update_val recomputes it at every step).
        if self.name == 'Micro_in_Macro' and self.numconc is not None:
            self.massconcentration = (self.numconc * np.pi / 6.
                                      * self.diam * self.diam * self.diam * self.density)

        # Initial macrofloc diameter, D_F = N_c^(1/nf) D_P (row-0 diagnostic; FlocProcesses
        # computes it at every timestep)
        if self.name == 'Macroflocs':
            self.diam = ((self.coupled_Nt.numconc / self.numconc) ** (1 / self.nf_fractal_dim)
                         * self.coupled_Np.diam)

        # Handle TEP coupling: either dynamic (coupled_glue) or prescribed (from Setup)
        if self.prescribe_tep_from_setup:
            # A PrescribedTEP instance, updated at each timestep from Setup.TEP_array
            self.coupled_glue = PrescribedTEP()
        else:
            # Standard dynamic coupling
            self.coupled_glue = coupled_glue

        # The processes shared by the three pools, created by whichever of them is coupled
        # first and held by Microflocs
        master = self.coupled_Np
        if master._floc_processes is None:
            master._floc_processes = FlocProcesses(master, self.coupled_Nf, self.coupled_Nt)
        self._processes = master._floc_processes

    def update_val(self, numconc,
                   t=None,
                   t_idx=None,
                   debugverbose=False):
        self.numconc = numconc
        self._processes.stale = True

        if self.name == "Microflocs" or self.name == "Micro_in_Macro":
            # Flocculi: their diameter is D_P, constant
            self.volconcentration = numconc * np.pi / 6. * self.diam * self.diam * self.diam
            self.massconcentration = self.volconcentration * self.density

    def get_diagnostic_variables(self):
        self._processes.publish_diagnostics(self)
        return super().get_diagnostic_variables()

    def get_sources(self, t=None, t_idx=None):
        """Sources of this pool, from the processes shared by the three pools.

        The pool's source and sink terms are kept as attributes, as in the other
        components (and for the diagnostics).
        """
        processes = self._processes
        processes.update(t_idx)

        if self.name == 'Microflocs':
            self.source_breakage = processes.breakage_flocculi
            sources = self.source_breakage

        elif self.name == 'Macroflocs':
            self.diam = processes.D_F
            self.volconcentration = processes.volconc_F
            self.source_PP_collision = processes.PP_flocs
            self.source_breakage = processes.breakage
            self.source_resuspension = processes.resuspension
            self.sink_sedimentation = processes.sedimentation
            self.settling_loss = processes.settling_loss
            self.net_vertical_loss_rate = processes.net_vertical_loss_rate
            sources = self.source_PP_collision + self.source_breakage + self.source_resuspension

        else:  # Micro_in_Macro: the flocculi bound inside the macroflocs
            self.source_PP_collision = processes.PP_flocculi
            self.source_PF_collision = processes.PF
            self.source_resuspension = processes.resuspension_flocculi
            self.sink_sedimentation = processes.sedimentation_flocculi
            self.settling_loss = processes.settling_loss_flocculi
            sources = self.source_PP_collision + self.source_PF_collision + self.source_resuspension

        return np.array((sources,), dtype=self.dtype)

    def get_sinks(self, t=None, t_idx=None):
        processes = self._processes
        processes.update(t_idx)

        if self.name == 'Microflocs':
            self.sink_PP_collision = processes.PP_flocculi
            self.sink_PF_collision = processes.PF
            sinks = self.sink_PP_collision + self.sink_PF_collision

        elif self.name == 'Macroflocs':
            self.sink_FF_collision = processes.FF
            sinks = self.sink_FF_collision + processes.sedimentation

        else:  # Micro_in_Macro
            self.sink_breakage = processes.breakage_flocculi
            sinks = self.sink_breakage + processes.sedimentation_flocculi

        return np.array((sinks,), dtype=self.dtype)
