import numpy as np

from ..core.base import BaseStateVar
from src.config import varinfos
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
                 efficiency_break=2e-4 / np.sqrt(3600 * 24),  # 2e-4/np.sqrt(3600*24*100),
                 fyflocstrength=1e-10,
                 mu_viscosity=1e-6 / 3600 / 24,  # 1e-6/3600/24/10000,

                 d_p_microflocdiam=18e-6,
                 nf_fractal_dim=2,

                 density=2500,  # [kg/m3]

                 alpha_FF_ref = None,
                 alpha_PP_ref = None,
                 alpha_PF_ref = None,
                 K_glue = None,
                 deltaFymax=1e-9,

                 # eps_kd=2e-5 * varinfos.molmass_C,  # [m2 mgC-1] Diffuse attenuation cross section
                 # # value from
                 # eps_kd=2e-5,  # [m2 mg-1] Diffuse attenuation cross section of SPM with kd = ... + eps_kd * sqrt(SPM)
                 # # value from Tian et al., 2009
                 eps_kd = 0.066* 1e3, # m-1 (mg l-1)-1 Light attenuation due to SPM
                 # model is in kg m-3

                 sinking_leak=0,

                 dt2=True,
                 time_conversion_factor = 86400,   # Flocs run in s-1 while the rest of the model is in d-1
                 ):

        super().__init__()

        self.deltaFymax = deltaFymax
        self.fyflocstrength_ref = fyflocstrength
        self.K_glue = K_glue
        self.alpha_FF_ref = alpha_FF_ref
        self.alpha_PP_ref = alpha_PP_ref
        self.alpha_PF_ref = alpha_PF_ref
        self.SMS = None
        self.settling_vel = None
        self.breaksource = None
        self.ffloss = None
        self.ppsource = None
        self.settling_loss = None
        self.g_shear_rate_at_t = None
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
        self.diagnostics = None
        self.classname = 'Floc'
        self.name = name
        self.alpha_PP = alpha_PP
        self.alpha_PF = alpha_PF
        self.alpha_FF = alpha_FF
        self.p_exp = p_exp
        self.q_exp = q_exp
        self.f_frac_floc_break = f_frac_floc_break
        self.efficiency_break = efficiency_break
        self.fyflocstrength = fyflocstrength
        self.mu_viscosity = mu_viscosity
        self.sinking_leak = sinking_leak

        self.d_p_microflocdiam = d_p_microflocdiam
        self.diam = d_p_microflocdiam
        self.nf_fractal_dim = nf_fractal_dim
        self.density = density
        self.eps_kd = eps_kd
        self.dt2 = dt2
        self.time_conversion_factor = time_conversion_factor

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
            self.alpha_PP = self.coupled_Np.alpha_PP
            self.alpha_PF = self.coupled_Np.alpha_PF
            self.alpha_FF = self.coupled_Np.alpha_FF
            self.f_frac_floc_break = self.coupled_Np.f_frac_floc_break
            self.efficiency_break = self.coupled_Np.efficiency_break
            self.fyflocstrength = self.coupled_Np.fyflocstrength

            self.eps_kd = self.coupled_Np.eps_kd

            self.K_glue = self.coupled_Np.K_glue
            self.alpha_FF_ref = self.coupled_Np.alpha_FF_ref
            self.alpha_PP_ref = self.coupled_Np.alpha_PP_ref
            self.alpha_PF_ref = self.coupled_Np.alpha_PF_ref
            self.deltaFymax = self.coupled_Np.deltaFymax

        if self.name == 'Micro_in_Macro':
            self.sinking_leak = self.coupled_Nf.sinking_leak

        # Initialize macrofloc diameter if couplings are now established
        if self.name == 'Macroflocs':
            self.diam = self._calculate_macrofloc_diameter()


    def update_val(self, numconc,
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

        if self.coupled_glue and self.K_glue:
                self.alpha_FF = self.alpha_FF_ref * self.coupled_glue.C / (self.K_glue + self.coupled_glue.C)

                # Optionally modify alpha_PP if reference value exists
                if hasattr(self, 'alpha_PP_ref') and self.alpha_PP_ref is not None:
                    self.alpha_PP = self.alpha_PP_ref * self.coupled_glue.C / (self.K_glue + self.coupled_glue.C)

                # Optionally modify alpha_PF if reference value exists
                if hasattr(self, 'alpha_PF_ref') and self.alpha_PF_ref is not None:
                    self.alpha_PF = self.alpha_PF_ref * self.coupled_glue.C / (self.K_glue + self.coupled_glue.C)

                self.fyflocstrength = (self.fyflocstrength_ref + self.deltaFymax *
                                               self.coupled_glue.C / (self.K_glue + self.coupled_glue.C))

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
            # Calculate settling velocity for sinking loss
            self.settling_vel = 2000 / 9 * (1200 - 1000) * 9.81 * (self.diam * self.diam / 4) / 1e-3 * 1e-3 * 3600 * 24
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
                        # First very rough test - note that viscosity could (should) be reused from mu_viscosity.
                        self.settling_loss
                        )

        else:

            self.settling_vel = self.coupled_Nf.settling_vel

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
                        self.sinking_leak * self.settling_vel * self.numconc

                        )


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

