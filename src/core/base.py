import numpy as np
from ..utils import functions as fns

from src.config_model import varinfos


class BaseStateVar:
    def __init__(self, dtype=np.float64):
        self.dtype = dtype
        self.setup = None
        self.time_conversion_factor = 1
        self.default_diagnostics = True
        self.diagnostics = None
        self.spinup_days = 0

    def get_all_diagnostics(self):
        """Get all numeric non-private attributes as diagnostics."""
        if not hasattr(self, '_diagnostic_attrs'):
            self._diagnostic_attrs = [fns.get_nested_attr(self, attr) for attr in vars(self) if not attr.startswith('_')
                                      and not attr.startswith('coupled')
                                      and attr not in ('setup', 'diagnostic')
                                      and isinstance(getattr(self, attr), (int, float, np.number, type(None)))]
        return self._diagnostic_attrs

    def get_diagnostic_variables(self):
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics], dtype=self.dtype)

    def _precompute_temp_limitation(self, A_E, T_ref, boltz=False,
                                   bound_temp_to_1=True, suffix=''):
        """
        Pre-compute temperature limitation array for entire simulation.
        Optimization: Called once during set_coupling() to avoid repeated calculations.

        Args:
            A_E: Activation energy parameter
            T_ref: Reference temperature [K]
            boltz: Use Boltzmann constant if True
            bound_temp_to_1: Bound limitation between 0 and 1
            suffix: Identifier suffix (e.g., '_growth', '_grazing') for multiple limT arrays
        """
        if self.setup is None:
            raise ValueError("setup must be initialized before pre-computing temperature limitation")

        limT_array = fns.getlimT(
            self.setup.T_array,
            A_E=A_E, T_ref=T_ref, boltz=boltz,
            bound_temp_to_1=bound_temp_to_1,
            T_max=self.setup.T_max
        )
        setattr(self, f'limT{suffix}_array', limT_array)


class BaseOrg(BaseStateVar):
    def __init__(self, dtype=np.float64, vertical_coupling_alpha=0.0):

        super().__init__(dtype=dtype)

        self.C = None
        self.N = None
        self.Chl = None
        self.P = None
        self.Si = None

        self.QN = None
        self.QP = None
        self.QSi = None
        self.thetaC = None
        self.checkQmax = False

        self.ICs = None

        self.setup = None

        # Vertical coupling parameter (EWMA filter for resuspension ratios)
        # α=0: fixed ratio (Option A), α>0: adaptive smoothing (Option B1)
        self.vertical_coupling_alpha = vertical_coupling_alpha

        # Initialize smoothed ratios (set by initialize_vertical_coupling_ratios after coupling)
        self.smoothed_C_to_Nf_ratio = None
        self.smoothed_N_to_Nf_ratio = None
        self.smoothed_P_to_Nf_ratio = None
        self.smoothed_Si_to_Nf_ratio = None

        # SMS terms
        self.C_SMS = None
        self.N_SMS = None
        self.Chl_SMS = None
        self.P_SMS = None
        self.Si_SMS = None

        self.C_sources = None
        self.N_sources = None
        self.Chl_sources = None
        self.P_sources = None
        self.Si_sources = None

        self.C_sinks = None
        self.N_sinks = None
        self.Chl_sinks = None
        self.P_sinks = None
        self.Si_sinks = None

    def set_ICs(self, C,
                N=None,
                Chl=None,
                P=None,
                Si=None,
                Qmaxratio = 0.95,
                checkQmax = False):
        self.C = C
        self.N = N
        self.Chl = Chl
        self.P = P
        self.Si = Si

        if self.N is not None:
            if self.name == 'Phy':
                if self.QratioIC:
                    Qmaxratio = self.QratioIC
                    self.N = self.C * self.QN_max * Qmaxratio
            self.QN = self.N / self.C
            if self.name == 'Phy' and self.checkQmax  and self.QN > self.QN_max:
                self.N = self.C * self. QN_max * Qmaxratio
                print('Phytoplankton N initial pool too high compared to QN_max, changed from {} to {}'.format(N, self.N))
                self.QN = self.N / self.C
        if self.P is not None:
            if self.name == 'Phy':
                self.P = self.C * self.QP_max * Qmaxratio
            self.QP = self.P / self.C
            if self.name == 'Phy' and self.checkQmax  and self.QP > self.QP_max:
                self.P = self.C * self. QP_max * Qmaxratio
                print('Phytoplankton P initial pool too high compared to QP_max, changed from {} to {}'.format(P, self.P))
                self.QP = self.P / self.C
        if self.Si is not None:
            if self.name == 'Phy':
                self.Si = self.C * self.QSi_max * Qmaxratio
            self.QSi = self.Si / self.C
            if self.name == 'Phy' and self.checkQmax  and self.QSi > self.QSi_max:
                self.Si = self.C * self. QSi_max * Qmaxratio
                print('Phytoplankton Si initial pool too high compared to QSi_max, changed from {} to {}'.format(Si, self.Si))
                self.QSi = self.Si / self.C
        if self.P is not None and self.Si is not None:
            self.fnut = min(self.QN, self.QP, self.QSi)
        if self.Chl is not None:
            self.thetaC = self.Chl / self.C #

        self.ICs = np.array([pool for pool in [self.C, self.N, self.Chl, self.P, self.Si] if pool is not None], dtype=self.dtype)

    def initialize_vertical_coupling_ratios(self):
        """
        Initialize BGC/floc ratios for resuspension after ICs are set.
        Separates sedimentation (proportional to current C) from resuspension (absolute flux).

        With vertical_coupling_alpha:
        - α=0 (default): Fixed ratio from initial conditions (Option A)
        - α>0: Exponentially weighted moving average (Option B1)
        """
        if hasattr(self, 'coupled_aggregate') and self.coupled_aggregate is not None:
            Nf = self.coupled_aggregate.numconc
            if Nf > 0:
                self.smoothed_C_to_Nf_ratio = self.C / Nf if self.C is not None else None
                self.smoothed_N_to_Nf_ratio = self.N / Nf if self.N is not None else None
                self.smoothed_P_to_Nf_ratio = self.P / Nf if self.P is not None else None
                self.smoothed_Si_to_Nf_ratio = self.Si / Nf if self.Si is not None else None
            else:
                self.smoothed_C_to_Nf_ratio = None
                self.smoothed_N_to_Nf_ratio = None
                self.smoothed_P_to_Nf_ratio = None
                self.smoothed_Si_to_Nf_ratio = None
        else:
            self.smoothed_C_to_Nf_ratio = None
            self.smoothed_N_to_Nf_ratio = None
            self.smoothed_P_to_Nf_ratio = None
            self.smoothed_Si_to_Nf_ratio = None

    def update_val(self, C,
                   N=None,
                   Chl=None,
                   P=None,
                   Si=None,
                   t=None,
                   t_idx=None,
                   debugverbose=False):

        if debugverbose:
            print('Checking update_val for {} with values before: '.format(self.name),
                  self.C, self.N, self.Chl, self.P, self.Si)

        self.C = C
        self.N = N
        self.Chl = Chl
        self.P = P
        self.Si = Si

        # Diagnostic variables
        if self.N is not None:
            self.QN = self.N / self.C
        if self.P is not None:
            self.QP = self.P / self.C
        if self.Si is not None:
            self.QSi = self.Si / self.C
        if self.Chl is not None:
            self.thetaC = self.Chl / self.C

        if debugverbose:
            print('Values after update: ', self.C, self.N, self.Chl, self.P, self.Si)
            input()


class Elms:
    def __init__(self, dict=False):
        self.C = None if not dict else {}
        self.N = None if not dict else {}
        self.Chl = None if not dict else {}
        self.P = None if not dict else {}
        self.Si = None if not dict else {}
        self.NH4 = None if not dict else {}
        self.NO3 = None if not dict else {}
