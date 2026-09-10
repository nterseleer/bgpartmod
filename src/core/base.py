import numpy as np
from ..utils import functions as fns


class BaseStateVar:
    def __init__(self, dtype=np.float64):
        self.dtype = dtype
        self.setup = None
        self.debug_mode = False
        self.time_conversion_factor = 1
        self.diagnostics = None
        self.spinup_days = 0

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
    def __init__(self, dtype=np.float64):

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

        # Smoothed ratios for vertical coupling (set by initialize_vertical_coupling_ratios)
        # Controlled by coupled_aggregate.resusp_ewma_alpha
        self.smoothed_C_to_Nf_ratio = None
        self.smoothed_N_to_Nf_ratio = None
        self.smoothed_P_to_Nf_ratio = None
        self.smoothed_Si_to_Nf_ratio = None

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
                Qmaxratio = 0.95):
        """Set the initial pools and the derived quotas.

        ASYMETRIE CONNUE, A CORRIGER (Phy uniquement) : N n'est reecrit que si QratioIC
        est defini, alors que P et Si sont reecrits a C * Q_max * Qmaxratio dans TOUS les
        cas, ce qui ecarte silencieusement les IC de P et Si donnees en configuration
        (~ -4.5 % avec Qmaxratio = 0.95). Non corrige pour l'instant : la correction
        change les trajectoires et casserait la comparabilite avec les optimisations en
        cours. Correction = aligner P et Si sur N, c'est-a-dire encadrer les deux blocs
        ci-dessous par `if self.QratioIC:` comme pour N.
        """
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

    def set_vertical_coupling_state(self, state_dict):
        """
        Set smoothed BGC/floc ratios from extracted state (e.g., from a previous simulation).

        Args:
            state_dict: Dict with keys like 'smoothed_C_to_Nf_ratio', 'smoothed_N_to_Nf_ratio', etc.
        """
        for attr, value in state_dict.items():
            if hasattr(self, attr):
                setattr(self, attr, value)

    def initialize_vertical_coupling_ratios(self):
        """
        Initialize BGC/floc ratios for resuspension after ICs are set.
        Respects pre-set values from set_vertical_coupling_state().

        With resusp_ewma_alpha:
        - α=0 (default): Fixed ratio from initial conditions (Option A)
        - α>0: Exponentially weighted moving average (Option B1)
        """
        # Skip if values already set (e.g., from extracted state)
        if self.smoothed_C_to_Nf_ratio is not None:
            return

        if hasattr(self, 'coupled_aggregate') and self.coupled_aggregate is not None:
            Nf = self.coupled_aggregate.numconc
            if Nf > 0:
                self.smoothed_C_to_Nf_ratio = self.C / Nf if self.C is not None else None
                self.smoothed_N_to_Nf_ratio = self.N / Nf if self.N is not None else None
                self.smoothed_P_to_Nf_ratio = self.P / Nf if self.P is not None else None
                self.smoothed_Si_to_Nf_ratio = self.Si / Nf if self.Si is not None else None

    # -- Couplage vertical aux flocs mineraux ------------------------------------------
    # Mutualise entre DOM, Detritus et Heterotrophs : meme formulation, seul l'ensemble
    # des devises transportees change d'un composant a l'autre.

    def _attach_aggregate(self, coupled_aggregate):
        """Wire the mineral aggregate driving the vertical dynamics, and cache its two
        coupling parameters locally (read at every timestep otherwise).

        Returns the aggregate, real or prescribed from Setup (BGC-only runs).
        """
        if self.prescribe_aggregate_from_setup:
            from ..components.flocs import PrescribedFlocs
            self.coupled_aggregate = PrescribedFlocs(
                name="Macroflocs",
                resusp_ewma_alpha=self.prescribed_resusp_ewma_alpha,
                organomin_coupling_fraction=self.prescribed_organomin_coupling_fraction
            )
        else:
            self.coupled_aggregate = coupled_aggregate

        if self.coupled_aggregate is not None:
            self.resusp_ewma_alpha = self.coupled_aggregate.resusp_ewma_alpha
            self.organomin_coupling_fraction = self.coupled_aggregate.organomin_coupling_fraction
        else:
            self.resusp_ewma_alpha = 0.0
            self.organomin_coupling_fraction = 1.0
        return self.coupled_aggregate

    def _update_prescribed_aggregate(self, t_idx):
        """Refresh the prescribed aggregate from the Setup arrays (no-op otherwise)."""
        if self.prescribe_aggregate_from_setup and self.coupled_aggregate:
            self.coupled_aggregate.sink_sedimentation = self.setup.Macroflocs_sink_sed_array[t_idx]
            self.coupled_aggregate.source_resuspension = self.setup.Macroflocs_source_resusp_array[t_idx]
            self.coupled_aggregate.numconc = self.setup.Macroflocs_numconc_array[t_idx]

    def get_sink_vertical_loss(self):
        """Vertical loss coupled to mineral floc dynamics (sedimentation - resuspension).

        Separated formulation:
        - Sedimentation: proportional to current concentration in water column
        - Resuspension: absolute flux based on smoothed BGC/floc ratio (EWMA filter)

        Applies to the currencies the component actually carries; the others are set to 0.
        """
        currencies = ('C', 'N', 'P', 'Si')
        if self.coupled_aggregate is None:
            for cur in currencies:
                setattr(self.sink_vertical_loss, cur, 0.)
            return

        conv = self.coupled_aggregate.time_conversion_factor
        Nf = self.coupled_aggregate.numconc

        # Update smoothed ratios (EWMA filter: alpha=0 -> fixed, alpha>0 -> adaptive)
        # Ratios based on fraction forming organo-mineral aggregates
        if Nf > 0 and self.resusp_ewma_alpha > 0:
            for cur in currencies:
                pool = getattr(self, cur)
                ratio = getattr(self, f'smoothed_{cur}_to_Nf_ratio')
                if pool is not None and ratio is not None:
                    setattr(self, f'smoothed_{cur}_to_Nf_ratio',
                            self.resusp_ewma_alpha * (pool * self.organomin_coupling_fraction / Nf)
                            + (1 - self.resusp_ewma_alpha) * ratio)

        # Sedimentation rate [d-1]
        settling_rate = (self.coupled_aggregate.sink_sedimentation / Nf * conv) if Nf > 0 else 0.0

        # Net vertical loss (positive = loss from water column). Sedimentation applies
        # only to the fraction forming organo-mineral aggregates; resuspension is an
        # absolute flux [mmol m-3 d-1] rebuilt from the smoothed ratio.
        for cur in currencies:
            pool = getattr(self, cur)
            if pool is None:
                setattr(self.sink_vertical_loss, cur, 0.)
                continue
            ratio = getattr(self, f'smoothed_{cur}_to_Nf_ratio')
            resusp = (self.coupled_aggregate.source_resuspension * conv * ratio) \
                if ratio is not None else 0.0
            setattr(self.sink_vertical_loss, cur,
                    settling_rate * pool * self.organomin_coupling_fraction - resusp)

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
