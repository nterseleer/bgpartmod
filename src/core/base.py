import operator
import numpy as np
from ..utils import functions as fns


# Currencies that can follow the mineral flocs, and the attribute holding the smoothed
# <currency>/N_F ratio of each (public names: diagnostics and vertical_coupling_state)
VERTICAL_CURRENCIES = ('C', 'N', 'Chl', 'P', 'Si')
SMOOTHED_RATIO = {cur: f'smoothed_{cur}_to_Nf_ratio' for cur in VERTICAL_CURRENCIES}


class BaseStateVar:
    def __init__(self, dtype=np.float64):
        self.dtype = dtype
        self.setup = None
        self.debug_mode = False
        self.time_conversion_factor = 1
        self.diagnostics = None

    def get_diagnostic_variables(self):
        """Current values of the diagnostics, in the order of self.diagnostics.

        Read through an operator.attrgetter (it resolves dotted names like 'source_PP.C'),
        built once per diagnostics list: called at every timestep.
        """
        if getattr(self, '_diag_reader_for', None) is not self.diagnostics:
            self._diag_reader = operator.attrgetter(*self.diagnostics)
            self._diag_reader_for = self.diagnostics
        values = self._diag_reader(self)
        return np.array(values if len(self.diagnostics) > 1 else (values,), dtype=self.dtype)

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

        # Vertical coupling to the mineral flocs (see _attach_aggregate). The prescribed_*
        # values are overwritten by the subclasses that expose them as parameters.
        self.coupled_aggregate = None  # Link to Macroflocs for vertical dynamics
        self.prescribe_aggregate_from_setup = False
        self.prescribed_resusp_ewma_alpha = 0.0
        self.prescribed_organomin_coupling_fraction = 1.0
        self.resusp_ewma_alpha = 0.0  # [-] cached from coupled_aggregate
        self.organomin_coupling_fraction = 1.0  # [-] cached from coupled_aggregate
        self.sink_vertical_loss = Elms()  # [mmol m-3 d-1] None (NaN at t0) until first computed
        self._vertical_ratios = ()  # (currency, ratio attribute) of the pools carried
        self._vertical_absent = VERTICAL_CURRENCIES  # currencies not carried: loss set to 0
        self._smoothed_ratios_pending = False  # True while an aggregate ratio is still unset

        # Smoothed ratios for vertical coupling (set by initialize_vertical_coupling_ratios)
        # Controlled by coupled_aggregate.resusp_ewma_alpha
        self.smoothed_C_to_Nf_ratio = None
        self.smoothed_N_to_Nf_ratio = None
        self.smoothed_Chl_to_Nf_ratio = None
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
            if self.classname == 'Phy':
                if self.QratioIC:
                    Qmaxratio = self.QratioIC
                    self.N = self.C * self.QN_max * Qmaxratio
            self.QN = self.N / self.C
            if self.classname == 'Phy' and self.checkQmax  and self.QN > self.QN_max:
                self.N = self.C * self. QN_max * Qmaxratio
                print('Phytoplankton N initial pool too high compared to QN_max, changed from {} to {}'.format(N, self.N))
                self.QN = self.N / self.C
        if self.P is not None:
            if self.classname == 'Phy':
                self.P = self.C * self.QP_max * Qmaxratio
            self.QP = self.P / self.C
            if self.classname == 'Phy' and self.checkQmax  and self.QP > self.QP_max:
                self.P = self.C * self. QP_max * Qmaxratio
                print('Phytoplankton P initial pool too high compared to QP_max, changed from {} to {}'.format(P, self.P))
                self.QP = self.P / self.C
        if self.Si is not None:
            if self.classname == 'Phy':
                self.Si = self.C * self.QSi_max * Qmaxratio
            self.QSi = self.Si / self.C
            if self.classname == 'Phy' and self.checkQmax  and self.QSi > self.QSi_max:
                self.Si = self.C * self. QSi_max * Qmaxratio
                print('Phytoplankton Si initial pool too high compared to QSi_max, changed from {} to {}'.format(Si, self.Si))
                self.QSi = self.Si / self.C
        if self.P is not None and self.Si is not None:
            self.fnut = min(self.QN, self.QP, self.QSi)
        if self.Chl is not None:
            self.thetaC = self.Chl / self.C #

        self.ICs = np.array([pool for pool in [self.C, self.N, self.Chl, self.P, self.Si] if pool is not None], dtype=self.dtype)

    # -- Couplage vertical aux flocs mineraux ------------------------------------------
    # Mutualise entre Phyto, DOM, Detritus et Heterotrophs : meme formulation, seul
    # l'ensemble des devises transportees change d'un composant a l'autre. Protocole :
    # set_coupling -> _attach_aggregate ; Model -> initialize_vertical_coupling_ratios
    # (apres les ICs) ; get_sinks -> get_sink_vertical_loss(t_idx), somme dans *_sinks.

    def set_vertical_coupling_state(self, state_dict):
        """
        Set smoothed BGC/floc ratios from extracted state (e.g., from a previous simulation).

        Args:
            state_dict: Dict with keys like 'smoothed_C_to_Nf_ratio', 'smoothed_N_to_Nf_ratio', etc.
        """
        unknown = set(state_dict) - set(SMOOTHED_RATIO.values())
        if unknown:
            raise TypeError(f"{self.name}: unknown vertical_coupling_state keys {sorted(unknown)}, "
                            f"expected some of {list(SMOOTHED_RATIO.values())}")
        for attr, value in state_dict.items():
            setattr(self, attr, value)

    def initialize_vertical_coupling_ratios(self):
        """
        Initialize BGC/floc ratios for resuspension after ICs are set.
        Respects pre-set values from set_vertical_coupling_state().

        With resusp_ewma_alpha:
        - α=0 (default): Fixed ratio from initial conditions (Option A)
        - α>0: Exponentially weighted moving average (Option B1)

        A ratio that cannot be set yet (no flocs, N_F = 0) is set at the first timestep
        where N_F > 0 (see get_sink_vertical_loss).

        ECART CONNU, NON CORRIGE : le ratio initial vaut pool / N_F alors que le filtre
        EWMA le fait tendre vers pool * organomin_coupling_fraction / N_F. Avec une
        fraction < 1, la resuspension est donc surestimee au debut de la simulation (le
        temps que le filtre oublie le ratio initial), et pendant toute la simulation avec
        alpha = 0. Le corriger change les trajectoires et casserait la comparabilite avec
        les optimisations en cours (OPT639 : fraction 0.154).
        """
        if self.coupled_aggregate is None:
            return
        self._vertical_ratios = tuple((cur, SMOOTHED_RATIO[cur]) for cur in VERTICAL_CURRENCIES
                                      if getattr(self, cur) is not None)
        self._vertical_absent = tuple(cur for cur in VERTICAL_CURRENCIES if getattr(self, cur) is None)
        # A prescribed aggregate holds no value before its first refresh
        self._update_prescribed_aggregate(0)
        self._smoothed_ratios_pending = True
        self._set_initial_smoothed_ratios(self.coupled_aggregate.numconc)

    def _set_initial_smoothed_ratios(self, Nf):
        """Seed the smoothed ratios still unset from the current pools (see the known gap
        in initialize_vertical_coupling_ratios)."""
        if Nf <= 0:
            return
        for cur, ratio_name in self._vertical_ratios:
            if getattr(self, ratio_name) is None:
                setattr(self, ratio_name, getattr(self, cur) / Nf)
        self._smoothed_ratios_pending = False

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

    def get_sink_vertical_loss(self, t_idx):
        """Vertical loss coupled to mineral floc dynamics (sedimentation - resuspension).

        Separated formulation:
        - Sedimentation: proportional to current concentration in water column
        - Resuspension: absolute flux based on smoothed BGC/floc ratio (EWMA filter)

        Applies to the currencies the component actually carries; the others are set to 0,
        as are all of them without aggregate. Positive = loss from the water column.
        """
        loss = self.sink_vertical_loss
        for cur in self._vertical_absent:
            setattr(loss, cur, 0.)
        aggregate = self.coupled_aggregate
        if aggregate is None:
            return
        self._update_prescribed_aggregate(t_idx)

        Nf = aggregate.numconc
        if self._smoothed_ratios_pending:
            self._set_initial_smoothed_ratios(Nf)
        conv = aggregate.time_conversion_factor
        alpha = self.resusp_ewma_alpha
        fraction = self.organomin_coupling_fraction

        # Sedimentation applies only to the fraction forming organo-mineral aggregates
        settling_rate = (aggregate.sink_sedimentation / Nf * conv) if Nf > 0 else 0.0  # [d-1]
        resusp_rate = aggregate.source_resuspension * conv  # [# m-3 d-1]
        update_ratios = Nf > 0 and alpha > 0

        for cur, ratio_name in self._vertical_ratios:
            pool = getattr(self, cur)
            ratio = getattr(self, ratio_name)
            if ratio is None:  # N_F = 0 since the start: no deposit to resuspend yet
                setattr(loss, cur, settling_rate * pool * fraction)
                continue
            # EWMA filter of the ratio (alpha = 0 -> fixed, alpha > 0 -> adaptive)
            if update_ratios:
                ratio = alpha * (pool * fraction / Nf) + (1 - alpha) * ratio
                setattr(self, ratio_name, ratio)
            # Resuspension: absolute flux [mmol m-3 d-1] rebuilt from the smoothed ratio
            setattr(loss, cur, settling_rate * pool * fraction - resusp_rate * ratio)

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
