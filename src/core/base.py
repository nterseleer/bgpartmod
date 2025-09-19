import numpy as np
from ..utils import functions as fns

from src.config_model import varinfos


class BaseStateVar:
    def __init__(self):
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
        return np.array([fns.get_nested_attr(self, diag) for diag in self.diagnostics])


class BaseOrg(BaseStateVar):
    def __init__(self):

        super().__init__()

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

        self.ICs = [pool for pool in [self.C, self.N, self.Chl, self.P, self.Si] if pool is not None]

    def update_val(self, C,
                   N=None,
                   Chl=None,
                   P=None,
                   Si=None,
                   t=None,
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
