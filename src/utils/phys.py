"""
Physical setup configuration for BGC model Simulations.
Defines the physical and computational environment for model runs.
"""
import os.path
from dataclasses import dataclass

import numpy as np
import pandas as pd
from typing import Dict, Any, Optional, Set

from src.Config_system import path_config as path_cfg

DEFAULT_SETUPS = {
    'schartau07': {
        'T': 12.5,
        'I': 26.,  # 115*3600.*24
        'light_prop': 14. / 24.,
    },
    'onur22': {
        'T': 10.,
        'I': 0.5,  # 115*3600.*24
        'light_prop': 14. / 24.,

    },
    'youri_high': {
        'T': 18.,
        'PARfromfile': False,
        'I': 38.394 * 86400,
        'light_prop': 12. / 24.,
    },
    'youri_medium': {
        'T': 12.,
        'PARfromfile': False,
        'I': 12.798 * 86400,
        'light_prop': 12. / 24.,
    },

    'youri_medium_mean_strain': {
        'T': 15.,
        'PARfromfile': False,
        'I': 12.798 * 86400,
        'light_prop': 12. / 24.,
    },
    # auria medium light = 15 µmol photons.m-2.s-1
    'Odontella_medium': {
        'T': 19.,
        'PARfromfile': False,
        'I': 15. * 86400, #µmol photon m-2 d-1 of PAR
        'light_prop': 12. / 24.,
    },

    'youri_low': {
        'T': 4.,
        'PARfromfile': False,
        'I': 4.266 * 86400,
        'light_prop': 12. / 24.,
    }
}


@dataclass
class PhysicalConstants:
    """Physical constants used in the model"""
    degCtoK: float = 273.15  # Conversion from Celsius to Kelvin
    molmass_C: float = 12.0107  # Molar mass of carbon [g/mol]
    day_to_seconds: float = 86400  # Seconds in a day
    PAR_conversion: float = 0.5 / 0.0079 / 54  # PAR/Wm-2 conversion factor


class Setup:
    """
    Physical and computational setup for BGC model Simulations.

    Attributes are organized into categories:
    - Time settings: control simulation timespan and steps
    - Physical parameters: define environmental conditions
    - Light settings: control light/PAR handling
    - Computational settings: control model execution
    """

    # Core physical attributes that define the setup state
    PHYSICAL_ATTRIBUTES = {
        'time_settings': [
            'tmin',  # Start time in days
            'tmax',  # End time in days
            'dt',  # Primary timestep in days
            'dt2',  # Optional secondary timestep in days
            'dt_ratio',  # Ratio between primary and secondary timesteps
            'used_dt',  # Actually used timestep
        ],
        'physical_params': [
            'T',  # Temperature in Kelvin
            'z',  # Depth in meters
            'k_att',  # Light attenuation coefficient
            'kb',  # Background turbidity
            'pCO2',  # Partial pressure CO2 in µatm
        ],
        'light_settings': [
            'light_prop',  # Proportion of time with light
            'PARfromfile',  # Whether PAR is loaded from file
            'lightfirst',  # Whether light period is at start of day
        ],
        'shear_settings': [
            'g_shear_rate',  # Base shear rate
            'vary_g_shear',  # Whether shear rate varies
            'gshearfact',  # Factor for shear variation
            'gshearper',  # Period of shear variation
        ]
    }

    # Variables that are computed/updated during simulation
    DIAGNOSTIC_VARIABLES = [
        'Chl_tot',  # Total chlorophyll
        'Cphy_tot',  # Total phytoplankton carbon
    ]

    def __init__(self,
                 name: str = '',
                 tmin: float = 0,
                 tmax: float = 20,
                 dt: float = 0.001,
                 dt2: Optional[float] = None,
                 dt2_s_to_d_ratio: float = 3600 * 24,
                 start_date: str = '2023-02-01 10:00:00',
                 PARfromfile: bool = False,
                 I: float = 100. * 3600. * 24.,
                 light_prop: float = 0.5,
                 lightfirst: bool = True,
                 T: float = 18.,
                 k_att: float = 16.,
                 z: float = 0.25,
                 pCO2: float = 370,
                 g_shear_rate: float = 95,
                 vary_g_shear: bool = True,
                 gshearfact: float = 0.5,
                 gshearper: float = 0.5,
                 kb: float = 0.13,
                 plotPAR: bool = False,
                 verbose: bool = False):
        """Initialize physical setup for simulation."""
        # Store initialization parameters
        self.constants = PhysicalConstants()
        self.name = name
        self.verbose = verbose

        # Store time parameters
        self.tmin = tmin
        self.tmax = tmax
        self.dt = dt
        self.dt2 = dt2
        self.dt2_s_to_d_ratio = dt2_s_to_d_ratio
        self.start_date = pd.to_datetime(start_date)
        self.end_date = self.start_date + pd.Timedelta(days=tmax)
        self.t_span = [tmin, tmax]

        # Initialize basic timestep variables
        if self.dt2 is not None:
            self.used_dt = self.dt2
            self.dt_ratio = round(self.dt / self.dt2)
            if self.dt_ratio < 1:
                raise ValueError(f"dt2 ({self.dt2}) must be smaller than dt ({self.dt})")
            self.two_dt = True
        else:
            self.used_dt = self.dt
            self.dt_ratio = None
            self.two_dt = False

        # Create time arrays
        self.t_eval = np.arange(tmin, tmax, self.used_dt)
        self.dates = pd.date_range(start=self.start_date, end=self.end_date, periods=len(self.t_eval))

        # Now initialize two timestep specific arrays
        if self.two_dt:
            self.dates1 = self.dates[::self.dt_ratio]
            self.dates1_set = set(self.dates1)
        else:
            self.dates1 = None
            self.dates1_set = None

        # Store physical parameters
        self.T = T + self.constants.degCtoK  # Convert to Kelvin
        self.k_att = k_att
        self.z = z
        self.pCO2 = pCO2
        self.kb = kb

        # Light settings
        self.PARfromfile = PARfromfile
        self.I = I
        self.light_prop = light_prop
        self.lightfirst = lightfirst
        self.PAR = self._initialize_PAR(plotPAR)

        # Shear settings
        self.vary_g_shear = vary_g_shear
        self.gshearfact = gshearfact
        self.gshearper = gshearper
        self.g_shear_rate = pd.DataFrame(
            self._calculate_shear_rates(g_shear_rate),
            index=self.dates,
            columns=['ShearRate']
        )

        # Initialize diagnostic variables
        for var in self.DIAGNOSTIC_VARIABLES:
            setattr(self, var, 0)

        # For two timestep handling
        self.dates1_set: Optional[Set] = None
        self.dates1: Optional[pd.DatetimeIndex] = None
        self.two_dt: bool = dt2 is not None

    def _initialize_PAR(self, plotPAR: bool) -> pd.DataFrame:
        """Initialize Photosynthetically Active Radiation data."""
        if not self.PARfromfile:
            if self.lightfirst:
                light_starts = np.linspace(self.tmin, self.tmax - 1, self.tmax - self.tmin)
                dark_starts = np.linspace(
                    self.tmin + self.light_prop,
                    self.tmax - 1 + self.light_prop,
                    self.tmax - self.tmin
                )
            else:
                firstL = (1 - self.light_prop)
                light_starts = np.linspace(
                    self.tmin + firstL,
                    self.tmax + firstL - 1,
                    self.tmax - self.tmin
                )
                dark_starts = np.linspace(self.tmin + 1, self.tmax, self.tmax - self.tmin)

            light = np.array([
                np.any((light_starts <= t) & (t <= dark_starts))
                for t in self.t_eval
            ])
            return pd.DataFrame(self.I * light, index=self.dates, columns=['PAR'])  # Using self.I
        else:
            return self._load_PAR_from_file(plotPAR)

    def _load_PAR_from_file(self, plotPAR: bool) -> pd.DataFrame:
        """Load PAR data from file."""
        solrad_clim = pd.read_csv(os.path.join(path_cfg.DATA_DIR, 'solrad_clim.dat'), sep='\s+',
                                  header=None,
                                  names=['DOY', 'Hour', 'PAR1', 'PAR2', 'PAR3'])

        # Convert to PAR
        solrad_clim['PAR'] = (self.constants.PAR_conversion *
                              solrad_clim['PAR1'])  # µmol photon m-2 s-1 of PAR
        solrad_clim['PAR'] = (solrad_clim['PAR'] *
                              self.constants.day_to_seconds)  # µmol photon m-2 d-1 of PAR

        # Create datetime index
        solrad_clim['Date'] = (
                pd.to_datetime(self.start_date.year * 1000 + solrad_clim['DOY'],
                               format='%Y%j') +
                pd.to_timedelta(solrad_clim['Hour'], unit='hours')
        )
        solrad_clim.set_index('Date', inplace=True)

        # Filter and interpolate
        filtered_df = solrad_clim[self.start_date:self.end_date]
        newdf = pd.DataFrame(index=self.dates)
        combined_df = filtered_df.join(newdf, how='outer')
        combined_df['PAR'] = combined_df['PAR'].interpolate(method='time')

        if plotPAR:
            self._plot_PAR(filtered_df, combined_df)

        return combined_df.loc[self.dates, ['PAR']]

    def _plot_PAR(self, original_df: pd.DataFrame, interpolated_df: pd.DataFrame):
        """Plot PAR data for verification."""
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        plt.plot(interpolated_df.index, interpolated_df['PAR'],
                 label='PAR in model')
        plt.scatter(original_df.index, original_df['PAR'],
                    color='red', label='Original PAR', s=3)
        plt.xlabel('Date')
        plt.ylabel('PAR')
        plt.title('PAR Values Over Time')
        plt.legend()

    def _calculate_shear_rates(self, base_rate: float) -> np.ndarray:
        """Calculate shear rates for the simulation period."""
        if not self.vary_g_shear:
            return np.full_like(self.t_eval, base_rate)

        return np.array([
            base_rate + self.gshearfact * np.cos(t / self.gshearper * 2. * np.pi) * base_rate
            for t in self.t_eval
        ])

    def get_physical_state(self) -> Dict[str, Any]:
        """
        Get the current physical state as a dictionary.

        Returns:
            Dictionary of physical attributes and their values,
            organized by category.
        """
        state = {}
        for category, attrs in self.PHYSICAL_ATTRIBUTES.items():
            state[category] = {
                attr: getattr(self, attr)
                for attr in attrs
                if hasattr(self, attr)
            }
        return state

    def get_diagnostic_state(self) -> Dict[str, float]:
        """
        Get current values of diagnostic variables.

        Returns:
            Dictionary of diagnostic variables and their values.
        """
        return {
            var: getattr(self, var)
            for var in self.DIAGNOSTIC_VARIABLES
        }

    def summarize(self) -> str:
        """
        Create a human-readable summary of the setup.

        Returns:
            Multi-line string describing key setup parameters.
        """
        lines = [f"Setup: {self.name}", "-" * 40]

        for category, attrs in self.PHYSICAL_ATTRIBUTES.items():
            lines.append(f"\n{category.replace('_', ' ').title()}:")
            for attr in attrs:
                if hasattr(self, attr):
                    value = getattr(self, attr)
                    if isinstance(value, (np.ndarray, pd.Series, pd.DataFrame)):
                        lines.append(f"  {attr}: array[{len(value)}]")
                    else:
                        lines.append(f"  {attr}: {value}")

        return "\n".join(lines)
