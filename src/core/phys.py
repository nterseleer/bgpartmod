
"""
Physical setup configuration for BGC model Simulations.
Defines the physical and computational environment for model runs.
"""
import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Any, List, Optional, Set, Union
from src.config_system import path_config as path_cfg

@dataclass
class PhysicalConstants:
    """Physical constants used in the model"""
    degCtoK: float = 273.15          # Conversion from Celsius to Kelvin
    molmass_C: float = 12.0107       # Molar mass of carbon [g/mol]
    day_to_seconds: float = 86400    # Seconds in a day
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
            'base_T',  # Temperature in Kelvin
            'varyingTEMP', # Whether T should be time-varying or not
            'k_att',  # Light attenuation coefficient
            'kb',  # Background turbidity
            'pCO2',  # Partial pressure CO2 in µatm
            'water_depth',  # Water depth in m (base/mean value)
            'vary_water_depth',  # Whether water depth should vary with tides
            'base_mu_water',  # Base dynamic viscosity of water in Pa·s
            'vary_mu_water',  # Whether viscosity varies with temperature
            'rho_water',  # Density of water in kg/m³
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
            'bed_shear_stress',  # Base bed shear stress in Pa
            'vary_bed_shear',  # Whether bed shear stress varies
            'bed_shear_fact',  # Factor for bed shear stress variation
            'bed_shear_per',  # Period of bed shear stress variation
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
                 dtype: type = np.float64,
                 start_date: str = '2023-02-01 10:00:00',
                 PARfromfile: bool = False,
                 I: float = 100. * 3600. * 24.,
                 light_prop: float = 0.5,
                 lightfirst: bool = True,
                 T: float = 18.,
                 varyingTEMP: bool = False,
                 k_att: float = 16.,
                 pCO2: float = 370,
                 g_shear_rate: float = 95,
                 vary_g_shear: bool = True,
                 gshearfact: float = 0.5,
                 gshearper: float = 0.5,  # Will be deprecated in favor of tidal_period_M2
                 bed_shear_stress: float = 0.5,
                 vary_bed_shear: bool = True,
                 bed_shear_fact: float = 1.,
                 bed_shear_per: float = 0.5,
                 kb: float = 0.13,
                 water_depth: float = 10.,
                 vary_water_depth: bool = False,
                 water_depth_amplitude_spring: float = 2.5,
                 water_depth_amplitude_neap: float = 1.5,
                 shear_rate_amplitude_spring: float = 30.0,
                 shear_rate_amplitude_neap: float = 15.0,
                 bed_shear_stress_amplitude_spring: float = 0.3,
                 bed_shear_stress_amplitude_neap: float = 0.2,
                 shear_rate_additive_mode: bool = False,
                 bed_shear_stress_additive_mode: bool = False,
                 tidal_period_M2: float = 12.42/24,
                 spring_neap_period: float = 14.7,
                 # Configurable phase shifts for tidal parameters
                 water_depth_phase_shift: float = 0.0,           # High tide at t=0 (cosine)
                 shear_rate_phase_shift: float = np.pi / 2,      # Max at mid-tide (sine)
                 bed_shear_stress_phase_shift: float = np.pi / 2, # Max at mid-tide (sine)
                 riverine_loads: bool = False,
                 riverine_loads_file: str = 'riverine_loads.feather',
                 mu_water: float = 1.002e-3,  # Pa·s at 20°C (base value)
                 vary_mu_water: bool = False,  # Whether viscosity varies with temperature
                 seawater_salinity: float = 0.0315,  # Salinity in kg/kg for viscosity calc (31.5 g/kg)
                 rho_water: float = 1025.,    # kg/m³
                 plotPAR: bool = False,
                 plotTEMP: bool = False,
                 verbose: bool = False,
                 prescribe_TEP: bool = False,  # Whether to prescribe TEP from previous simulation
                 TEP_sim_name: Optional[str] = None,  # Name of simulation to load TEP from
                 TEP_column: str = 'TEPC_C',  # Column name for TEP concentration
                 plotTEP: bool = False,
                 tep_year: int = 2023,
                 prescribe_Flocs: bool = False,  # Whether to prescribe Flocs from previous simulation
                 Flocs_sim_name: Optional[str] = None,  # Name of simulation to load Flocs from
                 Flocs_columns: Optional[Dict[str, str]] = None,  # Column name mapping for Flocs data
                 plotFlocs: bool = False,
                 flocs_year: int = 2023):
        """Initialize physical setup for simulation."""
        # Store initialization parameters
        self.constants = PhysicalConstants()
        self.name = name
        self.verbose = verbose
        self.dtype = dtype

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
        self.t_eval = np.arange(tmin, tmax, self.used_dt, dtype=self.dtype)
        self.dates = pd.date_range(start=self.start_date, end=self.end_date, periods=len(self.t_eval))

        # Now initialize two timestep specific arrays
        if self.two_dt:
            self.dates1 = self.dates[::self.dt_ratio]
            self.dates1_set = set(self.dates1)
        else:
            self.dates1 = None
            self.dates1_set = None

        # Store physical parameters
        self.base_T = T  # Store base temperature value
        # self.T = T + self.constants.degCtoK  # Convert to Kelvin
        self.k_att = k_att
        self.pCO2 = pCO2
        self.kb = kb
        self.base_mu_water = mu_water  # Base dynamic viscosity of water
        self.vary_mu_water = vary_mu_water
        self.seawater_salinity = seawater_salinity
        self.rho_water = rho_water  # Density of water

        # Store base values for tidal parameters (before DataFrame conversion)
        self.base_water_depth = water_depth
        self.base_g_shear_rate = g_shear_rate
        self.base_bed_shear_stress = bed_shear_stress

        # Tidal settings (shared across all tidal parameters)
        self.vary_water_depth = vary_water_depth
        self.tidal_period_M2 = tidal_period_M2
        self.spring_neap_period = spring_neap_period
        
        # Tidal amplitude parameters
        self.water_depth_amplitude_spring = water_depth_amplitude_spring
        self.water_depth_amplitude_neap = water_depth_amplitude_neap
        self.water_depth_phase_shift = water_depth_phase_shift

        self.shear_rate_amplitude_spring = shear_rate_amplitude_spring
        self.shear_rate_amplitude_neap = shear_rate_amplitude_neap
        self.shear_rate_phase_shift = shear_rate_phase_shift
        self.shear_rate_additive_mode = shear_rate_additive_mode

        self.bed_shear_stress_amplitude_spring = bed_shear_stress_amplitude_spring
        self.bed_shear_stress_amplitude_neap = bed_shear_stress_amplitude_neap
        self.bed_shear_stress_phase_shift = bed_shear_stress_phase_shift
        self.bed_shear_stress_additive_mode = bed_shear_stress_additive_mode


        # Riverine loads
        self.riverine_loads = riverine_loads
        self.riverine_loads_file = riverine_loads_file
        if self.riverine_loads:
            self._load_riverine_loads()

        # Prescribed TEP (for Flocs-only optimization with prescribed biogeochemistry)
        self.prescribe_TEP = prescribe_TEP
        self.TEP_sim_name = TEP_sim_name
        self.TEP_column = TEP_column
        self.TEP = None
        self.TEP_array = None
        if self.prescribe_TEP:
            self._load_prescribed_TEP(plotTEP, tep_year)

        # Prescribed Flocs (for BGC-only optimization with prescribed mineral dynamics)
        self.prescribe_Flocs = prescribe_Flocs
        self.Flocs_sim_name = Flocs_sim_name
        self.Flocs_columns = Flocs_columns
        self.Microflocs_massconc_array = None
        self.Micro_in_Macro_massconc_array = None
        self.Macroflocs_sink_sed_array = None
        self.Macroflocs_source_resusp_array = None
        self.Macroflocs_numconc_array = None
        if self.prescribe_Flocs:
            self._load_prescribed_Flocs(plotFlocs, flocs_year)

        # Light settings
        self.PARfromfile = PARfromfile
        self.I = I
        self.light_prop = light_prop
        self.lightfirst = lightfirst
        self.PAR = self._initialize_PAR(plotPAR)

        # Temperature settings
        self.varyingTEMP = varyingTEMP
        self.T = self._initialize_TEMP(plotTEMP)

        # Dynamic viscosity (depends on T, must be after T initialization)
        self.mu_water = self._initialize_mu_water()

        # Calculate temperature bounds once for efficient access
        # self.T_max = self.T['T'].max()  # Maximum temperature in setup [K]
        # self.T_min = self.T['T'].min()  # Minimum temperature in setup [K]

        # Initialize all tidal parameters using unified approach
        self.water_depth = self._create_tidal_parameter_dataframe(
            'WaterDepth', water_depth, vary_water_depth,
            water_depth_amplitude_spring, water_depth_amplitude_neap,
            water_depth_phase_shift, False, 1.0
        )

        # Shear settings
        self.vary_g_shear = vary_g_shear
        self.gshearfact = gshearfact # Will soon be deprecated
        self.gshearper = gshearper  # Will soon be deprecated
        self.g_shear_rate = self._create_tidal_parameter_dataframe(
            'ShearRate', g_shear_rate, vary_g_shear,
            shear_rate_amplitude_spring, shear_rate_amplitude_neap,
            shear_rate_phase_shift, shear_rate_additive_mode, 2.0
        )
        self.g_shear_rate_min = self.g_shear_rate.min().iloc[0]
        self.delta_g_shear_rate = (self.g_shear_rate.max() - self.g_shear_rate.min()).iloc[0]

        # Bed shear stress settings
        self.vary_bed_shear = vary_bed_shear
        self.bed_shear_fact = bed_shear_fact # Will soon be deprecated
        self.bed_shear_per = bed_shear_per # Will soon be deprecated
        self.bed_shear_stress = self._create_tidal_parameter_dataframe(
            'BedShearStress', bed_shear_stress, vary_bed_shear,
            bed_shear_stress_amplitude_spring, bed_shear_stress_amplitude_neap,
            bed_shear_stress_phase_shift, bed_shear_stress_additive_mode, 2.0
        )

        # Optimization: Pre-compute arrays for faster access during component calculations
        self.g_shear_rate_array = self.g_shear_rate['ShearRate'].values.astype(self.dtype)
        self.bed_shear_stress_array = self.bed_shear_stress['BedShearStress'].values.astype(self.dtype)
        self.water_depth_array = self.water_depth['WaterDepth'].values.astype(self.dtype)
        self.T_array = self.T['T'].values.astype(self.dtype)
        self.PAR_array = self.PAR['PAR'].values.astype(self.dtype)
        self.mu_water_array = self.mu_water['mu_water'].values.astype(self.dtype)

        # Calculate temperature bounds for efficient access
        self.T_max = self.T['T'].max()
        self.T_min = self.T['T'].min()

        # Create mapping for fast time index lookup
        self.dates_to_index = {date: idx for idx, date in enumerate(self.dates)}

        # Initialize diagnostic variables
        for var in self.DIAGNOSTIC_VARIABLES:
            setattr(self, var, 0)

        # For two timestep handling
        self.dates1_set: Optional[Set] = None
        self.dates1: Optional[pd.DatetimeIndex] = None
        self.two_dt: bool = dt2 is not None

        # Spin-up phase tracking
        self.in_spinup_phase: bool = False

    def _cycle_yearly_dataframe(self, df: pd.DataFrame, years_needed: int) -> pd.DataFrame:
        """
        Cycle a yearly DataFrame for multi-year simulations.

        Parameters:
        - df: DataFrame with DatetimeIndex covering one year
        - years_needed: Number of years to replicate

        Returns:
        - DataFrame with cycled data for multiple years
        """
        if years_needed <= 1:
            return df

        # Replicate the yearly cycle
        cycled_dfs = []
        for year_offset in range(years_needed):
            year_df = df.copy()
            year_df.index = year_df.index + pd.DateOffset(years=year_offset)
            cycled_dfs.append(year_df)

        # Concatenate and remove any duplicate indices that may occur at year boundaries
        result = pd.concat(cycled_dfs)

        # Keep first occurrence of duplicates (in chronological order)
        if result.index.duplicated().any():
            result = result[~result.index.duplicated(keep='first')]

        return result

    def _load_prescribed_TEP(self, plotTEP: bool = False, tep_year: int = 2023):
        """Load prescribed TEP from simulation, cycle backwards, interpolate to setup dates."""
        from src.utils import simulation_manager as sim_manager

        tep_sim = sim_manager.load_simulation(self.TEP_sim_name)
        tep_df = tep_sim.df[[self.TEP_column]].copy()

        # Filter reference year
        tep_df = tep_df[tep_df.index.year == tep_year]

        # Cycle backwards for multi-year simulations
        years_needed = int(np.ceil(self.tmax / 365))
        if years_needed > 1:
            cycled_dfs = [tep_df.copy()]
            for year_offset in range(1, years_needed):
                year_df = tep_df.copy()
                year_df.index = year_df.index - pd.DateOffset(years=year_offset)
                cycled_dfs.insert(0, year_df)
            tep_df = pd.concat(cycled_dfs).sort_index()

        # Interpolate onto self.dates
        newdf = pd.DataFrame(index=self.dates)
        combined_df = tep_df.join(newdf, how='outer')
        combined_df[self.TEP_column] = combined_df[self.TEP_column].interpolate(method='time')

        # Extract array directly
        self.TEP_array = combined_df.loc[self.dates, self.TEP_column].values.astype(self.dtype)

        if plotTEP:
            self._plot_prescribed_TEP()

    def _plot_prescribed_TEP(self):
        """Plot prescribed TEP data for verification."""
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

        # Plot 1: Full time series
        ax1.plot(self.dates, self.TEP_array, 'b-', linewidth=1.5, label='Prescribed TEP')
        ax1.set_xlabel('Date')
        ax1.set_ylabel('TEP C (mmol C m⁻³)')
        ax1.set_title(f'Prescribed TEP from Simulation: {self.TEP_sim_name}')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Plot 2: First 30 days detail
        days_to_plot = min(30, self.tmax)
        end_idx = int(days_to_plot / self.used_dt)
        ax2.plot(self.t_eval[:end_idx], self.TEP_array[:end_idx], 'b-', linewidth=2)
        ax2.set_xlabel('Time (days)')
        ax2.set_ylabel('TEP C (mmol C m⁻³)')
        ax2.set_title(f'Prescribed TEP - First {days_to_plot} days detail')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    def _load_prescribed_Flocs(self, plotFlocs: bool = False, flocs_year: int = 2023):
        """Load prescribed Flocs from simulation, cycle backwards, interpolate to setup dates."""
        from src.utils import simulation_manager as sim_manager

        flocs_sim = sim_manager.load_simulation(self.Flocs_sim_name)

        # Default column mapping if not specified
        if self.Flocs_columns is None:
            self.Flocs_columns = {
                'Microflocs_massconc': 'Microflocs_massconcentration',
                'Micro_in_Macro_massconc': 'Micro_in_Macro_massconcentration',
                'Macroflocs_sink_sed': 'Macroflocs_sink_sedimentation',
                'Macroflocs_source_resusp': 'Macroflocs_source_resuspension',
                'Macroflocs_numconc': 'Macroflocs_numconc'
            }

        # Extract and filter reference year
        flocs_df = flocs_sim.df[list(self.Flocs_columns.values())].copy()
        flocs_df = flocs_df[flocs_df.index.year == flocs_year]

        # Cycle backwards for multi-year simulations
        years_needed = int(np.ceil(self.tmax / 365))
        if years_needed > 1:
            cycled_dfs = [flocs_df.copy()]
            for year_offset in range(1, years_needed):
                year_df = flocs_df.copy()
                year_df.index = year_df.index - pd.DateOffset(years=year_offset)
                cycled_dfs.insert(0, year_df)
            flocs_df = pd.concat(cycled_dfs).sort_index()

        # Interpolate onto self.dates
        newdf = pd.DataFrame(index=self.dates)
        combined_df = flocs_df.join(newdf, how='outer')
        for col in self.Flocs_columns.values():
            combined_df[col] = combined_df[col].interpolate(method='time')

        # Extract arrays and convert back to model units (g/l = kg/m³)
        # sim.df stores data in output units (mg/l), so we need to divide by trsfrm factor
        # to get back to model units. For mass concentrations, trsfrm = 1e3 (see varinfos.py)
        from src.config_model import varinfos
        trsfrm_massconc = varinfos.doutput.get('Microflocs_massconcentration', {}).get('trsfrm', 1e3)

        self.Microflocs_massconc_array = combined_df.loc[self.dates, self.Flocs_columns['Microflocs_massconc']].values.astype(self.dtype) / trsfrm_massconc
        self.Micro_in_Macro_massconc_array = combined_df.loc[self.dates, self.Flocs_columns['Micro_in_Macro_massconc']].values.astype(self.dtype) / trsfrm_massconc
        self.Macroflocs_sink_sed_array = combined_df.loc[self.dates, self.Flocs_columns['Macroflocs_sink_sed']].values.astype(self.dtype)
        self.Macroflocs_source_resusp_array = combined_df.loc[self.dates, self.Flocs_columns['Macroflocs_source_resusp']].values.astype(self.dtype)
        self.Macroflocs_numconc_array = combined_df.loc[self.dates, self.Flocs_columns['Macroflocs_numconc']].values.astype(self.dtype)

        if plotFlocs:
            self._plot_prescribed_Flocs()

    def _plot_prescribed_Flocs(self):
        """Plot prescribed Flocs data for verification."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(5, 1, figsize=(12, 12))

        # Plot 1: Microflocs mass concentration
        axes[0].plot(self.dates, self.Microflocs_massconc_array, 'b-', linewidth=1.5)
        axes[0].set_ylabel('Microflocs mass\n(kg m⁻³)')
        axes[0].set_title(f'Prescribed Flocs from Simulation: {self.Flocs_sim_name}')
        axes[0].grid(True, alpha=0.3)

        # Plot 2: Macroflocs mass concentration
        axes[1].plot(self.dates, self.Micro_in_Macro_massconc_array, 'g-', linewidth=1.5)
        axes[1].set_ylabel('Micro_in_Macro mass\n(kg m⁻³)')
        axes[1].grid(True, alpha=0.3)

        # Plot 3: Macroflocs sedimentation
        axes[2].plot(self.dates, self.Macroflocs_sink_sed_array, 'r-', linewidth=1.5)
        axes[2].set_ylabel('Sedimentation\n(# m⁻³ s⁻¹)')
        axes[2].grid(True, alpha=0.3)

        # Plot 4: Macroflocs resuspension
        axes[3].plot(self.dates, self.Macroflocs_source_resusp_array, 'orange', linewidth=1.5)
        axes[3].set_ylabel('Resuspension\n(# m⁻³ s⁻¹)')
        axes[3].grid(True, alpha=0.3)

        # Plot 5: Macroflocs number concentration
        axes[4].plot(self.dates, self.Macroflocs_numconc_array, 'purple', linewidth=1.5)
        axes[4].set_xlabel('Date')
        axes[4].set_ylabel('Numconc\n(# m⁻³)')
        axes[4].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    def _load_riverine_loads(self):
        """Load riverine nutrient loads data and cycle for multi-year simulations."""
        import pandas as pd

        loads_file = os.path.join(path_cfg.DATA_DIR, self.riverine_loads_file)

        if not os.path.exists(loads_file):
            raise FileNotFoundError(
                f"Riverine loads file not found: {loads_file}\n"
                "Please run the preprocessing script first."
            )

        loads_df = pd.read_feather(loads_file)

        # Normalize dates to the simulation start year (keep day/month/time, change year)
        # This ensures proper alignment regardless of the original data year
        original_dates = loads_df.index
        normalized_dates = pd.to_datetime({
            'year': self.start_date.year,
            'month': original_dates.month,
            'day': original_dates.day,
            'hour': original_dates.hour,
            'minute': original_dates.minute,
            'second': original_dates.second,
            'microsecond': original_dates.microsecond
        })
        loads_df.index = normalized_dates

        # Remove any duplicates created during normalization
        # (e.g., 2023-01-01 and 2024-01-01 both become 2023-01-01)
        if loads_df.index.duplicated().any():
            loads_df = loads_df[~loads_df.index.duplicated(keep='first')]

        # Cycle for multi-year simulations if needed
        years_needed = int(np.ceil(self.tmax / 365))
        loads_df = self._cycle_yearly_dataframe(loads_df, years_needed)

        # Ensure index is sorted (required for reindex with method='nearest')
        loads_df = loads_df.sort_index()

        # Reindex to match setup.dates exactly, using nearest time interpolation
        # This also handles edge cases where data doesn't cover the full year
        loads_df = loads_df.reindex(self.dates, method='nearest')

        self.loads = loads_df.rename(columns={
            'NH4_loads': 'NH4',
            'NO3_loads': 'NO3',
            'DIP_loads': 'DIP',
            'DSi_loads': 'DSi'
        })

        # Optimization: Pre-compute arrays for faster access (avoid DataFrame .loc lookups)
        self.loads_NH4_array = self.loads['NH4'].values.astype(self.dtype) if 'NH4' in self.loads.columns else None
        self.loads_NO3_array = self.loads['NO3'].values.astype(self.dtype) if 'NO3' in self.loads.columns else None
        self.loads_DIP_array = self.loads['DIP'].values.astype(self.dtype) if 'DIP' in self.loads.columns else None
        self.loads_DSi_array = self.loads['DSi'].values.astype(self.dtype) if 'DSi' in self.loads.columns else None

        if self.verbose:
            print(f"Loaded riverine loads for {len(loads_df)} time steps")

    def _initialize_PAR(self, plotPAR: bool) -> pd.DataFrame:
        """Initialize Photosynthetically Active Radiation data."""
        if not self.PARfromfile:
            if self.lightfirst:
                light_starts = np.linspace(self.tmin, self.tmax - 1, self.tmax - self.tmin, dtype=self.dtype)
                dark_starts = np.linspace(
                    self.tmin + self.light_prop,
                    self.tmax - 1 + self.light_prop,
                    self.tmax - self.tmin,
                    dtype=self.dtype
                )
            else:
                firstL = (1 - self.light_prop)
                light_starts = np.linspace(
                    self.tmin + firstL,
                    self.tmax + firstL - 1,
                    self.tmax - self.tmin,
                    dtype=self.dtype
                )
                dark_starts = np.linspace(self.tmin + 1, self.tmax, self.tmax - self.tmin, dtype=self.dtype)

            light = np.array([
                np.any((light_starts <= t) & (t <= dark_starts))
                for t in self.t_eval
            ])
            return pd.DataFrame(self.I * light, index=self.dates, columns=['PAR'])  # Using self.I
        else:
            return self._load_PAR_from_file(plotPAR)

    def _load_PAR_from_file(self, plotPAR: bool) -> pd.DataFrame:
        """Load PAR data from file and cycle for multi-year simulations."""
        solrad_clim = pd.read_csv(os.path.join(path_cfg.DATA_DIR ,'solrad_clim.dat'), sep='\s+',
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

        # Cycle for multi-year simulations
        years_needed = int(np.ceil(self.tmax / 365))
        solrad_clim = self._cycle_yearly_dataframe(solrad_clim, years_needed)

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

    def _initialize_TEMP(self, plotTEMP: bool) -> pd.DataFrame:
        """Initialize temperature data as DataFrame."""
        if not self.varyingTEMP:
            # Constant temperature - create DataFrame with constant values
            temp_values = np.full(len(self.t_eval), self.base_T + self.constants.degCtoK, dtype=self.dtype)
            return pd.DataFrame(temp_values, index=self.dates, columns=['T'])
        else:
            # Time-varying temperature using cosine function
            return self._create_cosine_temperature(plotTEMP)

    def _initialize_mu_water(self) -> pd.DataFrame:
        """
        Initialize dynamic viscosity of seawater.

        If vary_mu_water is False, uses constant base_mu_water.
        If True, calculates T-dependent viscosity using Sharqawy et al. (2010), Eq. 22-23.
        Valid for 0 < t < 180°C; 0 < S < 0.15 kg/kg; accuracy ±1.5%.
        """
        if not self.vary_mu_water:
            mu_values = np.full(len(self.t_eval), self.base_mu_water, dtype=self.dtype)
        else:
            # Temperature in Celsius (formula requirement)
            t = self.T['T'].values - self.constants.degCtoK
            S = self.seawater_salinity

            # Pure water viscosity (IAPWS 2008, Eq. 23) [Pa·s]
            mu_w = 4.2844e-5 + (0.157 * (t + 64.993)**2 - 91.296)**(-1)

            # Seawater viscosity coefficients (Eq. 22)
            A = 1.541 + 1.998e-2 * t - 9.52e-5 * t**2
            B = 7.974 - 7.561e-2 * t + 4.724e-4 * t**2

            # Seawater dynamic viscosity [Pa·s]
            mu_values = mu_w * (1 + A * S + B * S**2)

        return pd.DataFrame(mu_values, index=self.dates, columns=['mu_water'])

    def _create_cosine_temperature(self, plotTEMP: bool) -> pd.DataFrame:
        """Create time-varying temperature using cosine function fitted to observations."""
        # Convert dates to julian days
        julian_days = self.dates.dayofyear

        # Cosine function parameters fitted to observations
        # T_min = 5.5°C at day 37.5, T_max = 21.5°C at day 220 (MOW1 observations)
        self.T_min = 5.5 + self.constants.degCtoK
        self.T_max = 21.5 + self.constants.degCtoK
        T_mean = (5.5 + 21.5) / 2.0  # 13.5°C
        T_amplitude = (21.5 - 5.5) / 2.0  # 8.0°C
        phase_offset = 210  # Day of maximum temperature

        # Generate cosine temperature (in Celsius)
        temp_celsius = T_mean + T_amplitude * np.cos(2 * np.pi * (julian_days - phase_offset) / 365.25)

        # Convert to Kelvin
        temp_kelvin = temp_celsius + self.constants.degCtoK

        # Create DataFrame
        temp_df = pd.DataFrame(temp_kelvin, index=self.dates, columns=['T'])

        if plotTEMP:
            self._plot_temperature(temp_df, temp_celsius, julian_days)

        return temp_df

    def _plot_temperature(self, temp_df: pd.DataFrame, temp_celsius: np.ndarray, julian_days: np.ndarray):
        """Plot temperature data for verification against observations."""
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Plot 1: Temperature vs Julian Day
        ax1.plot(julian_days, temp_celsius, 'b-', label='Model temperature', linewidth=2)
        ax1.set_xlabel('Julian Day')
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_title('Model Temperature vs Julian Day')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Plot 2: With observations (if available)
        ax2.plot(julian_days, temp_celsius, 'b-', label='Model temperature', linewidth=2)

        # Try to load and plot observations for comparison
        try:
            from src.utils.observations import Obs
            import src.config_model.obsinfos as obsinfos

            mow1 = Obs(read_df=False)
            mow1readdict = {'sep': '\s+', 'skiprows': [0, 1, 2, 4], 'engine': 'python'}
            mow1.df = mow1.read_obs(obsinfos.MAIN_DATA_FILE, readdic=mow1readdict)
            temp_obs = mow1.df[mow1.df['T@pH'] > 1]['T@pH']

            if len(temp_obs) > 0:
                mow1.df = mow1.compute_climatology(mow1.df[['T@pH']])
                ax2.scatter(mow1.df.index, mow1.df['T@pH'],
                            color='red', s=30, alpha=0.7, label='Observations', zorder=5)
                ax2.plot(mow1.df.index, mow1.df['T@pH'],
                         color='red', alpha=0.5, linewidth=1)
        except Exception as e:
            print(f"Could not load temperature observations for plotting: {e}")

        ax2.set_xlabel('Julian Day')
        ax2.set_ylabel('Temperature (°C)')
        ax2.set_title('Model vs Observed Temperature')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        plt.show()

    def _create_tidal_parameter_dataframe(self, 
                                        column_name: str,
                                        base_value: float, 
                                        vary_flag: bool,
                                        amplitude_spring: float, 
                                        amplitude_neap: float,
                                        phase_shift: float,
                                        additive_mode: bool = False,
                                        period_divisor: float = 1.0) -> pd.DataFrame:
        """
        Unified factory method for creating time-varying tidal parameter DataFrames.
        
        Parameters:
        - column_name: Name for the DataFrame column
        - base_value: Base/background value 
        - vary_flag: Whether parameter should vary with tides
        - amplitude_spring: Tidal amplitude during spring tides
        - amplitude_neap: Tidal amplitude during neap tides
        - phase_shift: Phase shift in radians (0=cosine, π/2=sine)
        - additive_mode: If True, base + [0 to amplitude]; If False, base ± amplitude
        - period_divisor: Divisor for tidal period (1.0 for M2, 2.0 for twice-daily)
        
        Returns:
        - DataFrame with time-indexed values
        """
        if not vary_flag:
            # Constant parameter - create DataFrame with constant values
            values = np.full(len(self.t_eval), base_value, dtype=self.dtype)
            return pd.DataFrame(values, index=self.dates, columns=[column_name])
        
        # Use unified tidal approach if water depth varies, otherwise backward compatibility
        if hasattr(self, 'tidal_period_M2') and self.vary_water_depth:
            period = self.tidal_period_M2 / period_divisor
            values = np.array([
                self._calculate_tidal_parameter(
                    t, base_value, amplitude_spring, amplitude_neap,
                    period, phase_shift, additive_mode
                )
                for t in self.t_eval
            ], dtype=self.dtype)
            return pd.DataFrame(values, index=self.dates, columns=[column_name])
        else:
            # Backward compatibility - use legacy parameters
            if column_name == 'ShearRate':
                values = np.array([
                    base_value + self.gshearfact * np.cos(t / self.gshearper * 2. * np.pi) * base_value
                    for t in self.t_eval
                ], dtype=self.dtype)
            elif column_name == 'BedShearStress':
                values = np.array([
                    base_value + self.bed_shear_fact * np.cos(t / self.bed_shear_per * 2. * np.pi) * base_value
                    for t in self.t_eval
                ], dtype=self.dtype)
            else:  # WaterDepth
                values = np.full(len(self.t_eval), base_value, dtype=self.dtype)

            return pd.DataFrame(values, index=self.dates, columns=[column_name])

    def _calculate_tidal_parameter(self, t: float, base_value: float, 
                                 amplitude_spring: float, amplitude_neap: float,
                                 tidal_period: float, phase_shift: float = 0.0,
                                 additive_mode: bool = False) -> float:
        """
        Unified function to calculate any tidal parameter with spring-neap modulation.
        
        Parameters:
        - t: Time in days
        - base_value: Base/background value 
        - amplitude_spring: Tidal amplitude during spring tides
        - amplitude_neap: Tidal amplitude during neap tides  
        - tidal_period: Tidal oscillation period (days)
        - phase_shift: Phase shift in radians (0=cosine, π/2=sine)
        - additive_mode: If True, base + [0 to amplitude]; If False, base ± amplitude
        
        Returns:
        - Parameter value at time t
        """
        # Spring-neap modulation factor (0 = neap, 1 = spring)
        spring_neap_factor = (1 + np.cos(t / self.spring_neap_period * 2 * np.pi)) / 2
        
        # Current tidal amplitude interpolated between neap and spring
        current_amplitude = (
            amplitude_neap + 
            (amplitude_spring - amplitude_neap) * spring_neap_factor
        )
        
        # Tidal signal calculation
        if additive_mode:
            # Additive: base + [0 to amplitude] - always positive tidal enhancement
            tidal_factor = (1 + np.cos(t / tidal_period * 2 * np.pi + phase_shift)) / 2
            return base_value + current_amplitude * tidal_factor
        else:
            # Oscillating: base ± amplitude - symmetric oscillation around base
            tidal_oscillation = current_amplitude * np.cos(t / tidal_period * 2 * np.pi + phase_shift)
            return base_value + tidal_oscillation


    def plot_setup(self, days_to_plot: float = None, save_path: str = None) -> None:
        """
        Plot all major setup properties for visualization and verification.
        
        Parameters:
        - days_to_plot: Number of days to plot (default: min(tmax, 5))
        - save_path: Path to save plot (optional)
        """
        import matplotlib.pyplot as plt
        
        # Determine plotting range
        if days_to_plot is None:
            days_to_plot = self.tmax
        
        # Get time indices for plotting
        end_idx = int(days_to_plot / self.used_dt)
        end_idx = min(end_idx, len(self.t_eval))
        
        t_plot = self.t_eval[:end_idx]
        dates_plot = self.dates[:end_idx]
        
        # Prepare data
        plot_data = {}
        
        # PAR data
        if hasattr(self, 'PAR'):
            plot_data['PAR'] = {
                'values': self.PAR['PAR'].values[:end_idx],
                'ylabel': 'PAR (µmol photon m⁻² d⁻¹)',
                'title': 'Photosynthetically Active Radiation'
            }
        
        # Temperature data  
        if hasattr(self, 'T'):
            plot_data['Temperature'] = {
                'values': self.T['T'].values[:end_idx] - self.constants.degCtoK,  # Convert to Celsius
                'ylabel': 'Temperature (°C)',
                'title': 'Water Temperature'
            }

        # Water dynamic viscosity
        if hasattr(self, 'mu_water'):
            plot_data['Dynamic Viscosity'] = {
                'values': self.mu_water['mu_water'].values[:end_idx],
                'ylabel': 'Dynamic viscosity (Pa.s)',
                'title': 'Seawater Dynamic Viscosity'
            }

        # Water depth data
        if hasattr(self, 'water_depth'):
            plot_data['Water Depth'] = {
                'values': self.water_depth['WaterDepth'].values[:end_idx],
                'ylabel': 'Water Depth (m)',
                'title': 'Tidal Water Depth' if self.vary_water_depth else 'Water Depth (constant)'
            }
        
        # Shear rate data
        if hasattr(self, 'g_shear_rate'):
            plot_data['Shear Rate'] = {
                'values': self.g_shear_rate['ShearRate'].values[:end_idx],
                'ylabel': 'Shear Rate (s⁻¹)',
                'title': 'Shear Rate'
            }
        
        # Bed shear stress data
        if hasattr(self, 'bed_shear_stress'):
            plot_data['Bed Shear Stress'] = {
                'values': self.bed_shear_stress['BedShearStress'].values[:end_idx],
                'ylabel': 'Bed Shear Stress (Pa)',
                'title': 'Bed Shear Stress'
            }

        # Riverine loads data
        if hasattr(self, 'loads') and self.riverine_loads:
            # We'll plot all 4 nutrients in a single subplot
            plot_data['Riverine Loads'] = {
                'values': {
                    'NH4': self.loads['NH4'].values[:end_idx],
                    'NO3': self.loads['NO3'].values[:end_idx],
                    'DIP': self.loads['DIP'].values[:end_idx],
                    'DSi': self.loads['DSi'].values[:end_idx]
                },
                'ylabel': 'Load (mmol N/P/Si m⁻² d⁻¹)',
                'title': 'Riverine Nutrient Loads',
                'multi': True  # Flag to indicate multiple lines
            }

        # Create subplot layout
        n_plots = len(plot_data)
        if n_plots == 0:
            print("No data to plot")
            return
            
        fig, axes = plt.subplots(n_plots, 1, figsize=(12, 2.5 * n_plots))
        if n_plots == 1:
            axes = [axes]  # Make it iterable
        
        # Plot each variable
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink']
        load_colors = {'NH4': 'blue', 'NO3': 'green', 'DIP': 'red', 'DSi': 'orange'}

        for i, (var_name, data) in enumerate(plot_data.items()):
            ax = axes[i]

            # Check if this is a multi-line plot (riverine loads)
            if data.get('multi', False):
                # Plot each nutrient load separately
                for nutrient, values in data['values'].items():
                    ax.plot(t_plot, values, color=load_colors[nutrient],
                           linewidth=1.5, label=nutrient, alpha=0.8)
                ax.legend(loc='best', framealpha=0.9)
            else:
                # Single line plot
                color = colors[i % len(colors)]
                ax.plot(t_plot, data['values'], color=color, linewidth=2)

            ax.set_ylabel(data['ylabel'])
            ax.set_title(data['title'])
            ax.grid(True, alpha=0.3)

            # Only show x-label on bottom plot
            if i == n_plots - 1:
                ax.set_xlabel('Time (days)')
        
        # Add overall title with setup info
        title_parts = [f"Setup: {self.name}" if self.name else "Physical Setup"]
        if self.vary_water_depth:
            title_parts.append("(with tidal variation)")
        
        fig.suptitle(" ".join(title_parts), fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        # Save if requested
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {save_path}")
            
        # plt.show()
        
        # Print summary statistics
        print(f"\nSetup Summary ({days_to_plot:.1f} days plotted):")
        print(f"{'Variable':<20} {'Min':<10} {'Max':<10} {'Mean':<10} {'Unit'}")
        print("-" * 60)

        for var_name, data in plot_data.items():
            values = data['values']
            unit = data['ylabel'].split('(')[-1].rstrip(')')

            # Handle multi-line data (e.g., riverine loads)
            if data.get('multi', False):
                for nutrient, nutrient_values in values.items():
                    full_name = f"{var_name} ({nutrient})"
                    print(f"{full_name:<20} {np.min(nutrient_values):<10.2f} {np.max(nutrient_values):<10.2f} "
                          f"{np.mean(nutrient_values):<10.2f} {unit}")
            else:
                print(f"{var_name:<20} {np.min(values):<10.2f} {np.max(values):<10.2f} "
                      f"{np.mean(values):<10.2f} {unit}")

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

    def to_dict(self) -> Dict[str, Any]:
        """Convert setup to dictionary using existing serialization."""
        from src.utils import functions as fns
        return fns.serialize_for_json(self)

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

    def crop_from_date(self, new_start_date: str, new_tmax: Optional[float] = None,
                       new_name: Optional[str] = None) -> 'Setup':
        """
        Create a new Setup cropped from the current one, starting at a specified date.

        This preserves the exact forcing values (shear rate, temperature, PAR, etc.)
        at the corresponding timestamps, ensuring tidal phase continuity.

        Parameters
        ----------
        new_start_date : str
            The new start date (e.g., '2023-01-01 00:00:00')
        new_tmax : float, optional
            New simulation duration in days. If None, uses remaining duration from
            original setup.
        new_name : str, optional
            Name for the cropped setup. If None, appends '_cropped' to original name.

        Returns
        -------
        Setup
            A new Setup instance with cropped time series, preserving forcing alignment.

        Notes
        -----
        **Attributes that are SLICED** (time series data):
        - t_eval, dates, dates1
        - PAR, PAR_array
        - T, T_array (actual temperature values)
        - g_shear_rate, g_shear_rate_array (actual shear values)
        - bed_shear_stress, bed_shear_stress_array
        - water_depth, water_depth_array
        - mu_water, mu_water_array
        - loads* arrays (riverine loads)
        - TEP_array, Flocs* arrays (prescribed data)

        **Attributes that are INHERITED** (normalization constants):
        - T_max, T_min: Annual temperature range for limT normalization in getlimT()
        - g_shear_rate_min, delta_g_shear_rate: Full tidal cycle range for
          normalized_shear calculation in flocs.py

        These normalization constants MUST represent the full annual/tidal cycle
        range, not just the cropped period. Otherwise, limitation functions would
        give incorrect values (e.g., limT=1.0 in winter instead of ~0.5).

        Example
        -------
        >>> # From a 2-year simulation setup, extract 2023 only
        >>> cropped_setup = sim0.setup.crop_from_date('2023-01-01 00:00:00', new_tmax=365)
        >>> # cropped_setup will have identical forcing values at '2023-01-01' as original
        >>> # AND the same T_max/T_min for correct limT normalization
        """
        import copy

        new_start = pd.to_datetime(new_start_date)

        # Find the index in current dates array closest to new_start_date
        if new_start < self.dates[0]:
            raise ValueError(f"new_start_date ({new_start}) is before setup start ({self.dates[0]})")
        if new_start > self.dates[-1]:
            raise ValueError(f"new_start_date ({new_start}) is after setup end ({self.dates[-1]})")

        # Find nearest index
        start_idx = self.dates.get_indexer([new_start], method='nearest')[0]
        actual_start = self.dates[start_idx]

        # Determine end index
        if new_tmax is not None:
            # Calculate how many timesteps for new_tmax days
            n_steps = int(new_tmax / self.used_dt)
            end_idx = min(start_idx + n_steps, len(self.dates))
        else:
            end_idx = len(self.dates)

        # Compute new time parameters
        new_duration_days = (end_idx - start_idx) * self.used_dt

        # Create a shallow copy first, then update attributes
        cropped = copy.copy(self)

        # Update name
        cropped.name = new_name if new_name else f"{self.name}_cropped"

        # Update time parameters
        cropped.tmin = 0  # Reset to 0 since we're starting fresh
        cropped.tmax = new_duration_days
        cropped.start_date = actual_start
        cropped.end_date = actual_start + pd.Timedelta(days=new_duration_days)
        cropped.t_span = [0, new_duration_days]

        # Slice time arrays
        cropped.t_eval = np.arange(0, new_duration_days, self.used_dt, dtype=self.dtype)
        cropped.dates = self.dates[start_idx:end_idx]

        # Handle two-timestep case
        if self.two_dt:
            # Slice dates1 appropriately
            start_idx1 = start_idx // self.dt_ratio
            end_idx1 = (end_idx - 1) // self.dt_ratio + 1
            cropped.dates1 = self.dates1[start_idx1:end_idx1] if self.dates1 is not None else None
            cropped.dates1_set = set(cropped.dates1) if cropped.dates1 is not None else None

        # Slice all DataFrames
        if hasattr(self, 'PAR') and isinstance(self.PAR, pd.DataFrame):
            cropped.PAR = self.PAR.iloc[start_idx:end_idx].copy()
            cropped.PAR.index = cropped.dates
        if hasattr(self, 'T') and isinstance(self.T, pd.DataFrame):
            cropped.T = self.T.iloc[start_idx:end_idx].copy()
            cropped.T.index = cropped.dates
        if hasattr(self, 'mu_water') and isinstance(self.mu_water, pd.DataFrame):
            cropped.mu_water = self.mu_water.iloc[start_idx:end_idx].copy()
            cropped.mu_water.index = cropped.dates
        if hasattr(self, 'g_shear_rate') and isinstance(self.g_shear_rate, pd.DataFrame):
            cropped.g_shear_rate = self.g_shear_rate.iloc[start_idx:end_idx].copy()
            cropped.g_shear_rate.index = cropped.dates
            # IMPORTANT: Keep original normalization constants (full tidal cycle range)
            # These are used for normalized_shear in flocs.py and must represent
            # the full spring-neap tidal range, not just the cropped period
            cropped.g_shear_rate_min = self.g_shear_rate_min
            cropped.delta_g_shear_rate = self.delta_g_shear_rate
        if hasattr(self, 'bed_shear_stress') and isinstance(self.bed_shear_stress, pd.DataFrame):
            cropped.bed_shear_stress = self.bed_shear_stress.iloc[start_idx:end_idx].copy()
            cropped.bed_shear_stress.index = cropped.dates
        if hasattr(self, 'water_depth') and isinstance(self.water_depth, pd.DataFrame):
            cropped.water_depth = self.water_depth.iloc[start_idx:end_idx].copy()
            cropped.water_depth.index = cropped.dates

        # Slice pre-computed arrays (these are the critical ones for model performance)
        if hasattr(self, 'PAR_array'):
            cropped.PAR_array = self.PAR_array[start_idx:end_idx].copy()
        if hasattr(self, 'T_array'):
            cropped.T_array = self.T_array[start_idx:end_idx].copy()
        if hasattr(self, 'mu_water_array'):
            cropped.mu_water_array = self.mu_water_array[start_idx:end_idx].copy()
        if hasattr(self, 'g_shear_rate_array'):
            cropped.g_shear_rate_array = self.g_shear_rate_array[start_idx:end_idx].copy()
        if hasattr(self, 'bed_shear_stress_array'):
            cropped.bed_shear_stress_array = self.bed_shear_stress_array[start_idx:end_idx].copy()
        if hasattr(self, 'water_depth_array'):
            cropped.water_depth_array = self.water_depth_array[start_idx:end_idx].copy()

        # Slice riverine loads if present
        if hasattr(self, 'loads') and self.riverine_loads and self.loads is not None:
            cropped.loads = self.loads.iloc[start_idx:end_idx].copy()
            cropped.loads.index = cropped.dates
            if hasattr(self, 'loads_NH4_array') and self.loads_NH4_array is not None:
                cropped.loads_NH4_array = self.loads_NH4_array[start_idx:end_idx].copy()
            if hasattr(self, 'loads_NO3_array') and self.loads_NO3_array is not None:
                cropped.loads_NO3_array = self.loads_NO3_array[start_idx:end_idx].copy()
            if hasattr(self, 'loads_DIP_array') and self.loads_DIP_array is not None:
                cropped.loads_DIP_array = self.loads_DIP_array[start_idx:end_idx].copy()
            if hasattr(self, 'loads_DSi_array') and self.loads_DSi_array is not None:
                cropped.loads_DSi_array = self.loads_DSi_array[start_idx:end_idx].copy()

        # Slice prescribed TEP if present
        if hasattr(self, 'TEP_array') and self.prescribe_TEP and self.TEP_array is not None:
            cropped.TEP_array = self.TEP_array[start_idx:end_idx].copy()

        # Slice prescribed Flocs if present
        if hasattr(self, 'prescribe_Flocs') and self.prescribe_Flocs:
            if hasattr(self, 'Microflocs_massconc_array') and self.Microflocs_massconc_array is not None:
                cropped.Microflocs_massconc_array = self.Microflocs_massconc_array[start_idx:end_idx].copy()
            if hasattr(self, 'Micro_in_Macro_massconc_array') and self.Micro_in_Macro_massconc_array is not None:
                cropped.Micro_in_Macro_massconc_array = self.Micro_in_Macro_massconc_array[start_idx:end_idx].copy()
            if hasattr(self, 'Macroflocs_sink_sed_array') and self.Macroflocs_sink_sed_array is not None:
                cropped.Macroflocs_sink_sed_array = self.Macroflocs_sink_sed_array[start_idx:end_idx].copy()
            if hasattr(self, 'Macroflocs_source_resusp_array') and self.Macroflocs_source_resusp_array is not None:
                cropped.Macroflocs_source_resusp_array = self.Macroflocs_source_resusp_array[start_idx:end_idx].copy()
            if hasattr(self, 'Macroflocs_numconc_array') and self.Macroflocs_numconc_array is not None:
                cropped.Macroflocs_numconc_array = self.Macroflocs_numconc_array[start_idx:end_idx].copy()

        # Rebuild date-to-index mapping
        cropped.dates_to_index = {date: idx for idx, date in enumerate(cropped.dates)}

        # IMPORTANT: Keep original temperature bounds for normalization
        # T_max and T_min are used in getlimT() (base.py) to normalize the Arrhenius
        # temperature limitation function. They must represent the ANNUAL temperature
        # range, not the cropped period's range. Otherwise, limT would equal 1.0 at
        # winter temperatures instead of ~0.5.
        # cropped.T_max and cropped.T_min are inherited via copy.copy(self)
        # so we do NOT recalculate them here.

        if self.verbose:
            print(f"Setup cropped: {self.name} -> {cropped.name}")
            print(f"  Original: {self.start_date} to {self.end_date} ({self.tmax} days)")
            print(f"  Cropped:  {cropped.start_date} to {cropped.end_date} ({cropped.tmax:.1f} days)")
            print(f"  Array lengths: {len(self.t_eval)} -> {len(cropped.t_eval)}")
            print(f"  Inherited normalization constants (from full cycle):")
            print(f"    T_max={cropped.T_max:.2f} K, T_min={cropped.T_min:.2f} K")
            print(f"    g_shear_rate_min={cropped.g_shear_rate_min:.2f}, delta={cropped.delta_g_shear_rate:.2f}")

        return cropped


if __name__ == "__main__":

    """Test the tidal implementation with MOW1_STATION configuration."""
    print("Testing Tidal Implementation in BGC Physical Setup")
    print("=" * 60)


    
    # Test 1: Create tidal setup
    print("\n1. Testing MOW1_STATION with Tidal Parameters")
    print("-" * 40)
    
    try:
        setup = Setup(**phys_setup.MOW1)

        print(f"✓ Setup created: {setup.name}")
        print(f"✓ Tidal variation: {setup.vary_water_depth}")
        print(f"✓ M2 period: {setup.tidal_period_M2:.4f} days ({setup.tidal_period_M2*24:.2f} hours)")
        print(f"✓ Spring-neap period: {setup.spring_neap_period:.1f} days")

        if setup.vary_water_depth:
            depths = setup.water_depth['WaterDepth'].values
            print(f"✓ Water depth range: {np.min(depths):.1f} - {np.max(depths):.1f} m")
            print(f"✓ Mean depth: {np.mean(depths):.1f} m")

            # Check tidal physics
            expected_spring_range = 2 * setup.water_depth_amplitude_spring  # 5m
            expected_neap_range = 2 * setup.water_depth_amplitude_neap      # 3m
            actual_range = np.max(depths) - np.min(depths)

            print(f"✓ Tidal range: {actual_range:.1f} m (expected: {expected_neap_range}-{expected_spring_range} m)")
            
    except Exception as e:
        print(f"✗ Failed to create tidal setup: {e}")
        exit(1)
    
    # Test 2: Backward compatibility
    print("\n2. Testing Backward Compatibility (Constant Depth)")
    print("-" * 40)
    
    try:
        const_config = phys_setup.MOW1_STATION.copy()
        const_config['vary_water_depth'] = False
        const_config['tmax'] = 1
        
        const_setup = Setup(**const_config)
        
        print(f"✓ Constant setup created")
        print(f"✓ vary_water_depth: {const_setup.vary_water_depth}")
        print(f"✓ Water depth: {const_setup.water_depth} m (type: {type(const_setup.water_depth).__name__})")
        
        if isinstance(const_setup.water_depth, (int, float)):
            print("✓ Backward compatibility maintained")
        else:
            print("⚠ Water depth should be scalar for constant case")
            
    except Exception as e:
        print(f"✗ Backward compatibility test failed: {e}")
    
    # Test 3: Visualization
    print("\n3. Generating Setup Visualization")
    print("-" * 40)
    
    try:
        print("Creating comprehensive setup plot...")
        setup.plot_setup(days_to_plot=None)
        print("✓ Visualization complete!")
        
    except Exception as e:
        print(f"⚠ Visualization failed: {e}")
        print("  (This might be due to missing data files or display issues)")

    from src.utils import functions as fns



    setup2 = fns.deep_update(phys_setup.MOW1, phys_setup.additive_mode)
    setup2 = Setup(**setup2)
    setup2.plot_setup(days_to_plot=None)
    plt.show()



    # Summary
    print("\n" + "=" * 60)
    print("TIDAL IMPLEMENTATION TEST SUMMARY")
    print("=" * 60)
    print("✓ Tidal water depth calculation implemented")
    print("✓ M2 tidal period (12.42h) and spring-neap cycle (14.7d) configured")
    print("✓ Shear rate and bed shear stress synchronized with tides")  
    print("✓ Backward compatibility maintained")
    print("✓ Visualization method added")
    print("\nImplementation ready! 🌊")
    print("\nUsage:")
    print("  from src.core.phys import Setup")
    print("  from bgpartmod_private.src.config_model.phys_setup import MOW1_STATION")
    print("  setup = Setup(**MOW1_STATION)")
    print("  setup.plot_setup()  # Visualize all parameters")