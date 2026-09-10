# Marine Biogeochemical Model Framework

A Python framework for implementing 0D marine biogeochemical models with comprehensive simulation management and optimization facilities. This modular code provides a flexible platform for marine ecosystem modeling, making it particularly suitable for research applications, MSc and PhD theses, and model development studies.

## What This Code Is

This is a **Python implementation framework** that enables users to:

- Build and run 0D marine biogeochemical models through simple configuration
- Manage complex simulation workflows with automatic logging and version control
- Optimize model parameters against observational data using advanced algorithms
- Analyze and visualize results with comprehensive plotting utilities

## Note about the license
This software is currently **under embargo** until scientific publication (submission planned for autumn 2025). 
After publication, it will be released under the European Union Public Licence (EUPL). 
For any inquiry, please contact @nterseleer.

## What It Was Developed For

The framework was specifically developed to implement a **pioneering coupled biogeochemical-mineral flocculation model**. This coupling addresses the traditionally separate treatment of biological processes and physical particle dynamics, which are fundamentally interconnected in marine systems. The coupled approach enables simulation of:

- Phytoplankton exudation and transparent exopolymer particle (TEP) formation
- Bacterial processing of dissolved and particulate organic matter  
- Bimodal flocculation dynamics (microflocs and macroflocs)
- Feedback mechanisms between biological activity and particle aggregation

## Framework Philosophy

The code follows the **Framework for Aquatic Biogeochemical Models (FABM; Bruggeman and Bolding, 2014)** approach, using runtime composition where model structure is defined through dictionary-based configuration. This FABM-inspired design makes it exceptionally well-suited for biogeochemical modeling applications, particularly in academic settings where flexibility and ease of use are paramount.

**Key advantages for BGC modeling:**
- **Runtime model composition**: Add/remove processes without code modification
- **Parameter experimentation**: Easy sensitivity analysis and calibration
- **Reproducible configurations**: Version-controlled model setups
- **Educational value**: Clear separation between model and implementation

## Model Components

The framework is built around the concept of **components** - self-contained modules that represent specific biological, chemical, or physical processes. Each component can be independently configured and coupled with others to create complex ecosystem models.

### Available Components

#### Biological Components
- **Phytoplankton (`Phyto`)**: Primary producers with multi-nutrient (N, P, Si) limitation, variable stoichiometry, and growth-dependent DOM exudation
- **Heterotrophs**: Multiple bacterial types (free-living bacteria, particle-attached bacteria), heterotrophic flagellates, and ciliates with distinct feeding preferences and metabolic pathways

#### Organic Components  
- **Dissolved Inorganic Matter (`DIM`)**: Individual nutrient pools (NH4, NO3, DIP, DSi) with competitive uptake dynamics
- **Dissolved Organic Matter (`DOM`)**: Size-structured pools with different bioavailability and aggregation properties
- **Detritus**: Particulate organic matter with size-dependent bacterial processing and settling

#### Mineral Components
- **Flocculation (`Flocs`)**: Bimodal population balance model representing microflocs and macroflocs with TEP-mediated aggregation processes

### Key Processes

- **DOM Dynamics**: Phytoplankton exudation shifts from small to large molecules under nutrient stress
- **TEP Formation**: DOC coagulation creates sticky particles that enhance flocculation
- **Bimodal Flocculation**: Efficient representation of natural particle size distributions
- **Organo-mineral interactions**: TEP of biological origin affect mineral flocculation dynamics
- **Bacterial Processing**: Distinct communities handle dissolved vs. particulate substrates

## Installation

### Prerequisites

- Python 3.8+
- Required packages: numpy, pandas, matplotlib, scipy, dill

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/nterseleer/bgpartmod.git
   cd bgpartmod
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Create necessary directories:
   ```bash
   mkdir -p Figs Observations Simulations/{Model_runs,Reference_simulations,Optimizations}
   ```

4. Verify installation:
   ```bash
   python src/main_example.py
   ```

## Quick Start

### Basic Simulation

```python
from src.config_model import base_config
from src.core import model, phys
from src.utils import plotting

# Reference biogeochemical configuration (Kerimoglu et al., 2022)
setup = phys.Setup(tmax=30., dt=1e-2, dt2=1e-3)
simulation = model.Model(base_config.Kerimoglu2022, setup=setup, name="basic_simulation")

plotting.plot_results(simulation, ['Phy_C', 'Phy_Chl', 'NO3_concentration'], observations=None)
```

See `src/main_example.py` for a runnable version, and `src/utils/simulation_manager.py`
(`run_or_load_simulation`) to save and reload runs instead of recomputing them.

## Model Configuration

Models are configured using nested dictionaries that define components, parameters, coupling relationships, and initial conditions:

```python
from src.components import flocs
from src.config_model import base_config
from src.utils import functions as fns

# Add the flocculation module to the reference biogeochemistry, and couple the two ways:
# TEP feeds floc aggregation, mineral flocs attenuate light for the phytoplankton.
coupled_config = fns.deep_update(base_config.Kerimoglu2022, {
    'Phy': {
        'coupling': {'coupled_SPM': ['Microflocs', 'Micro_in_Macro']},
    },
    'Microflocs': {
        'class': flocs.Flocs,
        'parameters': {
            'alpha_PP_base': 0.10,   # [-] mineral collision efficiency
            'K_glue': 15.0,          # [mmolC m-3] half-saturation of the TEP effect
            'delta_alpha_PP': 0.03,  # [-] increment of the collision efficiency at saturating TEP
        },
        'coupling': {'coupled_Nf': 'Macroflocs', 'coupled_Nt': 'Micro_in_Macro',
                     'coupled_glue': 'TEPC'},
        'initialization': {'numconc': 1.0e12},   # [# m-3]
    },
    # ... Macroflocs and Micro_in_Macro, same class, mirrored couplings
})
```

Parameters not listed keep the class defaults, which reproduce the reference publications
(Kerimoglu et al. 2022 for the biogeochemistry, Lee et al. 2011 for the flocculation).

## Configuration Modules You Provide

Some things are your own choices rather than part of the framework. Write them as modules
and drop them into `src/config_model/`; if a module is absent, the library falls back to
the minimal version in `src/config_model/_defaults/`, so nothing here is required to get
started.

| Module | Holds | Fallback |
|---|---|---|
| `vars_to_plot.py` | which variable sets you plot | `_defaults/` |
| `config_diagnostics.py` | which diagnostics you store (the main lever on memory use) | `_defaults/` |
| `plot_config.py` | figure geometry, styles, default observation dataset | `_defaults/` |
| `src/utils/observations.py` | your observation loader | none — plots simply omit observations |

Your own model configurations go in `src/config_model/config.py`, built on top of
`base_config.Kerimoglu2022` with `config_tools.deep_update` (see Model Configuration above).

Site-specific configurations, observation data and analysis notebooks are kept in a
separate private repository, mounted as `_private/` alongside `src/`. It is not part of
this distribution, and nothing here depends on it.

## Model Output and Analysis

The model produces comprehensive output including:

- **State variables**: All component concentrations over time
- **Diagnostic variables**: Growth rates, limitation factors, flocculation rates
- **Aggregate variables**: Total C, N, P, Si pools and fluxes
- **Size distributions**: Particle size spectra and settling velocities

Output formats: pickle (full data), feather (fast I/O), CSV (external analysis)


## File Structure

```
bgpartmod/
├── src/
│   ├── components/           # Model components
│   │   ├── phytoplankton.py    # Primary producers
│   │   ├── heterotrophs.py     # Bacteria, flagellates, ciliates
│   │   ├── dim.py             # Dissolved inorganic nutrients
│   │   ├── dom.py             # Dissolved organic matter
│   │   ├── detritus.py        # Particulate organic matter
│   │   └── flocs.py           # Flocculation processes
│   ├── core/                # Core model functionality
│   │   ├── base.py             # Base classes
│   │   ├── phys.py             # Physical setup and forcings
│   │   ├── legacy.py           # Backward compatibility of configuration dictionaries
│   │   └── model.py            # Main model orchestration
│   ├── config_model/        # Model configurations
│   │   ├── base_config.py      # Reference configuration (Kerimoglu et al., 2022)
│   │   ├── varinfos.py         # Output variables: units, labels, derived expressions
│   │   └── _defaults/          # Fallbacks for the modules you are meant to provide
│   ├── config_system/       # System configurations
│   │   └── path_config.py      # Directory paths
│   ├── utils/               # Utilities and analysis tools
│   │   ├── simulation_manager.py # Simulation workflow management
│   │   ├── optimization.py     # Parameter optimization
│   │   ├── plotting.py         # Visualization utilities
│   │   ├── plotted_variables_sets.py # Variable sets with their display preferences
│   │   ├── evaluation.py       # Model-data comparison
│   │   ├── desolver.py         # Differential evolution solver
│   │   ├── flux_network.py     # Annual flux network (Sankey analysis)
│   │   ├── config_tools.py     # Configuration dictionaries: merge, compare, transform
│   │   ├── observations_example.py # Minimal observation loader
│   │   └── functions.py        # General utilities
│   ├── main_example.py      # Basic usage examples
│   └── optim_main_example.py  # Optimization examples
├── Figs/                    # Generated figures
├── Observations/            # Observational data
├── Simulations/            # Model output storage
│   ├── Model_runs/           # Regular simulations
│   ├── Reference_simulations/ # Reference runs
│   └── Optimizations/        # Parameter optimization results
└── data/                   # Input data files
```

## Troubleshooting

**Import errors**: Ensure you're running Python from the repository root directory
**Missing data**: Check that required directories exist (created in setup step 3)
**Optimization issues**: Reduce population size for testing; check parameter bounds

## Acknowledgements

This model builds upon several key works:

- **Biogeochemical model**: Kerimoglu, O., Hintz, N. H., Lücken, L., Blasius, B., Böttcher, L., Bunse, C., ... & Simon, M. (2022). Growth, organic matter release, aggregation and recycling during a diatom bloom: a model-based analysis of a mesocosm experiment. bioRxiv, 2022-05.
- **Flocculation approach**: Lee, B. J., Toorman, E., Molz, F. J., & Wang, J. (2011). A two-class population balance equation yielding bimodal flocculation of marine or estuarine sediments. Water research, 45(5), 2131-2145.
- **FABM framework**: Bruggeman, J., & Bolding, K. (2014). A general framework for aquatic biogeochemical models. Environmental modelling & software, 61, 249-265.
- **Differential evolution**: Storn, R., & Price, K. (1997). Differential evolution–a simple and efficient heuristic for global optimization over continuous spaces. Journal of global optimization, 11, 341-359.


## Citation

When using this model, please cite the forthcoming publication. 
