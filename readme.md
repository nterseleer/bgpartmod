# Marine Biogeochemical Model Framework

A Python framework for implementing 0D marine biogeochemical models with comprehensive simulation management and optimization facilities. This modular code provides a flexible platform for marine ecosystem modeling, making it particularly suitable for research applications, MSc and PhD theses, and model development studies.

## What This Code Is

This is a **Python implementation framework** that enables users to:

- Build and run 0D marine biogeochemical models through simple configuration
- Manage complex simulation workflows with automatic logging and version control
- Optimize model parameters against observational data using advanced algorithms
- Analyze and visualize results with comprehensive plotting utilities

## Note about the license
This software is currently **under embargo** until scientific publication (submission planned for summer 2025). 
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
   mkdir -p Figs Observations simulations/{Simulations,References,Optimizations}
   ```

4. Verify installation:
   ```bash
   python src/main_example.py
   ```

## Quick Start

### Basic Simulation

```python
from src.Config_model import base_config
from src.utils import simulation_manager, plotting
from core import phys

# Use base configuration (available in repository)
model_config = base_config.Onur

# Set up physical environment
setup = phys.Setup(
    **phys.DEFAULT_SETUPS['onur22'],
    tmax=20,  # Simulation duration (days)
    PARfromfile=True  # Use light data if available
)

# Run simulation
simulation = simulation_manager.run_or_load_simulation(
    model_config, setup,
    name="basic_simulation"
)

# Plot key results
plotting.plot_results(simulation, plotting.phy_nuts)
```

## Model Configuration

Models are configured using nested dictionaries that define components, parameters, coupling relationships, and initial conditions:

```python
coupled_config = {
    'formulation': 'Coupled_BGC_Floc',
    'Phy': {
        'class': phyto.Phyto,
        'parameters': {
            'mu_max': 4.5,        # Maximum growth rate (d⁻¹)
            'QN_max': 0.15,       # Max N:C ratio (mol:mol)
            'exud_rate': 0.1      # DOM exudation rate
        },
        'coupling': {
            'coupled_NH4': 'NH4',
            'coupled_TEP_source': 'TEP'
        },
        # ... initialization and aggregation rules
    },
    'Microflocs': {
        'class': flocs.Flocs,
        'parameters': {
            'K_glue': 10.0,       # TEP stickiness factor
            'alpha_FF_ref': 0.05  # Base collision efficiency
        },
        'coupling': {
            'coupled_TEP': 'TEP',
            'coupled_detritus': 'Detritus'
        }
    }
    # ... other components
}
```

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
│   │   └── model.py           # Main model orchestration
│   ├── Config_model/        # Model configurations
│   │   └── base_config.py      # Standard model setups
│   ├── Config_system/       # System configurations
│   │   └── path_config.py      # Directory paths
│   ├── utils/               # Utilities and analysis tools
│   │   ├── simulation_manager.py # Simulation workflow management
│   │   ├── optimization.py     # Parameter optimization
│   │   ├── plotting.py         # Visualization utilities
│   │   ├── evaluation.py       # Model-data comparison
│   │   ├── desolver.py         # Differential evolution solver
│   │   ├── phys.py            # Physical setup
│   │   └── functions.py        # General utilities
│   ├── main_example.py      # Basic usage examples
│   └── optim_main_example.py  # Optimization examples
├── Figs/                    # Generated figures
├── Observations/            # Observational data
├── Simulations/            # Model output storage
│   ├── Model_runs/           # Regular simulations  
│   ├── References_simulations/ # Reference runs
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
