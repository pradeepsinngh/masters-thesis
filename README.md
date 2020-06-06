# Codebase for my Masters Thesis: 
## Models for propagating facilitation in the insect visual system.



### Structure: 
```bash
├── model
│   ├── SD_rates.mat
│   ├── genrate.m
│   ├── simulate.m
│   ├── plot_Ca.m
│   ├── wavefront_speed.m
└── 
```

### Usage:
- Step 0: Run ```kinetic_rates_SD.m``
  - This will create ```SD_rates.mat``` file, make sure that is present in your current directory.
- Step 1: Run ```genrate_str.m```.
  - This will create ```Data_files.mat```, that will store positions (cordinates) for all compartments in all cells.
- Step 2: Run ```simulate.m``` 
  - This will create ```Ca_store.mat```, that will hold Ca history for all cells for all time setps.
- Step 3: Run ```plot_Ca.m``` 
  - This script will visualize (plot) the calcium waves in network. You can also store this visualization as video.
- Step 4: Run ```wavefront_speed.m``` 
  - Computes average wavefront speed in a network.
  
NOTE: All files should be in the same directory to work probably. 

### Dependencies:
1. MOLE Library: https://github.com/jcorbino/mole
   - Download and install MOLE (Refer MOLE documentation.
   - We use MOLE for numerical simulations.
2. MATLAB (any stable version, we suggest R2017A)
