# Masters Thesis: 
## Models for propagating facilitation in the insect visual system.

Flying insect species like dragonflies are capable of predicting the path or location of their target even if the target has occluded by some object for some period of time. This ability to predict the path is supported by a processing mechanism which is called *Response Facilitation*. This facilitation is known to increases sensitivity to small objects that move along continuous paths and is thought to increase reliability of small target detection. Because the locus of facilitation that is induced by a moving target propagates in visual space even after a small target stimulus ceases, we have proposed that it could be supported by traveling wave phenomena in retinotopically-organized regions of the visual system. 

Accordingly, we have proposed two models that could be the biological substrate for this sort of mechanism: 
- A network of astrocyte-like glia cells, and 
- A network of neurons; 

both in which calcium waves are initiated by target stimuli and propagate via diffusion with the participation of regenerative mechanisms

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
