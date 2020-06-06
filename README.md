# Masters Thesis: 
## Models for propagating facilitation in the insect visual system.

### Abstract: 
Detecting and tracking moving targets within a visual scene is a complex task. Over thousands of years, many species of animals like flying insects have evolved neural mechanisms for tracking path or location of flying targets that move against the visually cluttered background. Small target motion-detecting (STMD) neurons found in the visual pathway of flying insects display remarkable selectivity and sensitivity to small moving targets and are thought to support these mechanisms. Contributing to their sensitivity is a form of facilitation, in which the responsiveness of an STMD is enhanced by prior exposure to a small target moving along a continuous path in visual space. The locus of facilitation in the receptive field is found to be local to the area of the target and to continue propagating in the direction of target motion even after a stimulus ceases. 

In this thesis, we are modeling this phenomenon with the propagation of traveling waves in densely interconnected, retinotopic layers of cells. We hypothesize that waves are initiated and reinforced by the presence of a moving target stimulus, and the network, in turn, interacts with STMDs to modulate their excitability. Membrane potentials travel too fast to play this role, so we have studied and modeled propagating calcium waves as a possible mechanism. Accordingly, we have proposed two models that could be the biological substrate for this mechanism: 1) a network of astrocyte-like glia cells, and 2) a network of neurons; both in which calcium waves are initiated by target stimuli and propagate via diffusion with the participation of regenerative mechanisms. Finally we show qualitatively and quantitatively how facilitation would look in any such network and discuss ranges of parameters that would support the facilitation mechanism. 


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
