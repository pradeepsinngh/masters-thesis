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
- Step 0: Make sure you have ```SD_rates.mat``` is in your current directory.
- Step 1: Run ```genrate_str.m```.
  - This will create 
- Step 2: Run ```simulate.m``` 
- Step 3: Run ```plot_Ca.m``` 
- Step 4: Run ```wavefront_speed.m``` 

### Dependencies:
1. MOLE Library: https://github.com/jcorbino/mole
   - Download and install MOLE (Refer MOLE documentation.
   - We use MOLE for numerical simulations.
2. MATLAB (any stable version, we suggest R2017A)
