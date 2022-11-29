# Naming conventions for plots/simulations

## Dump of things to take into account
- name of figure file
- name of folder
- different inputs: SM data
- optimization: best options

model = irri/wcm
timestr = time.strftime("%y%m%d-%H%M%S")
optim_choice = 

name = '{model}+'_'+timestr+f'_{n_particles}_{n_step}_{optim_choice}_{norma}'


## User input
Inputs are to be stored in a separate variable or if-else, named with
`opt` (=option) followed by parameters.

```opt_{parameter} = input(...)```

Filenames are to be stored in a separate variable, named with
`filename` followed by brief description

ex.
```
filename_plot_{...}
filename_df_{...}

```

