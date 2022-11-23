# ToDo

## General
- [ ] temporal resolution of your data IN EVERY SINGLE CODE
- [x] [biblio research] bounds for WCM with CR
- [ ] orbit normalization (simple mean bias elimination)


## GEE
- [ ]     [2] implement S2 data download
- [ ]     [2] routine for NDVI calculation

## WCM
- [x] [1] check input data (to solve divergence)
    - [ ] [2] comparsion between different satellite products
- [ ] [1] backscattering normalization for different acquisition geom and angles: try both cos^2 and cdf normalization 
- [ ]         [3] run with different vegetation indexes (try Cross Ratio)
- [x] [1] check SM normalization
- [ ]     [2] weighted mean SM in 3 hours around hour of passage of s-1
- [ ]     [2] add plot input IRR, RAIN
- [ ]     [2] plot RT1 over retrieved SM

To run this code in a significant way, follow:
- 1. check spatial mean: in linear or db scale? --> linear (check Reading_summaries)
    More on spatial distribution of backscattering intensity:
    https://developers.google.com/earth-engine/tutorials/community/detecting-changes-in-sentinel-1-imagery-pt-1
- 2. check if one-line WCM (single equation, not more) is different from standard
    --> it is, standard performs better but can present wider divergence
- 3. normalize backscattering for acquisition geometry (angle) HOW??? options: cosine, distribution bias elimination

## IRR
- [ ]     [2] implement calibration with hardcoded s_fc and s_w
- [x] [1] run with different satellite products (THEIA, RT1)
- [ ]         [5] bridge gap in SM data by using as benchmark similar irrigation+rain event
- [ ] [1] work on normalization of SM from different sources
- [x]     [2] compare retrieved SM with satellite measurements (scatterplots)
- [ ] 

## Optimizer (PSO) performance check
- [ ]     [2] produce plots of distribution of parameters over multiple runs, fit them
- [ ]     [2] produce animation of particles' trajectories
- [ ]     [2] implement montecarlo for finding best options for optimizer
- [ ]     [2] implement montecarlo for finding best bounds (you will need to work with some fixed params) (not sure if this makes sense)
- [ ]     [2] study options of optimizer: https://pyswarms.readthedocs.io/en/latest/examples/tutorials/options_handler.html 
- [ ]     [2] study also: https://pyswarms.readthedocs.io/en/latest/examples/tutorials/custom_optimization_loop.html

------------------------------------------------------------------------------
# Questions

- how do I normalize SM products from different satellites (always max-min norm?)