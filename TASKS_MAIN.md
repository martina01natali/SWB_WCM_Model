# ToDo

GEE
- [ ]     [2] implement S2 data download
- [ ]     [2] routine for NDVI calculation

WCM
- [x] [1] check input data (to solve divergence)
    - [ ] [2] comparsion between different satellite products
- [ ] [1] backscattering normalization for different acquisition geom and angles: try both cos^2 and cdf normalization 
- [ ]         [3] run with different vegetation indexes
- [ ] 

IRR
- [ ]     [2] implement calibration with hardcoded s_fc and s_w
- [ ] [1] run with different satellite products (THEIA, RT1)
- [ ]         [5] bridge gap in SM data by using as benchmark similar irrigation+rain event
- [ ] 
- [ ] 
- [ ] 

Optimizer (PSO) performance check
- [ ] produce plots of distribution of parameters over multiple runs, fit them
- [ ] produce animation of particles' trajectories
- [ ] implement montecarlo for finding best options for optimizer
- [ ] implement montecarlo for finding best bounds (you will need to work with some fixed params) (not sure if this makes sense)
- [ ] 
- [ ] 