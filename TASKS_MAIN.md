# ToDo

[#] #=priority in linspace(0,5) (decreasing)

## General
- [ ] temporal resolution of your data IN EVERY SINGLE CODE
- [x] [biblio research] bounds for WCM with CR --> check Notion under Master Thesis/WCM
- [x] orbit normalization (simple mean bias elimination) --> ref. Mladenova 2013, code: Sigma_norm
- [ ] build database (.csv) of satellite SM products (to aggregate sources in a single place)
- [x] check data on IRRINET (ask Matteo): missing data in 2020, 2018 (possible to have only 2 irrigation events?), must check diario di campo with CER people
- [ ] check data in diario di campo (ask Matteo, CER)
- [ ] tillage periods? (ask Matteo, CER)
- [x] extend period in sigma0 normalization (download all data from starting date of satellite up until now)
- [ ] Comparison with literature
- [ ] Consider different vegetation indices (e.g. cr, Vegetation Water Content, LAI from S2 (see mail Domenico & TerraScope))
- [ ] Do the analysis over a longer time period (2017-2020) and so collect precipitation, EPOT, irrigation and crop data for 2018 and 2019 in Budrio...TEST

## Other products
- [ ] Comparison THEIA, RT1 with OBSERVED and compare goodness with our model
- [ ] Compare with PLANET!!! (also VWC)!
- [ ] RT1 old VS RT1 new


## GEE
- [x]     [2] implement S2 data download
- [x]     [2] routine for NDVI calculation

## WCM
- [x] [1] check input data (to solve divergence)
    - [x] [2] comparsion between different satellite products (vd IRRI_WCM)
- [x] [1] backscattering normalization for different acquisition geom and angles: try both cos^2 and cdf normalization 
- [x]         [3] run with different vegetation indexes (try Cross Ratio)
- [x] [1] check SM normalization
- [ ]     [2] weighted mean SM in 3 hours around hour of passage of s-1
- [x]     [2] add plot input IRR, RAIN

To run this code in a significant way, follow:
- 1. check spatial mean: in linear or db scale? --> linear (check Reading_summaries)
    More on spatial distribution of backscattering intensity:
    https://developers.google.com/earth-engine/tutorials/community/detecting-changes-in-sentinel-1-imagery-pt-1
- 2. check if one-line WCM (single equation, not more) is different from standard
    --> it is, standard performs better but can present wider divergence
- 3. normalize backscattering for acquisition geometry (angle) HOW??? options: cosine, **distribution bias elimination**

## IRR
- [x]     [2] implement calibration with hardcoded s_fc and s_w
- [x] [1] run with different satellite products (THEIA, RT1)
- [ ]         [5] bridge gap in SM data by using as benchmark similar irrigation+rain event
- [x] [1] work on normalization of SM from different sources
- [x]     [2] compare retrieved SM with satellite measurements (scatterplots)
- [x] [1] add BIAS function to statistical metrics in title of plot and also in aux functions
- [ ] Check input soil moisture and solve scaling issue between observed and modeled
- [ ] Running IRRmodel on hourly dataset and then calibrating 𝜎^0 on an hourly basis: this could give better results



## Optimizer (PSO) performance check
- [ ]     [2] produce plots of distribution of parameters over multiple runs, fit them
- [ ]     [2] produce animation of particles' trajectories
- [ ]     [2] implement montecarlo for finding best options for optimizer
- [ ]     [2] implement montecarlo for finding best bounds (you will need to work with some fixed params) (not sure if this makes sense)
- [ ]     [2] study options of optimizer: https://pyswarms.readthedocs.io/en/latest/examples/tutorials/options_handler.html 
- [ ]     [2] study also: https://pyswarms.readthedocs.io/en/latest/examples/tutorials/custom_optimization_loop.html
- [ ] Fix seed of optimizer for stability
- [ ] Build a best-fit to retrieve parameters that maximize R coefficient on SM, IRR

## Templates
- [ ] [5] make it a class (at the end!)
- [ ] statistics retrieval: function that produces a dict with R, RMSE, NS, KGE etc
- [x] hist fitted with gaussian
- [x] fitting functions (gauss, skew gauss, others)

------------------------------------------------------------------------------
# Questions

- how do I normalize SM products from different satellites (always max-min norm?) --> YES
- Is NDVI the best vegetation descriptor to be used in WCM? Using just 1 vegetation descriptor is the best choice or should one use 2 different products in the model?
- Can we use the crop specific parameters in IRRmodel to build a vegetation descriptor and thus eliminate the need for external input?  should look at parameters timeseries (to begin with…)!
- Check if daily aggregation of rain, irrigation (sum) over the whole 24 hours could affect 𝜎^0 calibration (we have 30% of total satellite points with rain in the same day and 5% with irrigation)
- Model is VERY UNSTABLE in parameters calibration: why? Optimizer's fault? Few data?
