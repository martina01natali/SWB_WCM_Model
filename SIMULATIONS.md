# 15-02-23
working on IRRI_WCM_v8_Copy1, following local minima of r_st and Kc0, adding parameters instead of starting simulation from the beginning, 20 runs for each configuration of bounds

## 22.49
run with new bounds for w_max, ww_w, ww_fc, @ 3 sigma

# 16-02-23
## 12.02
working on v8\_copy1 e \_copy2, 2 parallel simulations with 100 run with respectively large and narrow boundaries

# 17-02-23
## 15.11
- running simulation on \_v9 with threshold on R of SM @ 0.75: takes forever --> change strategy
- running simulation on \_v9 with threshold on KGE of SM @ 0.75 (too high), @ 0.7 (too high), @ 0.6 (too high, takes forever) --> change strategy
- running simulation on \_v9_Copy1 with new cost function with sum of KGE on sigma0 and SM: takes a standard time

# 19-02-23
## 9.33
- \_narrow-bounds: 2017, bounds all narrow but for Kc in [0.1, 1.5]
- \_large-bounds: 2017, bounds all narrow but for Kc in [0.5,1.5]
- \_large-bounds-Copy1: 2020, large Wmax, Wfc, Ww, narrow rho_st, large Kc0 (+- 50%)

- \_large-bounds-Copy1: 2020, large Wmax, narrower Wfc [0.4,0.5], large Ww, narrow rho_st, large Kc0 (+- 50%)