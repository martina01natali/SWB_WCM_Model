#############################################################################
# FIRST BUILD FUNCTIONS IN NOTEBOOK THEN MOVE THEM HERE WHEN PERFECT
#############################################################################
#############################################################################
#############################################################################
#############################################################################


#############################################################################
# Triple plot
#############################################################################


def triple_plot(s0:list, ww:list, ):
    """
    
    """
    
    #----------------------------------------------------------------------
    # Plot of sigma0
    
    obs = s0[0]; obs_label=r'$\sigma^0_{obs}$'
    sim = sigma0; sim_label=r'$\sigma^0_{sim}$'
    labely = r'$\sigma^0$[dB]'
    times = t_sat
    marker='o'; linestyle='-'
    
    RMSE=np.mean((sim-obs)**2)**0.5; print('RMSE =', RMSE)
    R=np.corrcoef(sim,obs)[0][1]; print('R=', R)
    BIAS=bias(sim,obs); print('bias =', BIAS)
    print('KGE=', KGE[0])
    
    title=f'{sim_label} VS {obs_label} - RMSE={RMSE:.2f}, R={R:.2f}, bias={BIAS:.2f}, KGE={KGE[0]:.2f}'
    
    ax[0].set_xlim(xmin=times[0], xmax=times[len(times)-1])
    ax[0].plot(times, sim, c='tab:red', label=sim_label,
               linestyle=linestyle, marker=marker, )#alpha=.4, zorder=-1)
    ax[0].plot(times, obs, c='tab:blue', label=obs_label,
               linestyle=linestyle, marker=marker, alpha=.4, zorder=-1)
    ax[0].legend(loc='upper left')
    ax[0].set_title(title)
    ax[0].set_ylabel(labely)
    
    ax1 = ax[0].twinx()
    ax1.plot(times, veg, label=opt_veg, color='tab:green')
    ax1.legend(loc='upper right')
    ax1.set_ylabel(opt_veg)
    
    #-----------------------------------------------------------------------
    # Plot of SM
    
    obs = WW_obs; obs_label='SM_obs'
    sim = WW; sim_label='SM_sim'
    times = t
    
    # RMSE, NS, R, bias calculation
    RMSE=np.nanmean((sim-obs)**2)**0.5; print('RMSE =', RMSE)
    NS=1-np.nansum((sim-obs)**2)/np.nansum((obs-np.nanmean(obs))**2); print('NS =', NS)
    NS_radQ=1-np.nansum((np.sqrt(sim+0.00001)-np.sqrt(obs+0.00001))**2)/np.nansum((np.sqrt(obs+0.00001)-np.nanmean(np.sqrt(obs+0.00001)))**2)
    NS_lnQ=1-np.nansum((np.log(sim+0.0001)-np.log(obs+0.0001))**2)/np.nansum((np.log(obs+0.0001)-np.nanmean(np.log(obs+0.0001)))**2)
    NS_lnQ=NS_lnQ.real; # print(NS_lnQ) 
    NS_radQ=NS_radQ.real; # print(NS_radQ)
    
    simmatrix = np.array( [ [sim[i], obs[i]] for i in range(len(sim)) if not np.isnan(obs[i]) ] )
    R=np.corrcoef(simmatrix,rowvar=False)[0][1]; print('R (sim vs obs) =', R)
    BIAS=bias(np.array([e[0] for e in simmatrix]), np.array([e[1] for e in simmatrix]))
    
    if irri:
        IRRmatrix = np.array( [ [IRR[i], IRR_obs[i]] for i in range(len(IRR)) if not np.isnan(IRR_obs[i]) ] )
        R_IRR=np.corrcoef(IRRmatrix,rowvar=False)[0][1]; print('R_IRR (IRR vs IRR_obs)=', R_IRR)
        B_IRR=bias(np.array([e[0] for e in IRRmatrix]), np.array([e[1] for e in IRRmatrix]))
        irri_title = f'sumIRR_obs={np.sum(IRR_obs):.2f}, '+\
                     f'sumIRR_sim={np.sum(IRR):.2f}, '+\
                     f'R_IRR={R_IRR:.2f}, '+\
                     f'bias_IRR={B_IRR:.2f}, '
    else: irri_title=''
    
    title=f'SM_obs VS SM_sim - '+f'{irri_title}'+\
        f'R_SM={R:.2f}, bias_SM={BIAS:.2f}'
    
    ax[1].set_xlim(xmin=times[0], xmax=times[-1])
    ax[1].plot(times, sim, c='tab:red', label=sim_label)
    ax[1].plot(times, obs, c='tab:blue', label=obs_label, linestyle='-', alpha=.8, zorder=-1)
    ax[1].legend(loc='upper left')
    ax[1].set_title(title)
    ax[1].set_ylabel('Relative SM [-]')
    
    #-----------------------------------------------------------------------
    # Plot of inputs P, IRR, veg
    
    times = t
    
    ax[2].bar(times, P, color='tab:gray', label=r'rainfall')
    ax[2].bar(times, IRR_obs, color='tab:blue', label=r'$IRR_{obs}$', zorder=2)
    ax[2].legend(loc='upper left')
    ax[2].set_ylabel('Irrigation and rain [mm]')
    
    ax2 = ax[2].twinx()
    ax2.plot(times, EPOT, label='ET0 [mm]', color='tab:green')
    ax2.legend(loc='upper right')
    ax2.set_ylabel('ET0 [mm]')
    
    name=''
    if opt_save:
        optim_choice = 'glo' if (optim=='')or(optim=='global') else 'local'
        name = timestr+filename # +f'_{n_particles}_{n_step}_{optim_choice}_{norma}'
        plt.savefig('Plot\\'+name+'.png')
        
    
#############################################################################
# Scatter plot
#############################################################################

    if automate: opt_save = True
    else: opt_save = True if input('Save scatterplot of SM? [y/n]')=='y' else False
    
    import matplotlib.gridspec as gridspec
    
    def linear(x,a,b):
        return a+b*x
    
    quantity = 'SM'# r'$\sigma^0$' # SM
    sim = WW # WW, sigma0
    obs = WW_obs # WW_d, VV
    
    title = f'Observed {quantity} VS simulated - ' # y VS x
    xlabel = f'{quantity}_sim [-]'
    ylabel = f'{quantity}_obs [-]'
    filename = f'scatter_'+'sm_'+units+PAR_str_add # 'SM'
    
    data = pd.DataFrame({'sim': sim,'obs': obs})
    data.dropna(inplace=True)
    x = data.sim.values
    y = data.obs.values
    
    fig = plt.figure(figsize=(6, 6), dpi=200)
    gs = gridspec.GridSpec(nrows=1, ncols=1, width_ratios=[1], height_ratios=[1])
    ax = plt.subplot(gs[0])
    ax.plot(x, y, marker='o', linestyle='', color='tab:blue')
    ax.set_xlim(np.min([x,y])-0.1*abs(np.mean([x,y])), np.max([x,y])+0.1*abs(np.mean([x,y])))
    ax.set_ylim(np.min([x,y])-0.1*abs(np.mean([x,y])), np.max([x,y])+0.1*abs(np.mean([x,y])))
    lin_grid = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100); ax.plot(lin_grid, lin_grid, color='k')
    
    # Fit
    popt, pcov = curve_fit(linear, x, y)
    ax.plot(lin_grid, linear(np.array(lin_grid),*popt), color='tab:orange')
    R=np.corrcoef(x,y)[0][1]; print('R=', R, 'R^2=', R**2)
    B=bias(x,y)
    
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    xtext=0.8*(np.max(x)-np.min(x))+np.min(x)  ; ytext=0.1*(np.max(y)-np.min(y))+np.min(y)
    ax.text(xtext, ytext,
            f'y={popt[0]:.2f}+{popt[1]:.2f}x\n'+
            r'$R^2$'+f'={R**2:.2f}',
            ha="center", va="center", size=15,
            bbox=dict(boxstyle="round,pad=0.3", fc="tab:orange", ec="k", lw=2, alpha=.5))
    
    ax.set_title(title+f'R={R:.2f},'+r' $R^2$'+f'={R**2:.2f}, bias={B:.2f}')
    ax.set_aspect('equal', adjustable='box')
    
    if opt_save: plt.savefig('Plot\\'+timestr+'_'+filename+'.png')