function [output,R,R_IRR,NS] = IRRmodel(DPEIS,PAR,FIG,namefig)

% IRRmodel function
% INPUT: TEST_SITE_country_city.nc
% OUTPUT: [WW, IRR] = [soil moisture, irrigation needed]
% Ref. paper: https://doi.org/10.1371/journal.pone.0250979

% Loading input data
[M]=size(DPEIS,1); % dimension of date vector
D=DPEIS(:,1); % datetime, vector
P=DPEIS(:,2); % precipitation, vector
EPOT=DPEIS(:,3); % evapotrans
IRRobs=DPEIS(:,4); % irrigation obs
WWobs=DPEIS(:,5); % data_SM obs and cleaned from nan and interpolated

% Model parameter
W_0=      PAR(1); % [-]     % water min
W_max=    PAR(2); % [mm]    % water max
S_fc=     PAR(3); % [-]     % field capacity
S_w=      PAR(4); % [-]     % wilting point
rho_st=   PAR(5); % [-]     % time-dependent depletion fraction
Kc=       PAR(6); % [-]     % constant for EVO

W_fc = S_fc*W_max;
W_w  = S_w*W_max;

W=NaN(M,1);
IRR=zeros(M,1);
PS=zeros(M,1);
Ks=zeros(M,1);
rho=zeros(M,1);
W(1)=W_0*W_max;

% Main ROUTINE
for t=2:M
    rho(t)=rho_st+0.04.*(5-Kc.*EPOT(t));
    
    if W(t-1)>=(1-rho(t))*W_fc
        Ks(t)=1;
    elseif (W(t-1)>W_w)&&(W(t-1)<(1-rho(t))*W_fc)
        Ks(t)=(W(t-1)-W_w)./((1-rho(t)).*(W_fc-W_w));
    else
        Ks(t)=0;
    end
    DOY=D(t)-datenum(year(D),1,1);
    if (DOY>134)&(DOY<230) % summer season
        if W(t-1)<=(1-rho(t))*W_fc
            IRR(t)=W_fc-W(t-1);
        end
    end
    
    W(t)=W(t-1)+P(t)+IRR(t)-EPOT(t).*Kc.*Ks(t);
    if W(t)>W_fc
        PS(t)=W(t)-W_fc;
        W(t)=W_fc;
    end
end
WW=W./W_max;

% Calculation of model performance
RMSE=nanmean((WW-WWobs).^2).^0.5;
NS=1-nansum((WW-WWobs).^2)./nansum((WWobs-nanmean(WWobs)).^2);
NS_radQ=1-nansum((sqrt(WW+0.00001)-sqrt(WWobs+0.00001)).^2)./...
    nansum((sqrt(WWobs+0.00001)-nanmean(sqrt(WWobs+0.00001))).^2);
NS_lnQ=1-nansum((log(WW+0.0001)-log(WWobs+0.0001)).^2)...
    ./nansum((log(WWobs+0.0001)-nanmean(log(WWobs+0.0001))).^2);
NS_lnQ=real(NS_lnQ);
NS_radQ=real(NS_radQ);

R=corr(WW,WWobs,'rows','complete');
R_IRR=corr(IRR,IRRobs,'rows','complete');

% Figure
if FIG==1
    set(gcf,'paperpositionmode','manual','paperposition',[1 1 20 16],'Color','white')
    set(gcf,'position',[50 50 1000 500])
    
    h(1) = axes('Position',[0.1 0.5 0.8 0.40]);
    set(gca,'Fontsize',12)
    s=(['sumIRRobs= ',num2str(nansum(IRRobs),'%4.1f'),...
        ' sumIRRsim= ',num2str(nansum(IRR),'%4.1f'),...
        ' R-SM= ',num2str(R,'%4.3f'),...
        ' R-IRR= ',num2str(R_IRR,'%4.3f')]);
    title(['\bf',s]);
    hold on
    plot(D,WWobs,'g','Linewidth',3)
    plot(D,WW,'r','Linewidth',2);
    
    legend('\theta_r_e_f','\theta_s_i_m');
    datetick('x',20)
    set(gca,'Xticklabel','')
    ylabel('relative soil moisture [-]')
    grid on, box on
    axis([D(1) D(end) min(WWobs)-0.05 max(WWobs)+0.05])
    
    h(2) = axes('Position',[0.1 0.1 0.8 0.40]);
    set(gca,'Fontsize',12)
    hold on
    area(D,P,'edgecolor','none','facecolor',[.5 .5 .5])
    plot(D,IRRobs,'color','r','Linewidth',3)
    plot(D,IRR,'color','b','Linewidth',1.5)
    grid on, box on
    legend('rainfall','irrigation_o_b_s','irrigation_s_i_m');
    ylabel('rain and irrigation (mm)')
    datetick('x',20)
    grid on, box on
    axis([D(1) D(end) 0 1.05.*max([P;IRR])])
    
    export_fig(['IRRmodel_',namefig], '-png','-q60','-r150')
end

output=[WW,IRR];
