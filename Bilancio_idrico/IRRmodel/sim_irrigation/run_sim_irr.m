% Main code for running simulation

clc,clear,close all
namesite='GERMANY_BRANDENBURG';
siteID='4';

% Data reading of specific namesie and siteID
tot=ncread(['TEST_SITE_',namesite,'.nc'],:);
D=ncread(['TEST_SITE_',namesite,'.nc'],'Time_days')+datenum('1-Jan-2000'); % 
IRR=ncread(['TEST_SITE_',namesite,'.nc'],['Irrigation_',(siteID)]);
P=ncread(['TEST_SITE_',namesite,'.nc'],['Rainfall_',(siteID)]);
ET=ncread(['TEST_SITE_',namesite,'.nc'],['ET_',(siteID)]);
PET=ncread(['TEST_SITE_',namesite,'.nc'],['PET_',(siteID)]);
SSM_CCI_active=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_CCI_active_',(siteID)]);
SSM_SMAP=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_SMAP_',(siteID)]);
SSM_THEIA=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_THEIA_',(siteID)]);
SSM_RT1=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_RT1_',(siteID)]);
prod={'CCI_ACT','SMAP','THEIA','RT1'};
SSM=[SSM_CCI_active,SSM_SMAP,SSM_THEIA,SSM_RT1];
IRR(isnan(IRR))=0;

dataSM=[D,SSM(:,3)]; % matrix date-data
dataSM=rem_nan(dataSM); % clean from NaN
dataSM_int=interp1gap(dataSM(:,1),dataSM(:,2),D,30); % interpolate

W_0=      0.5; % [-]
W_max=    80; % [mm]
S_fc=     0.8; % [-]
S_w=      0.1; % [-]
rho_st=   0.2; % [-]
Kc=       0.4; % [-]

DPEIS=[D,P,PET,IRR,dataSM_int];
% DPEIS=DPEIS(find(D==datenum('1-apr-2017')):end,:);
% DPEIS=DPEIS(1:find(D==datenum('1-apr-2017')),:);
PAR=[W_0,W_max,S_fc,S_w,rho_st,Kc];
FIG=1;
namefig=[namesite,'_',siteID];

[output] = IRRmodel(DPEIS,PAR,FIG,namefig);

nansum(IRR)
sum(output(:,2))