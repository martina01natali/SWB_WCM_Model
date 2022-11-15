clc,clear,close all
namesite='ITALY_BUDRIO';
siteID='5';

dirfobs='D:\WORK\PROJECTS\02_EU\ESA\IRRIGATION+\elab\myelab\Budrio\';
fobs='Irrigation_Test_Site_Budrio.nc';
SMobs=ncread([dirfobs,fobs],'Soil_Moisture_first_layer_5');
Dobs=ncread([dirfobs,fobs],'Time_hours_5')+datenum('2000-01-01');
SMobs(1:13)=[];
Dobs(1:13)=[];
SMobs=time_res_mean(SMobs,12);
Dobs=Dobs(1:24:end);
SM_OBS=[Dobs,norma(SMobs)];SM_OBS=rem_nan(SM_OBS);


D=ncread(['TEST_SITE_',namesite,'.nc'],'Time_days')+datenum('1-Jan-2000');
IRR=ncread(['TEST_SITE_',namesite,'.nc'],['Irrigation_',(siteID)]);
P=ncread(['TEST_SITE_',namesite,'.nc'],['Rainfall_',(siteID)]);
ET=ncread(['TEST_SITE_',namesite,'.nc'],['ET_',(siteID)]);
PET=ncread(['TEST_SITE_',namesite,'.nc'],['PET_',(siteID)]);

% SSM_CCI_active=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_CCI_active_',(siteID)]);
% SSM_SMAP=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_SMAP_',(siteID)]);
% SSM_THEIA=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_THEIA_',(siteID)]);
% SSM_RT1=ncread(['TEST_SITE_',namesite,'.nc'],['SSM_RT1_',(siteID)]);
% prod={'CCI_ACT','SMAP','THEIA','RT1'};
% SSM=[SSM_CCI_active,SSM_SMAP,SSM_THEIA,SSM_RT1];
% IRR(isnan(IRR))=0;

dataSM=[SM_OBS];
dataSM=rem_nan(dataSM);
dataSM_int=interp1gap(dataSM(:,1),dataSM(:,2),D,30);
dataSM_int=dataSM_int.*.45+.55;

W_0=      0.7; % [-]
W_max=    80; % [mm]
S_fc=     0.9; % [-]
S_w=      0.4; % [-]
rho_st=   0.2; % [-]
Kc=       0.4; % [-]

DPEIS=[D,P,PET,IRR,dataSM_int];
DPEIS=DPEIS(find(D==datenum('1-apr-2017')):end,:);
PAR=[W_0,W_max,S_fc,S_w,rho_st,Kc];
FIG=1;
namefig=[namesite,'_',siteID,'_OK'];

cal_IRRmodel(DPEIS,[.5,.5,.5,.5,.5,.5]',1,namefig)
%%
[output] = IRRmodel(DPEIS,PAR,FIG,namefig);

nansum(IRR)
sum(output(:,2))