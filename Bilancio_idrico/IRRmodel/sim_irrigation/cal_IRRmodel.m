function [X]=cal_IRRmodel(DPEIS,X_ini,FIG,namefig)
if nargin==2,X_ini=[0.1,0.05,0.1,0.5,0.1,0.1]';end
if nargin<4,FIG=0;end
[RES]=fmincon(@calibOK,X_ini,[],[],[],[],...
     zeros(6,1),ones(6,1),[],...
     optimset('Display','iter','MaxIter',300,'MaxFunEvals',5000,...
     'TolFun',1E-20,'TolCon',4,'Largescale','off','Algorithm','active-set'),...
     DPEIS);
X=convert_adim(RES);
[output,R,R_IRR]=IRRmodel(DPEIS,X,1,namefig);

%---------------------------------------------------------------------------------
function [err]=calibOK(X_0,DPEIS)

X=convert_adim(X_0);
[output,R,R_IRR,NS]=IRRmodel(DPEIS,X,0,'');
% err=RMSEcum+2*(1-R^2);
err=1-R;
% save X_PAR

%---------------------------------------------------------------------------------
function X=convert_adim(X_0)
% LOW=[   1,   0.0,  1,  0.0, 0.000]';
% UP =[ 500, 800.0, 50, 20.0, 0.040]';
LOW=[  0.5,  50, .80, .01, .1, .3]';
UP =[  0.9, 120, .95, .45, .7, .5]';
X=LOW+(UP-LOW).*X_0;
