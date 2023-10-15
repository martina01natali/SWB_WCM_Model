function [BF,QS,W,W2,E]=SWB(PIO,EPOT,SWE,W,W2,Ks,m2,Ks2,m22,W_max,W_max2,alpha)    

%Vol_Ini=W+W2+PIO+SWE;

IE=(PIO+SWE).*((W/W_max).^alpha);

E1=EPOT.*W/W_max; 
E2=0;
E=E1+E2;

PERC=Ks.*(W/W_max).^(m2).*(W2<W_max2); 
PERC_amm=W+PIO+SWE-IE-E1; PERC_amm=PERC_amm.*(PERC_amm>0);
PERC=PERC.*(PERC<=PERC_amm)+PERC_amm.*(PERC>PERC_amm);   
W=max(0,W+PIO+SWE-IE-PERC-E1);

PERC2=Ks2.*(W2/W_max2).^(m22);
PERC2_amm=W2+PERC-E2; PERC2_amm=PERC2_amm.*(PERC2_amm>0);
PERC2=PERC2.*(PERC2<=PERC2_amm)+PERC2_amm.*(PERC2>PERC2_amm);
W2=max(0,W2+PERC-PERC2-E2);  

SE=(W-W_max).*(W>W_max);
SE2=(W2-W_max2).*(W2>W_max2);
W(W>W_max)=W_max;
W2(W2>W_max2)=W_max2;
W(W<0)=0;
W2(W2<0)=0;

BF=PERC2;
QS=IE+SE+SE2;   

%Vol_Fin=W+W2+E1+E2+BF+QS;

# -------------------
depth = depth(S, C, W[i-1])... [mm]
static theta_sat
W_max = depth*theta_sat [mm]
W[i] = W[i-1] + P + I - ET - PERC [mm]
WW[i] = W[i]/W_max [sat degree] = W[i]/(depth*theta_sat)
sat_degree -> vol? sat_degree*theta_sat = vol

theta[i] = W[i]/depth [m3m-3]





end
