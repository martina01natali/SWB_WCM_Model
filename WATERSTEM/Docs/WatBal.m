function [Theta2]=WatBal(Theta1,IntPio,ET0mod)
	global TheSat The_FC The_WP f_min Coeff_p Ksat Coe_c SpeStr PasTemSec %stress   
    SatEff=Theta1/TheSat;
    CapInf=f_min+Coeff_p*(1-SatEff);
    Inf=IntPio*(IntPio<=CapInf)+CapInf*(IntPio>CapInf);
    Ric=Ksat*SatEff^Coe_c;
    stress=The_FC;
    AET=ET0mod*(Theta1>=stress)+...
        ET0mod*(Theta1-The_WP)/(stress-The_WP)*(Theta1<stress)*(Theta1>The_WP);
    Theta2=Theta1+((Inf-Ric-AET)/SpeStr*PasTemSec);
    Theta2=Theta2*(Theta2<=TheSat)*(Theta2>=0)+TheSat*(Theta2>TheSat);
end