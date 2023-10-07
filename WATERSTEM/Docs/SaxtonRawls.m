function [lambda,SMsat,Ksat]=SaxtonRawls(Sand,Clay,OM)
Sand=Sand/100;
Clay=Clay/100;
SM33t=-0.251*Sand+0.195*Clay+0.011*OM+0.006*Sand*OM-0.027*Clay*OM+0.452*Sand*Clay+0.299;
SM33=SM33t+1.283*SM33t^2-0.374*SM33t-0.015;
SM1500t=-0.024*Sand+0.487*Clay+0.006*OM+0.005*Sand*OM-0.013*Clay*OM+0.068*Sand*Clay+0.031;
SM1500=SM1500t+0.14*SM1500t-0.02;
SMsat33t=0.278*Sand+0.034*Clay+0.022*OM-0.018*Sand*OM-0.027*Clay*OM-0.584*Sand*Clay+0.078;
SMsat33=SMsat33t+0.636*SMsat33t-0.107;
lambda=(log(SM33)-log(SM1500))/(log(1500)-log(33));
SMsat=SM33+SMsat33-0.097*Sand+0.043;
Ksat=1930*(SMsat-SM33)^(3-lambda);%Ksat in mm/hr
Ksat=Ksat/1000/3600;%Ksat in m/s
end