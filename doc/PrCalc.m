%% SNR Calculation Function
% EE725 Project Function
% Rich Mullin 04-28-07

function SNR = PrCalc(R,theta)
 
    Pt = 40e3;
    G = 316.2*cos(20*pi/180 - theta);
    lambda = .015;
    sigma = .1;
    k = 1.38e-23;
    T = 273;
    B = (1/13)*(1/1e-6);

    SNR = (Pt*G^2*lambda^2*sigma)/((4*pi)^3*R^4*k*T*B);