%% EE725 Project Main Program %%
%% Written by Rich Mullin, 4-28-07 %%
 
 
clc;
clear all;
close all;
d2r = pi/180;
r2d = 180/pi;
c = 3e8;
 
%% Initialization Parameters
v = [250 0];
RinitCh1 = [10000 500];
RinitCh2 = [10000 500-.0117];
freq = 20e9;
lambda = c/freq;
t2 = 0:0.5:40;
sampFreq = 26e6;
prf = 50e3;
numPulses = 128;
BC = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
%BC = [1 1 -1 1];%For Debug Purposes
%BC = [1 1 1 1];%For Debug Purposes
 
%% Generate Target Truth Data Set for Both Radar Channels
index1 = 1;
for t = 0:0.5:40;
trueRtCh1(index1,:) = RinitCh1 - v.*t;
trueRtCh2(index1,:) = RinitCh2 - v.*t;
trueLosRCh1(index1,:) = sqrt(trueRtCh1(index1,1)^2 + trueRtCh1(index1,2)^2);
trueLosRCh2(index1,:) = sqrt(trueRtCh2(index1,1)^2 + trueRtCh2(index1,2)^2);
trueElAngCh1(index1,:) = atan(trueRtCh1(index1,2)/trueRtCh1(index1,1));
trueElAngCh2(index1,:) = atan(trueRtCh2(index1,2)/trueRtCh2(index1,1));
trueLosRdotCh1(index1,:) = cos(trueElAngCh1(index1,1))*v;
trueLosRdotCh2(index1,:) = cos(trueElAngCh2(index1,1))*v;
fdCh1(index1,:) = 2*trueLosRdotCh1(index1,:)/lambda;
fdCh2(index1,:) = 2*trueLosRdotCh2(index1,:)/lambda;
PrCh1(index1,:) = PrCalc(trueLosRCh1(index1,1),trueElAngCh1(index1,1));
PrCh2(index1,:) = PrCalc(trueLosRCh2(index1,1),trueElAngCh2(index1,1));
index1 = index1+1;
end;
 
%% Generate Range-Pulse and Range-Doppler Matrices for Both Channels for
%% All Ranges
for index4 = 1:length(t2)
%index4 = 1;%For Debug Purposes
    for index2 = 1:numPulses
        priCutCh1(index2,:) = matfiltv2(trueLosRCh1(index4,1),trueLosRdotCh1(index4,1),BC,index2);
        priCutCh2(index2,:) = matfiltv2(trueLosRCh2(index4,1),trueLosRdotCh2(index4,1),BC,index2);
    end
 
    for index3 = 1:length(priCutCh1)
        priCutWinCh1(:,index3) = priCutCh1(:,index3).*chebwin(length(priCutCh1(:,index3)),50);
        priCutWinCh2(:,index3) = priCutCh2(:,index3).*chebwin(length(priCutCh2(:,index3)),50);
    end
 
    rngPulseMatCh1 = priCutWinCh1;%cat(1,priCut,zeros(128,532));
    rngPulseMatCh2 = priCutWinCh2;%cat(1,priCut,zeros(128,532));
    
    for index3 = 1:size(priCutCh1,2)
        rngDopMatTransposeCh1(index3,:) = fft(rngPulseMatCh1(:,index3));
        rngDopMatTransposeCh2(index3,:) = fft(rngPulseMatCh2(:,index3));
    end
 
    rngDopMatCh1 = rngDopMatTransposeCh1';
    rngDopMatCh2 = rngDopMatTransposeCh2';
    rngDopMatCh1dB = 10*log10(abs(rngDopMatCh1));
    rngDopMatCh2dB = 10*log10(abs(rngDopMatCh2));
    maxResponseCh1 = max(max(rngDopMatCh1dB));
    [indexY,indexX] = find(rngDopMatCh1dB == maxResponseCh1);
    measRng(index4) = (536-indexX)*1/sampFreq*c/2;
    measDopFreq(index4) = indexY*(1/((1/prf)*numPulses));
    measPhaseDiff(index4) = angle(rngDopMatCh1(indexY,indexX)-rngDopMatCh2(indexY,indexX));
    measElAng(index4) = asin(measPhaseDiff(index4)/(2*pi*0.78));
end
 
%% Set Sample Time
Ts = 0.5;
t3 = 0:Ts:40-Ts*2;
 
%% Generate Smoothed Position Data with alpha-beta filter
for index2 = 2:length(measRng)
    if index2 == 2
    xHatNow = measRng(index2-1);
    xDotBarNowMinus = (measRng(index2-1)-measRng(index2))/Ts;
    else
        measXNow = measRng(index2);
        
        %% Set Gains
        alpha = 0.5;
        beta = 0.3;
 
        %% Smoothed Position
        xBarNow = xHatNow + alpha*(measXNow - xHatNow);
 
        %% Smoothed Velocity
        xDotBarNow = xDotBarNowMinus + beta*((measXNow - xHatNow)/Ts);
 
        %% Next Predicted Position
        xHatNowPlus = xHatNow + xDotBarNow*Ts;
        
        smoothPos(index2-2) = xBarNow;
        smoothVel(index2-2) = xDotBarNow;
        xHatNow = xHatNowPlus;
    end
end
 
%% Generate Smoothed Elevation Angle Data with alpha-beta filter
for index3 = 2:length(measElAng)
    if index3 == 2
    xHatNow1 = measElAng(index3-1);
    xDotBarNowMinus1 = (measElAng(index3-1)-measElAng(index3))/Ts;
    else
        measXNow1 = measElAng(index3);
        
        %% Set Gains
        alpha1 = 0.55;
        beta1 = 0.3;
 
        %% Smoothed Position
        xBarNow1 = xHatNow1 + alpha1*(measXNow1 - xHatNow1);
 
        %% Smoothed Velocity
        xDotBarNow1 = xDotBarNowMinus1 + beta1*((measXNow1 - xHatNow1)/Ts);
 
        %% Next Predicted Position
        xHatNowPlus1 = xHatNow1 + xDotBarNow1*Ts;
        
        smoothElAng(index3-2) = xBarNow1;
        smoothElAngRate(index3-2) = xDotBarNow1;
        xHatNow1 = xHatNowPlus1;
    end
end
 
%% Plotting for Debug and Report Generation
% figure()
% image(rngDopMatCh1dB);
% xlabel('range bin');ylabel('doppler bin');
% figure()
% image(rngDopMatCh2dB);
% xlabel('range bin');ylabel('doppler bin');
% figure()
% mesh(real(rngPulseMatCh1));
% figure()
% mesh(imag(rngPulseMatCh1));
% figure()
% mesh(abs(rngPulseMatCh1));
% xlabel('time samples');ylabel('pulse number');zlabel('amplitude')
 
 
figure()
plot(t2,trueLosRCh1,'LineWidth',1.5);
xlabel('time (sec)'); ylabel('Line of Sight Range (m)');
hold on;
plot(t2,measRng,'r','LineWidth',1.5);
hold on;
plot(t3,smoothPos,'k','LineWidth',1.5);
legend('True Line of Sight Range','Measured Line of Sight Range','\alpha-\beta Smoothed Line of Sight Range');
% 
figure()
plot(t2,trueElAngCh1.*r2d,'LineWidth',1.5);
hold on;
plot(t2,measElAng.*r2d,'r','LineWidth',1.5);
xlabel('time (sec)'); ylabel('Elevation Angle (deg)');
hold on;
plot(t3,smoothElAng.*r2d,'k','LineWidth',1.5);
legend('True Elevation Angle','Measured Elevation Angle','\alpha-\beta Smoothed Elevation Angle');
% 
% figure()
% plot(t2,fdCh1(:,1)./1000,'LineWidth',1.5);
% xlabel('time (sec)'); ylabel('Doppler Frequency (kHz)');
% hold on;
% plot(t2,measDopFreq./1000,'r','LineWidth',1.5);
% legend('True Doppler Frequency','Measured Doppler Frequency');
% 
% figure()
% plot(t2,10*log10(PrCh1),'LineWidth',1.5);
% xlabel('time (sec)'); ylabel('P_R (dBw)');





%%Errored Truth Model
%% EE725 Project Target Errored Truth Data %%
%% Written by Rich Mullin, 4-28-07 %%
 
 
clc;
clear all;
% close all;
d2r = pi/180;
r2d = 180/pi;
 
v = [250 0];
Rinit = [10000 500];
freq = 20e9;
c = 3e8;
lambda = c/freq;
Ts = 0.2;
t2 = 0:Ts:40;
 
index1 = 1;
for t = 0:Ts:40;
trueRt(index1,:) = Rinit - v.*t;
trueLosR(index1,:) = sqrt(trueRt(index1,1)^2 + trueRt(index1,2)^2);
trueElAng(index1,:) = atan(trueRt(index1,2)/trueRt(index1,1));
trueRdot(index1,:) = cos(trueElAng(index1,:))*v;
fd(index1,:) = 2*trueRdot(index1,:)/lambda;
Pr(index1,:) = PrCalc(trueRt(index1,1),trueElAng(index1,1));
index1 = index1+1;
end;
 
%% Generate Errored Range Data
sigmaRng = (c*1e-6/2)./sqrt(8.*Pr(:,1));
rngErr = 100.*randn(1,length(t2))';%sigmaRng.*randn(1,81)';
errLosR = trueLosR(:,1) + rngErr;
LosRerror = trueLosR - errLosR;
 
%% Set Sample Time
t3 = 0:Ts:40-Ts*2;
 
%% Generate Smoothed Position Data with alpha-beta filter
for index2 = 2:length(errLosR)
    if index2 == 2
    xHatNow = errLosR(index2-1);
    xDotBarNowMinus = (errLosR(index2-1)-errLosR(index2))/Ts;
    else
        measXNow = errLosR(index2);
        
        %% Set Gains
        alpha = 0.5;
        beta = 0.3;
 
        %% Smoothed Position
        xBarNow = xHatNow + alpha*(measXNow - xHatNow);
 
        %% Smoothed Velocity
        xDotBarNow = xDotBarNowMinus + beta*((measXNow - xHatNow)/Ts);
 
        %% Next Predicted Position
        xHatNowPlus = xHatNow + xDotBarNow*Ts;
        
        smoothPos(index2-2) = xBarNow;
        smoothVel(index2-2) = xDotBarNow;
        xHatNow = xHatNowPlus;
    end
end
 
%% Generate Errored Elevation Angle Data
sigmaElAng = 1./sqrt(4.9.*Pr(:,1));
ElAngErr = sigmaElAng.*randn(1,length(t2))';%100.*randn(1,length(t2))';%
errElAng = trueElAng(:,1) + ElAngErr;
ElAngError = trueElAng - errElAng;
 
 
 
%% Generate Smoothed Elevation Angle Data with alpha-beta filter
for index3 = 2:length(errElAng)
    if index3 == 2
    xHatNow1 = errElAng(index3-1);
    xDotBarNowMinus1 = (errElAng(index3-1)-errElAng(index3))/Ts;
    else
        measXNow1 = errElAng(index3);
        
        %% Set Gains
        alpha1 = 0.5;
        beta1 = 0.3;
 
        %% Smoothed Position
        xBarNow1 = xHatNow1 + alpha1*(measXNow1 - xHatNow1);
 
        %% Smoothed Velocity
        xDotBarNow1 = xDotBarNowMinus1 + beta1*((measXNow1 - xHatNow1)/Ts);
 
        %% Next Predicted Position
        xHatNowPlus1 = xHatNow1 + xDotBarNow1*Ts;
        
        smoothElAng(index3-2) = xBarNow1;
        smoothElAngRate(index3-2) = xDotBarNow1;
        xHatNow1 = xHatNowPlus1;
    end
end
 
 
figure()
plot(t2,trueLosR(:,1),'LineWidth',1.5);
hold on;
plot(t2,errLosR(:,1),'r','LineWidth',1.5);
xlabel('time (sec)'); ylabel('distance (m)');
hold on;
plot(t3,smoothPos,'k','LineWidth',1.5);
legend('True Line of Sight Range','Measured Line of Sight Range','\alpha-\beta Smoothed Line of Sight Range');
% 
% figure()
% plot(t2,LosRerror);
 
 
 
figure()
plot(t2,trueElAng.*r2d,'LineWidth',1.5);
xlabel('time (sec)'); ylabel('Elevation Angle (deg)');
hold on;
plot(t2,errElAng.*r2d,'r','LineWidth',1.5);
hold on;
plot(t3,smoothElAng.*r2d,'k','LineWidth',1.5);
legend('True Elevation Angle','Measured Elevation Angle','\alpha-\beta Smoothed Elevation Angle');
% 
% figure()
% plot(t2,fd(:,1)./1000,'LineWidth',1.5);
% xlabel('time (sec)'); ylabel('Doppler Frequency (kHz)');
% 
% figure()
% plot(t2,10*log10(Pr),'LineWidth',1.5);
% xlabel('time (sec)'); ylabel('SNR (dBw)');


