%% Matched Filter Function
% EE 725 Barker Code Ambiguity Response
% Written by Rich Mullin 04-21-07

function [priCut] = matfiltv2(Rng,v,BC,priNum)

    c = 3e8;
    k = 1.38e-23;
    Tn = 273;

    % v = 250;
    tau = 1e-6;
    prf = 50e3;
    T = 1/prf;
    freq = 20e9;
    sampFreq = 26e6;
    lambda = c/freq;
    fd = 2*v/lambda;
    %BC = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    tauPrime = tau/length(BC);
    B = 1/tauPrime;
    %Rng = 10000;
    A = sqrt(2*PrCalc(Rng,atan(500/Rng)));
    t = 0:1/(1*sampFreq):tau;
    td = 2*Rng/c;

    for index1 = 1:length(BC)
        t2 = 1;%tauPrime/2:1/sampFreq:tauPrime;
        %t2 = tauPrime:1/sampFreq:tauPrime;
        if BC(index1) == 1
            %xmtChip = exp(-j*2*pi*t2);
            xmtChip = 1;
            %recChip = exp(-j*2*pi*fd*t2);
            %recChip = xmtChip*exp(-j*2*pi*t2*fd);
            recChip = xmtChip*exp(-j*4*pi*T*(priNum-1)*v*t2/lambda);
        else
            %xmtChip = j*exp(-j*2*pi*t2);
            xmtChip = -1;
            %recChip = j*exp(-j*2*pi*fd*t2);
            %recChip = xmtChip*exp(-j*2*pi*fd*t2);
            recChip = xmtChip*exp(-j*4*pi*T*(priNum-1)*v*t2/lambda);
        end
        if index1 == 1
           xmtPulse = xmtChip;
           recPulse = recChip;
        else
           xmtPulse = cat(2,xmtPulse,xmtChip);
           recPulse = cat(2,recPulse,recChip);
        end
    end



    %recPulseNoisy = recPulse;

    sampTime = 1/sampFreq;
    numSamps = T/sampTime;
    recPulse = [A.*recPulse zeros(1,numSamps-length(recPulse))];
    numShifts = ceil(td/sampTime);
    %noiseSamp = exp(-j*randn(1,numSamps-length(recPulse)));
    qnoise = randn(1,numSamps);
    inoise = randn(1,numSamps);
    noiseSamp = sqrt(qnoise.^2+inoise.^2).*exp(-j.*atan2(qnoise,inoise));

    recNoisySamp = recPulse + noiseSamp;

    tempShiftSampNoisy = recNoisySamp';
    tempRngSampNoisy = circshift(tempShiftSampNoisy,numShifts);
    rngSampNoisy = tempRngSampNoisy';

    priCut = conv(xmtPulse,fliplr(conj(rngSampNoisy)));

    % figure()
    % stem(1:length(xmtPulse),xmtPulse);
    % hold on;
    % stem(1:length(xmtPulse),imag(xmtPulse),'r');
    % title('Transmit Pulse');
    % legend('Real Part','Imaginary Part');
    % 
    % figure()
    % stem(1:length(recPulse),recPulse,'g');
    % hold on;
    % stem(1:length(recPulse),imag(recPulse),'k');
    % title('Received Pulse Amplitude Scaled without Noise');
    % legend('Real Part','Imaginary Part');
    % 
    % figure()
    % stem(1:length(recNoisySamp),real(recNoisySamp),'g');
    % hold on;
    % stem(1:length(recNoisySamp),imag(recNoisySamp),'k');
    % title('Received Pulse Amplitude Scaled, with Noise');
    % legend('Real Part','Imaginary Part');
    % 
    % figure()
    % plot(1:length(xmtPulse),real(xmtPulse));
    % hold on;
    % plot(1:length(rngSampNoisy),real(rngSampNoisy),'r');
    % title('Transmit Pulse and Received Pulse with Noise Samples Real Part Only');
    % legend('Transmit Pulse','Received Noisy Sample');
    % 
    % figure()
    % plot(1:length(priCut),abs(priCut));
    % title('Matched Filter Output');