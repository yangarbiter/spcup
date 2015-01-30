% Program: srps_v2.m
% Author: Jerry Chou
% Last Modified Date: 2015/1/23
%
% 1. Introduction:
% In this program, I implement the SRPS framework to the PPG treconstruct 
% ECG signal processing problem. SRPS stands for Signal Reconstruction and
% Peak Selection. The algorithm will be divided by two parts, one is for
% SR, the other is for PS.
%
% 2. Instruction
%   To run the program, please follow the instructions below:
% 
% a. Load data
%   User can set the dataset you want to examine in the Set dataset section
%
% b. Algorithm
%   The algorithm is online, which maens that it does the estimation
%   sequencially. So, there's a loop run through every window in order, and
%   doing estimation inside the loop.
%
% c. Remarks
%   In current algorithm, I only use the first PPG channel.
%
% d. Results
%   The best result is in dataset9:
%      error1 = 0.5085 and error2 = 0.45%.
%   The worst resilt is in dataset12:
%      error1 = 7.9370 and error2 = 5.58%.
%   The average error and standard deviation of the twelve dataset are:
%      error1 = 2.66910    std = 2.592
%      error2 = 2.1503%    std = 1.89%

clc;
clear;
tic
%% Universal initialization of the program

% Parameters setting
peaktol = 0.05; % peak detection tolerance for peakdet()
ftol = 0.1; % I'm going to remove frquencies in (MAf-ftol, MAf+ftol)
removerate = 0.3;
doubleremoverate = 0.3;
isdoubleMA = 0;
avgwin = 4; % Window size of average previous BPM
mawin = 3;
MINBPM = 65;
MAXBPM = 180;
ssatype = 2;

% features
prequessBPM = 120;
quessBPMavg = 0;

% Set dataset index
sdataidx = 2;
edataidx = 2;
error1 = zeros(12,1);
error2 = zeros(12,1);

% datatype:
%   1: data with answer
%   2: test data
datatype = 1;

% ffttype:
%   1: rectangular
%   2: gaussian window
%   3: multi-window
%   4: two-channel
ffttype = 1;
fftsigma = 0.01;
fftwindow = 4;

% display result figure or not
isdisplay = 1;

%***********************************************************************
for dataidx = sdataidx:edataidx
%% Initialization of each dataset
% load data
if datatype == 1
    if dataidx ~= 1
        dataname = sprintf('./dataset/DATA_%02d_TYPE02.mat',dataidx);
        load(dataname);
        dataname = sprintf('./dataset/DATA_%02d_TYPE02_BPMtrace.mat',dataidx);
        load(dataname);
    else
        dataname = sprintf('./dataset/DATA_%02d_TYPE01.mat',dataidx);
        load(dataname);
        dataname = sprintf('./dataset/DATA_%02d_TYPE01_BPMtrace.mat',dataidx);
        load(dataname);
    end;
else
    dataname = './dataset/TEST_S01_T01.mat';
    load(dataname);
end;

% The original signal

if datatype == 1
    ECG = sig(1,:);
    PPG1 = sig(2,:);
    PPG2 = sig(3,:);
    MA1 = sig(4,:);
    MA2 = sig(5,:);
    MA3 = sig(6,:);
    bpm = BPM0;
elseif datatype == 2
    PPG1 = sig(1,:);
    PPG2 = sig(2,:);
    MA1 = sig(3,:);
    MA2 = sig(4,:);
    MA3 = sig(5,:);
end;

sollen = floor(length(PPG1)/250 - 3);
% the BPM calculated by SRPS framwork, our estimate BPM will be stored in
% this vector
BPM = [];
ii = 1;
under = 0;
over = 0;
t0 = 0; % The start time, default value is 0
% Status of signal: currently unused
%  1: stationary low frequency
%  2: stationary high frequency
%  3: accelerating
status = 0;

% frequency data
fdata = [];
Adata = [];

% MA power init
MAmaxp1 = [];
MAmaxp2 = [];
MAmaxp3 = [];
maA1 = zeros(1,sollen);
maA2 = zeros(1,sollen);
maA3 = zeros(1,sollen);
maA = zeros(1,sollen);
maAgrad = zeros(1,sollen);

% Algorithm Begins
fprintf('*** SRPS Algorithm Begins ***\n');
% debug figure init
bpmh = figure(1);
clf(1);
title('BPM');
peakh = figure(2);
clf(2);
title('Peak');

for t = t0:2:sollen*2-2

% Previous Information
%   prequessBPM: the estimated BPM in the last window
%   prequessavgBPM: the average estimated BPM in the last avgwin windows
niter = (t+2)/2;
if niter > 1
    prequessBPM = BPM(end);
    if niter > avgwin
        prequessavgBPM = sum(BPM(niter-avgwin:niter-1))/avgwin;
    end;
end;
    
% Initial setting of this window
wl = 8; % window length regulated by dataset
ts = t; % start time
te = ts+wl; % end time
% Signals in this window
x1 = PPG1(125*ts+1:125*te);
x2 = PPG2(125*ts+1:125*te);

m1 = MA1(125*ts+1:125*te);
m2 = MA2(125*ts+1:125*te);
m3 = MA3(125*ts+1:125*te);

% Band Pass
% x1 = filter(bp,x1);
% x2 = filter(bp,x2);

%%  MA processing
% Calculate dominant frequency of MAs, MAfi records the dominant frequency
% of each MA, and MApi records the power of the dominant frquency.
[MAf1, MAp1] = MADominantFrequency(m1);
[MAf2, MAp2] = MADominantFrequency(m2);
[MAf3, MAp3] = MADominantFrequency(m3);
% Collect the frequencies
MAf = unique([MAf1;MAf2;MAf3]);
maA1(t/2+1) = MAp1;
maA2(t/2+1) = MAp2;
maA3(t/2+1) = MAp3;
maA(t/2+1) = (maA1(t/2+1)+maA2(t/2+1)+maA3(t/2+1))/3;

if niter >= mawin
     maA(niter-mawin+1:niter) = sum(maA(niter-mawin+1:niter))/mawin*ones(1,mawin);
end;

%% Stage 2
%[xica1, newPrev] = icaForSSA(sig, t, prev);
%prev = newPrev;
% psd frquency
%x = linspace(0, 62.5, 1000);
% Singular Spectrum Analysis
% The component parts I am going to select after the decomposition
per = 0.8;
L = 40; % Size of trajectory matrix
if ssatype == 1
    L = 80;
    xssa1 = HLssa(x1,sig,t,prequessBPM/60,L,per); % Doing SSA on the first PPG channel
    xssa2 = HLssa(x2,sig,t,prequessBPM/60,L,per); % Doing SSA on the second PPG channel
    ffttype = 1;
else
    I = [1,2];
    xssa1 = ssaSignal(x1,L,I);
    xssa2 = ssaSignal(x2,L,I);
end;
x1 = xssa1;
x2 = xssa2;

%% Using FFT to look at the spectrum of this window after SSA
% Preprocess before FFT
if ffttype == 3
    % FFT
    T = 1/125; % sampling rate
    Fs = 1/T;
    L = 5000;
    tt = (0:L-1)*T;
    Y11 = fft(x1(1:500),L)/L;
    f11 = Fs*linspace(0,1,L);
    Y12 = fft(x1(501:1000),L)/L;
    f12 = Fs*linspace(0,1,L);
    Y21 = fft(x2(1:500),L)/L;
    f21 = Fs*linspace(0,1,L);
    Y22 = fft(x2(501:1000),L)/L;
    f22 = Fs*linspace(0,1,L);

    f11 = f11(1:L/10); % frequency data
    A11 = 2*abs(Y11(1:L/10)); % amplitude data 
    f12 = f12(1:L/10); % frequency data
    A12 = 2*abs(Y12(1:L/10)); % amplitude data 
    f21 = f21(1:L/10); % frequency data
    A21 = 2*abs(Y21(1:L/10)); % amplitude data 
    f22 = f22(1:L/10); % frequency data
    A22 = 2*abs(Y22(1:L/10)); % amplitude data     
    
    f = zeros(4,length(f11));
    A = zeros(4,length(f11));
    f(1,:) = f11;
    f(2,:) = f12;
    f(3,:) = f21;
    f(4,:) = f22;
    A(1,:) = A11;
    A(2,:) = A12;
    A(3,:) = A21;
    A(4,:) = A22;
    
else
    if ffttype == 2
        x = x1;
        for fftt = 1:1000
            x1(fftt) = (fftsigma)^(1/4)*x1(fftt)*exp(-fftsigma*pi*((fftt - 500)/1000)^2);  
        end;    
    end;
    % FFT
    T = 1/125; % sampling rate
    Fs = 1/T;
    L = 10000;
    tt = (0:L-1)*T;
    Y1 = fft(x1,L)/L;
    f1 = Fs*linspace(0,1,L);
    Y2 = fft(x2,L)/L;
    f2 = Fs*linspace(0,1,L);

    if ffttype == 1
        f = zeros(1,length(f1(1:L/10)));
        A = zeros(1,length(f1(1:L/10)));    
        f(1,:) = f1(1:L/10); % frequency data
        A(1,:) = 2*abs(Y1(1:L/10))'; % amplitude data
        %A(1,:) = 0.5*( 2*abs(Y1(1:L/10))' + 2*abs(Y2(1:L/10))' );
    elseif ffttype == 4
        f = zeros(2,length(f1(1:L/10)));
        A = zeros(2,length(f1(1:L/10)));    
        f(1,:) = f1(1:L/10); % frequency data
        A(1,:) = 2*abs(Y1(1:L/10))'; % amplitude data
        f(2,:) = f2(1:L/10); % frequency data
        A(2,:) = 2*abs(Y2(1:L/10)); % amplitude data 
    end;
end;
%% temp test
%{
if t>5 
location=find(f>Np);
location=location(1);
filtert=zeros(1,size(A,2));
filtert( 1,max(1,location-10):location+10 )=1;
A(1,:)=A(1,:).*filtert;
%A(2,:)=A(2,:).*filtert;
%YY2=[YY2;A];
end
%}
%% Stage 3
% Finding peak
numf = size(f,1);
peakf = zeros(numf,20);
peakA = zeros(numf,20);
% Finding peaks in the spectrum
for jj = 1:numf
    [maxfpeak,minfpeak] = peakdet(A(jj,:),peaktol); % peak detection
    peakf(jj,1:length(f(maxfpeak(:,1)))) = f(jj,maxfpeak(:,1))';
    peakA(jj,1:length(f(maxfpeak(:,1)))) = maxfpeak(:,2)';
end;

% Peak filter with MA dominant frquencies
%
if prequessBPM > 100
matol = 0.2;
maidx = [];
for jj = 1:length(MAf)
    for kk = 1:numf
        % remove peak
        maidx = find(abs(peakf(kk,:)-MAf(jj)) < matol);
        peakA(kk,maidx) = peakA(kk,maidx)*removerate;
        % recover peaks that are closed
        if niter > 1
            maidxclose = find(abs(peakf(kk,maidx) - prequessBPM/60) < matol*1);
            for maidxclosei = 1:length(maidxclose)
                if abs(peakf(kk,maidx(maidxclose(maidxclosei))) - MAf(jj)) < matol
                    peakA(kk,maidx(maidxclose(maidxclosei))) = peakA(kk,maidx(maidxclose(maidxclosei)))*2;
                end;
            end;
        end;
        % remove double frequency
        if isdoubleMA == 1
            maidx = find(abs(peakf(kk,:)-MAf(jj)*2) < matol);
            peakA(kk,maidx) = peakA(kk,maidx)*doubleremoverate;
            % recover peaks that are closed
            if niter > 1
                maidxclose = find(abs(peakf(kk,maidx) - prequessBPM/60) < matol*1);
                for maidxclosei = 1:length(maidxclose)
                    if abs(peakf(kk,maidx(maidxclose(maidxclosei))) - MAf(jj)*2) < matol
                        peakA(kk,maidx(maidxclose(maidxclosei))) = peakA(kk,maidx(maidxclose(maidxclosei)))*2;
                    end;
                end;
            end;   
        end;
%         maidx = find(abs(peakf(kk,:)-MAf(jj)) < matol);
%         peakA(kk,maidx) = peakA(kk,maidx)*removerate;
%         maidx = find(abs(f(kk,:)-MAf(jj)) < matol);
%         if niter > avgwin && abs(MAf(jj)-prequessavgBPM/60) < 0.1
%             A(kk,maidx) = A(kk,maidx);
%         else
%             A(kk,maidx) = A(kk,maidx)*removerate;
%         end;
    end;
end;

end
%}
% Remove dominant frequency of MA from PPG
% Currently unused

%% Peak Selection
if niter <= avgwin % Initial Guess
    guessf = peakselectioninit2(peakf,peakA,f,A,0);
    guessBPM = guessf(1) * 60;
else
    maAgrad(niter) = (maA(niter-1)-maA(niter-mawin+1-1))/mawin;
    if niter > 78 && niter < 145 && false
        guessf = peakselection2(niter,peakf,peakA,f,A,prequessBPM,prequessavgBPM,maA,bpm(niter),MAf,1,peakh);
        
    else
        guessf = peakselection2(niter,peakf,peakA,f,A,prequessBPM,prequessavgBPM,maA,bpm(niter),MAf,0,peakh);
    end;
%     if abs(guessf(1) - prequessBPM) < abs(guessf(2) - prequessBPM)
%         guessBPM = guessf(1) * 60;
%     else
%         guessBPM = guessf(2) * 60;
%     end;
    guessBPM = guessf(1) * 60;
end;


% Optimization on guessed BPM
if niter > avgwin
    BPMtol = 6;
    preBPM = BPM(end);
    if guessBPM > preBPM + BPMtol
        guessBPM = preBPM + BPMtol;
        over = over + 1;
        under = 0;
        if over >= 5
            BPMtol = 6;
        else
            BPMtol = 4;
        end;
    elseif guessBPM < preBPM - BPMtol
        guessBPM = preBPM - BPMtol;
        over = 0;
        under = under + 1;
        if under >= 3
            BPMtol = 6;
        else
            BPMtol = 4;
        end;
    end;
    
end;

% Display banner
if mod(niter,10) == 1 && isdisplay == 1
    fprintf('#win  diff  diff val  estimated   real\n');
end;

% guess BPM
if niter > 1
    alpha = 0.75;
    guessBPM = (1-alpha)*BPM(end) + alpha*guessBPM;
end;

BPM = [BPM;guessBPM];
ii = ii+1;
if isdisplay == 1
    fprintf('%4d ', ii);
    fprintf('%f\n',BPM(end));
end;

figure(bpmh); hold on;
plot(niter,BPM(niter),'bo');
plot(niter,bpm(niter),'ro')
% legend('guess','real');
hold off;

% % store f and A
% fdata = [fdata,f1'];
% Adata = [Adata,A1];

end;

%% Error Calculation and analysis
if datatype == 1
    fprintf('err1 = %f\n',err1(BPM,bpm(1:length(BPM))));
    fprintf('err2 = %f\n',err2(BPM,bpm(1:length(BPM))));
    error1(dataidx) = err1(BPM,bpm(1:length(BPM)));
    error2(dataidx) = err2(BPM,bpm(1:length(BPM)));
end;

%% Sowing the real BPM and estimated BPM
if isdisplay == 1
    if datatype == 2
        plot(1:length(BPM), BPM);
    else
        showestimatedBPM(bpm,BPM,dataidx);
    %   Showing the spectrum, real BPM, and estimated BPM of every window, and
    %   display in movie
    %     showfandAdata(fdata,Adata,length(BPM),bpm,BPM,0.5);
    end;
end;
    
end;
toc

%% Error calculation
if datatype == 1
    fprintf('avg error1: %f\n', sum(error1(sdataidx:edataidx)/(edataidx-sdataidx+1)));
    fprintf('avg error2: %f\n', sum(error2(sdataidx:edataidx)/(edataidx-sdataidx+1)));
end;
TEMP=[BPM0 BPM];