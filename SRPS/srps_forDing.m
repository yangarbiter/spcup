% Program: srps.m
% Author: Jerry Chou
% Last Modified Date: 2015/1/5
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

% Parameters setting
peaktol = 0.05; % peak detection tolerance for peakdet()
ftol = 0.1; % I'm going to remove frquencies in (MAf-ftol, MAf+ftol)
removerate = 0.3;
avgwin = 4; % Window size of average previous BPM
mawin = 3;
MINBPM = 65;
MAXBPM = 180;

% features
prequessBPM = 0;
quessBPMavg = 0;

%% Set dataset
sdataidx = 1;
edataidx = 12;
error1 = zeros(12,1);
error2 = zeros(12,1);

for dataidx = sdataidx:edataidx
%% Initialization
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
bpm = BPM0; % the answer BPM
% The original signal
ECG = sig(1,:);
PPG1 = sig(2,:);
PPG2 = sig(3,:);
MA1 = sig(4,:);
MA2 = sig(5,:);
MA3 = sig(6,:);
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

% MA power
MAmaxp1 = [];
MAmaxp2 = [];
MAmaxp3 = [];
maA1 = zeros(1,length(bpm));
maA2 = zeros(1,length(bpm));
maA3 = zeros(1,length(bpm));
maA = zeros(1,length(bpm));
maAgrad = zeros(1,length(bpm));

%% Algorithm Begins
fprintf('*** SRPS Algorithm Begins ***\n');
for t = t0:2:2*length(BPM0)-2

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
y1 = ECG(125*ts+1:125*te);
m1 = MA1(125*ts+1:125*te);
m2 = MA2(125*ts+1:125*te);
m3 = MA3(125*ts+1:125*te);

% Band Pass
% x1 = filter(bp,x1);
% x2 = filter(bp,x2);

%%  MA
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

%% Singular Spectrum Analysis
% The component parts I am going to select after the decomposition
I = [1,2];
L = 40; % Size of trajectory matrix
xssa1 = ssaSignal(x1,L,I); % Doing SSA on the first PPG channel
x1 = xssa1;
xssa2 = ssaSignal(x2,L,I); % Doing SSA on the second PPG channel
x2 = xssa2;

%% Using FFT to look at the spectrum of this window after SSA
% FFT
T = 1/125; % sampling rate
Fs = 1/T;
L = 10000;
tt = (0:L-1)*T;
Y1 = fft(x1,L)/L;
f1 = Fs*linspace(0,1,L);
Y2 = fft(x2,L)/L;
f2 = Fs*linspace(0,1,L);



f1 = f1(1:L/10); % frequency data
A1 = 2*abs(Y1(1:L/10)); % amplitude data 
f2 = f2(1:L/10); % frequency data
A2 = 2*abs(Y2(1:L/10)); % amplitude data 

% Finding peaks in the spectrum
[maxfpeak1,minfpeak1] = peakdet(A1,peaktol); % peak detection
peakf1 = f1(maxfpeak1(:,1));
peakA1 = maxfpeak1(:,2);
[maxfpeak2,minfpeak2] = peakdet(A2,peaktol); % peak detection
peakf2 = f2(maxfpeak2(:,1));
peakA2 = maxfpeak2(:,2);

% Peak filter with MA dominant frquencies
matol = 0.2;
maidx = [];
for jj = 1:length(MAf)
    maidx = find(abs(peakf1-MAf(jj)) < matol);
    peakA1(maidx) = peakA1(maidx)*removerate;
    maidx = find(abs(f1-MAf(jj)) < matol);
    if niter > avgwin && abs(MAf(jj)-prequessavgBPM/60) < 0.1
        A1(maidx) = A1(maidx);
    else
        A1(maidx) = A1(maidx)*removerate;
    end;
end;

% Remove dominant frequency of MA from PPG
% Currently unused

%% Peak Selection
if niter <= avgwin % Initial Guess
    guessf1 = peakselectioninit1(peakf1,peakA1,f1',A1,bpm((t+2)/2),0);
    % guessf2 = peakselectioninit1(peakf2,peakA2,f2',A2,bpm((t+2)/2),0);
else
    maAgrad(niter) = (maA(niter-1)-maA(niter-mawin+1-1))/mawin;
    if mod(niter,5) == 3 && niter > 80 && niter < 110
        guessf1 = peakselection1(niter,peakf1,peakA1,f1',A1,bpm((t+2)/2),prequessBPM,prequessavgBPM,maA,0);
        % guessf2 = peakselection1(niter,peakf2,peakA2,f2',A2,bpm((t+2)/2),prequessBPM,prequessavgBPM,maA,0);
    else
        guessf1 = peakselection1(niter,peakf1,peakA1,f1',A1,bpm((t+2)/2),prequessBPM,prequessavgBPM,maA,0);
        % guessf2 = peakselection1(niter,peakf2,peakA2,f2',A2,bpm((t+2)/2),prequessBPM,prequessavgBPM,maA,0);
    end;
end;
% guessf = (guessf1+guessf2)/2;
guessf = guessf1;
guessBPM = guessf * 60;

%% Optimization on guessed BPM
if niter > avgwin
    BPMtol = 4;
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

% store f and A
fdata = [fdata,f1'];
Adata = [Adata,A1];

% Display banner
if mod(niter,10) == 1
    fprintf('#win  diff  diff val  estimated   real\n');
end;

% BPM
BPM = [BPM;guessBPM];
fprintf('%4d ', ii);
if BPM(end) > bpm(ii)
    fprintf(' over  ');
else
    fprintf(' under ');
end;
fprintf('%f  %f  %f\n',abs(BPM(end)-bpm(ii)),BPM(end),bpm(ii));
ii = ii+1;
    
end;



%% Error Calculation and analysis
fprintf('err1 = %f\n',err1(BPM,bpm(1:length(BPM))));
fprintf('err2 = %f\n',err2(BPM,bpm(1:length(BPM))));
error1(dataidx) = err1(BPM,bpm(1:length(BPM)));
error2(dataidx) = err2(BPM,bpm(1:length(BPM)));

% Sowing the real BPM and estimated BPM
showestimatedBPM(bpm,BPM,dataidx);

% Showing the spectrum, real BPM, and estimated BPM of every window, and
% display in movie
% showfandAdata(fdata,Adata,length(BPM),bpm,BPM,0.5);

end;

%% Error calculation
fprintf('avg error1: %f\n', sum(error1(sdataidx:edataidx)/(edataidx-sdataidx+1)));
fprintf('avg error2: %f\n', sum(error2(sdataidx:edataidx)/(edataidx-sdataidx+1)));
