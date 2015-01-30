% Program: peakselectioninit2.m
function [guessf] = peakselectioninit2(peakf,peakA,f,A, isplot)

% test:
% peakselection1(peakf,peakA,f',A,bpm((t+2)/2),prequessBPM,prequessavgBPM)

% input:
%   peakf
%   peakA
%   f
%   A
%   bpm: real BPM, uesd for plotting compared to estimated BPM
%   maxA
%   isplot

f = f(1,:);
A = A(1,:);
peakf = peakf(1,:);
peakA = peakA(1,:);


if isplot ~= 0
    figure();
    stem(f,A,'o'); hold on;
end;
    
%% Peak Selection
heartftol = 0.5;
heartf = 1.5;
maxA = max(A((f < heartf + heartftol) & (f > heartf - heartftol)));
heartA = maxA;
heartAtol = maxA/2;

% Select Peak

heartindex = find((peakA > heartA - heartAtol));
fcandidate = peakf(heartindex);


% Approach 2: Select the peak with highest power in the heart frequency
% zone (heartf - hearttol, heartf + hearttol) near the previous estimated
% BPM
if ~isempty(fcandidate)
    delta = psapp2(fcandidate(1),A(f == fcandidate(1)));
    findex = 1;
    for j = 1:length(fcandidate)
        if psapp2(fcandidate(j),A(f == fcandidate(j))) < delta
            delta = psapp2(fcandidate(j),A(f == fcandidate(j)));
            findex = j;
        end;
    end;
    guessf = fcandidate(findex);
else
    guessf = f(A == maxA);
end;


if isplot ~= 0
    plot(peakf(heartindex),peakA(heartindex),'r*');
    plot(guessf,A(f==guessf),'g*');
    plot(bpm/60,0,'ro');
end;

end

% Approach 2 cost function
function [err] = psapp2(f,A)
lambdaA = 10;
sigmaf = 0.5;
meanf = 1.5;
lambdaf = 20;
errf = lambdaf / (exp(-(f-meanf)^2/(2*sigmaf)));
err = lambdaA/A + errf;
end