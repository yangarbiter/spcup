% Program: peakselection1.m
% Author: Jerry Chou
% Last Modified: 2015/1/5
%
% This program select an estimated BPM from the peak candidates given from
% the SRPS algorithm with peak informations, previous guess, motion
% artifect informations
%
function [guessf] = peakselection1(niter,peakf,peakA,f,A, bpm,preBPM,preavgBPM,maA,isplot)

% test:
% peakselection1(peakf,peakA,f',A,bpm((t+2)/2),prequessBPM,prequessavgBPM)

% input:
%   peakf
%   peakA
%   f
%   A
%   bpm: real BPM, uesd for plotting compared to estimated BPM
%   preBPM
%   preavgBPM
%   isplot

if isplot ~= 0
    figure();
    stem(f,A,'o'); hold on;
end;
    
mawin = 4;
magrad = (maA(niter-1)-maA(niter-mawin+1-1))/mawin;
%% Peak Selection
maxA = max(A);
heartA = 2.2;
heartAtol = 0.8;
heartftol = 0.2;

% Select Peak

% heartindex1 = find((peakA > heartA - heartAtol) & (peakA < heartA + heartAtol));
heartindex1 = find(peakA > maxA/2);
fcandidate1 = peakf(heartindex1);
heartindex2 = find((peakf > preBPM/60 - heartftol) & (peakf < preBPM/60 + heartftol));
fcandidate2 = peakf(heartindex2);
fcandidate = [];
heartindex = [];
% combine fcandidate1 & 2
[n,m] = size(fcandidate1);
if n < m
    fcandidate1 = fcandidate1';
    hearindex1 = heartindex1';
end;
[n,m] = size(fcandidate2);
if n < m
    fcandidate2 = fcandidate2';
    heartindex2 = heartindex2';
end;
if ~isempty(fcandidate1) || ~isempty(fcandidate2)
    fcandidate = [fcandidate1;fcandidate2];
    fcandidate = unique(fcandidate);
    heartindex = [heartindex1;heartindex2];
    heartindex = unique(heartindex);
end;

% add local max peak
localrange = [1,1.3; 1.3,1.6; 1.6,1.9;1.9,2.2; 2.2,2.5;2.5,2.8;2.8,3.1;3.1,3.4];
for i = 1:7
    localmaxA = max(peakA((peakf > localrange(i,1)) & (peakf < localrange(i,2))));
    if isempty(localmaxA)
        continue;
    end;
    localmaxindex = find(peakA == localmaxA);
    fcandidate = [fcandidate;peakf(localmaxindex)];
    fcandidate = unique(fcandidate);
    heartindex = [heartindex;localmaxindex];
    heartindex = unique(heartindex);
end;

if ~isempty(fcandidate)
    fcost = zeros(1,length(fcandidate));
    delta = psapp2(fcandidate(1),preBPM,preavgBPM,A(f == fcandidate(1)),magrad);
    findex = 1;
    for j = 1:length(fcandidate)
        fcost(j) = psapp2(fcandidate(j),preBPM,preavgBPM,A(f == fcandidate(j)),magrad);
        if psapp2(fcandidate(j),preBPM,preavgBPM,A(f == fcandidate(j)),magrad) < delta
            delta = psapp2(fcandidate(j),preBPM,preavgBPM,A(f == fcandidate(j)),magrad);
            findex = j;
        end;
    end;
    guessf = fcandidate(findex);
else
    guessf = f(A == maxA);
end;


if isplot ~= 0
    if ~isempty(heartindex)
        plot(peakf(heartindex),peakA(heartindex),'r*');
    end;
    plot(guessf,A(f==guessf),'g*');
    plot(bpm/60,0,'ro');
    plot(preBPM/60,0,'go');
    plot(preavgBPM/60,0,'yo');
    niter
    magrad
    fcost
end;

end

% Approach 2 cost function
function [err] = psapp2(f,preBPM,preavgBPM, A,magrad)
if magrad > 0.0015 % increasing
    lambda = 12;
    lambdaf = 8;
    lambdafgrad = 6;
    err = lambdaf*abs(f-preBPM/60) + lambdaf*abs(f-preavgBPM) + lambda/A;
    if f < preBPM/60
        err = err + lambdafgrad*abs(f-preBPM/60);
    end;
    if f < preavgBPM/60
        err = err + lambdafgrad*abs(f-preavgBPM/60);
    end;
elseif magrad < -0.0015 % decreasing
    lambda = 12;
    lambdaf = 8;
    lambdafgrad = 6;
    err = lambdaf*abs(f-preBPM/60) + lambdaf*abs(f-preavgBPM) + lambda/A;
    if f > preBPM/60
        err = err + lambdafgrad*abs(f-preBPM/60);
    end;
    if f > preavgBPM/60
        err = err + lambdafgrad*abs(f-preavgBPM/60);
    end;
else % normal
    lambda = 10;
    err = abs(f-preBPM/60) + abs(f-preavgBPM) + lambda/A;    
end;

end