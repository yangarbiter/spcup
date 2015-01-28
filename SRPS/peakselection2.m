% Program: peakselection2.m
% Author: Jerry Chou
% Last Modified: 2015/1/23
%
% This program select an estimated BPM from the peak candidates given from
% the SRPS algorithm with peak informations, previous guess, motion
% artifect informations
%
function [guessf] = peakselection2(niter,peakf,peakA,f,A, preBPM,preavgBPM,maA,realbpm,isplot,peakh)

% test:
% peakselection1(peakf,peakA,f',A,bpm((t+2)/2),prequessBPM,prequessavgBPM)

% input:
%   peakf: each row is a channel
%   peakA: each row is a channel
%   f: each row is a channel
%   A: each row is a channel
%   bpm: real BPM, uesd for plotting compared to estimated BPM
%   preBPM
%   preavgBPM
%   isplot

% peakA = peakA(1,:);
% peakf = peakf(1,:);
% f = f(1,:);
% A = A(1,:);
numchannel = size(A,1);

% Motion artifect preprocessing
mawin = 4;
magrad = (maA(niter-1)-maA(niter-mawin+1-1))/mawin;

%% Start Peak Selection
heartftol = 0.2;
% Add local max peak
localrange = [1,1.3; 1.3,1.6; 1.6,1.9;1.9,2.2; 2.2,2.5;2.5,2.8;2.8,3.1;3.1,3.4];
% Select Peak Candidate
fcandidate = zeros(numchannel,20);
heartindex = zeros(numchannel,20);
for c = 1:numchannel
    maxA = max(A(c,:));
    heartindex1 = find(peakA(c,:) > maxA/2);
    fcandidate1 = peakf(c,heartindex1);
    heartindex2 = find((peakf(c,:) > preBPM/60 - heartftol) & (peakf(c,:) < preBPM/60 + heartftol));
    fcandidate2 = peakf(c,heartindex2);
    
    % Add local max peak
    fcandidate3 = [];
    heartindex3 = [];
    for i = 1:7
        localmaxA = max(peakA(c,(peakf(c,:) > localrange(i,1)) & (peakf(c,:) < localrange(i,2))));
        if isempty(localmaxA)
            continue;
        end;
        localmaxindex = find(peakA(c,:) == localmaxA);
        fcandidate3 = [fcandidate3,peakf(c,localmaxindex)];
        fcandidate3 = unique(fcandidate3);
        heartindex3 = [heartindex3,localmaxindex];
        heartindex3 = unique(heartindex3);
    end;
    
    tmpfcandidate = [fcandidate1,fcandidate2,fcandidate3];
    tmpheartindex = [heartindex1,heartindex2,heartindex3];
    tmpfcandidate = unique(tmpfcandidate);
    tmpheartindex = unique(tmpheartindex);
    fcandidate(c,1:length(tmpfcandidate)) = tmpfcandidate;
    heartindex(c,1:length(tmpheartindex)) = tmpheartindex;
end;


guessf = zeros(numchannel,1);
fcost = zeros(numchannel,20);
for c = 1:numchannel
    delta = psapp2(fcandidate(c,1),preBPM,preavgBPM,A(c,f(c,:) == fcandidate(c,1)),magrad);
    findex = 1;
    guessf(c,1) = fcandidate(c,1);
    if ~isempty(fcandidate(c,:))
        for j = 1:length(fcandidate(c,:))
            if fcandidate(c,j) == 0
                continue;
            end;
            fcost(c,j) = psapp2(fcandidate(c,j),preBPM,preavgBPM,A(c,f(c,:) == fcandidate(c,j)),magrad);
            if psapp2(fcandidate(c,j),preBPM,preavgBPM,A(c,f(c,:) == fcandidate(c,j)),magrad) < delta
                delta = psapp2(fcandidate(c,j),preBPM,preavgBPM,A(c,f(c,:) == fcandidate(c,j)),magrad);
                findex = j;
                guessf(c,1) = fcandidate(c,j);
            end;
        end;
    else
        guessf(c,:) = f(c,A(c,:) == maxA);
    end;
end;
    
    



if isplot ~= 0
    clf(peakh);
    figure(peakh);
    for c = 1:numchannel
        subplot(numchannel,1,c); hold on;
        stem(f(c,1:300),A(c,1:300),'o'); 
        hidx = find(heartindex(c,:)~=0);
        if ~isempty(heartindex(c,hidx))
            plot(peakf(c,heartindex(c,hidx)),peakA(c,heartindex(c,hidx)),'r*');
        end;
        plot(guessf(c,1),A(c,f(c,:)==guessf(c,1)),'g*');
        plot(preBPM/60,0,'go');
        plot(preavgBPM/60,0,'yo');
        plot(realbpm/60,0,'ro'); hold off;
    end;
    pause();
    hold off;
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
    if abs(f-preBPM/60) > 1/2
        err = 999999;
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
    if abs(f-preBPM/60) > 1/2
        err = 999999;
    end;
else % normal
    lambda = 10;
    err = abs(f-preBPM/60) + abs(f-preavgBPM) + lambda/A;    
end;

end