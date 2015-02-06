% Program: peakselection2.m
% Author: Jerry Chou
% Last Modified: 2015/2/5
%
% This program select an estimated BPM from the peak candidates given from
% the SRPS algorithm with peak informations, previous guess, motion
% artifect informations
%
function [guessf] = peakselection2_0(niter,peakf,peakA,f,A, preBPM,preavgBPM,maA,realbpm,MAf,MAp,isplot,peakh)

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

% cost f
costftype = 1;
% cost A
costAtype = 2;

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
fcost = ones(numchannel,50) * 9999;
Acost = ones(numchannel,numchannel,50) * 9999;
cost = ones(numchannel,50) * 9999;
for c = 1:numchannel
    if ~isempty(fcandidate(c,:))
        for j = 1:length(fcandidate(c,:))
            if fcandidate(c,j) == 0
                continue;
            end;
            fcost(c,j) = costf(fcandidate(c,j),preBPM,preavgBPM,magrad,costftype);
            for k = 1:numchannel
                Acost(c,k,j) = costA(peakA(c,heartindex(c,j))/max(peakA(:)),magrad,costAtype);
            end;
            cost(c,j) = fcost(c,j) + 1.5*((Acost(c,c,j))+sum(Acost(c,:,j)))/(numchannel+1);
        end;
        guessf(c,1) = fcandidate(c,find(cost(c,:) == min(cost(c,:))));
    else
        guessf(c,1) = f(c,A(c,:) == maxA);
    end;
end;

% plot
if isplot ~= 0
    clf(peakh);
    figure(peakh);
    for c = 1:numchannel
        subplot(numchannel,1,c); hold on;
        stem(f(c,1:300),A(c,1:300),'o');
        stem(MAf,ones(size(MAf)),'go');
        hidx = find(heartindex(c,:)~=0);
        if ~isempty(heartindex(c,hidx))
            plot(peakf(c,heartindex(c,hidx)),peakA(c,heartindex(c,hidx)),'r*');
        end;
        plot(guessf(c,1),A(c,f(c,:) == guessf(c,1)),'g*');
        plot(preBPM/60,0,'go');
        plot(preavgBPM/60,0,'yo');
        plot(realbpm/60,0,'ro'); hold off;
    end;
    pause();
    hold off; 
end;

end % function end

% cost: frequency
% type:
%    1: linear
%    2: square
function [cost] = costf(f,preBPM,preavgBPM,magrad,type)

if type == 1
    cost = costflinear(f,preBPM,preavgBPM,magrad);
elseif type == 2
    cost = costfsquare(f,preBPM,preavgBPM,magrad);
end;

end

% cost f: linear model
function [cost] = costflinear(f,preBPM,preavgBPM,magrad)
lambdaf = 8;
lambdafavg = 10;
lambdafgrad = 6;
lambdafavggrad = 8;
BPMdiff = 1/2;
if magrad > 0.0015 % increasing
    preBPM = preBPM + 5;
    preavgBPM = preavgBPM + 5;
    df = abs(f-preBPM/60) / BPMdiff;
    dfavg = abs(f-preavgBPM/60) / BPMdiff;
    if f > preBPM/60
        costf = lambdafgrad*df / 2;
    else
        costf = lambdafgrad*df;
    end;
    if f > preavgBPM/60
        costfavg = lambdafavggrad*dfavg / 2;
    else
        costfavg = lambdafavggrad*dfavg;
    end;
elseif magrad < -0.0015 % decreasing
    df = abs(f-preBPM/60) / BPMdiff;
    dfavg = abs(f-preavgBPM/60) / BPMdiff;
    costf = lambdafgrad*df;
    costfavg = lambdafavggrad*dfavg;
else % normal
    df = abs(f-preBPM/60) / BPMdiff;
    dfavg = abs(f-preavgBPM/60) / BPMdiff;
    costf = lambdaf*df;
    costfavg = lambdafavg*dfavg;    
end;

cost = costf + costfavg;
% exceed threshold
if abs(f-preBPM/60) > BPMdiff
    cost = 9999;
end;

end

% cost f: square model
function [cost] = costfsquare(f,preBPM,preavgBPM,magrad)
magradupper = 0.0015;
magradlower = -0.0015;
lambdaf = 8;
lambdafavg = 10;
lambdafgrad = 6;
lambdafavggrad = 8;
BPMdiff = 1/2; % When the difference between fcandidate and previous guess is larger than this threshold, we set the cost to be very large 

df = abs(f-preBPM/60) / BPMdiff;
dfavg = abs(f-preavgBPM/60) / BPMdiff;

if magrad > magradupper
    costf = lambdaf * df.^2;
    costfavg = lambdafavg * dfavg.^2;
    if f < preBPM/60
        costf = costf + lambdafgrad*df;
    end;
    if f < preavgBPM/60
        costfavg = costfavg + lambdafavggrad*dfavg;
    end;
elseif magrad < magradlower
    costf = lambdaf*df.^2;
    costfavg = lambdafavg*dfavg.^2;
    if f > preBPM/60
        costf = costf + lambdafgrad*df;
    end;
    if f > preavgBPM/60
        costfavg = costfavg + lambdafavggrad*dfavg;
    end;
else
    costf = lambdaf*df;
    costfavg = lambdafavg*dfavg;    
end;

cost = costf + costfavg;
% exceed threshold
if abs(f-preBPM/60) > BPMdiff
    cost = 9999;
end;

end

% amplitude cost
function [cost] = costA(A,magrad,type)
if type == 1
    cost = costAexponential(A,magrad);
elseif type == 2
    cost = costAlinear(A,magrad);
end;

end

% A cost: exponential model
function [cost] = costAexponential(A,magrad)
magradupper = 0.0015;
magradlower = -0.0015;
lambda = 12;
if magrad > magradupper
    cost = 1 - exp(-lambda/A);
elseif magrad < magradlower
    cost = 1 - exp(-lambda/A);
else
    lambda = 10;
    cost = 1 - exp(-lambda/A);    
end;

end

% A cost: linear model
function [cost] = costAlinear(A,magrad)
lambda = 10;
cost = lambda/A;
end

%A cost:IDK
function [cost] = costAPML(A,magrad)

end