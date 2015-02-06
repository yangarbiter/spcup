% Program: peakselectioninit2.m
function [guessf] = peakselectioninit2(peakf,peakA,f,A,bpm,MAf,MAp,initcosttype,isplot,peakh)

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

% f = f(1,:);
% A = A(1,:);
% peakf = peakf(1,:);
% peakA = peakA(1,:);

numchannel = size(f,1);
    
%% Peak Selection
heartftol = 0.8;
heartf = 2;

for c = 1:numchannel
    maxA(c) = max(peakA(c,(peakf(c,:) < heartf + heartftol) & (peakf(c,:) > heartf - heartftol)));
    heartA(c) = maxA(c);
    heartAtol(c) = maxA(c)/2;
end

% Select Peak
fcandidate = zeros(numchannel,20);
heartindex = zeros(numchannel,20);
for c = 1:numchannel
    tmpidx = find((peakA(c,:) > heartA(c) - heartAtol(c)));
    heartindex(c,1:length(tmpidx)) = tmpidx;
    fcandidate(c,1:length(tmpidx)) = peakf(c,tmpidx);
end

% Approach 2: Select the peak with highest power in the heart frequency
% zone (heartf - hearttol, heartf + hearttol) near the previous estimated
% BPM
for c = 1:numchannel
    if ~isempty(fcandidate(c,:))
        delta = initcost(fcandidate(c,1),A(c,f(c,:) == fcandidate(c,1)),initcosttype);
        findex = 1;
        for j = 1:length(fcandidate(c,:))
            if initcost(fcandidate(c,j),A(c,f(c,:) == fcandidate(c,j)),initcosttype) < delta
                delta = initcost(fcandidate(c,j),A(c,f(c,:) == fcandidate(c,j)),initcosttype);
                findex = j;
            end;
        end;
        guessf(c) = fcandidate(c,findex);
    else
        guessf(c) = f(c,A(c,:) == maxA);
    end;
end

if isplot ~= 0
    clf(peakh);
    figure(peakh); hold on;
    for c = 1:numchannel
        subplot(2,1,c); hold on;
        stem(f(c,1:300),A(c,1:300),'o');
        stem(MAf,MAp/max(MAp),'go');
        hidx = find(heartindex(c,:)~=0);
        if ~isempty(heartindex(c,hidx))
            plot(peakf(c,heartindex(c,hidx)),peakA(c,heartindex(c,hidx)),'r*');
        end;
        plot(guessf,A(c,f(c,:)==guessf(c)),'g*');
        plot(bpm/60,0,'ro'); hold off;
    end
%     pause();
    hold off;
end;

end

% Approach 2 cost function
function [err] = initcost(f,A,initcosttype)
if initcosttype == 1
    err = initcost1(f,A);
elseif initcosttype == 2
    err = initcost2(f,A);
end

end

function [err] = initcost1(f,A)
lambdaA = 10;
sigmaf = 0.5;
meanf = 1.5;
lambdaf = 20;
errf = lambdaf / (exp(-(f-meanf)^2/(2*sigmaf)));
err = lambdaA/A + errf;
end

function [err] = initcost2(f,A)
lambdaA = 10;
meanf = 1.8;
widthf = 0.9;
err = lambdaA/A;
if abs(f-meanf) > widthf
    err = 9999;
end
end