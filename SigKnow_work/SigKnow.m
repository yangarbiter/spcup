% Program: SigKnow.m
% Team ID: 20219
% Last Modified Time: 2015/02/01

function [guessBPM, maA, over, under, f, A, MAf] = SigKnow(varargin)

% =============== Parsing Input ===============
p = inputParser;
addOptional(p, 'curSig', @(x) ~isempty(x));
addOptional(p, 'BPM', @(x) ~isempty(x));
addOptional(p, 'maA', @(x) ~isempty(x));
addOptional(p, 'niter', @isnumeric);
addOptional(p, 'avgwin', @isnumeric);
addOptional(p, 'mawin', @isnumeric);
addOptional(p, 'frequencytype', @isnumeric);
addOptional(p, 'peaktol', @isnumeric);
addOptional(p, 'removerate', @isnumeric);
addOptional(p, 'doubleremoverate', @isnumeric);
addOptional(p, 'isdoubleMA', @isnumeric);
addOptional(p, 'isdisplay', @isnumeric);
addOptional(p, 'realbpm', @isnumeric);
addOptional(p, 'over', @isnumeric);
addOptional(p, 'under', @isnumeric);
addOptional(p, 'bpmh', @(x) 1);
addOptional(p, 'peakh', @(x) 1);
addOptional(p, 'isHHT', @isnumeric);
addOptional(p, 'fftsigma', @isnumeric);
addOptional(p, 'isICA', @isnumeric);
addOptional(p, 'isSSA', @isnumeric);
addOptional(p, 'icaPrev', @isnumeric);
addOptional(p, 'recoveryrate', @isnumeric);
addOptional(p, 'peakdettype', @isnumeric);
addOptional(p, 'isneighboring', @isnumeric);
addOptional(p, 'displayiter', @isnumeric);
addOptional(p, 'initcosttype', @isnumeric);
% addOptional(p, '', @isnumeric);
% addOptional(p, '', @isnumeric);
% addOptional(p, '', @isnumeric);


parse(p, varargin{:});

curSig              = p.Results.curSig;
BPM                 = p.Results.BPM;
maA                 = p.Results.maA;
niter               = p.Results.niter;
avgwin              = p.Results.avgwin;
mawin               = p.Results.mawin;
frequencytype             = p.Results.frequencytype;
peaktol             = p.Results.peaktol;
removerate          = p.Results.removerate;
doubleremoverate    = p.Results.doubleremoverate;
isdoubleMA          = p.Results.isdoubleMA;
isdisplay           = p.Results.isdisplay;
realbpm             = p.Results.realbpm;
over                = p.Results.over;
under               = p.Results.under;
bpmh                = p.Results.bpmh;
peakh               = p.Results.peakh;
isHHT               = p.Results.isHHT;
fftsigma            = p.Results.fftsigma;
isICA               = p.Results.isICA;
isSSA               = p.Results.isSSA;
icaPrev             = p.Results.icaPrev;
recoveryrate        = p.Results.recoveryrate;
peakdettype         = p.Results.peakdettype;
isneighboring       = p.Results.isneighboring;
displayiter         = p.Results.displayiter;
initcosttype        = p.Results.initcosttype;

if niter > 1
    preguessBPM = BPM(niter-1);
    if niter > avgwin
        preguessavgBPM = sum(BPM(niter-avgwin:niter-1))/avgwin;
    end;
end;

Fs=125;
% Signals in this window
if frequencytype ~= 5
    x1 = curSig(1,:); x2 = curSig(2,:);
    m1 = curSig(3,:); m2 = curSig(4,:); m3 = curSig(5,:);
else
    lwa=8;
    lw = (8+lwa) * Fs;   lw1 = (4+lwa) * Fs;  lw2 = 4 * Fs;
    gw1 = exp(-15*(([-lw1:-1]+0.5)/lw1).^2);
    gw2 = exp(-4*(([1:lw2]-0.5)/lw2).^2);
    gw=[gw1,gw2]; 
    curSig1 =curSig;   curSig = curSig(:,lwa * Fs+1: (8+lwa) * Fs);
    x1 = curSig1(1,:).*gw; x2 = curSig1(2,:).*gw;
    m1 = curSig1(3,:).*gw; m2 = curSig1(4,:).*gw; m3 = curSig1(5,:).*gw;
%     x1c = curSig(1,:); x2c = curSig(2,:);
%     m1c = curSig(3,:); m2c = curSig(4,:); m3c = curSig(5,:);
end;
%% Stage 1
% ======================= HHT ===========================
if isHHT == 1
    % temp=x1;
    % x1=x2;
    % x2=temp;
    % Hilbert-Huang Transform
    %Sorting range of IMF
    if(niter == 1)
        mid = 12;
    else
        mid = round(preguessBPM*8/60);
    end
    lb = mid-5;
    ub = round(mid*2.2+1);
    energyx1 = abs(fft(x1)).^2;
    energyx2 = abs(fft(x2)).^2;

    %wx=((sum(abs(energyx1(1,lb:ub)))/sum(abs(energyx1(1,:))))/(sum(abs(energyx2(1,lb:ub)))/sum(abs(energyx2(1,:)))));
    inData1=(x1)';
    % if ((sum(abs(energyx1(1,lb:lb+10)))/sum(abs(energyx1(1,:))))>(sum(abs(energyx2(1,lb:lb+10)))/sum(abs(energyx2(1,:)))))
    % inData=x1';
    % else
    % inData=x2';
    % end
    rslt = eemd(inData1,0.3,5);
    l = length(rslt(1,:));
    rslt = rslt(:,2:l)';
    l = l-1;
    rslt = rslt';
    energy = (abs(fft(rslt))).^2;
    [b,idx] = sort(sum(abs(energy(lb:ub,:)))./sum(abs(energy)),2,'descend');
    rslt = rslt';
    for i = 1:l-5
        if(idx(i) == l||idx(i) == 1)
           if(idx(i+1) == l||idx(i+1) == 1)
            k = i:l-2;
            idx(k) = idx(k+2);
           else
            k = i:l-1;
            idx(k) = idx(k+1);
           end
        end
    end
    %Weighted Coefficient
    % w=max(rslt(idx(4),:))/max(max(rslt(idx(1:3),:)));
    % if(w>1)
    %     w=1/w;
    % end
    %z1=(sum(abs(energy(lb:ub,idx(1))))/sum(abs(energy(:,idx(1)))))*rslt(idx(1),:)+(sum(abs(energy(lb:ub,idx(2))))/sum(abs(energy(:,idx(2)))))*rslt(idx(2),:)+(sum(abs(energy(lb:ub,idx(3))))/sum(abs(energy(:,idx(3)))))*rslt(idx(3),:)+(sum(abs(energy(lb:ub,idx(4))))/sum(abs(energy(:,idx(4)))))*rslt(idx(4),:);
    z1=rslt(idx(1),:)+rslt(idx(2),:)+rslt(idx(3),:)+rslt(idx(4),:);
    x1=z1;
    % Band Pass
    % x1 = filter(bp,x1);
    % x2 = filter(bp,x2);
end;
    
    
% ========================= MA =========================
% Calculate dominant frequency of MAs, MAfi records the dominant frequency
% of each MA, and MApi records the power of the dominant frquency.
[MAf1, MApmax1, MAp1] = MADominantFrequency(m1);
[MAf2, MApmax2, MAp2] = MADominantFrequency(m2);
[MAf3, MApmax3, MAp3] = MADominantFrequency(m3);
% Collect the frequencies
MAf = unique([MAf1;MAf2;MAf3]);
MAp = [];
for fi = 1:length(MAf)
    count = 0;
    for mi = 1:length(MAf1)
        if MAf(fi) == MAf1(mi)
            count = count + MAp1(mi)^2;
        end
    end
    for mi = 1:length(MAf2)
        if MAf(fi) == MAf2(mi)
            count = count + MAp2(mi)^2;
        end
    end
    for mi = 1:length(MAf3)
        if MAf(fi) == MAf3(mi)
            count = count + MAp3(mi)^2;
        end
    end
    MAp = [MAp;sqrt(count)];
end

maA1 = MApmax1;
maA2 = MApmax2;
maA3 = MApmax3;
maA(niter) = (maA1 + maA2 + maA3)/3;

if niter >= mawin
     maA(niter-mawin+1:niter) = sum(maA(niter-mawin+1:niter))/mawin*ones(1,mawin);
end;

%% Stage 2
% ============== Singular Spectrum Analysis ===================
% The component parts I am going to select after the decomposition
I = [1,2];
L = 40; % Size of trajectory matrix
xssa1 = ssaSignal(x1,L,I); % Doing SSA on the first PPG channel
x1 = xssa1;
xssa2 = ssaSignal(x2,L,I); % Doing SSA on the second PPG channel
x2 = xssa2;

% ======================== FFT ==========================
% Using FFT to look at the spectrum of this window after SSA
% Preprocess before FFT

if frequencytype == 2
    x = x1;
    for fftt = 1:1000
        x1(fftt) = (fftsigma)^(1/4)*x1(fftt)*exp(-fftsigma*pi*((fftt - 500)/1000)^2);  
    end;    
end;
% FFT
if frequencytype == 5
    L = 15000;
    Y1 = fft(x1,L)/L; f1 = [0:L-1]/L*Fs;
    Y2 = fft(x2,L)/L; f2 = f1;
else
    L = 10000;
    Y1 = fft(x1,L)/L; f1 = Fs*linspace(0,1,L);
    Y2 = fft(x2,L)/L; f2 = Fs*linspace(0,1,L);
end;
    
if frequencytype == 1
    f = zeros(1,length(f1(1:L/10)));
    A = zeros(1,length(f1(1:L/10)));    
    f(1,:) = f1(1:L/10);            % frequency data
    A(1,:) = 2*abs(Y1(1:L/10))';    % amplitude data
elseif frequencytype == 3
    f = zeros(2,length(f1(1:L/10)));
    A = zeros(2,length(f1(1:L/10)));    
    f(1,:) = f1(1:L/10);            % frequency data
    A(1,:) = 2*abs(Y1(1:L/10))';    % amplitude data
    f(2,:) = f2(1:L/10);            % frequency data
    A(2,:) = 2*abs(Y2(1:L/10));     % amplitude data 
elseif frequencytype == 4
    [Y1,f1] = periodogram(x1,[],10000,125);
    [Y2,f2] = periodogram(x2,[],10000,125);
    A(1,:) = Y1; f(1,:) = f1;
    A(2,:) = Y2; f(2,:) = f2;
elseif frequencytype == 5
    f = zeros(2,length(f1(1:L/10)));
    A = zeros(2,length(f1(1:L/10)));    
    f(1,:) = f1(1:L/10);            % frequency data
    A(1,:) = 2*abs(Y1(1:L/10))';    % amplitude data
    f(2,:) = f2(1:L/10);            % frequency data
    A(2,:) = 2*abs(Y2(1:L/10));     % amplitude data    
    Am=(A(1,:).^2+A(2,:).^2).^.5;  
end;


%% Stage 3: Peak Selection
% Finding peak
numf = size(f,1);
peakf = zeros(numf,50);
peakA = zeros(numf,50);

% Finding peaks in the spectrum
for jj = 1:numf
    if peakdettype == 1
        [maxfpeak,minfpeak] = peakdet(A(jj,:),peaktol); % peak detection
        peakf(jj,1:length(f(maxfpeak(:,1)))) = f(jj,maxfpeak(:,1))';
        peakA(jj,1:length(f(maxfpeak(:,1)))) = maxfpeak(:,2)';
    elseif peakdettype == 2
        [maxfpeak] = peakdetectarbiter(A(jj,:),10^(-8)*10); % peak detection
        peakf(jj,1:length(f(jj,maxfpeak(1,:)))) = f(jj,maxfpeak(1,:))';
        peakA(jj,1:length(f(jj,maxfpeak(1,:)))) = A(jj,maxfpeak(1,:))';
    end
end;


% Peak filter with MA dominant frquencies
matol = 0.1;
maidx = [];
for jj = 1:length(MAf)
    for kk = 1:numf
        % remove peak
        maidx = find(abs(peakf(kk,:)-MAf(jj)) < matol);
        peakA(kk,maidx) = peakA(kk,maidx)*removerate;
%         peakA(kk,maidx) = peakA(kk,maidx)*removerate.*(abs(peakA(kk,maidx)-MAf(jj))+0.1)*3;
        % recover peaks that are closed
        if niter > 1
            maidxclose = find(abs(peakf(kk,maidx) - preguessBPM/60) < matol);
            for maidxclosei = 1:length(maidxclose)
                if abs(peakf(kk,maidx(maidxclose(maidxclosei))) - MAf(jj)) < matol
                    peakA(kk,maidx(maidxclose(maidxclosei))) = peakA(kk,maidx(maidxclose(maidxclosei)))*recoveryrate;
                end;
            end;
        end;
        % remove double frequency
        if isdoubleMA == 1
            maidx = find(abs(peakf(kk,:)-MAf(jj)*2) < matol);
            peakA(kk,maidx) = peakA(kk,maidx)*doubleremoverate;
            % recover peaks that are closed
            if niter > 1
                maidxclose = find(abs(peakf(kk,maidx) - preguessBPM/60) < matol);
                for maidxclosei = 1:length(maidxclose)
                    if abs(peakf(kk,maidx(maidxclose(maidxclosei))) - MAf(jj)*2) < matol
                        peakA(kk,maidx(maidxclose(maidxclosei))) = peakA(kk,maidx(maidxclose(maidxclosei)))*recoveryrate;
                    end;
                end;
            end;   
        end;
    end;
end;
% Remove dominant frequency of MA from PPG
% Currently unused

%% Peak Selection
if niter <= avgwin % Initial Guess
    if isdisplay == 1 && niter >= displayiter
        guessf = peakselectioninit2(peakf,peakA,f,A,realbpm,MAf,MAp,initcosttype,1,peakh);
    else
        guessf = peakselectioninit2(peakf,peakA,f,A,realbpm,MAf,MAp,initcosttype,0,peakh);
    end
%     guessBPM = guessf(1) * 60;
    guessBPM = (guessf(1)+guessf(2))*30;
else
%     maAgrad(niter) = (maA(niter-1)-maA(niter-mawin+1-1))/mawin;
    if isdisplay == 1 && niter >= displayiter
        guessf = peakselection2(niter,peakf,peakA,f,A,preguessBPM,preguessavgBPM,maA,realbpm,MAf,MAp,isneighboring,1,peakh);
    else
        guessf = peakselection2(niter,peakf,peakA,f,A,preguessBPM,preguessavgBPM,maA,realbpm,MAf,MAp,isneighboring,0);
    end;
    
    if frequencytype ~= 1
        if abs(guessf(1) - preguessBPM/60) < abs(guessf(2) - preguessBPM/60)
            guessBPM = guessf(1) * 60;
        else
            guessBPM = guessf(2) * 60;
        end;
    elseif frequencytype == 1
        guessBPM = guessf(1)*60;
    
    end;
end;


% Optimization on guessed BPM
if niter > avgwin
    BPMtol = 10;
    preBPM = BPM(niter-1);
    if guessBPM > preBPM + BPMtol
        guessBPM = preBPM + BPMtol;
        over = over + 1;
        under = 0;
        if over >= 2
            BPMtol = 20;
        else
            BPMtol = 10;
        end;
    elseif guessBPM < preBPM - BPMtol
        guessBPM = preBPM - BPMtol;
        over = 0;
        under = under + 1;
        if under >= 2
            BPMtol = 20;
        else
            BPMtol = 10;
        end;
    else
        over = 0;
        under = 0;
        BPMtol = 10;
    end;
    
end;


% guess BPM
if niter > 1
    alpha = 0.75;
    guessBPM = (1-alpha)*BPM(niter-1) + alpha*guessBPM;
end;



%% Print information
% BPM = [BPM;guessBPM];

if isdisplay == 1
    fprintf('%4d ', niter);
    fprintf('%f\n', guessBPM);
    figure(bpmh); hold on;
    subplot(2,1,1); hold on;
    plot(niter, guessBPM,'bo');
    plot(niter,realbpm,'ro'); hold off;
    if niter == 1
        legend('guess','real');
    end;
    subplot(2,1,2); hold on;
    plot(niter,maA(niter),'go'); hold off;
    if niter == 1
        legend('MA A');
    end;
end;



end
