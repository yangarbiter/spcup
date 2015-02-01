% Program: SigKnow.m
% Team ID: 20219
% Last Modified Time: 2015/02/01

function [guessBPM, maA, over, under] = SigKnow(varargin)

% =============== Parsing Input ===============
p = inputParser;
addOptional(p, 'curSig', @(x) ~isempty(x));
addOptional(p, 'BPM', @(x) ~isempty(x));
addOptional(p, 'maA', @(x) ~isempty(x));
addOptional(p, 'niter', @isnumeric);
addOptional(p, 'avgwin', @isnumeric);
addOptional(p, 'mawin', @isnumeric);
addOptional(p, 'ffttype', @isnumeric);
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
% addOptional(p, '', @isnumeric);


parse(p, varargin{:});

curSig              = p.Results.curSig;
BPM                 = p.Results.BPM;
maA                 = p.Results.maA;
niter               = p.Results.niter;
avgwin              = p.Results.avgwin;
mawin               = p.Results.mawin;
ffttype             = p.Results.ffttype;
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

if niter > 1
    preguessBPM = BPM(niter-1);
    if niter > avgwin
        preguessavgBPM = sum(BPM(niter-avgwin:niter-1))/avgwin;
    end;
end;
    
% Signals in this window
x1 = curSig(1,:); x2 = curSig(2,:);
m1 = curSig(3,:); m2 = curSig(4,:); m3 = curSig(5,:);

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
% MAp = [];
% for fi = 1:length(MAf)
%     count = 0;
%     if MAf(fi)
%     
% end;
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

if ffttype == 2
    x = x1;
    for fftt = 1:1000
        x1(fftt) = (fftsigma)^(1/4)*x1(fftt)*exp(-fftsigma*pi*((fftt - 500)/1000)^2);  
    end;    
end;
% FFT
T = 1/125;                          % sampling rate
Fs = 1/T;
L = 10000;
tt = (0:L-1)*T;
Y1 = fft(x1,L)/L; f1 = Fs*linspace(0,1,L);
Y2 = fft(x2,L)/L; f2 = Fs*linspace(0,1,L);

if ffttype == 1
    f = zeros(1,length(f1(1:L/10)));
    A = zeros(1,length(f1(1:L/10)));    
    f(1,:) = f1(1:L/10);            % frequency data
    A(1,:) = 2*abs(Y1(1:L/10))';    % amplitude data
elseif ffttype == 3
    f = zeros(2,length(f1(1:L/10)));
    A = zeros(2,length(f1(1:L/10)));    
    f(1,:) = f1(1:L/10);            % frequency data
    A(1,:) = 2*abs(Y1(1:L/10))';    % amplitude data
    f(2,:) = f2(1:L/10);            % frequency data
    A(2,:) = 2*abs(Y2(1:L/10));     % amplitude data 
end;


%% Stage 3: Peak Selection
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
matol = 0.1;
maidx = [];
for jj = 1:length(MAf)
    for kk = 1:numf
        % remove peak
        maidx = find(abs(peakf(kk,:)-MAf(jj)) < matol);
        peakA(kk,maidx) = peakA(kk,maidx)*removerate;
        % recover peaks that are closed
        if niter > 1
            maidxclose = find(abs(peakf(kk,maidx) - preguessBPM/60) < matol*1);
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
                maidxclose = find(abs(peakf(kk,maidx) - preguessBPM/60) < matol*1);
                for maidxclosei = 1:length(maidxclose)
                    if abs(peakf(kk,maidx(maidxclose(maidxclosei))) - MAf(jj)*2) < matol
                        peakA(kk,maidx(maidxclose(maidxclosei))) = peakA(kk,maidx(maidxclose(maidxclosei)))*2;
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
    guessf = peakselectioninit2(peakf,peakA,f,A,0);
    guessBPM = guessf(1) * 60;
else
%     maAgrad(niter) = (maA(niter-1)-maA(niter-mawin+1-1))/mawin;
    if isdisplay == 1 && niter > 82
        guessf = peakselection2(niter,peakf,peakA,f,A,preguessBPM,preguessavgBPM,maA,realbpm,MAf,0,peakh);
    else
        guessf = peakselection2(niter,peakf,peakA,f,A,preguessBPM,preguessavgBPM,maA,realbpm,MAf,0);
    end;
    
    if abs(guessf(1) - preguessBPM/60) < abs(guessf(2) - preguessBPM/60)
        guessBPM = guessf(1) * 60;
    else
        guessBPM = guessf(2) * 60;
    end;
end;


% Optimization on guessed BPM
if niter > avgwin
    BPMtol = 8;
    preBPM = BPM(niter-1);
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
    plot(niter, guessBPM,'bo');
    plot(niter,realbpm,'ro');
if niter == 1
    legend('guess','real');
end;
hold off;
end;


end
