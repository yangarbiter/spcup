% This program is for 2015 IEEE Signal Processing Cup
% Program: SigKnow.m
% Team ID: 20219
% Last Modified Time: 2015/02/01
% 
% 

clc; clear; close all;
tic;

%%  ==================== Program Initial Settings =======================

addpath(genpath('TrainData'));
addpath(genpath('TestData'));
addpath(genpath('lib'));

% Train Dataset IDs
TRID = { 'DATA_01_TYPE01', 'DATA_02_TYPE02', 'DATA_03_TYPE02', 'DATA_04_TYPE02', ...
         'DATA_05_TYPE02', 'DATA_06_TYPE02', 'DATA_07_TYPE02', 'DATA_08_TYPE02', ...
         'DATA_09_TYPE02', 'DATA_10_TYPE02', 'DATA_11_TYPE02', 'DATA_12_TYPE02', ...
         'DATA_13_TYPE02' };
% Test Dataset IDs
ID = { 'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
       'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
       'TEST_S07_T02', 'TEST_S08_T01'}; 
REID = { 'RESULT_S01_T01', 'RESULT_S02_T01', 'RESULT_S02_T02', 'RESULT_S03_T02', ...
       'RESULT_S04_T02', 'RESULT_S05_T02', 'RESULT_S06_T01', 'RESULT_S06_T02',...
       'RESULT_S07_T02', 'RESULT_S08_T01'};
   
% Parameters setting
srate = 125;                        % 125 Hz
window   = 8 * srate;               % window length is 8 seconds
lwa = 8;                            % left extension
lw = (8+lwa) * srate;
step     = 2 * srate;               % step size is 2 seconds
peaktol = 0.05;                     % peak detection tolerance for peakdet()
ftol = 0.1;                         % remove frquencies in (MAf-ftol, MAf+ftol)
removerate = 0.3;   
doubleremoverate = 0.3;
isdoubleMA = 0;
avgwin = 4;                         % Window size of average previous BPM
mawin = 3;
MINBPM = 65; MAXBPM = 180;
prequessBPM = 80; quessBPMavg = 80;
recoveryrate = 2.5;
isBP = 0;
isHHT = 0;
isICA = 0;
isSSA = 1;
isneighboring = 0;

% cost
initcosttype = 2;

% peakdettype
%   1: peakdet
%   2: critical point detect
peakdettype = 2;

% For testing
guessBPMA = zeros(10,2,200);        % guessBPM amplitude
error1 = zeros(13,1); error2 = zeros(13,1);
elapsedtime = zeros(13,1);

% Data id range
sdataidx = 10; edataidx = 10;
% datatype:
%   1: data with answer
%   2: test data
datatype = 1;
% frequencytype:
%   1: rectangular
%   2: gaussian window
%   3: two-channel
%   4: periodogram
%   5: gaussian
frequencytype = 3;
fftsigma = 0.01;
fftwindow = 4;

% display result figure or not
isdisplay = 1;
displayiter = 14;
isrecord = 0;

%% ===================== Start ==========================
for dataidx = sdataidx:edataidx
tdata = tic;
% ============ load data ==============
if datatype == 1
    load(TRID{dataidx});
    tracename = sprintf('%s_BPMtrace', TRID{dataidx});
    load(tracename);
else
    load(ID{dataidx});
end;
if datatype == 1
    ECG = sig(1,:); PPG1 = sig(2,:); PPG2 = sig(3,:);
    MA1 = sig(4,:); MA2 = sig(5,:); MA3 = sig(6,:);
    bpm = BPM0';
elseif datatype == 2
    PPG1 = sig(1,:); PPG2 = sig(2,:);
    MA1 = sig(3,:); MA2 = sig(4,:); MA3 = sig(5,:);
    bpm = zeros(1,200);
end;

% ======================= BP ===========================
if isBP == 1
    Lf=size(sig,2);
    ff1=round(0.2/srate*Lf);  % 0.2: lower cutoff frequency
    ff2=round(4/srate*Lf);    % 4:   higher cutoff frequency  
    for i=1:size(sig,1)
        sf1=fft(sig(i,:));
        sf1([1:ff1,ff2+2:Lf-ff2,Lf+2-ff1:Lf])=0;
        sig(i,:)=real(ifft(sf1));
        % sig(i,:)=filter(bp,sig(i,:));
    end
%     subplot(211);plot(sig1(2,1:20000))
%     subplot(212);plot(sig(2,1:20000))
end

% ========== Initial Setting for Dataset =============
windowNb = floor(length(PPG1)/250 - 3);  % total number of windows(estimates)
BPM = zeros(1,windowNb);
ii = 1;
under = 0;
over = 0;
t0 = 0;                             % The start time, default value is 0
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
maA1 = zeros(1,windowNb);
maA2 = zeros(1,windowNb);
maA3 = zeros(1,windowNb);
maA = zeros(1,windowNb);
maAgrad = zeros(1,windowNb);


% ============= Algorithm Begins ============
fprintf('\n*** Algorithm Begins on Data %02d ***\n', dataidx);

% Initialize Debug Information Figure
if isdisplay == 1
    bpmh = figure();
    BPMtitle = sprintf('Data: %d BPM', dataidx);
    title(BPMtitle);
    peakh = figure();
    title('Peak');
end;
icaPrev = zeros(2,1000);

if frequencytype == 5
    if datatype == 1
        sig1=[zeros(5,srate*lwa),sig(2:6,:)];
    else
        sig1=[zeros(5,srate*lwa),sig];
    end
end
% =================== Iteration Start =====================
    for i = 1 : windowNb
        
%         if i < 5 
%             frequencytype = 5;
%         else
%             frequencytype = 3;
%         end;
        
        if frequencytype == 5
            curSegment = (i-1)*step+1 : (i-1)*step+lw;
            inputSig = sig1(:,curSegment);
        else
            curSegment = (i-1)*step+1 : (i-1)*step+window;
            if datatype == 1
                inputSig = sig(2:6,curSegment);
            else
                inputSig = sig(:,curSegment);
            end;
        end
        
        if isdisplay == 1
            [BPM(i), maA, over, under, f, A, MAf] = ... 
                SigKnow('curSig', inputSig, 'BPM', BPM, 'maA', maA, 'niter', i, ...
                    'avgwin', avgwin, 'mawin', mawin, 'frequencytype', frequencytype, ...
                    'peaktol', peaktol, 'removerate', removerate, ...
                    'doubleremoverate', doubleremoverate, ...
                    'isdoubleMA', isdoubleMA, 'isdisplay', isdisplay, ...
                    'realbpm',bpm(i), 'over', over, 'under', under, ...
                    'bpmh', bpmh, 'peakh', peakh,'isHHT', isHHT,...
                    'isICA', isICA, 'isSSA', isSSA, 'icaPrev', icaPrev, ...
                    'recoveryrate', recoveryrate, 'peakdettype', peakdettype, ...
                    'isneighboring', isneighboring, 'displayiter', displayiter, ...
                    'initcosttype', initcosttype); 
        else
            [BPM(i), maA, over, under, f, A, MAf] = ... 
                SigKnow('curSig', inputSig, 'BPM', BPM, 'maA', maA, 'niter', i, ...
                    'avgwin', avgwin, 'mawin', mawin, 'frequencytype', frequencytype, ...
                    'peaktol', peaktol, 'removerate', removerate, ...
                    'doubleremoverate', doubleremoverate, ...
                    'isdoubleMA', isdoubleMA, 'isdisplay', isdisplay, ...
                    'realbpm',bpm(i), 'over', over, 'under', under, ...
                    'isHHT', isHHT, 'isICA', isICA, 'isSSA', isSSA, ...
                    'icaPrev', icaPrev, 'recoveryrate', recoveryrate, ...
                    'peakdettype', peakdettype, 'isneighboring', isneighboring, ...
                    'displayiter', displayiter, 'initcosttype', initcosttype);   
        end;
                
                
        if isdisplay == 0
            fprintf('.');
        end;
        

        if isrecord == 1
            figure(peakh);
            titlename = sprintf('Iteration: %d', i);
            title(titlename);
            Mpeak(i) = getframe(peakh);
        end
        
        if datatype == 2
            filename = sprintf('BPM_SSA_GAU_%s',ID{dataidx});
            save(filename,'BPM');
        end
        
    end;
    
    % Error Calculation and analysis
    if datatype == 1
        fprintf('\nerr1 = %f\n',err1(BPM,bpm(1:length(BPM))));
        fprintf('err2 = %f\n',err2(BPM,bpm(1:length(BPM))));
        error1(dataidx) = err1(BPM,bpm(1:length(BPM)));
        error2(dataidx) = err2(BPM,bpm(1:length(BPM)));
    end;
    
    elapsedtime(dataidx) = toc(tdata);

    if isrecord == 1
        filename = sprintf('TestVideo/TestData%02d.avi', dataidx);
        movie2avi(Mpeak, filename,'Compression', 'None', 'fps', 1);
        clear('Mpeak');
    end;
    
end;

totalelapsedtime = toc;

%% ==================== Results ==========================

% Elapsed Time Analysis

fprintf('\n============== Elapsed Time =============\n');
for i = sdataidx:edataidx
    fprintf('Data %02d Elapsed Time: %f\n', i, elapsedtime(i));
end;
fprintf('Total Elapsed Time: %f\n', totalelapsedtime);


% Error Calculation and analysis
if datatype == 1
    fprintf('\n================ Error ==================\n');
    for i = sdataidx:edataidx
        fprintf('Data %02d Error1: %f\n', i, error1(i));
    end;
    fprintf('Avg Error1: %f\n', sum(error1(sdataidx:edataidx)/(edataidx-sdataidx+1)));
    fprintf('Avg Error2: %f\n', sum(error2(sdataidx:edataidx)/(edataidx-sdataidx+1)));
end;




%% =================== Close ======================


