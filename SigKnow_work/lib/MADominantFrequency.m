% Program: MADominantFrequency.m
% Author: Jerry Chou
% Last Modified Date: 2014/12/30
%
% Introduction:
% This program calculate the dominant frquency of MA with periodogram
% provided by built-in function of Matlab. The program return all the
% peak frequencies whose power are greater than half of the dominant
% frquency, and also return the power of them.
function [f, maxp,p] = MADominantFrequency(MA)
% type: 1: periodogram; 2: FFT
type = 1;
if type == 1
    [f,maxp,p] = MADominantPeriodogram(MA);
elseif type == 2
    [f,maxp,p] = MADominantFFT(MA);
end;

end

function [f,maxp,p] = MADominantPeriodogram(MA)
MAmax = 1;

[MApxx,MAf] = periodogram(MA,[],[],125); % periodogram
maxp = max(MApxx(4:end));
maxp = max(maxp,MAmax);
f = MAf(MApxx > maxp/2); % frequency that has the spectrum power exeed the half of the maximun spectrum power
p = MApxx(MApxx > maxp/2);

end

function [f,maxp,p] = MADominantFFT(MA)
T = 1/125;
Fs = 1/T;
L = 1000;
MAA = abs(fft(MA,L)/L)';
MAf = Fs*linspace(0,1,L)';
maxp = max(MAA(4:end));
f = MAf(MAA > maxp/2); % frequency that has the spectrum power exeed the half of the maximun spectrum power
p = MAA(MAA > maxp/2);
f = f(f<4);
p = p(f<4);

end