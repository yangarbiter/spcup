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
[MApxx,MAf] = periodogram(MA,[],[],125); % periodogram
maxp = max(MApxx(4:end));
f = MAf(MApxx > maxp/2); % frequency that has the spectrum power exeed the half of the maximun spectrum power
p = MApxx(MApxx > maxp/2);
if length(f) > 6
    f = [];
end;
end