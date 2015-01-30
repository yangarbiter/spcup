function [err] = err1(BPM,realBPM)
n = length(BPM);
err = abs(BPM-realBPM);
err = sum(err)/n;
end