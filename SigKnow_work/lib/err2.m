function [err] = err2(BPM,realBPM)
n = length(BPM);
err = abs(BPM-realBPM);
err = err ./ realBPM;
err = sum(err)/n;
end