function [idx] = peakdetectarbiter(A, thrs)
en = 501;
dA = A(2:en) - A(1:en-1);
dAd = dA(2:en-1) .* dA(1:en-2);
fi = find(dAd < thrs);
idx = [];
lasti = 0;
for i = 1: length(fi)-1
    if fi(i) ~= fi(i+1)-1
        idx = [idx, 1+floor(mean(fi(lasti+1: i)))];
        lasti = i;
    end
end
end