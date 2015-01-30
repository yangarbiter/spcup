function y = myParse(sig, t0, n, m)

% parse the 8-seconds signal from channel n-m of "sig" starting at time "t0" 

first = t0*125+1;
last = t0*125+1000;
[s1,s2]=size(sig);
y(1:m-n+1,:) = sig(n:m,first:last) ;

end


