function y = myParse(sig, n, m)

% parse the 8-seconds signal from channel n-m of "sig" starting at time "t0" 
y(1:m-n+1,:) = sig(n:m,:) ;

end


