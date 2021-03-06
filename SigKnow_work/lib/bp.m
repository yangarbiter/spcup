function Hd = bp
%BP Returns a discrete-time filter object.

%
% MATLAB Code
% Generated by MATLAB(R) 8.0 and the Signal Processing Toolbox 6.18.
%
% Generated on: 21-Nov-2014 17:02:08
%

% FIR Window Bandpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 125;  % Sampling Frequency

Fstop1 = 0.5;               % First Stopband Frequency
Fpass1 = 1;             % First Passband Frequency
Fpass2 = 5;               % Second Passband Frequency
Fstop2 = 5.5;             % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
flag   = 'scale';         % Sampling Flag

% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                             1 0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dfilt.dffir(b);

% [EOF]
