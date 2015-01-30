function result = icaForSSA(sig, t0, prev)
% INPUT SIGNAL ============================================================

% parse the signal
raw_xyz = myParse(sig, t0, 4, 6);
raw_ppg = myParse(sig, t0, 2, 3);

% STAGE 1 =================================================================

% bandpass filter
Hd = PPGBandpass ;
raw_xyz(1, :) = filter(Hd, raw_xyz(1, :)) ;
raw_xyz(2, :) = filter(Hd, raw_xyz(2, :)) ;
raw_xyz(3, :) = filter(Hd, raw_xyz(3, :)) ;
raw_ppg(1, :) = filter(Hd, raw_ppg(1, :)) ;
raw_ppg(2, :) = filter(Hd, raw_ppg(2, :)) ;

% do ica
[ica_xyz] = fastica (raw_xyz, 'numOfIC', 3);
[ica_ppg] = fastica (raw_ppg, 'numOfIC', 2);

% prevention for diverge ICA result
xyz_size = size(ica_xyz);
ppg_size = size(ica_ppg);
if ( (xyz_size(1) ~= 3) || (ppg_size(1) ~= 2) )
    result = prev ;
else
    result(1, :) = ica_ppg(1, :) ;
    result(2, :) = ica_ppg(2, :) ;
end

end