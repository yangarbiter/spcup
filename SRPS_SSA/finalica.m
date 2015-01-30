function result = finalica(sig, t0, prev)
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
    % normalize the signal components
    xyz = ica_xyz(1,:).^2 + ica_xyz(2,:).^2 + ica_xyz(3,:).^2 ;
    xyz = xyz/max(xyz) ;
    ppg1 = myNormalize(ica_ppg(1, :)) ;
    ppg2 = myNormalize(ica_ppg(2, :)) ;
    ppg3 = myNormalize(raw_ppg(1, :) + raw_ppg(1, :)) ;

    % computing auto correlation
    acor_xyz = xcorr(xyz, xyz) ;
    acor_ppg(1,:) = xcorr(ppg1, ppg1) ;
    acor_ppg(2,:) = xcorr(ppg2, ppg2) ;
    acor_ppg(3,:) = xcorr(ppg3, ppg3) ;

    % computing the coefficients by cross correlation functions
    c(1) = max(xcorr(xyz, ppg1)) ;
    c(2) = max(xcorr(xyz, ppg2)) ;
    c(3) = max(xcorr(xyz, ppg3)) ;

    %STAGE 2 ==================================================================
    psd_xyz = abs(fft(acor_xyz)) ;
    psd_ppg1 = abs(fft(acor_ppg(1, :))) ;
    psd_ppg2 = abs(fft(acor_ppg(2, :))) ;
    psd_ppg3 = abs(fft(acor_ppg(3, :))) ;
    result_xyz = psd_xyz(8:64) ;
    result_ppg(1, :) = psd_ppg1(8:64) ;
    result_ppg(2, :) = psd_ppg2(8:64) ;
    result_ppg(3, :) = psd_ppg3(8:64) ;

    result(1, :) = result_xyz ;
    result(2, :) = result_ppg(1, :) ;
    result(3, :) = result_ppg(2, :) ;
    result(4, :) = result_ppg(3, :) ;
end

% % %STAGE 3 ==================================================================
% a = zeros(1,3) ;
% a(1) = max(result_ppg(1, :)) ;
% a(2) = max(result_ppg(2, :)) ;
% a(3) = max(result_ppg(3, :)) ;
% 
% [m1, index1] = max(psd_xyz(8:34)) ;
% if index1 == 1
%     result_psd = result_ppg(3, :) ;
% else
%     [m2, index2] = min(a) ;
%     [m3, index3] = min(c) ;
%     if index2 ~= index3
%         result_psd = result_ppg(index3, :) ;
%     else
%        [m4, index4] = max(a) ;
%        result_psd = result_ppg(index4, :) ;
%     end
% end
% result = result_psd/max(result_psd) ;

% % OUTPUT RESULT ===========================================================
% 
% % plot the results
x2 = linspace(0, 62.5, 1000);
x_axis2 = x2(1, 8:64);
% 
% figure(1) ;
% subplot(3,1,1) ;
% plot(ppg1) ;
% subplot(3,1,2) ;
% plot(ppg2) ;
% subplot(3,1,3) ;
% plot(ppg3) ;
% 
% figure(2) ;
% subplot(4,1,1) ;
% plot(x_axis2, result_xyz) ;
% subplot(4,1,2) ;
% plot(x_axis2, result_ppg(1, :)) ;
% subplot(4,1,3) ;
% plot(x_axis2, result_ppg(2, :)) ;
% subplot(4,1,4) ;
% plot(x_axis2, result_ppg(3, :)) ;
% 
% figure(3) ;
% plot(x_axis2, result) ;

end