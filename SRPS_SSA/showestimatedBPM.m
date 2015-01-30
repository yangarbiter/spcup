% Program: showestimatedBPM.m
% Author: Jerry Chou
% Last Modified: 2015/1/5
%
% This program show the comparison between real BPM and estimated BPM

function [] = showestimatedBPM(bpm,BPM,dataidx)
figure();
plot(1:length(bpm),bpm)
grid on
hold on
plot(1:length(bpm),BPM,'r')
error1 = err1(BPM,bpm(1:length(BPM)));
error2 = err2(BPM,bpm(1:length(BPM)));
titlestr = sprintf('Dataset: %d Comparison between real BPM and estimated BPM, err1: %f err2: %f', dataidx,error1,error2);
title(titlestr);
legend('Real', 'Estimates');

end

