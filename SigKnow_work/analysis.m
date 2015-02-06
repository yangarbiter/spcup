clc; clear;

addpath('TestData');
load 'DATA_10_TYPE02'

m1 = sig(4,1:1000);
m2 = sig(5,1:1000);
m3 = sig(6,1:1000);


[MAf1, MApmax1, MAp1] = MADominantFrequency(m1);
[MAf2, MApmax2, MAp2] = MADominantFrequency(m2);
[MAf3, MApmax3, MAp3] = MADominantFrequency(m3);


MAf = unique([MAf1;MAf2;MAf3]);
MAp = [];
for fi = 1:length(MAf)
    count = 0;
    for mi = 1:length(MAf1)
        if MAf(fi) == MAf1(mi)
            count = count + MAp1(mi)^2;
        end
    end
    for mi = 1:length(MAf2)
        if MAf(fi) == MAf2(mi)
            count = count + MAp2(mi)^2;
        end
    end
    for mi = 1:length(MAf3)
        if MAf(fi) == MAf3(mi)
            count = count + MAp3(mi)^2;
        end
    end
    MAp = [MAp;count];
end



