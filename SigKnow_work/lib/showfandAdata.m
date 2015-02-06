function [] = showfandAdata(fdata,Adata,i,bpm,BPM,period)
figure();
% stem(fdata(:,i),Adata(:,i)); hold on;
% plot(bpm(i)/60,0,'ro');

% h = plot(fdata(:,1),Adata(:,1), 'EraseMode', 'xor');
% h2 = plot(bpm(1)/60,0,'ro');

for j = 2:i
%     set(h,'xdata', fdata(:,j), 'ydata', Adata(:,j));
%     set(h2,'xdata', bpm(j)/60);
    stem(fdata(:,j),Adata(:,j)); hold on;
    plot(bpm(j)/60,0,'ro');
    plot(BPM(j)/60,0,'go');
    plot(BPM(j-1)/60,0,'yo');
    legend('energy','real','result','past');
    axis([0 5 0 8]); hold off;
    str = sprintf('error: %f  guess: %f  real: %f\n',abs(BPM(j)-bpm(j)),BPM(j),bpm(j));
    xlabel(str);
    str = sprintf('peroid: %d\n', j);
    title(str);
    pause(period);
    hold off;
end;

end