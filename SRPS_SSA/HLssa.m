function [y,r,vr]=HLssa(x1,sig,t0,Np,L,per)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------                           
%    Author: Francisco Javier Alonso Sanchez    e-mail:fjas@unex.es
%    Departament of Electronics and Electromecanical Engineering
%    Industrial Engineering School
%    University of Extremadura
%    Badajoz
%    Spain
% -----------------------------------------------------------------
%
% SSA generates a trayectory matrix X from the original series x1
% by sliding a window of length L. The trayectory matrix is aproximated 
% using Singular Value Decomposition. The last step reconstructs
% the series from the aproximated trayectory matrix. The SSA applications
% include smoothing, filtering, and trend extraction.
% The algorithm used is described in detail in: Golyandina, N., Nekrutkin, 
% V., Zhigljavsky, A., 2001. Analisys of Time Series Structure - SSA and 
% Related Techniques. Chapman & Hall/CR.

% x1 Original time series (column vector form)
% L  Window length
% y  Reconstructed time series
% r  Residual time series r=x1-y
% vr Relative value of the norm of the approximated trajectory matrix with respect
%	  to the original trajectory matrix

% The program output is the Singular Spectrum of x1 (must be a column vector),
% using a window length L. You must choose the components be used to reconstruct 
%the series in the form [i1,i2:ik,...,iL], based on the Singular Spectrum appearance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin<5)
    L=400;
    per=0.8;
end


% Step1 : Build trayectory matrix
   ma=myParse(sig,t0,4,6);
   N=length(x1); 
   if L>N/2;L=N-L;end
	K=N-L+1; 
   X=zeros(L,K);  
	for i=1:K
	  X(1:L,i)=x1(i:L+i-1); 
	end
    
% Step 2: SVD

   S=X*X'; 
	[U,autoval]=eig(S);
	[d,i]=sort(-diag(autoval));  
   d=-d;
   U=U(:,i);sev=sum(d); 
	%plot((d./sev)*100),hold on,plot((d./sev)*100,'rx');
	%title('Singular Spectrum');xlabel('Eigenvalue Number');ylabel('Eigenvalue (% Norm of trajectory matrix retained)')
   V=(X')*U; 
   rc=U*V';

 % periodogram
   [Uu,Ss,Vd]=svd(ma);
   z=Uu'*ma;
   [pxx,f] = periodogram(z(1,:),[],4096,125);

   pxx(1)=0;
   maxpeak=max(pxx);
   [pxx_row pxx_col]=find(pxx>(maxpeak*0.5));

   Delta=10;
   [locationx locationy]=find(f>Np);
   location=locationx(1);   
   harm_freq=f( max(1,location-Delta) :location+Delta);
   noise_freq=f(pxx_row);
   for i=1:size(noise_freq,1)
       if ~isempty(find(harm_freq==noise_freq(i)))
        noise_freq(i)=0;
       end
   end
   [tempx tempy]=find(noise_freq==0);
   noise_freq(tempx)=[];
   
   I=[1:L];
   Vt=V';
%   fprintf('Eliminating group:');
   for i=1:30
       rca=U(:,i)*Vt(i,:);
       y=[rca(1,1:size(rca,1)) rca(end,:)];
       [pyy,fy] = periodogram(y,[],4096,125);
       [maxypeak,maxyidx]=max(pyy);
       [pyy_row pyy_col]=find(pyy>(maxypeak*per));%variable
       for j=1:size(pyy_row,1)
       if ~isempty(find(noise_freq==fy(pyy_row(j))))
%           fprintf(' %d',i);
           [tc tp]=find(I==i);
           I(tp)=0;
           break
       end
       end
   end
%   fprintf('\n');
   [tempx tempy]=find(I==0);
   I(tempy)=[];
% Step 3: Grouping
   %I=input('Choose the agrupation of components to reconstruct the series in the form I=[i1,i2:ik,...,iL]  ')
   Vt=V';
   rca=U(:,I)*Vt(I,:);
   
   
% Step 4: Reconstruction

   y=zeros(N,1);  
   Lp=min(L,K);
   Kp=max(L,K);

   for k=0:Lp-2
     for m=1:k+1;
      y(k+1)=y(k+1)+(1/(k+1))*rca(m,k-m+2);
     end
   end

   for k=Lp-1:Kp-1
     for m=1:Lp;
      y(k+1)=y(k+1)+(1/(Lp))*rca(m,k-m+2);
     end
   end

   for k=Kp:N
      for m=k-Kp+2:N-Kp+1;
       y(k+1)=y(k+1)+(1/(N-k))*rca(m,k-m+2);
     end
   end
   
   r=x1-y';
   vr=(sum(d(I))/sev)*100;
   y=y';
   
   %Y=abs(fft(y));
   %Y=Y/max(Y);
   
  %{
   f=linspace(0, 62.5, 500);
   figure(1);
   stem(f(1:100),2*abs(Y(1:100)));
    %}
%    figure;subplot(2,1,1);hold on;xlabel('Data poit');ylabel('Original and reconstructed series')
%    plot(x1);grid on;plot(y,'r')
% 
%    subplot(2,1,2);plot(r,'g');xlabel('Data poit');ylabel('Residual series');grid on
end
