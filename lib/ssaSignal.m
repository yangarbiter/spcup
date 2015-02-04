% Program: ssaSignal.m
% Author: Jerry Chou
% Last Modified Date: 2014/12/30
%
% Introduction:
% Most of the code is copied from the internet. The output is the
% reconstruction signal after SSA with the component specified in input I.

function [y] = ssaSignal(x,L,I)

%% Step1 : Build trayectory matrix

N=length(x); 
if L>N/2
    L=N-L;
end;
K=N-L+1; 
X=zeros(L,K);  
for i=1:K
    X(1:L,i)=x(i:L+i-1); 
end
    
%% Step 2: SVD

S=X*X'; 
[U,autoval]=eig(S);
[d,i]=sort(-diag(autoval));  
d=-d;
U=U(:,i);sev=sum(d);
eigenP = (d./sev)*100; % eigenvalue percentage
% plot(eigenP),hold on,plot(eigenP,'rx');
% title('Singular Spectrum');xlabel('Eigenvalue Number');ylabel('Eigenvalue (% Norm of trajectory matrix retained)')
V=(X')*U; 
rc=U*V';
Vt=V';


%% Step 3: Grouping
% I = [1;2];
rca=U(:,I)*Vt(I,:);

%% Step 4: Reconstruction

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

end