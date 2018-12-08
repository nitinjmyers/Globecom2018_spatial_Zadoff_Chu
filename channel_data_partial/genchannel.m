clear all;
clc;
DirPDPInfo=dlmread('DirPDPInfo.txt');
M=DirPDPInfo;
T=10; %ns
W=1/T; %GHz
Na=16;
Ne=32;
Ua=dftmtx(Na)/sqrt(Na);
Ue=dftmtx(Ne)/sqrt(Ne);
L=13; % make even
d=[];
Nit=100;
tic
ls=zeros(Na,Ne*L,Nit);
Nam=0;
s=[];
del=[];
for i=1:1:Nit
  % H=zeros(Na,Ne,L);
   ind=find(M(:,1)==i);
   G=M(ind,:);
   P=M(ind,4);
   d=M(ind,3);
   P=10.^(P/10);
   dmean=sum(P.*d)/sum(P);
   diffsq=(d-dmean).^2;
   drms=sum(P.*diffsq)/sum(P);
   del=[del,sqrt(drms)];

   tm=([0:1:L-1]*T);
   tm=tm-mean(tm)+dmean-(T/2);  
   Z=zeros(Na*Ne,L);
   for k=1:1:length(ind)
        pulse=sinc(W*(tm-G(k,3)));
        pulse=pulse/norm(pulse);
        theta_a=pi*G(k,8)/180; % in rad
        theta_e=pi*G(k,9)/180;
        Tmp=exp(1i*pi*sin(theta_a)*[1:1:Na].')*exp(1i*pi*sin(theta_e)*[1:1:Ne]);
        Tmp=Tmp*sqrt(10^(G(k,4)/10))*exp(1i*G(k,5));
        Z=Z+(vec(Tmp)*pulse);
   end
   s=[s,norm(Z,'fro')/sqrt(Na*Ne)];
   H=Z*sqrt(Na*Ne)/norm(Z,'fro');
   dlmwrite(strcat('it_',num2str(i),'time_dom_channel.txt'),H);

end
%dlmwrite('channelscalings.txt',s)
toc