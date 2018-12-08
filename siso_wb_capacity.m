function Rate=siso_wb_capacity(h,Nsub,SNR)

hfreq=fft(h,Nsub);   % effectively Nsub*norm(h)

P=10^(SNR/10);

noise=1;

%Nsub transmissions
% frequency domain
eigval=abs(hfreq).^2;
noisefreq=noise./eigval;

Pfill=wfill(noisefreq,P*Nsub);

Rate=0;

for i=1:1:Nsub
    SNRsub= Pfill(i)/noisefreq(i);
    Rate=Rate+log2(1+SNRsub);
end

Rate=Rate/Nsub;
end
