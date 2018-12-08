%% Source code for
%% [1] N.J.Myers, A. Mezghani and R. W. Heath Jr.,"Spatial Zadoff-Chu modulation for rapid beam alignment in mmWave phased arrays",
%% in proc. of the IEEE Globecom 2018
%% If you use this code, please cite our Globecom paper [1] (or)
%% [2] N.J.Myers, A. Mezghani and R. W. Heath Jr., "Swift-Link: A compressive beam alignment algorithm for practical mmWave radios" IEEE Trans. on Signal Process. 2018
%% Note that the dictionary used for sparse representation of mmWave channel is standard 2D-DFT and not any oversampled version of 2D-DFT.
%% Using oversampled dictionaries introduces additional complexity to the problem and introduces the frame coherence problem.
%% The whole point of this exercise is to show that the ZC-based CS matrices achieve better beam alignment than random phase shifts under the same settings
clear all;
close all;
clc;
rng(24);
addpath('../channel_data_partial')         % adds path containing NYU channel data for wideband; parameters in paper

% System parameters
Nrx=16;
Ntx=32;

Nantsq=Nrx*Ntx;
Nant=sqrt(Nantsq);              % The normalization factor used for fft2 and ifft2 operations

qbrx=2;                         % resolution of phase shifters at RX
qbtx=3;                         % resolution of phase shifters at TX
% Algorithm paramters
Stand=1;                        %  0 for just the proposed method, 1 to also evaluate Standard IID random phase shifts
result=1;                       % 1- Performance with SNR, 2- Performance with measurements
Nrealz=20;                   % number of channels-100 available in folder conf_chan; 100 used for [1] using DirPDPInfo.txt
Nruns=10;                    % random training realizations; 50 used for [1]
Ntaps=13;                    % taps in wideband channel
Codelen=13;                  % length of aperiodic autocorrelation sequence, barker seq. is used here
Nomp=100;                   % Max. number of iterations for OMP ; Max itern. and stopping criterion is maintained same for both CS matrices

Baserx=quantz(exp(1i*(9*pi/16)*([0:1:15].^2)),qbrx).'/sqrt(Nrx);        % quantized core ZC sequence at RX
Basetx=quantz(exp(1i*(11*pi/32)*([0:1:31].^2)),qbtx)/sqrt(Ntx);          % quantized core ZC sequence at RX
% It is important that abs(fft(Baserx)) must be uniform- choose root
% accordingly; Our upcoming work will address this issue in a different way
Specmkrx=fft(Baserx);                                             % Spectral mask-Refer Swift-Link
Specmktx=fft(Basetx);

if(result==1)          %% For rate vs SNR
    Measvec=[20];      % 20 channel measurements corresponds to 5 micro seconds per parameters in [1]
    SNRvec= [-15:2.5:5];
else
    Measvec=[12:4:48];  %%For Rate vs channel measurements
    SNRvec=[0];                          % Note that this is ~ "Omnidirectional" SNR at chip rate
end
Ratenew=zeros(length(Measvec),length(SNRvec));      % Placeholder for rate using proposed ZC-based training in OMP
Rateold=Ratenew;                                    % Placeholder for rate using standard IID random phase shifts in OMP
Ratemax=zeros(length(SNRvec),1);                    % Placeholder for perfect CSI case
Cmat=match_filterop(0);                             % Barker code correlation at zero CFO - used as spread sequence/ LAAC sequence
Timenew=zeros(length(Measvec),length(SNRvec));      % Placeholder for computation time-ZC training
Timeold=Timenew;                                    % Placeholder for computation time- random training
for real=1:1:Nrealz         % Over channel realizations
    % read channel
    Hv=dlmread(strcat('it_',num2str(real),'time_dom_channel.txt'));
    % pick max energy tap
    [~,tp]=max(norms(Hv));
    H=reshape(Hv(:,tp),[Nrx,Ntx]);          % Perform beamforming according to max energy tap
    [Uh,~,Vh]=svd(H);
    Uh=quantz(Uh(:,1),qbrx);
    Vh=quantz(Vh(:,1),qbtx);
    Uh=Uh/norm(Uh);Vh=Vh/norm(Vh);
    maxch=kron(Vh.',Uh')*Hv;                %  wideband SISO channel seen when Uh, Vh are applied as beamformers
    
    for Measind=1:1:length(Measvec)
        Meas=Measvec(Measind);                      % Number of spatial measurements to be used
        for run=1:1:Nruns                                    % randomize training and noise
            [Rw,Cl,indx]=gentrajectory(Nrx,Ntx,Meas);       % random trajectory/ random partial 2D-DFT sampling
            Sampler=zeros(Nrx,Ntx);                            % Sampler is 1 if cooresponding partial 2D-DFT entry is sampled;
            Ypure=zeros(Meas,Ntaps);
            for tps=1:1:Ntaps
                Allsamp=fft2(ifft2(reshape(Hv(:,tps),[Nrx,Ntx])).*(Specmkrx*Specmktx));
                Ypure(:,tps)=Allsamp(indx);                           % make this as a matrix-Meas x Ntaps  - this is noiseless Y in [1]
            end
            Sampler(indx)=1;
            Noise=(randn(Meas,Ntaps)+1i*randn(Meas,Ntaps))/sqrt(2);  % make this as a matrix-Meas x Ntaps
            if(Stand==1)
                Ast=zeros(Meas,Nantsq);                              % Ast is CS matrix obtained with standard IID random phase shifts
                Yst=zeros(Meas,Ntaps);
                PSrx=exp(1i*(randi(2^qbrx,[Nrx,Meas])*2*pi/(2^qbrx)))/sqrt(Nrx);     %  Meas # of random IID phase shift vectors for spatial channel measurements
                PStx=exp(1i*(randi(2^qbtx,[Ntx,Meas])*2*pi/(2^qbtx)))/sqrt(Ntx);
                for mm=1:1:Meas
                    PSmatmm=PSrx(:,mm)*PStx(:,mm).';
                    Yst(mm,:)=(conj(PSmatmm(:)).')*Hv;
                    Ast(mm,:)=vec(fft2(conj(PSmatmm))).'/Nant;                      % scaling because fft2 multiplies by Nant
                end
            end
            
            for SNRind=1:1:length(SNRvec)
                SNR=SNRvec(SNRind);
                sig=10^(-SNR/20);
                Yall=(Ypure*Cmat)+(sqrt(Codelen)*sig*Noise);              % Meas x Ntaps matrix - include noise enhancement due to spread seq.
                % Cmat incorporates spreading gain of barker_code here
                [~,bindx]=max(norms(Yall));                               % pick the best index corresponding to the "best" tap
                Y=Yall(:,bindx);
                Nsz=size(Sampler);
                ts1=tic;
                xest=FastFourOMP(Sampler, Y, sig,Nomp);                   % run fast partial 2D-DFT OMP  over the masked beamspace
                Timenew(Measind,SNRind)=Timenew(Measind,SNRind)+toc(ts1);
                Hest=fft2(reshape(xest,[Nrx,Ntx]).*conj((Specmkrx*Specmktx))); % invert the mask and find Hest
                [NUh,~,NVh]=svd(Hest);
                NUh=quantz(NUh(:,1),qbrx);
                NVh=quantz(NVh(:,1),qbtx);
                NUh=NUh/norm(NUh);NVh=NVh/norm(NVh);
                chestnew=kron(NVh.',NUh')*Hv;                            % wideband SISO channel seen after beamforming with the proposed training + OMP
                Ratenew(Measind,SNRind)=Ratenew(Measind,SNRind)+siso_wb_capacity(chestnew,128,SNR); % Use 128 subcarriers and waterfilling
                
                if(Stand==1)
                    Yallstd=(Yst*Cmat)+(sqrt(Codelen)*sig*Noise);
                    % incorporate spreading gain of barker_code here
                    [~,bindxst]=max(norms(Yallstd));
                    Ysns=Yallstd(:,bindxst);
                    ts1=tic;
                    xstest=OMPf(Ysns,Ast,sig,Nomp);             % standard OMP with measurements, CS matrix, noise, Niter
                    Timeold(Measind,SNRind)=Timeold(Measind,SNRind)+toc(ts1);
                    Hstest=fft2(reshape(xstest,[Nrx,Ntx]));
                    [OUh,~,OVh]=svd(Hstest);
                    OUh=quantz(OUh(:,1),qbrx);
                    OVh=quantz(OVh(:,1),qbtx);
                    OUh=OUh/norm(OUh);OVh=OVh/norm(OVh);
                    chestold=kron(OVh.',OUh')*Hv;                       % wideband SISO channel seen after beamforming with standard raandom training + OMP
                    Rateold(Measind,SNRind)=Rateold(Measind,SNRind)+siso_wb_capacity(chestold,128,SNR);
                end
                
                if(run==1 && Measind==1)
                    Ratemax(SNRind)=Ratemax(SNRind)+siso_wb_capacity(maxch,128,SNR);        % Max rate needs to be computed for a given channel realization , no need to randomize over training/ Meas
                end
            end
        end
        
    end
end

Ratenew=100*Ratenew/(Nrealz*Nruns);         % Scale by 100 as Bandwidth is 100 MHz, get rate in Mbps
Ratemax=100*Ratemax.'/Nrealz;
Rateold=100*Rateold/(Nrealz*Nruns);
Timenew=Timenew/(Nrealz*Nruns);
Timeold=Timeold/(Nrealz*Nruns);

% Plot results here

% if looping over SNRs
if(length(SNRvec)>1 && length(Measvec)==1)
    figure(1)
    plot(SNRvec,Ratemax,'r-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on;
    plot(SNRvec,Ratenew,'k^-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r')
    plot(SNRvec,Rateold,'bs-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
    set(gca,'fontsize',16);
    legend1=legend('Perfect CSI','OMP with shifted ZC training','OMP with IID random phase shifts')
    set(legend1,'Interpreter','Latex','FontSize',14)
    xlabel('SNR (dB)','Interpreter', 'latex')
    ylabel('Achievable rate (Mbps)','Interpreter', 'latex')
end

%looping over Meas for given SNR
if(length(SNRvec)==1 && length(Measvec)>1)
    figure(1)
    plot(Measvec,Ratemax*ones(size(Rateold)),'r-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on;
    plot(Measvec,Ratenew,'k^-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r')
    plot(Measvec,Rateold,'bs-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
    set(gca,'fontsize',16);
    legend1=legend('Perfect CSI','OMP with shifted ZC training','OMP with IID random phase shifts')%,'DFT-beam with shifted ZC','DFT-beam with IID Random')
    set(legend1,'Position',[0.429572677612305 0.121666668710255 0.463284465244837 0.131904759861174],'Interpreter','Latex','FontSize',14)
    xlabel('Number of channel measurements $N_{\mathrm{p}}$','Interpreter', 'Latex')
    ylabel('Achievable rate (Mbps)','Interpreter', 'latex')
    
end
