%% If you use this code, please cite our Globecom paper [1] N.J.Myers, A. Mezghani and R. W. Heath Jr.,"Spatial Zadoff-Chu modulation for rapid beam alignment in mmWave phased arrays",
%% in proc. of the IEEE Globecom 2018 (or)
%% [2] N.J.Myers, A. Mezghani and R. W. Heath Jr., "Swift-Link: A compressive beam alignment algorithm for practical mmWave radios" IEEE Trans. on Signal Process. 2018

% fast algorithm for partial DFT CS
% y= P(F(X))+N; measurements are partial 2D-DFT of sparse X 
% OMP- inputs [mask, y, sig,Niter]
% res=y-A*x, append index of maximizer(A'(y-Ax)), solve LS problem
function xe= FastFourierOMP(mask, y, sig,Niter)
    order=numel(size(mask));
    xe=zeros(size(mask));
    Nsz=size(mask);
    Nel=prod(size(mask));
    cind=[];
    suppv=[];
    tic;
    thresh=sig*sqrt(sum(mask(:)));
    for it=1:1:Niter
        Fxe=fft2(xe);
        Axe=Fxe(mask==1);
        resd=y-Axe;
            if(norm(resd)<thresh)       % stopping criteria
                break;
            end
        xdual=zeros(Nsz);
        xdual(mask==1)=resd;
        xdual=ifft2(xdual)*Nel;   % computes xdual=A'*(y-A*x) 
        [vl,ind]=max(abs(xdual(:)));
        cind=[cind,ind];
        [c1,c2]=ind2sub(Nsz,ind);
        vand1=exp(-1i*2*pi*(c1-1)*[0:1:Nsz(1)-1]/Nsz(1)); % Fourier vector corresponding to non-zero support locn
         vand2=exp(-1i*2*pi*(c2-1)*[0:1:Nsz(2)-1]/Nsz(2)); % Fourier vector corresponding to non-zero support locn
         [v11, v22] = ndgrid(1:Nsz(1), 1:Nsz(2));
         M = vand1(v11) .* vand2(v22);                   % outer product 
         suppv=[suppv,vec(M(mask==1))];                % suupv= A(:,suppindx) required for LS step in OMP
        % subsample M at locations of mask;
        xatsupp=pinv(suppv)*y;
        xe(cind)=xatsupp;
        % find maxindx of this tensor    
    end
end

