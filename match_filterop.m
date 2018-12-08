function B=match_filterop(we)
    N=13;
    taps=13;
    %C=zeros(1,13);
    %C(1)=1;
    C=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];  % length 13 barker code
    Cr=C.*exp(1i*we*[0:1:N-1]);

    MF=conv(Cr,C(N:-1:1));
    Len=length(MF);
    B=[];
    for j=0:1:taps-1
        B=[B;exp(1i*we)*MF(N-j:Len-j)];
    end
end