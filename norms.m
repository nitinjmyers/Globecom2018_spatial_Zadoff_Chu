function N=norms(A)
    N=[];
    Acol=size(A,2);
    for vv=1:1:Acol
        N=[N,norm(A(:,vv))];
    end
end
