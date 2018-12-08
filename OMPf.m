function x=OMPf(Y,A,sig,it)
dum=zeros(size(A,2),1);
x=dum;
iter=0;
epsilon=sig*sqrt(size(A,1));
cind=[];
            while(1)
                    PH=[];
                    res=Y-A*x;                   
                    ind=retin(abs(A'*res),1);
                    D=union(cind,ind);
                    for ii=1:1:length(D)
                        PH=[PH,A(:,D(ii))];
                    end
                    tmpv=pinv(PH)*Y;%(inv(PH'*PH))*PH'*Y;
                    x=dum;
                    x(D)=tmpv;
                    cind=D;
                    iter=iter+1;
                    if(norm(res)<=epsilon | iter==it)
                        break;
                    end
            end  
         
end