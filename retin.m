function v=retin(p,lm)
v=[]; 
for k=1:1:lm
    [a,b]=max(p);
    v=[v;b];
    p(b)=-1;
end
end