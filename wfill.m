function palloc=wfill(Noise_Power,Trans_Power)

Number_Channel= length(Noise_Power) ; 
[S_Number dt]=sort(Noise_Power);
sum(Noise_Power);
for p=length(S_Number):-1:1
    T_P=(Trans_Power+sum(S_Number(1:p)))/p;
    Input_Power=T_P-S_Number;
    Pt=Input_Power(1:p);
    if(Pt(:)>=0),
        break
    end
end
Allocated_Power=zeros(1,Number_Channel);
Allocated_Power(dt(1:p))=Pt;

palloc=Allocated_Power;
end