function nvx=quantz(vx,qb)
 ovx=vx;
 vx= wrapTo2Pi(angle(vx));
 ql=2^qb;
 vx=round(ql*vx/(2*pi))*(2*pi/ql);
 nvx=exp(1i*vx);
% nvx=nvx*norm(ovx,'fro')/norm(nvx,'fro');
end