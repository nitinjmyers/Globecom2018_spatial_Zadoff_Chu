function [R,C,indx]=gentrajectory(Nx,Ny,Meas)
        indx=randperm(Nx*Ny,Meas);
        [R,C]=ind2sub([Nx,Ny],indx);
        indx=sort(indx.');
end