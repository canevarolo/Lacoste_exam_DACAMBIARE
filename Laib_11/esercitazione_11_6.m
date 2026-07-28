% Esercitazione
% Simone Canevarolo
% S269893
% 10/05/2024

clear all
close all
clc

ww = 6e-3; % m
ll = 48e-3; % m

kk = 50; % W/m/K
hh = 500; % W/m^2/K

Taria = 303; % K
Tbase = 373; % K

dx = 1e-5;
dy = dx;

xx = (0:dx:ww)';
yy = (0:dy:ll)';

Nx = length(xx);
Ny = length(yy);

Ntot = Nx*Ny;

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (ii-1)*Nx+jj;
        AA(kk,kk-Nx) = 1/dx^2;
        AA(kk,kk-1) = 1/dy^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dy^2;
        AA(kk,kk+Nx) = 1/dx^2;

    end
end

ii = 1;
for jj = 2:Ny-1

    kk = (ii-1)*Nx+jj;
    AA(kk,kk-Nx) = 1/2/dy^2;
    AA(kk,kk) = -(1/dx^2+1/dy^2+hh/kk/dx);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 1/2/dy^2;
    bb(kk) = -Taria*hh/kk/dx;

end

jj = Ny;
for ii = 2:Nx-1

    kk = (ii-1)*Nx+jj;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/2/dx^2;
    bb(kk) = 0;

end

ii = Nx;
for jj = 2:Ny-1

    kk = (ii-1)*Nx+jj;
    AA(kk,kk-Nx) = 1/2/dy^2;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -(1/dx^2+1/dy^2+hh/kk/dx);
    AA(kk,kk+Nx) = 1/2/dy^2;
    bb(kk) = -Taria*hh/kk/dx;

end

jj = 1;
for ii = 1:Nx

    kk = (ii-1)*Nx+jj;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tbase;

end



