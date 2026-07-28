% Esercitazione 11, esercizio 1
% Simone Canevarolo
% S269893
% 8/04/2024

clear all
close all
clc

lato = 1;
Tfissa = 500; % temperatura imposta su tre lati, K
Taria = 300; % temperatura aria, K

ro = 1920; % densità, kg/m^3
cond = 0.72; % conducibilità termica, W/K/m
cp = 835; % calore specifico, J/kg/K
hh = 10; % coefficiente di scambio termico, W/m^2/K

dx = 1e-2;
dy = dx;

xx = (0:dx:lato)';
yy = xx;

Nx = length(xx);
Ny = Nx;
Ntot = Nx*Ny;

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);
bb = zeros(Ntot,1);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = Nx*(jj-1)+ii;
        AA(kk,kk-Nx) = 1/dx^2;
        AA(kk,kk-1) = 1/dy^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dy^2;
        AA(kk,kk+Nx) = 1/dx^2;

    end
end

% bordo ovest - Dirichlet
kk = 1:Nx;
% AA(kk,kk-1) = 0; superfluo
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
% AA(kk,kk+1) = 0;
% AA(kk,kk+Nx) = 0; superfluo
bb(kk) = Tfissa;

% bordo sud - Dirichlet
kk = Nx:Nx:Ntot;
% AA(kk,kk-1) = 0;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tfissa;

% bordo est - Dirichlet
kk = (jj-1)*Nx+1:Ntot;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
 % AA(kk,kk-Nx) = 0;
bb(kk) = Tfissa;

% % bordo Nord - Dirichlet
% kk = 1:Nx:Ntot-Nx+1;
% AA(kk,:) = 0;
% % AA(kk,kk+1) = 0;
% AA(kk,kk) = eye(length(kk));
% bb(kk) = Tfissa;

% bordo nord - Robin
ii=1;
for jj=2:Ny-1

    kk = Nx*(jj-1)+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk) = -2/dy^2-2/dx^2-2*hh/dx;
    AA(kk,kk+1) = 2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -2*hh/cond/dx*Taria;

end

TT = AA\bb;

TTT = reshape(TT,Nx,Ny);

figure(1)
surf(xmat',ymat',TTT)

