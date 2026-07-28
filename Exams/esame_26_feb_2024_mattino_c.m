% Esame del 26 febbraio 2024 - turno MATTINO
% Simone Canevarolo
% S269893
% 22 aprile 2024

clear all
close all
clc

ll = 0.1; % m
Tsud = 278; % K
Tnord = 296; % K

Taria = 283; % K
haria = 20; % W/M^2/K

cond1 = 0.1; % W/m/K
cond2 = 10; % W/m/K

dx = 1e-3;
dy = dx;

xx = (0:dx:ll)';
yy = (0:dy:ll)';

Nx = length(xx);
Ny = length(yy);

Ntot = Nx*Ny;

yy1 = (0:dy:ll/2)';
Ny1 = length(yy1);
Ntot1 = Ny1*Nx;

% yy2 = (ll/2:dy:ll)';
% Ny2 = length(yy2);
% Ntot2 = Ny2*Nx;

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

bb = zeros(Ntot,1);

% OVEST

jj = 1;
for ii = 2:Nx-1

    kk = ii;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-haria/dy;
    AA(kk,kk+1) = 1/2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -Taria*haria/dy;

end

% NORD

ii = 1;
for jj = 1:Ny

kk = (jj-1)*Nx+ii;
AA(1,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tnord;

end

% EST

jj = Ny;
for ii = Ntot-Nx+1:Ntot-1 % vedere se deve essere +2 e non +1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-haria/dy;
    AA(kk,kk+1) = 1/2/dx^2;
    bb(kk) = -Taria*haria/dy;

end

% SUD

ii = Nx;
for jj = 1:Ny

kk = (jj-1)*Nx+ii;
AA(end,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tsud;

end

TT = AA\bb;
TTT = reshape(TT,Nx,Ny);

surf(xmat,ymat,TTT')