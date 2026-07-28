% Esercitazione 12, esercizio 3
% Simone Canevarolo
% S269893
% 06/02/2024

clear all
close all
clc

ll1 = 0.25; % lato orizzontale, m
ll2 = 0.15; % lato verticale, m
kk = 2; % conducibilità termica, W/m/K
ro = 8e3; % densità, kg/m^3
cs = 500; % calore specifico, J/kg/K

dx = 5e-3;
dy = dx;

xx = (0:dx:ll2)';
yy = (0:dy:ll1)';

Nx = length(xx);
Ny = length(yy);
Ntot = Nx*Ny;

%%

[xmat,ymat] = meshgrid(xx,yy);

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

% Condizioni al contorno

% Est
jj = 1;
for ii = 2:Nx-1

    kk = ii;
    AA(:,kk) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = 0;

end

% Nord

ii = 1;
for jj = 2:Ny-1

    

end


