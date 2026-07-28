% Esercitazione 11, esercizio 4
% Simone Canevarolo
% S269893
% 13/02/2024

clear all
close all
clc

kchip = 50; % conducibilità termica chip, W/m/K
ksub = 5; % conducibilità termica substrato, W/m/K
Tout = 20+273; % temperatura esterna lato nord, K
hout = 500; % coefficiente di scambio termico, W/m^2/K
qvol = 1e7; % generazione interna di calore, W/m^3

lungh = 13.5e-3; % lunghezza, m
alt = 12e-3; % altezza, m

dx = 1e-4;
dy = dx;

yy = (0:dy:lungh)';
xx = (0:dx:alt)';

Nx = length(xx);
Ny = length(yy);

Ntot = Nx*Ny;

lungh1 = 9e-3; % lunghezza primo rettangolo, m
lungh2 = lungh-lungh1; % lunghezza secondo e terzo rettangolo, m
alt1 = 3e-3; % altezza secondo rettangolo, m
alt2 = alt-alt1; % altezza terzo rettangolo, m

yy1 = (0:dx:lungh1)';
yy2 = (lungh1+1:dx:lungh)';

xx1 = (0:dx:alt1)';
xx2 = (alt1+1:dx:alt)';

Nx1 = length(xx1);
Nx2 = length(xx2);
Ny1 = length(yy1);
Ny2 = length(yy2);

Ntot1 = Nx*Ny1;
Ntot2 = Ntot1 + Nx1*Ny2;
% Ntot3 = xx2*yy2;

%%

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

[xmat,ymat] = meshgrid(xx,yy);

% griglia interna rettangolo 1

for ii = 2:Nx-1
    for jj = 2:Ny1-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;
        bb(kk) = 0;

    end
end

% lato ovest rettangolo 1, Neumann adiabatico

jj = 1;
for ii = 2:Nx-1

    kk = ii;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 2/dy^2;
    bb(kk) = 0;

end

% lato sud rettangolo 1, Neumann adiabatico

ii = Nx;
for jj = 2:Ny-1

    kk = ii*jj;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 2/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = 0;

end

% lato nord rettangolo 1, Robin

ii = 1;
for jj = 2:Ny-1
    
    kk = (jj-1)+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk) = -2*(1/dy^2+1/dx^2+hout/ksub/dx);
    AA(kk,kk+1) = 2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -2*hout*Tout/ksub/dx;

end

% Griglia interna rettangolo 2 "chip"

for ii = 2:Nx1-1
    for jj = yy1+1:yy-1

        kk = Ntot1+(jj-1)*xx1+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;
        bb(kk) = -qvol/kchip;

    end
end

% lato ovest rettangolo 2 "chip", Neumann adiabatico

jj = yy1;
for ii = 2:xx1-1

    kk = Ntot1+ii;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+hout/kchip/dx);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx1) = 2/dy^2;
    bb(kk) = 0;

end


% lato nord rettangolo 2 "chip", Robin

ii = 1;
for jj = xx1+2:xx-1

    kk = Ntot1+(jj-Ny1-1)*Nx1+ii;
    AA(kk,kk-Nx1) = 1/dy^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+hout/kchip/dy);
    AA(kk,kk+1) = 2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -2*hout*Tout/kchip/dy;

end


% lato est rettangolo 2 "chip", Neumann adiabatico

jj = yy;
for ii = 2:xx1-1

    kk = Ntot1+(Ny2-1)*ii;
    AA(kk,kk-Nx1) = 2*(1/dy^2);
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/dx^2;
    bb(kk) = 0;

end

% lato sud rettangolo 2 "chip", Robin

ii = xx1;
for jj = yy1+1:yy-1

    kk = Ntot1+(yy-yy1+ii-1)*Nx1+jj;
    AA(kk,kk-Nx1) = 1/dy^2;
    AA(kk,kk-1) = 2/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+Nx1) = 1/dy^2;
    bb(kk) = -2*hout*Tout/kchip/dx;

end


% Griglia interna rettangolo 3

for ii = xx1+1:xx-1
    for jj = yy1+1:yy-1 

        kk = Ntot2 + (jj-yy1-1)*Nx2 + ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;
        bb(kk) = 0;

    end
end


% lato ovest rettangolo 3, Neumann adiabatico



% lato nord rettangolo 3, Neumann adiabatico



% lato est rettangolo 3, Neumann adiabatico



% lato sud rettangolo 3, Neumann adiabatico



