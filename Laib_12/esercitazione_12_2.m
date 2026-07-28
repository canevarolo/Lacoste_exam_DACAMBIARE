% Esercitazione 12, esercizio 2
% Simone Canevarolo
% S269893
% 31/01/2024

clear all
close all
clc

%% DATI

ro = 170; % densità, kg/m^3
qs = 210; % J/kg/K
cond = 2.5; % conducibilità, W/m/K
Th2o = 5; % temperatura acqua, °C
Taria = 30; % temperatura aria, °C
hh = 10; % coefficiente di scambio termico, W/m^2/K
qq = 1e3; % flusso termico, W/m^2

llx = 18e-2; % m
lly = 25e-2; % m
llh2o = 3e-2; % altezza pelo libero acqua considerato perfettamente laminare, m
llaria = llx-llh2o; % m
llrisc = 5e-2; % lunghezza riscaldatore, m

%% STAZIONARIO

dx = 1e-2;
dy = dx;

xx = (0:dx:llx)';
yy = (0:dy:lly)';
xxh2o = (llaria:dx:llh2o)';
xxaria = (0:dx:llaria)';

Nx = length(xx);
Ny = length(yy);
Nxh2o = length(xxh2o);
Nxaria = length(xxaria);

Ntot = Nx*Ny;
Ntotaria = Ny*Nxaria;
Ntoth2o = Ny*Nxh2o;

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond*1/dy^2;
        AA(kk,kk-1) = cond*1/dx^2;
        AA(kk,kk) = -2*cond*(1/dy^2+1/dx^2);
        AA(kk,kk+1) = cond*1/dx^2;
        AA(kk,kk+Nx) = cond*1/dy^2;

    end
end

bb = zeros(Ntot,1);

% Condizioni al contorno
% Nord ovest - Robin aria
jj = 1;
for ii = 2:Nxaria-1

    kk = ii;
    AA(kk,kk-1) = cond/2/dx^2;
    AA(kk,kk) = -(cond/dx^2+cond/dy^2+2*hh/kk/dy);
    AA(kk,kk+1) = cond/2/dx^2;
    AA(kk,kk+Nx) = cond/dy^2;
    bb(kk) = -2*Taria*hh/kk/dy;

end

% Nord - Robin aria

ii = 1;
for jj = 2:Ny-1

   kk = (jj-1)*Nx+ii;
   AA(kk,kk-Nx) = cond/2/dy^2;
   AA(kk,kk) = -(cond/dx^2+cond/dy^2+2*hh/kk/dy);
   AA(kk,kk+1) = cond/dx^2;
   AA(kk,kk+Nx) = cond/2/dy^2;
   bb(kk) = -2*Taria*hh/kk/dx;

end

% Nord est - Robin aria

jj = Ny;
for ii = 2:Nxaria-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = cond/dy^2;
    AA(kk,kk-1) = cond/2/dx^2;
    AA(kk,kk) = -(cond/dx^2+cond/dy^2+2*hh/kk/dy);
    AA(kk,kk+1) = cond/2/dx^2;
    bb(kk) = -2*Taria*hh/dy;

end

% Sud est - Dirichlet acqua

jj = Ny;
for ii = Nxaria:Nx

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = 1;
    bb(kk) = Th2o;

end

% Sud - Dirichlet aria

ii = Nx;
for jj = 1:Ny

    kk = jj*Nx;
    AA(kk,:) = 0;
    AA(kk,kk) = 1;
    bb(kk) = Th2o;

end

% Sud Ovest - Dirichlet aria

jj = 1;
for ii = Nxaria:Nx

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk) = 1;
    bb(kk) = Th2o;

end

% Vertice nord ovest

% AA(1,1) = -cond/2/dy^2-cond/2/dx^2-hh/2/dx-hh/2/dy;
% AA(1,2) = cond/2/dx^2;
% AA(1,Nx+1) = cond/2/dy^2;
% bb(1) = hh/2/dx+hh/2/dy;

kk=1;
AA(kk,kk)=-(1/dx^2+1/dy^2+((hh/cond)*(1/dx+1/dy)));
AA(kk,kk+1)=1/dx^2;
AA(kk,kk+Nx)=1/dy^2;
bb(kk)=-(hh/cond)*(1/dx+1/dy)*Taria;

% Vertice nord est

% jj = Ny;
% AA(1,(jj-1)*Nx+1) = -cond/2/dy^2-cond/2/dx^2-hh/2/dx-hh/2/dy;
% AA(1,(jj-1)*Nx+2) = cond/2/dx^2;
% AA(1,(jj-2)*Nx+1) = cond/2/dy^2;
% bb((jj-1)*Nx+1) = hh/2/dx+hh/2/dy;

kk=Ntot-Nx+1;
AA(kk,kk)=-(1/dx^2+1/dy^2+((hh/cond)*(1/dx+1/dy)));
AA(kk,kk-Nx)=1/dy^2;
AA(kk,kk+1)=1/dx^2;
bb(kk)=-(hh/cond)*(1/dx+1/dy)*Taria;

TT = AA\bb;

TTT=reshape(TT,Nx,Ny);

figure(1)
surf(xmat',ymat',TTT)




