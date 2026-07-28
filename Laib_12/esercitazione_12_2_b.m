% Esercitazione 12, esercizio 2 b
% Simone Canevarolo
% S269893
% 10 aprile 2024

clear all
close all
clc

ro = 170; % densità, kg/m^3
altezza = 18e-2; % m
base = 25e-2; % m
alt1 = 15e-2; % altezza dall'alto del punto del pelo dell'acqua, m

ro = 170; % densità, kg/m^3
cs = 210; % calore specifico, J/kg/K
cond = 2.5; % conducibilità termica, W/m/K

Taria = 30+273; % temperatrura dell'aria, K
hh = 10; % coefficiente di scambio termico dell'aria, W/m^2/K
Th2o = 5+273; % temperatura acqua, K

dxvett = [1e-3 5e-3 1e-2];
zz = 1;
err = zeros(length(dxvett),1);

for zz = 1:length(dxvett)

dx = dxvett(zz);
dy = dx;

xx = (0:dx:altezza)';
xx1 = (0:dx:alt1)';
yy = (0:dy:base)';

Nx = length(xx);
Nx1 = length(xx1);
Ny = length(yy);

Ntot = Nx*Ny;

Tm = Taria*ones(Ntot,1);

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

% bordo nord ovest - Robin aria

jj = 1;
for ii = 2:Nx1-1

    kk = ii;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+hh/(dy*cond));
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 2/dy^2;
    bb(kk) = -Taria*hh*2/(dy*cond);

end

% bordo sud ovest - Dirichlet acqua

kk = Nx1:Nx;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Th2o;

% bordo sud - Dirichlet acqua

kk = Nx:Nx:Ntot; %% Verificare che la sovrapposizione non causi problemi
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Th2o;

% Bordo sud est - Dirichlet acqua

kk = Ntot-Nx+Nx1:Ntot;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Th2o;

% Bordo nord est - Robin aria

jj = Ny;
for ii = 2:Nx1-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-hh/dx;
    AA(kk,kk+1) = 1/2/dx^2;
    bb(kk) = -Taria*hh/dy;

end

% Bordo nord - Robin aria

ii = 1;
for jj = 2:Ny-1

    kk = Nx*(jj-1)+ii;
    AA(kk,kk-Nx) = 1/2/dy^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-hh/dy;
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 1/2/dy^2;
    bb(kk) = -Taria*hh/dx;

end

% Studio vertice nord west

kk = 1;
AA(kk,kk) = -hh/dx-hh/dy-1/dx^2-1/dy^2;
AA(kk,kk+1) = 1/dx^2;
AA(kk,kk+Nx) = 1/dy^2;
bb(kk) = -Taria*(hh/dx+hh/dy);

% Studio vertice nord est

kk = Ntot-Nx+1;
AA(kk,kk-Nx) = 1/dy^2;
AA(kk,kk) = -hh/dx-hh/dy-1/dx^2-1/dy^2;
AA(kk,kk+1) = 1/dx^2;
bb(kk) = -Taria*(hh/dx+hh/dy);

TT = AA\bb;

err(zz) = norm(TT(Ntot-Nx+1)-Tm(Ntot-Nx+1))/norm(Tm(Ntot-Nx+1)-Taria);

Tm = TT;

% TTT = reshape(TT,Nx,Ny);

end

figure(1)
loglog(dxvett,err,'linewidth',2)
grid on
hold on

% figure(1)
% surf(xmat',ymat',TTT)
