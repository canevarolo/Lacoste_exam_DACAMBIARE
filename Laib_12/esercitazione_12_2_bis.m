% Esercizio 2, esercitazione 12
% Simone Canevarolo
% S269893
% 16/02/2024

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
haria = 10; % coefficiente di scambio termico dell'aria, W/m^2/K
Th2o = 5+273; % temperatura acqua, K

dx = 5e-3;
dy = dx;

xx = (0:dx:altezza)';
xx1 = (0:dx:alt1)';
yy = (0:dy:base)';

Nx = length(xx);
Nx1 = length(xx1);
Ny = length(yy);

Ntot = Nx*Ny;

%%

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

bb = zeros(Ntot,1);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;
        bb(kk) = 0;

    end
end

% bordo nord ovest, Robin

jj = 1;
for ii = 2:Nx1-1

    kk = ii;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+haria/cond/dy);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 2/dy^2;
    bb(kk) = -2*Taria*haria/dy/cond;

end

% bordo sud ovest, Dirichlet

jj = 1;
for ii = Nx1:Nx

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;

end

% bordo sud, Dirichlet

ii = Nx;
for jj = 1:Ny

    kk = ii*jj;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;

end

% bordo sud est, Dirichlet

jj = Ny;
for ii = Nx1:Nx

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;

end

% bordo nord est, Robin

jj = Ny;
for ii = 2:Nx1-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 2/dy^2;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+haria/cond/dx);
    AA(kk,kk+1) = 1/dx^2;
    bb(kk) = -2*Taria*haria/dy/cond;

end

% bordo nord, Robin

ii = 1;
for jj = 2:Ny-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+haria/cond/dx);
    AA(kk,kk+1) = 2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -2*Taria*haria/cond/dx;

end

% Vertice nord ovest

kk = 1;
AA(kk,kk) = -1/dx^2-1/dy^2-haria/cond/dy-haria/cond/dx;
AA(kk,kk+1) = 1/dx^2;
AA(kk,kk+Nx) = 1/dy^2;
bb(kk) = -Taria*haria/cond/dx-Taria*haria/cond/dy;

% Vertice nord est

kk = Ntot-Nx+1;
AA(kk,kk-Nx) = 1/dy^2;
AA(kk,kk) = -1/dx^2-1/dy^2-haria/cond/dy-haria/cond/dx;
AA(kk,kk+1) = 1/dx^2;
bb(kk) = -Taria*haria/cond/dx-Taria*haria/cond/dy;

% Soluzione

TT = AA\bb;

TTT = reshape(TT,Nx,Ny);

figure(1)
surf(xmat',ymat',TTT-273)
title('Andamento bidimensionale della temperatura sulla superficie')
% xlabel('Lunghezza')
% ylabel('Altezza')
% zlabel('Temperatura')

% T1 = TTT(1,1)-273

%% STUDIO INGEGNERISTICO

% Vettore nodi

Nodi = [50 100 200 500 1000];

dd = zeros(length(Nodi),1);
Tsave = zeros(length(Nodi),1);

for ss = 1:length(Nodi)

    Nxx = Nodi(ss);
    Nyy = Nxx;

    xx=linspace(0,base/2,Nxx);
    yy=linspace(0,altezza,Nyy);

dx=xx(2)-xx(1);
dy=yy(2)-yy(1);

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

bb = zeros(Ntot,1);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;
        bb(kk) = 0;

    end
end

% bordo nord ovest, Robin

jj = 1;
for ii = 2:Nx1-1

    kk = ii;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+haria/cond/dy);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 2/dy^2;
    bb(kk) = -2*Taria*haria/dy/cond;

end

% bordo sud ovest, Dirichlet

jj = 1;
for ii = Nx1:Nx

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;

end

% bordo sud, Dirichlet

ii = Nx;
for jj = 1:Ny

    kk = ii*jj;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;

end

% bordo sud est, Dirichlet

jj = Ny;
for ii = Nx1:Nx

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;

end

% bordo nord est, Robin

jj = Ny;
for ii = 2:Nx1-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 2/dy^2;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+haria/cond/dx);
    AA(kk,kk+1) = 1/dx^2;
    bb(kk) = -2*Taria*haria/dy/cond;

end

% bordo nord, Robin

ii = 1;
for jj = 2:Ny-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+haria/cond/dx);
    AA(kk,kk+1) = 2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -2*Taria*haria/cond/dx;

end

% Vertice nord ovest

kk = 1;
AA(kk,kk) = -1/dx^2-1/dy^2-haria/cond/dy-haria/cond/dx;
AA(kk,kk+1) = 1/dx^2;
AA(kk,kk+Nx) = 1/dy^2;
bb(kk) = -Taria*haria/cond/dx-Taria*haria/cond/dy;

% Vertice nord est

kk = Ntot-Nx+1;
AA(kk,kk-Nx) = 1/dy^2;
AA(kk,kk) = -1/dx^2-1/dy^2-haria/cond/dy-haria/cond/dx;
AA(kk,kk+1) = 1/dx^2;
bb(kk) = -Taria*haria/cond/dx-Taria*haria/cond/dy;

% Soluzione

TT = AA\bb;

TTT = reshape(TT,Nx,Ny);

Tsave(ss) = TTT(round(Nx/2),round(Ny,2));
dd(ss) = dx;

end

figure(2)
semilogx(dd,Tsave,'linewidth',2);