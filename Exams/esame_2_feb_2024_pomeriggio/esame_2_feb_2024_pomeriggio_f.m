% Esame del 2 febbraio 2024 pomeriggio, punto e
% Simone Canevarolo
% S269893
% 21 aprile 2024

clear all
close all
clc

raggio = 10e-2; % m
alt = 40e-2; % m
altcaduta = 5e3; % m
uu = 5; % m/s

ro = 950; % kg/m^3
hh = 100; % W/m^2/K
Taria = 280; % K
T0 = 233; % k

K = @(T) 4.5*(T/T0);
cpfunz = @(T) 670+1.2*(T-233);

dx = 1e-3;
xx = (0:dx:alt)';
Nx = length(xx);

tmax = altcaduta/uu;
dt = 1;
tt = (0:dt:tmax);
Nt = length(tt);

Tm = T0*ones(Nx,1);
kk = K(T0);
cp = cpfunz(T0);

Tmedia = T0*ones(Nx,1);

Potenza1D(1) = abs(hh*(T0-Taria)*pi*raggio^2)+abs(hh*(T0-Taria)*pi*raggio^2)+abs(trapz(xx,hh*(Taria-Tm)))*2*pi*raggio*alt;

for ii = 2:Nt

    kk = K(Tm);
    cp = cpfunz(Tm);

    aa = dt*kk./ro./cp./dx^2;

    sup_diag = -1/2*aa;
    main_diag = 1+aa;
    sub_diag = -1/2*aa;

    % Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    Band = [[sub_diag(2:end);0], main_diag, [sup_diag(1:end-1);0]];

    AA = spdiags(Band,-1:1,Nx,Nx);

    bb = Tm;

    AA(1,1) = +kk(1)/dx+hh;
    AA(1,2) = -kk(1)/dx;
    bb(1) = +hh*Taria;

    AA(end,end-1) = -kk(end)/dx;
    AA(end,end) = +kk(end)/dx+hh;
    bb(end) = +hh*Taria;

    TT = AA\bb;

    potenza1D(ii) = abs(hh*(TT(1)-Taria))*pi*raggio^2+abs(hh*(TT(end)-Taria))*pi*raggio^2+abs(trapz(xx,hh*(TT-Taria)))*2*pi*raggio*alt;

    Tm = TT;

end

figure(1)
plot(tt,potenza1D,'LineWidth',2)
hold on

Tm0 = T0;
TTzerod(1) = T0;
ii = 1;

Potenza0D(1)=abs(hh*(2*pi*raggio*alt+2*pi*raggio^2)*(Taria-Tm0));

for ii = 2:Nt

    cpzero = cpfunz(Tm0);

    cc = Tm0+dt*hh/ro/cpzero*2/raggio*Taria;

    Tzerod = (1+dt*hh*2/raggio/ro/cpzero)\cc;

%     TTzerod(ii) = Tzerod;

    Potenza0D(ii) = abs(hh*(Taria-Tzerod)*(2*pi*raggio*alt+2*pi*raggio^2));

    Tm0 = Tzerod;

end

figure(1)
plot(tt,Potenza0D,'LineWidth',2)
legend('Potenza 1D','Potenza 0D')

