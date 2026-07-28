% Esame del 26 febbraio 2024 - pomeriggio c EULERO ALL'INDIETRO
% Simone Canevarolo
% S269893
% 22 aprile 2024

clear all
close all
clc

alt = 30e-2; % m
raggio = 10e-2; % m
altcaduta = 7.5e3; % m

T0 = 233; % K
hamb = 100; % W/m^2/K
Tamb = 270; % K

uu = 5; % m/s
ro = 950; % kg/m^3

K = @(T) 4.5*(T/T0);
cpfunz = @(T) 670+1.2*(T-233);

dx = 5e-3;

xx = (0:dx:alt/2)';
Nx = length(xx);

dt = 1;
tmax = altcaduta/uu;
tt = (0:dt:tmax);
Nt = length(tt);

Tm = T0*ones(Nx,1);

kk = K(T0);
cp = cpfunz(T0);

Tcentro = T0*ones(Nt,1);
Tsuperficie = T0*ones(Nt,1);

for zz = 2:Nt

        kk = K(Tm);
        cp = cpfunz(Tm);

        aa = dt.*kk./ro./cp./dx^2; 

        sub_diag_1 = -aa;
        main_diag_1 = 1+2*aa;
        sup_diag_1 = sub_diag_1;

        Band = [sub_diag_1,main_diag_1,sup_diag_1];

        AA_1 = spdiags(Band,-1:1,Nx,Nx);

        AA_1(1,1) = +kk(end)/dx+hamb;
        AA_1(1,2) = -kk(end)/dx;

        AA_1(end,end-1) = -1;
        AA_1(end,end) = 1;

        BB = Tm;

        BB(1) = +hamb*Tamb;
        BB(end) = 0;

        TT = AA_1\BB;

        Tsuperficie(zz) = TT(1);
        Tcentro(zz) = TT(end);

        Tm = TT;

end

    figure(1)
    plot(tt,Tcentro,'LineWidth',2)
    hold on
    plot(tt,Tsuperficie,'LineWidth',2)
    legend('T al centro','T superficie', 'Location','Best')
    grid on
