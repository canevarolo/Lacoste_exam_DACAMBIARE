% % Esame del 26 febbraio 2024 - pomeriggio e
% % Simone Canevarolo
% % S269893
% % 22 aprile 2024
% 
% clear all
% close all
% clc
% 
% alt = 30e-2; % m
% raggio = 10e-2; % m
% altcaduta = 7.5e3; % m
% 
% T0 = 233; % K
% hamb = 100; % W/m^2/K
% Tamb = 270; % K
% 
% uu = 5; % m/s
% ro = 950; % kg/m^3
% 
% K = @(T) 4.5*(T/T0);
% cpfunz = @(T) 670+1.2*(T-233);
% 
% kk = K(T0);
% cp = cpfunz(T0);
% 
% dtvett = [1 2 5 10 20 50 100];
% dxvett = [1e-3 5e-3 1e-2];
% 
% dx = dxvett(1);
% xx = (0:dx:alt/2)';
% Nx = length(xx);
% 
% Tm = T0*ones(Nx,1);
% 
% for zz = 1:length(dtvett)
% 
%     dt = dtvett(zz);
%     tmax = altcaduta/uu;
%     tt = (0:dt:tmax);
%     Nt = length(tt);
% 
% for ii = 2:length(dxvett)
% 
%         dx = dxvett(ii);
%         xx = (0:dx:alt/2)';
%         Nx = length(xx);
% 
%         kk = K(Tm);
%         cp = cpfunz(Tm);
% 
%         aa = dt.*kk./ro./cp./dx^2; 
% 
%         sub_diag_1 = -aa;
%         main_diag_1 = 1+2*aa;
%         sup_diag_1 = sub_diag_1;
% 
%         Band = [sub_diag_1,main_diag_1,sup_diag_1];
% 
%         AA_1 = spdiags(Band,-1:1,Nx,Nx);
% 
%         AA_1(1,1) = +kk(end)/dx+hamb;
%         AA_1(1,2) = -kk(end)/dx;
% 
%         AA_1(end,end-1) = -1;
%         AA_1(end,end) = 1;
% 
%         BB = Tm;
% 
%         BB(1) = +hamb*Tamb;
%         BB(end) = 0;
% 
%         TT = AA_1\BB;
% 
%         err = norm(TT-Tm)/norm(Tm-T0);
% 
%         Tm = TT;
% 
%     end
% 
%     figure(1)
%     loglog(dxvett,err,'LineWidth',2)
%     hold on
% 
% end
% 
% 



% Esame del 26 febbraio 2024 - pomeriggio e
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

kk = K(T0);
cp = cpfunz(T0);

dtvett = [0.1 0.5 1 2 5 10 20];
% dxvett = [1e-3 5e-3 1e-2];

dx = 1e-3;
xx = (0:dx:alt/2)';
Nx = length(xx);

Tm = T0*ones(Nx,1);
tmax = altcaduta/uu;

for zz = 1:length(dtvett)

    dt = dtvett(zz);
    tt = (0:dt:tmax);
    Nt = length(tt);

    for ii = 2:Nt

        kk = K(Tm);
        cp = cpfunz(Tm);

        aa = dt*kk./ro./cp./dx^2; 

        sub_diag = -aa;
        main_diag = 1+2*aa;
        sup_diag = sub_diag;

        Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

        AA = spdiags(Band,-1:1,Nx,Nx);

        AA(1,1) = +kk(1)/dx+hamb;
        AA(1,2) = -kk(1)/dx;

        AA(end,end-1) = -1;
        AA(end,end) = 1;

        bb = Tm;

        bb(1) = +hamb*Tamb;
        bb(end) = 0;

        TT = AA\bb;

        err(zz) = norm(TT-Tm)/norm(Tm-T0);

        Tm = TT;

    end

end

    figure(1)
    loglog(dtvett,err,'LineWidth',2)
   grid on

