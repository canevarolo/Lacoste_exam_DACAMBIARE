% Esame del 2 febbraio 2024, CRANK NICOLSON 
% Simone Canevarolo
% S269893
% 21 aprile 2024

clear all
close all


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
xx = (0:dx:alt/2)';
Nx = length(xx);

tmax = altcaduta/uu;
dt = 1;
tt = (0:dt:tmax);
Nt = length(tt);

Tm = T0*ones(Nx,1);
kk = K(T0);
cp = cpfunz(T0);

Tcentro(1) = T0;
Tsuperficie(1) = T0;

tic

for ii = 2:Nt
    

    kk = K(Tm);
    cp = cpfunz(Tm);

    aa = dt*kk./ro./cp/dx^2;

    % Matrice per Tm+1

    sup_diag_1 = -1/2*aa;
    main_diag_1 = (1+aa);
    sub_diag_1 = -1/2*aa;

    Band_1 = [[sub_diag_1(2:end);0], main_diag_1, [0;sup_diag_1(1:end-1)]];
  

    AA_1 = spdiags(Band_1,-1:1,Nx,Nx);

    bb_1 = 0;

    AA_1(1,1) = -kk(1)/dx;
    AA_1(1,2) = +kk(1)/dx+hh;


    AA_1(end,end-1) = -1;
    AA_1(end,end) = 1;


    % Matrice per Tm

    sup_diag_2 = -1/2*aa;
    main_diag_2 = (1+aa);
    sub_diag_2 = -1/2*aa;

    Band_2 = [[sub_diag_2(2:end);0], main_diag_2, [0;sup_diag_2(1:end-1)]];

    AA_2 = spdiags(Band_1,-1:1,Nx,Nx);

    bb_2 = 0;

    AA_2(1,1) = -kk(1)/dx;
    AA_2(1,2) = +kk(1)/dx+hh;

    AA_2(end,end-1) = -1;
    AA_2(end,end) = 1;

    % BB=(1/2*(bb_1+bb_2))*zeros(Nx,1)+Tm*AA_2;
BB = AA_2*Tm;
    BB(1) = +hh*Taria;
    BB(end) = 0;


    TT = AA_1\BB;

    Tcentro(ii) = TT(end);
    Tsuperficie(ii) = TT(1);
    Tm = TT;

end

toc

figure(1)
plot(tt,Tcentro,'LineWidth',2)
hold on
plot(tt,Tsuperficie,'linewidth',2)
% plot(xx,TT)





% % Esame del 2 febbraio 2024 pomeriggio, punto b EULERO ALL'INDIETRO
% % Simone Canevarolo
% % S269893
% % 21 aprile 2024
% 
% clear all
% close all
% clc
% 
% raggio = 10e-2; % m
% alt = 40e-2; % m
% altcaduta = 5e3; % m
% uu = 5; % m/s
% 
% ro = 950; % kg/m^3
% hh = 100; % W/m^2/K
% Taria = 280; % K
% T0 = 233; % k
% 
% K = @(T) 4.5*(T/T0);
% cpfunz = @(T) 670+1.2*(T-233);
% 
% dx = 1e-3;
% xx = (0:dx:alt/2)';
% Nx = length(xx);
% 
% tmax = altcaduta/uu;
% dt = 1;
% tt = (0:dt:tmax);
% Nt = length(tt);
% 
% Tm = T0*ones(Nx,1);
% kk = K(T0);
% cp = cpfunz(T0);
% 
% Tsuperficie = T0*ones(Nx,1);
% Tcentro = T0*ones(Nx,1);
% 
% for ii = 2:Nt
% 
%     kk = K(Tm);
%     cp = cpfunz(Tm);
% 
%     aa = dt*kk./ro./cp./dx^2;
% 
%     sup_diag1 = -1/2*aa;
%     main_diag1 = 1+aa;
%     sub_diag1 = sup_diag1;
% 
%     % Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
%     Band1 = [sub_diag1, main_diag1, sup_diag1];
% 
%     AA1 = spdiags(Band1,-1:1,Nx,Nx);
% 
%     AA1(1,1) = +kk(1)/dx+hh;
%     AA1(1,2) = -kk(1)/dx;
% 
%     AA1(end,end-1) = -1;
%     AA1(end,end) = 1;
% 
% 
% 
%     sup_diag2 = 1/2*aa;
%     main_diag2 = 1-aa;
%     sub_diag2 = sup_diag2;
% 
%     % Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
%     Band2 = [sub_diag2, main_diag2, sup_diag2];
% 
%     AA2 = spdiags(Band2,-1:1,Nx,Nx);
% 
%     AA2(1,1) = kk(1)/dx+hh;
%     AA2(1,2) = -kk(1)/dx;
% 
%     AA2(end,end-1) = -1;
%     AA2(end,end) = 1;
% 
%     bb = Tm;
% 
%     bb(1) = -hh*Taria;
%     bb(end) = 0;
% 
%     TT = AA1\(AA2*bb);
% 
%     Tcentro(ii) = TT(end);
%     Tsuperficie(ii) = TT(1);
% 
%     Tm = TT;
% 
% end
% 
% figure(1)
% plot(tt,Tcentro,'LineWidth',2)
% hold on
% plot(tt,Tsuperficie,'linewidth',2)
% % plot(xx,TT)
% 
