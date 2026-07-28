% Esame del 2 febbraio 2024 punto b - turno mattino 
% Simone Canevarolo
% S269893

clear all
close all
clc

raggio = 20e-2; % m
alt = 7500; % m
vv = 5; % m/s
T0 = -40+273; % K
hh = 100; % W/m^2/K
Tair = 270; % K

rovol = 1e3; % kg/m^3

kk = @(T) 5.*(T./T0); % W/m/K

cp = @(T) 670 + 1.2*(T-233); % J/kg/K


dr = 1e-3; % m
rr = (0:dr:raggio/2)';
Nr = length(rr);

dt = 1; % s
tmax = alt/vv; % s
tt = (0:dt:tmax);
Nt = length(tt);

Tm = T0*ones(Nr,1);
Tcentro = T0*ones(Nt,1);
Tsup = T0*ones(Nt,1);


for ii = 2:Nt

    aa = kk(Tm)*dt/dr^2/rovol./cp(Tm);
    
    sub_diag = -aa.*(1-dr./rr);
    main_diag = 1+2.*aa;
    sup_diag = -aa.*(1+dr./rr);
    
    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    
    AA = spdiags(Band,-1:1,Nr,Nr);
    
    bb = Tm;

    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = 0;

    AA(end,end-1) = -kk(Tm(end))/dr;
    AA(end,end) = hh+kk(Tm(end))/dr;
    bb(end) = hh*Tair;
    
    TT = AA\bb;

    Tcentro(ii) = TT(1);
    Tsup(ii) = TT(end);

    Tm = TT;

end


figure(1)
plot(tt,Tcentro,'LineWidth',2)
title('Temperature at the center')
xlabel('time [s]')
ylabel('Temperature [K]')
hold on
plot(tt,Tsup,'LineWidth',2)
title('Temperature on surface')
xlabel('time [s]')
ylabel('Temperature [K]')
legend('T centro','T superficie','Location','best')
