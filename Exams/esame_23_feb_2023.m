% Esame 23 febbraio 2023
% Simone Canevarolo
% S269893
% 1/07/2024

clear all
close all
clc

ssacc = 1e-2;
roacc = 7800;
cpacc = 500;
kkacc = 10;

ssis = 1.2e-2;
rois = 2300;
cpis = 2100;
kkis = 0.045;

ss = ssis+ssacc;

uu = 2;
Th2o = 373;
hhh2o = 250;
hharia = 50;

Taria = @(t) 40-25*sin(pi*t/250);
% tt = linspace(0,1000)
% plot(tt,Taria(tt))

dr = 1e-4; % m
rracc = (0:dr:ssacc)';
Nracc = length(rracc);

rris = (0:dr:ssis)';
Nris = length(rris);

rr = (0:dr:ss)';
Nr = length(rr);

dt = 1; % s
tt = (0:dt:40000);
Nt = length(tt);

Tm = Th2o*ones(Nr,1);
TTA = Th2o*ones(Nt,1);

for ii = 2:Nt

    kk = kkacc*ones(Nr,1);
    kk(Nracc:Nr) = kkis;

    ro = roacc*ones(Nr,1);
    ro(Nracc:Nr) = rois;

    cp = cpacc*ones(Nr,1);
    cp(Nracc:Nr) = cpis;

    aa = kk*dt./ro./cp/dr^2;

    sub_diag = -aa.*(1-dr./rr);
    main_diag = (1+2.*aa);
    sup_diag = -aa.*(1+dr./rr);

    Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

    AA = spdiags(Band,-1:1,Nr,Nr);

    bb = Tm;

    AA(1,1) = hhh2o+kkacc/dr;
    AA(1,2) = -kkacc/dr;
    bb(1) = hhh2o*Th2o;

    AA(Nracc,Nracc-1) = -1;
    AA(Nracc,Nracc) = 1+kkis/kkacc;
    AA(Nracc,Nracc+1) = -kkis/kkacc;
    bb(Nracc) = 0;

    AA(end,end-1) = -kkis/dr;
    AA(end,end) = hharia+kkis/dr;
    bb(end) = hharia*(Taria(ii)+273);

    TT = AA\bb;

    TTA(ii) = TT(Nracc);

    Tm = TT;

end

figure(1)
plot(rr,TT-273)

figure(2)
plot(tt,TTA-273)




