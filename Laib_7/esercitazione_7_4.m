% Esercitazione 7, esercizio 4
% Simone Canevarolo
% S269893
% 1/12/2023

%%

ll = 5e-3; % spessore m
Dout = 2e-2; % diametro esterno m
Din = 1e-2; % diametro interno m
kk = 0.1; % conduciblità termica W/(m*K)
qv = 2e6; % generazione interna di calore W/m^3
Text = 300; % temperatura °C
hext = 10; % coefficiente di scambio termico W/(m^2*K)

dx = 5e-6; % m
xx = (0:dx:ll)';
NN = length(xx);

%%

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = -qv*dx^2/kk*ones(NN,1);

%%

% Robin

AA(1,1) = kk/dx+hext;
AA(1,2) = -kk/dx;
bb(1) = hext*Text;

% Dirichlet

AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = 0;

TT = AA\bb;
Tmax = max(TT)

%% Conservazione dell'energia

Flusso_generato = qv*pi*(Dout^2-Din^2)/4;
Flusso_uscente = abs(hh*(Text-TT(end)))*pi*Dout

%%

figure(1)
plot(xx,TT,'linewidth',2)
xlabel('spessore [m]')
ylabel('Temperatura [°C]')
title('Distribuzione di temperatura')
