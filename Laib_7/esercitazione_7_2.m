% Esercitazione 7, esercizio 2
% Simone Canevarolo
% S269893
% 28/11/2023

clear all
close all
clc

%%

DD = 1e-2; % diametro m
lltot = 4; % lunghezza totale pezzo m
ll = 2; % lunghezza del pezzo considerato m
II = 1e3; % corrente A
Tb = 77; % temperatura imposta pareti K
hh = 500; % coefficiente di scambio termico W/(m^2*K)

roel = 1.75e-8; % resistività elettrica ohm*m
kk = 350; % conducibilità elettrica W/(m*K)
As = ll*2*pi*(DD/2); % Area laterale m^2
Asez = pi*DD^2/4; % Area sezione m^2
Vol = ll*pi*(DD/2)^2; % Volume m

qv = roel*II^2/Asez^2;

%%

dx = 1e-3; % m
xx = (0:dx:ll/2)';
NN = length(xx);

sub_diag = ones(NN,1);
main_diag = -(2+hh*As/Vol*dx^2/kk)*ones(NN,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = ((-qv-hh*As*Tb/Vol)*dx^2/kk)*ones(NN,1);

%%

% Dirichlet

AA(1,1) = 1;
AA(1,2) = 0;
bb(1) = Tb;

% Neumann, adiabatica

AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = 0;

TT = AA\bb;

%%

figure(1)
plot(xx,TT,'linewidth',2)
xlabel('lunghezza [m]')
ylabel('temperatura [K]')