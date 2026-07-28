% Esercitazione 7, esercizio 3
% Simone Canevarolo
% S269893
% 29/11/2023

clear all
close all
clc

%%

kk = 0.55; % conducibilità termica, W/(m*K)
hext = 25; % coefficiente di scambio termico, W/(m^2*K)
Text = -5; % temperatura, °C
Tint = 20; % temperatura, °C
qv = 1e3; % generazione interna di calore, W/m^2
ll = 25e-2; % spessore, m

dx = 5e-4; % m
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

AA(1,1) = 1+kk/(dx*hext);
AA(1,2) = -kk/(hext*dx);
bb(1) = Text;

% Dirichlet

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tint;

TT = AA\bb;

Tmax = max(TT)

%%

figure(1)
plot(xx,TT,'LineWidth',2)
xlabel('Spessore [m]');
ylabel('Temperatura [°C]');