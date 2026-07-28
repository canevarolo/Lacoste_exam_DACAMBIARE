% Esercitazione 7, esercizio 5
% Simone Canevarolo
% S269893
% 8/12/2023

clear all
close all
clc

%%

Dout = 30e-2; % diametro esterno [m]
Din = 20e-2; % diametro interno [m]
Taria = 5; % temperatura [°C]
haria = 100; % coefficiente di scambio termico [W/m^2*K]
kk = 0.2; % conducibilità termica [W/m*K]
qv = 50; % generazione interna di calore [W/m^2]

rr = Dout/2-Din/2; % m
dx = 1e-4;
xx = (Din/2:dx:Dout/2)';
NN = length(xx);

%%

sup_diag = (1+dx./xx);
main_diag = -2*ones(NN,1);
sub_diag = (1-dx./xx);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = zeros(NN,1);

%%

% Robin

AA(1,1) = kk/dx;
AA(1,2) = -kk/dx;
bb(1) = qv;

% Robin

AA(end,end-1) = -kk/dx;
AA(end,end) = kk/dx+haria;
bb(end) = haria*Taria;

TT = AA\bb;

%%

figure(1)
plot(xx,TT,'LineWidth',2)
xlabel('spessore [m]')
ylabel('Temperatura [°C]')
title('Andamento temperatura nel guscio della zucca')
