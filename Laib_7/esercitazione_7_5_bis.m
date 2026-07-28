% Esercitazione 7, esercizio 5
% Simone Canevarolo
% S269893
% 06/04/2024

clear all
close all
clc

din = 20e-2; % diametro interno, m
dout = 30e-2; % diametro esterno, m

kk = 0.2; % conducibilità termica, W/m/K
haria = 100; % coefficiente di scambio termico, W/K/m^2
Taria = 5; % temperatura aria, K
qsup = 50; % calore irradiato in superficie interna, W/m^2

raggio = (dout/2-din/2);
dr = 1e-4;
rr = (din/2:dr:dout/2)';
Nr = length(rr);

sup_diag = (1+dr./rr);
main_diag = -2*ones(Nr,1);
sub_diag = (1-dr./rr);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nr,Nr);

bb = zeros(Nr,1);

% Neumann interno
AA(1,1) = 1;
AA(1,2) = -1;
bb(1) = qsup*dr/kk;

% Robin esterno
AA(end,end-1) = -kk/dr;
AA(end,end) = haria+kk/dr;
bb(end) = haria*Taria;

TT = AA\bb;

figure(1)
plot(rr,TT,'LineWidth',2,'Color','r')
xlabel('spessore [m]')
ylabel('Temperatura [°C]')
title('Andamento temperatura nel guscio della zucca')