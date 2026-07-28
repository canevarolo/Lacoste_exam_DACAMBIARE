% Esercitazione 11, esercizio 2
% Simone Canevarolo
% S269893
% 25/01/2024

clear all
close all
clc

%%

ll = 0.2; % larghezza pezzo, m
lbase = 0.1; % larghezza della base, m
% lcentro = 0.2; % larghezza del pezzo nella zona centrale, m dato poco utile
altezza = 0.1; % altezza complessiva del pezzo, m
altezzabase = 0.05; % altezza delle parti che sorreggono la base, m

roacc = 7800; % densità acciaio, kg/m^3
cond = 1; % conducibilità termica, W/m/K ----> non lo posso indicare con kk
cs = 450; % calore specifico, J/kg/K

Test = 5; % temperatura imposta esterna, K
Tin = 20; % temperatura imposta interna, K

%%

dx = 1e-2; % m
dy = dx;

xx = (0:dx:altezza)';
xx1 = (0:dx:altezzabase)';

yy = (0:dy:ll)';
yy1 = (0:dy:lbase)';

Nx = length(xx);
Nx1 = length(xx1);

Ny = length(yy);
Ny1 = length(yy1);

Ntot = Nx*Ny1+Nx1*(Ny-Ny1);
Ntot1 = Nx*Ny1;
Ntot2 = Nx1*(Ny-Ny1);

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

% escludo il cornicione della matrice

bb = zeros(Ntot,1);

% Studio il dominio quadrato a sinistra

for ii = 2:Nx-1
    for jj = 2:Ny1-1

    kk = Nx*(jj-1)+ii;
    AA(kk,kk-Nx) = 1/dy^2;    % posso anche non moltiplicare per cond, dato che si semplifica con lo zero all'altro membro nell'equazione generale
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;

    end
end

% Studio il dominio rettangolare a destra

for ii = 2:(Nx1-1)
    for jj = Ny1:Ny-1

        if jj == Ny1    % Neumann adiabatico sul lato

            kk = Ntot1-Nx+ii;
            AA(kk,kk-Nx) = 1/dy^2;
            AA(kk,kk-1) = 1/dx^2;
            AA(kk,kk) = -2*(1/dx^2+1/dy^2);
            AA(kk,kk+1) = 1/dx^2;
            AA(kk,kk+Nx) = 1/dy^2;

        elseif  jj == Ny1+1
               
            kk=Ntot1+ii;
            AA(kk,kk-Nx)=1/dy^2;
            AA(kk,kk-1)=1/dx^2;
            AA(kk,kk)=-2*(1/dx^2+1/dy^2);
            AA(kk,kk+1)=1/dx^2;
            AA(kk,kk+Nx1)=1/dy^2;
                   
        else

            kk = Ntot1+(jj-Ny1-1)*Nx1+ii;       % uso Nx1 perchè devo "scendere" o "salire" sulla colonna di un altro valore
            AA(kk,kk-Nx1) = 1/dy^2;  
            AA(kk,kk-1) = 1/dx^2;
            AA(kk,kk) = -2*(1/dx^2+1/dy^2);
            AA(kk,kk+1) = 1/dx^2;
            AA(kk,kk+Nx1) = 1/dy^2;

        end
    end
end


%% Inizio con le condizioni al contorno

% Dirichlet, quadrato sx a ovest
jj = 1;
for ii = 2:Nx

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Test;

end

% Dirichlet, quadrato sx a nord

ii = 1;
for jj = 1:Ny1

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Test;

end

% Dirichlet, quadrato a sx a sudest

jj = Ny1;
for ii = Nx1:Nx

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tin;

end

% Neumann adiabatico, quadrato a sx a sud

ii = Nx;
for jj = 2:Ny1-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 2*(1/dx^2);
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+Nx) = 1/dy^2;
    % bb non lo metto perchè è già dichiarato in precedenza che è zeros

end

% Dirichlet, rettangolo dx a nord

ii=1;
for jj = Ny1+1:Ny

    kk = Ntot1+(jj-Ny1-1)*Nx1+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Test;

end

% Neumann adiabatico, rettangolo dx a est

jj = Ny;
for ii = 2:Nx1-1

    kk = Ntot-Nx1+ii;
    AA(kk,kk-Nx1) = 2*(1/dy^2);
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/dx^2;
    % bb non lo metto perchè è già dichiarato in precedenza che è zeros

end

% Dirichlet, rettangolo dx a sud

ii = Nx1;
for jj = Ny1+1:Ny

     kk = Ntot1+(jj-Ny1-1)*Nx1+ii;
     AA(kk,kk) = eye(length(kk));
     bb(kk) = Tin;

end

%%

TT = AA\bb;

TTT1 = reshape(TT(1:Ntot1),Nx,Ny1);
TTT2 = reshape(TT(Ntot1+1:end),Nx1,Ny-Ny1);
TTT3 = NaN*ones(Nx-Nx1,Ny-Ny1);

TTT = [TTT1,[TTT2;TTT3]];

figure (1)
surf(xmat,ymat,TTT')

xlabel('x','fontsize',16)
ylabel('y','fontsize',16)
zlabel('Temperature','fontsize',16)
