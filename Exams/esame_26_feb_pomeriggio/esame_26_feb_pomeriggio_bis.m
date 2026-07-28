% Esame 26 febbraio 2024
% Simone Canevarolo
% S269893
% 1/07/2024

clear all
close all
clc

llx = 10e-2; % m
lly = llx;
lly1 = lly/2;
cond1 = 0.1; % W/m/K
cond2 = 1; % W/m/K

condave = (cond1+cond2)/2;

Taria = 293; % K
haria = 20; % W/m^2/K
halto = 3e5; % W/m^2/K

Tnord = 296; % K
Tsud = 278; % K

dx = 5e-3;
xx = (0:dx:llx)';
Nx = length(xx);

dy = dx;
yy = (0:dy:lly)';
yy1 = (0:dy:lly1)';
Ny = length(yy);
Ny1 = length(yy1);

Ntot1 = Ny1*Nx;
Ntot = Nx*Ny;

% caso vettori "in piedi"
% cond = cond1*ones(Nx,1)*(xx,yy<=lly1)+cond2*ones(Nx,1)*(xx,yy>lly1);

cond = cond1*ones(Nx,Ny);
cond(Ny1:Ny,:) = cond2;

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);
[xmat,ymat] = meshgrid(xx,yy);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,:) = 0;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

bb = zeros(Ntot,1);

% Ovest

% jj = 1;
for ii = 2:Nx-1

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-haria/dy/cond1;
    AA(kk,kk+1) = 1/2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -Taria*haria/dy/cond1;

end

% Est

% jj = Ny;
for ii = 2:Nx-1

    kk = Ntot-Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-haria/dy/cond2;
    AA(kk,kk+1) = 1/2/dx^2;
    bb(kk) = -Taria*haria/dy/cond2;

end

% % Nord
% 
% ii = 1;
% for jj = 2:Ny1-1
% 
%     kk = (jj-1)*Nx+ii;
%     AA(kk,:) = 0;
%     AA(kk,kk-Nx) = 1/2/dy^2;
%     AA(kk,kk) = -1/dx^2-1/dy^2-halto/dx/cond1;
%     AA(kk,kk+1) = 1/dx^2;
%     AA(kk,kk+Nx) = 1/2/dy^2;
%     bb(kk) = Taria*halto/dx/cond1;
% 
% end
% 
% for jj = Ny1:Ny-1
% 
%     kk = (jj-1)*Nx+ii;
%     AA(kk,:) = 0;
%     AA(kk,kk-Nx) = 1/2/dy^2;
%     AA(kk,kk) = -1/dx^2-1/dy^2-halto/dx/cond2;
%     AA(kk,kk+1) = 1/dx^2;
%     AA(kk,kk+Nx) = 1/2/dy^2;
%     bb(kk) = Taria*halto/dx/cond2;
% 
% end

% % Sud
% 
% ii = Nx;
% for jj = 2:Ny1-1
% 
%     kk = (jj-1)*Nx+ii;
%     AA(kk,:) = 0;
%     AA(kk,kk-Nx) = 1/2/dy^2;
%     AA(kk,kk-1) = 1/dx^2;
%     AA(kk,kk) = -1/dx^2-1/dy^2-halto/dx/cond1;
%     AA(kk,kk+Nx) = 1/2/dy^2;
%     bb(kk) = Taria*halto/dx/cond1;
% 
% end
% 
% for jj = Ny1+1:Ny-1
% 
%     kk = (jj-1)*Nx+ii;
%     AA(kk,:) = 0;
%     AA(kk,kk-Nx) = 1/2/dy^2;
%     AA(kk,kk-1) = 1/dx^2;
%     AA(kk,kk) = -1/dx^2-1/dy^2-halto/dx/cond2;
%     AA(kk,kk+Nx) = 1/2/dy^2;
%     bb(kk) = Taria*halto/dx/cond2;
% 
% end



kk = 1:Nx:Ntot-Nx+1;
AA(1,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tnord;

kk = Nx:Nx:Ntot;
AA(end,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tsud;

% interfaccia

for kk = Ntot1-Nx+2:Ntot1-1

AA(kk,kk-Nx) = cond1/dy^2;
AA(kk,kk-1) = condave/dx^2;
AA(kk,kk) = -cond1/dy^2-2*condave/dx^2-cond2/dy^2;
AA(kk,kk+1) = condave/dx^2;
AA(kk,kk+Nx) = cond2/dy^2;
bb(kk) = 0;

end


% AA(1, :) = 0;
% AA(1, 1) = 1;
% bb(1) = Tnord;
% 
% AA(1,Ny) = 1;
% bb(Ntot-Nx+1) = Tnord;
% 
% AA(Nx, :) = 0;
% AA(Nx, 1) = 1;
% bb(Nx) = Tsud;
% 
% AA(Nx, Ny) = 1;
% bb(Ntot) = Tsud;

TT = AA\bb;

TTT = reshape(TT,Nx,Ny);

figure(1)
surf(xmat',ymat',TTT)




xx1 = (0:dx:llx/2)';
Nx1 = length(xx1);

% caso vettori "sdraiati"

cond = cond1*ones(Nx,Ny);
cond(Nx1:Nx,:) = cond2;

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);
[xmat,ymat] = meshgrid(xx,yy);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,:) = 0;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

bb = zeros(Ntot,1);

% Ovest

jj = 1;
for ii = 2:Nx-1

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-haria/dy/cond(ii,jj);
    AA(kk,kk+1) = 1/2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;
    bb(kk) = -Taria*haria/dy/cond(ii,jj);

end

% Est

jj = Ny;
for ii = 2:Nx-1

    kk = Ntot-Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-haria/dy/cond(ii,jj);
    AA(kk,kk+1) = 1/2/dx^2;
    bb(kk) = -Taria*haria/dy/cond(ii,jj);

end

kk = 1:Nx:Ntot-Nx+1;
AA(1,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tnord;

kk = Nx:Nx:Ntot;
AA(end,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tsud;

% interfaccia

for kk = Nx1+Nx:Nx:Ntot-2*Nx+Nx1

AA(kk,kk-Nx) = condave/dy^2;
AA(kk,kk-1) = cond1/dx^2;
AA(kk,kk) = -cond1/dx^2-2*condave/dy^2-cond2/dx^2;
AA(kk,kk+1) = cond2/dx^2;
AA(kk,kk+Nx) = condave/dy^2;
bb(kk) = 0;

end


% AA(1, :) = 0;
% AA(1, 1) = 1;
% bb(1) = Tnord;
% 
% AA(1,Ny) = 1;
% bb(Ntot-Nx+1) = Tnord;
% 
% AA(Nx, :) = 0;
% AA(Nx, 1) = 1;
% bb(Nx) = Tsud;
% 
% AA(Nx, Ny) = 1;
% bb(Ntot) = Tsud;

TT = AA\bb;

TTT = reshape(TT,Nx,Ny);

figure(2)
surf(xmat',ymat',TTT)


