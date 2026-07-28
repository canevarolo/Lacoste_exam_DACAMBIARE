clear all
clear vars
clc

%Dati
Ly=25e-2;
Lx=15e-2;
Lx1=Lx/3;
Lx2=Lx-Lx1;
Ta=20;
ha=20;
Tc=5;
cond=2;

%% Mappa di temperatura in condizioni stazionarie
%Discretizzazione
dx=5e-3;
dy=dx;

x1=(0:dx:Lx1)';
nx1=length(x1);
x2=(0:dx:Lx2)';
nx2=length(x2);

xx=(0:dx:Lx)';
nx=length(xx);
yy=(0:dy:Ly)';
ny=length(yy);

ntot=nx*ny;

[xmat,ymat]=meshgrid(xx,yy);

%Matrice
AA=sparse([],[],[],ntot,ntot,5*ntot);
bb=zeros(ntot,1);

for ii=2:nx-1
    for jj=2:ny-1
        kk=ii+(jj-1)*nx;
        kn=kk-1;
        ks=kk+1;
        ko=kk-nx;
        ke=kk+nx;
        AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
        AA(kk,kn)=cond/dx^2;
        AA(kk,ks)=cond/dx^2;
        AA(kk,ko)=cond/dy^2;
        AA(kk,ke)=cond/dy^2;
    end
end

%Bordi
%Nord condizione di neumann adiabatica
% kk=1:nx:ntot;
% AA(kk,:)=0;
ii=1;
for jj=2:ny-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,ks)=2*cond/dx^2;
    AA(kk,ke)=cond/dy^2;
    AA(kk,ko)=cond/dy^2;
    bb(kk)=0;
end

%Sud condizione di Robin
ii=nx;
for jj=2:ny-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2)-ha/dx;
    AA(kk,kn)=2*cond/dx^2;
    AA(kk,ko)=cond/dy^2;
    AA(kk,ke)=cond/dy^2;
    bb(kk)=-ha/dx*Ta;
end

%Ovest Neumann adiabatico
jj=1;
for ii=2:nx-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,ks)=cond/dx^2;
    AA(kk,kn)=cond/dx^2;
    AA(kk,ke)=2*cond/dy^2;
    bb(kk)=0;
end

%Est divisa tra adiabatica e temperatura imposta
jj=ny;
for ii=2:nx1-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,kn)=cond/dx^2;
    AA(kk,ks)=cond/dx^2;
    AA(kk,ko)=2*cond/dy^2;
    bb(kk)=0;
end
for ii=nx1:nx2
    kk=ii+(jj-1)*nx;
    AA(kk,kk)=1;
    bb(kk)=Tc;
end
for ii=nx2+1:nx-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,kn)=cond/dx^2;
    AA(kk,ks)=cond/dx^2;
    AA(kk,ko)=2*cond/dy^2;
    bb(kk)=0;
end

%Nodi
kk=nx;
AA(kk,:)=0;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2)-ha/dx;
AA(kk,kk-1)=cond/dx^2;
AA(kk,kk+nx)=cond/dy^2;
bb(kk)=-ha/dx*Ta;
kk=ntot;
AA(kk,:)=0;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2)-ha/dx;
AA(kk,kk-1)=cond/dx^2;
AA(kk,kk-nx)=cond/dy^2;
bb(kk)=-ha/dx*Ta;
kk=1;
AA(kk,:)=0;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2);
AA(kk,kk+1)=cond/dx^2;
AA(kk,kk+nx)=cond/dy^2;
bb(kk)=0;
kk=ntot-nx+1;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2);
AA(kk,kk+1)=cond/dx^2;
AA(kk,kk-nx)=cond/dy^2;
bb(kk)=0;

TT=AA\bb;

Tplot=reshape(TT,nx,ny);
surf(xmat',ymat',Tplot)

%% Calcolo della potenza parassita
deltay=[dy/2 dy*ones(1,ny-2) dy/2];
potenza=(cond*(Tplot(:,end)-Tc))/dx;
qq=sum(potenza.*deltay);

