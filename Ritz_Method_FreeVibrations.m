%RITZ METHOD FREEE VIBRATIONS

clear;clc; 

%Plate (dimensions and elastic properties)

a=1000;              %Dimensions a,b,t in millimeters
b=1000;
E=72000;
nu=0.3;
t=10;
J0=2.7*1e-9;
rho=2696;         %Density in kg/m^3
I0=rho*t/1e3;   %Mass/Area, Areal density
bcs='CCCC';

d=E.*t.^3./(12*(1-nu.^2));
D=[d,d*nu,0;d*nu,d,0;0,0,d*(1-nu)./2];

%Ritz (degree of freedom number and definition of c)

R=5;  %R number of trial functions in x and y

C=zeros(R*R,1);

%Trial definitions

%You can make different combinations changing x-y direction depending on
%what constraints you have. You should put into the code also the
%derivatives.

%You need define also a,i: length and index.

%I've already put b where it should be

switch  bcs
    case 'SSSS'

        phi=@(x,i) sin(i*pi*x/a);
        dphi=@(x,i) (i*pi/a)*cos(i*pi*x/a);
        ddphi=@(x,i) -(i*pi/a)^2*sin(i*pi*x/a);

        psi=@(x,i) sin(i*pi*x/b);
        dpsi=@(x,i) (i*pi/b)*cos(i*pi*x/b);
        ddpsi=@(x,i) -(i*pi/b)^2*sin(i*pi*x/b);

    case 'CCCC'
        phi=@(x,i) (x./a).^(i+3)-2.*(x./a).^(i+2)+(x./a).^(i+1);
        dphi=@(x,i) ((i+3)./a^(i+3)).*x.^(i+2)-2.*((i+2)/a^(i+2)).*x.^(i+1)+((i+1)/a.^(i+1)).*x.^(i);
        ddphi=@(x,i) ((i+3).*(i+2)./a^(i+3)).*x.^(i+1)-2.*((i+2).*(i+1)./a^(i+2)).*x.^(i)+...
            ((i+1).*i./a^(i+1)).*x.^(i-1);

        psi=@(x,i) (x./b).^(i+3)-2.*(x./b).^(i+2)+(x./b).^(i+1);
        dpsi=@(x,i) ((i+3)./b^(i+3)).*x.^(i+2)-2.*((i+2)/b^(i+2)).*x.^(i+1)+((i+1)/b.^(i+1)).*x.^(i);
        ddpsi=@(x,i) ((i+3).*(i+2)./b^(i+3)).*x.^(i+1)-2.*((i+2).*(i+1)./b^(i+2)).*x.^(i)+...
            ((i+1).*i./b^(i+1)).*x.^(i-1);
        
    case 'CFCF'

        phi=@(x,i) (a^(-i-1)).*x.^(i+1);
        dphi=@(x,i) (a^(-i-1)).*(i+1).*x.^i;
        ddphi=@(x,i) (a^(-i-1)).*i.*(i+1).*x.^(i-1);

        psi=@(x,i) (b^(-i-1)).*x.^(i+1);
        dpsi=@(x,i) (b^(-i-1)).*(i+1).*x.^i;
        ddpsi=@(x,i) (b^(-i-1)).*i.*(i+1).*x.^(i-1);


    case 'CFCC'

        phi=@(x,i) (a^(-i-1)).*x.^(i+1);
        dphi=@(x,i) (a^(-i-1)).*(i+1).*x.^i;
        ddphi=@(x,i) (a^(-i-1)).*i.*(i+1).*x.^(i-1);

        psi=@(x,i) (x./b).^(i+3)-2.*(x./b).^(i+2)+(x./b).^(i+1);
        dpsi=@(x,i) ((i+3)./b^(i+3)).*x.^(i+2)-2.*((i+2)/b^(i+2)).*x.^(i+1)+((i+1)/b.^(i+1)).*x.^(i);
        ddpsi=@(x,i) ((i+3).*(i+2)./b^(i+3)).*x.^(i+1)-2.*((i+2).*(i+1)./b^(i+2)).*x.^(i)+...
            ((i+1).*i./b^(i+1)).*x.^(i-1);
        
        %Ti puoi divertire a mischiare essential BCs
end

%IRITZ (definitions on Idirdir_kk)

%These values change only in dimensions with R,S

Ixx_00=zeros(R,R);
Ixx_02=zeros(R,R);
Ixx_11=zeros(R,R);
Ixx_20=zeros(R,R);
Ixx_22=zeros(R,R);
Iyy_00=zeros(R,R);
Iyy_02=zeros(R,R);
Iyy_11=zeros(R,R);
Iyy_20=zeros(R,R);
Iyy_22=zeros(R,R);
Jx=zeros(R,1);
Jy=zeros(R,1);     

%Evaluate integrals

for i=1:R
    for j=1:R
        Ixx_00(i,j)=integral(@(x) phi(x,i).*phi(x,j),0,a);
        Ixx_02(i,j)=integral(@(x) phi(x,i).*ddphi(x,j),0,a);
        Ixx_11(i,j)=integral(@(x) dphi(x,i).*dphi(x,j),0,a);
        Ixx_20(i,j)=integral(@(x) ddphi(x,i).*phi(x,j),0,a);
        Ixx_22(i,j)=integral(@(x) ddphi(x,i).*ddphi(x,j),0,a);
        Iyy_00(i,j)=integral(@(x) psi(x,i).*psi(x,j),0,b);
        Iyy_02(i,j)=integral(@(x) psi(x,i).*ddpsi(x,j),0,b);
        Iyy_11(i,j)=integral(@(x) dpsi(x,i).*dpsi(x,j),0,b);
        Iyy_20(i,j)=integral(@(x) ddpsi(x,i).*psi(x,j),0,b);
        Iyy_22(i,j)=integral(@(x) ddpsi(x,i).*ddpsi(x,j),0,b);
    end
end

%Build stiffness ans mass matrices

D11=D(1,1);
D12=D(1,2);
D22=D(2,2);
D66=D(3,3);

K=D11*kron(Ixx_22,Iyy_00)+D12*(kron(Ixx_20,Iyy_02)+kron(Ixx_02,Iyy_20))+...
    D22*kron(Ixx_00,Iyy_22)+...
    4*D66*kron(Ixx_11,Iyy_11);

M=I0*kron(Ixx_00,Iyy_00);

%Eigenvalue/Eigenvectors problem

[D,V]=eig(K,M);    %ATT: D EIGENVECTORS AND V EIGENVALUE

%Correction of V and D with only correct angular frequencies (real ones)

V=sqrt(V);
coeffcorr=a^2*sqrt(I0/d);
V=V.*coeffcorr;

%Plot of the solutions

vecx=linspace(0,a,100);
vecy=linspace(0,b,100);

[X,Y]=meshgrid(vecx,vecy);

%Modal shapes

w0=zeros(length(b),length(a));
w=0;
ind=1;


for h=1:R^2
    C=D(:,h);
for k=1:length(vecx)
    for z=1:length(vecy)
       posx=vecx(k);
       posy=vecy(z);

       for i=1:R
           for j=1:R
               w=w+C(ind).*phi(posx,i).*psi(posy,j);
               ind=ind+1;
           end
       end

       w0(z,k)=w;
       w=0;
       ind=1;
    end
end
omega=V(h,h);
figure(h);
surf(X,Y,w0);
title(sprintf('Modal shape of w= %d rad/s as natural angular frequency', round(omega)));
zlabel('[mm]');

end








