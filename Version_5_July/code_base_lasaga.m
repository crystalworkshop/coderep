function code_base_lasaga(varargin)
clear global; clear functions;
param;
mtr;
global t y x Xs V xm is_stop Nfile dt  Cm Cs xs dt0 Tfin  ...
     name_variables T Fo DT Tliq Tscale Clen VT Tmin n

is_stop=true;
InitCond;  %initial conditions
V=Gr_rate(y(:,1));
T00=Tliq;
iru=0;
dt=dt0;
%
dt=1d-6;
while(t<=Tfin && is_stop); %main loop in
    T=max(Tmin,T00- VT*t); 
    OneStep;
    [xx]=Cxx(y(:,1));
    %%%%%
    set(f,'CurrentAxes',ha2)
    hold on
    plot(t*Tscale/3600,Fo/100,'.b', 'MarkerSize', 10)
    plot(t*Tscale/3600,xx(2),'.g', 'MarkerSize', 10)
    plot(t*Tscale/3600,xx(3),'.r','MarkerSize', 10)
    legend('Fo/100', 'Fe', 'Mg',  'Location','NorthEast');

    %%%%%%%
    %%%%%
    set(f,'CurrentAxes',ha3)
    hold on
    plot(t*Tscale/3600,V*Clen/Tscale*1d8,'.b','MarkerSize', 10)
set(f,'CurrentAxes',ha4)
        hold on
    plot(t,Tliq,'.b','MarkerSize', 10)
set(f,'CurrentAxes',ha5)
        hold on
    plot(t*Tscale/3600,DT,'.b' ,'MarkerSize', 10)
    %%%%%%%
    Diff_trace;
    Xs=Xs+V*dt;
    Tmp=Cxx(y(:,1));
    fprintf(fid1, '%6d %6d %6d %6d %6d\r\n', t,V,Xs,Tmp(1),Tmp(2)); 
    xm=Xs+x*(1-Xs);    
    xs=Xs*x;    
    dt=dt0;
    t=t+dt
    set(a2, 'String',dt)
    set(a3, 'String',t*100/Tfin)
    if(rem(iru,Nfile)==0);
       set(f,'CurrentAxes',ha)
       plot(ha,xm(1,:),y(1,:),xm(1,:),y(2,:),xm(1,:),y(3,:),xm(1,:),y(4,:),xm(1,:),y(5,:),xm(1,:),1-sum(y)); 
       legend(name_variables{3:end}, 'Location','NorthEast');
       hold on
       set(f,'CurrentAxes',ha1)
       plot(ha1,xs(1,1:n),Cs(1:n),xm(1,1:n),fliplr(Cm(1:n)));
       xlabel(ha1,'Normalised Distance')
       ylabel(ha1,'Trace Elements')
       hold on
       for i=1:n
            fprintf(fid, ['%6d %6d ' repmat('%6d ',1,k) '\r\n'], xm(1,i),t,y(:,i)); 
            fprintf(fid2, '%6d %6d %6d\r\n',t,xs(1,i),Cs(i));   
       end;
       for i=1:n
            fprintf(fid2, '%6d  %6d %6d\r\n',t,xm(i),Cm(n-i+1));
       end;
       fprintf(fid, '\r\n');
       fprintf(fid2, '\r\n');
    end;
    iru=iru+1;
    pause(0.001)
end;
fclose(fid);
fclose(fid1);
%plot 3D view of the results for each component
%load result file fid -> space, time, k * components
% n is number of space point
% k is the number of component -> k figures
result=importdata(sprintf('%s_Majore%s',name,ext));
% we don't know the total numer of time step
% nTime correspond to the number of time step
nTime=size(result,1)/n;
matTime = reshape(result(:,2),n,nTime);
matSpace = reshape(result(:,1),n,nTime);
for iC=1:k
    %2+iC to take into account time and space
   matC{iC} = reshape(result(:,2+iC),n,nTime);
   figure(100+iC)
   clf
   surf(matSpace,matTime,matC{iC})
   shading interp
   xlabel(name_variables{1})
   ylabel(name_variables{2})
   zlabel(name_variables{2+iC})
   view([20 20])
end

hold off   
end %program main


function [s]=RHS(z1,z2,z3,z0,ix)

global dt x Xs V 

h12=x(ix)-x(ix-1);
h23=x(ix+1)-x(ix);

Di12=D_MTR((z1+z2)/2);
Di23=D_MTR((z3+z2)/2);
dth=dt/h23;

J=dth*(Di23*(z3-z2)/h23-Di12*(z2-z1)/h12);

Adv=dt*V*(1-x(ix))/(1-Xs)*(z2-z1)/h23;

s =z2-z0+Adv-J/(1-Xs)/(1-Xs);

end %subroutine RHS

function Di=D_MTR(Z) %Diffusion matrix
global k Dii charge T D0
% Si, Mg, Fe, Ca, Na-K, TI+Al
Do  =  [-12.79 -10.71 -4.76 -11.07 -9.32 -10.08];              
Ea  =  [20119 21057 31936 20684 19729 24600];
Dii=exp(Do).*exp(-Ea/(T+273))/D0;

zn=0d0;      
C00=1;
for i=1:k;
    zn=zn+charge(i)*max(0,Z(i))*Dii(i);
    C00=C00-max(0,Z(i));
end;
zn=zn+Dii(k+1)*charge(k+1)^2*C00;
for i=1:k;
    for j=1:k;
        krd=0;
        if(i==j) krd=1; end;
        Di(i,j)=Dii(i)*krd-Dii(i)*charge(i)*charge(j)*max(0,Z(i))*(Dii(j)-Dii(k+1))/zn;
    end;
end;
if(zn~=zn) 
    Di
    error('zero mtr'); 
end;
end

function Grr=Gr_rate(Z)%Crystal growth rate
global Nw Tscale Clen T DT Tliq
AlTi=Nw(3)/(Nw(2)+Nw(3));
Fe3Fe=Nw(4)/(Nw(4)+Nw(5));
NaK=Nw(8)/(Nw(8)+Nw(9));
Al=1-sum(Z);
a=[1180.13	-330.83	-285.45 ...
    1938.42	-403.98	1439.62...	
    -538.63	-1165.35	-1143.04 ...
    35.41	1265.78];

%a=[1022  88 -953 149 62 1208 345 -872 -274 50 1346];
rec=[Z(1) (1-AlTi)*Al AlTi*Al  Fe3Fe*Z(2) (1-Fe3Fe)*Z(2) ...
    Z(3) Z(4) NaK*Z(5) (1-NaK)*Z(5)];
srec=sum(rec);
rec=rec/srec;
st=((rec(6) + rec(5))^2*rec(1))^(1./3.);
%T=a(1)+a(2)*Al+a(3)*Ti+a(4)*Fe3 +a(5)*Fe2+a(6)*Mg+a(7)*Ca+a(8)*Na+a(9)*K+a(10)*log(st)+a(11)*sqrt*(Al*(Na+K));
%(SiO2; TiO2; Al2O3; Fe2O3; FeO; MgO; CaO; Na2O; K2O)
Tliq=a(1)+a(2)*rec(3)+a(3)*rec(2)+a(4)*rec(4) +a(5)*rec(5)...
    +a(6)*rec(6)+a(7)*rec(7)+a(8)*rec(8)+a(9)*rec(9)...
    +a(10)*log(st)+a(11)*sqrt(abs(rec(3)*(rec(8)+rec(9)))) ;

% Growth rate polynomial fit valid for DT up to 35 degrees, not higher!
% data are from Jambon 1992 melt inclusion olivine growth experiments
% and Schiano melt migration experiments (very low DT)
DT= min(35,Tliq-T);
Grr=(4.4d-4*DT^3-8.91d-3*DT^2+1.21d-1*DT)*1d-8*Tscale/Clen ;  %DT is in deg, V units are *10E-8 m/s
if(DT<0) Grr=0; end;

%Grr=0.1; %50*(y(3,1)-0.2)^3;
end

function InitCond % initial conditios

global x  n  Xs  Cs Cm Kd  xs
h=1d0/n;
hs=Xs/n;
xs(1)=0d0;
for i=1:n;
x(i)=1.-cos(pi*double(i-1)/2/(n-1));
end;
for i=1:n;
    Cs(i)=1d-5;
end;

for i=1:n;
    Cm(i)=Cs(n)/Kd;
end;

end %subroutine InitCond

function bound %boundary conditions - now fixed concentrations. it'll be updated later
global n cx fx bx ax y k
persistent  z2 z3
ax(:,:,1)=0d0;
ax(:,:,n)=0d0;
bx(:,:,1)=0d0;
bx(:,:,n)=0d0;
cx(:,:,1)=0d0;
cx(:,:,n)=0d0;

% if isempty(z2), z2=zeros(k); end;
% if isempty(z3), z3=zeros(k); end;

z2=y(:,1);
z3=y(:,2);
    dz=1.0d-7;
    
    [s]=RHSB(z2,z3);
    fx(:,1)=s;
    
    for j = 1:k;
        zt=z2;
        zt(j)=z2(j)+dz;
        [s]=RHSB(zt,z3);
        cx(:,j,1)=-(s(:)-fx(:,1))./dz;
    end;
    for j = 1:k;
        zt=z3;
        zt(j)=z3(j)+dz;
        [s]=RHSB(z2,zt);
        bx(:,j,1)=(s(:)-fx(:,1))./dz;
    end;

fx(:,n)=0d0;   %fixed values
for i=1:k;
    cx(i,i,n)=1d0;
    ax(i,i,n)=1d0;
end;
end %subroutine bound

function [s]=RHSB(z2,z3)

global x Xs V
persistent Cxtl
V=Gr_rate(z2);


Di23=D_MTR((z2+z3)/2);

J=-Di23*(z3-z2)/(x(2)-x(1));
Cxtl=Cxx(z2);
s = V*(z2-Cxtl')-J/(1-Xs);

end %subroutine RHS

function [xx]=Cxx(z2)
global Nw Fo
Kd1=1/0.31;
Fe3Fe=Nw(4)/(Nw(4)+Nw(5));
XmgFe=Kd1*z2(3)/z2(2)/(1-Fe3Fe);
xx(1)=1d0/3d0;
xx(2)=2d0/3d0/(1+XmgFe);
xx(3)=2d0/3d0-xx(2);
xx(4)=0d0;
xx(5)=0d0;
Fo=100*xx(3)/(xx(3)+xx(2));
end

function OneStep    %1 step in time
global tol1 y y0

y0=max(0,y);
toler=1d-8;
tol1=1d10;
it=0;
while(tol1>toler && it<=500);
    OneIter;
    it=it+1; 
end;
end %subroutine OneStep


function OneIter %One itteration
global n tol1 dy k y fx bx cx ax y0  
persistent  s z0 z1 z2 z3 zt

if isempty(z2), z2=zeros(k); end;
if isempty(z1), z1=zeros(k); end;
if isempty(z3), z3=zeros(k); end;
if isempty(z0), z0=zeros(k); end;
if isempty(s),  s =zeros(k); end;
if isempty(zt), zt=zeros(k); end;


tol1=0.0d0;
for ix = 2:n-1;
    i1=ix-1;
    i3=ix+1;
    
    z1=y(:,i1);
    z3=y(:,i3);
    z2=y(:,ix);
    z0=y0(:,ix);
    dz=1.0d-7;
    
    [s]=RHS(z1,z2,z3,z0,ix);
    fx(:,ix)=s;
    for j = 1:k;
        zt=z1;
        zt(j)=z1(j)+dz;
        [s]=RHS(zt,z2,z3,z0,ix);
        ax(:,j,ix)=(s(:)-fx(:,ix))./dz;
    end;
    
    for j = 1:k;
        zt=z2;
        zt(j)=z2(j)+dz;
        [s]=RHS(z1,zt,z3,z0,ix);
        cx(:,j,ix)=-(s(:)-fx(:,ix))./dz;
    end;
    for j = 1:k;
        zt=z3;
        zt(j)=z3(j)+dz;
        [s]=RHS(z1,z2,zt,z0,ix);
        bx(:,j,ix)=(s(:)-fx(:,ix))./dz;
    end;
    %a, B, C, F
    
end;

%boundary conditions
bound;

%Solution of linear equations
progon;

for ix =1 :n;
    y(:,ix)=y(:,ix)+dy(:,ix);
    for l= 1:k;       
        tol1=tol1+dy(l,ix).*dy(l,ix)/max(1d-5,y(l,ix)*y(l,ix));
    end;
end;
tol1=sqrt(tol1)./n;
end %subroutine OneIter




function progon
global tetax dy Psix n ax fx cx bx 
persistent tmp
tmp=inv(cx(:,:,1));

Psix(:,:,1)=tmp*bx(:,:,1);
tetax(:,1)=tmp*fx(:,1);

for k1=1:n-1;
    tmp=inv(cx(:,:,k1)-ax(:,:,k1)*Psix(:,:,k1));
    Psix(:,:,k1+1)=tmp*bx(:,:,k1);
    tetax(:,k1+1)=tmp*(ax(:,:,k1)*tetax(:,k1)+fx(:,k1));
end;
tmp=inv(cx(:,:,n)-ax(:,:,n)*Psix(:,:,n));

dy(:,n)=tmp*(fx(:,n)+ax(:,:,n)*tetax(:,n));
for i=1:n-1;
    j=n-i;
    dy(:,j)=(Psix(:,:,j+1)*dy(:,j+1))+tetax(:,j+1);
end;
end %subroutine progon

function Diff_trace
global dt n x Cm Cs V Kd Xs Ds Dm D0 T
persistent als alm bets betm as bs cs fs  Cs0 Cm0 am bm cm fm  

if isempty(as), as=zeros(n); end;
if isempty(bs), bs=zeros(n); end;
if isempty(cs), cs=zeros(n); end;
if isempty(fs), fs=zeros(n); end;
%if isempty(Cs0), Cs0=zeros(n); end;

if isempty(am), am=zeros(n); end;
if isempty(bm), bm=zeros(n); end;
if isempty(cm), cm=zeros(n); end;
if isempty(fm), fm=zeros(n); end;
%if isempty(Cm0), Cm0=zeros(n); end;

if isempty(als), als=zeros(n); end;
if isempty(bets),bets=zeros(n); end;
if isempty(alm), alm=zeros(n); end;
if isempty(betm),betm=zeros(n); end;
Cs0=Cs;
Cm0=Cm;
Ds=1d-5;
Dm=0.03*exp(-13)*exp(-19194/(T+273))/D0;
Kd = 0.100; % P - experimental data from Brunet & Chazot, 2001; Kd = 0.1;
as(1)=0;
bs(1)=1;
cs(1)=1;
fs(1)=0;

for i=2:n-1
    as(i)=Ds*dt/(x(i)-x(i-1))^2/Xs/Xs;
    bs(i)=Ds*dt/(x(i+1)-x(i))^2/Xs/Xs+V*x(i)*dt/(x(i+1)-x(i))/Xs;
    cs(i)=1+as(i)+bs(i);     
    fs(i)=Cs0(i);
end;

als(1)=bs(1)/cs(1); %in the crystal
bets(1)=fs(1)/cs(1);

for i=1:n-1;
    zn=1.0d0/(cs(i)-as(i)*als(i));
    als(i+1)=bs(i)*zn;
    bets(i+1)=(as(i)*bets(i)+fs(i))*zn;
end;

am(1)=0;
bm(1)=1;
cm(1)=1;
fm(1)=0 ;%Cm0(1);

for i=2:n-1
    Adv1=dt*V*(1-x(i))/(1-Xs)/(x(i)-x(i-1));
    am(i)=Dm*dt/(x(i)-x(i-1))^2/(1-Xs)/(1-Xs)+Adv1;
    bm(i)=Dm*dt/(x(i+1)-x(i))^2/(1-Xs)/(1-Xs);
    cm(i)=1+am(i)+bm(i);
    fm(i)=Cm0(i);
end;


alm(1)=bm(1)/cm(1); % in the melt
betm(1)=fm(1)/cm(1);
for i=1:n-1;
    zn=1.0d0/(cm(i)-am(i)*alm(i));
    alm(i+1)=bm(i)*zn;
    betm(i+1)=(am(i)*betm(i)+fm(i))*zn;
end;
hs=x(n)-x(n-1);
hm=-(x(n)-x(n-1));
t4 = (Ds * bets(n) * hm);
t12 = (Xs ^ 2);
t15 = (V * hs);
t16 = als(n) - 1;
Cs(n) = Kd * ((Dm * betm(n) * hs - t4) * Xs + t4) / (V * hm * hs * (Kd - 1) * t12 + (((Ds * t16 - t15) * Kd + t15) * hm - Dm * hs * (alm(n) - 1)) * Xs - Ds * Kd * hm * t16);
 
for i=1:n-1;
    j=n-i;
    Cs(j)=bets(j+1)+als(j+1).*Cs(j+1);
end;

Cm(n)=Cs(n)/Kd;

for i=1:n-1;
    j=n-i;
    Cm(j)=betm(j+1)+alm(j+1).*Cm(j+1);
end;


end %subroutine progon
   function surfbutton_Callback(~,~) 
   global is_stop
   is_stop=false;
   end
function [] = callb1(H,~)
global Nfile
Nfile = str2num(get(H,'string'));
disp(Nfile)
end
function [] = callb2(H,~)
global dt0
dt0 = str2double(get(H,'string'));
disp(dt0)
end
