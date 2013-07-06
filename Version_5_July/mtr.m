%%%--- module Mtr;
global n
global k; if isempty(k), k=5; end;
 global x; if isempty(x), x=zeros(1,n); end;
 global xm; if isempty(xm), xm=zeros(1,n); end;
 global xs; if isempty(xs), xs=zeros(1,n); end;
 global y; if isempty(y), y=zeros(k,n); end;
 global y0; if isempty(y0), y0=zeros(k,n); end;
 global dy; if isempty(dy), dy=zeros(k,n); end;
 global cx; if isempty(cx), cx=zeros(k,k,n); end;
 global ax; if isempty(ax), ax=zeros(k,k,n); end;
 global bx; if isempty(bx), bx=zeros(k,k,n); end;
 global Psix; if isempty(Psix), Psix=zeros(k,k,n); end;
 global tetax; if isempty(tetax), tetax=zeros(k,n); end;
 global Cs; if isempty(Cs), Cs=zeros(n); end;
 global Cm; if isempty(Cm), Cs=zeros(n); end;
 
 global Dii; if isempty(Dii), Dii=zeros(k+1); end;
 global charge; if isempty(charge), charge=zeros(k+1); end;
 
 global name_variables
 global tol1; if isempty(tol1), tol1=0; end;
 global dt dt0 t Nw
 
 
 %this will be from input file (SiO2; TiO2; Al2O3; Fe2O3; FeO; MgO; CaO; Na2O; K2O)
 %Ox = [49.3 1.4 13.5 0 10.4 8.2 12.4 2.5 0.3];
 Mw = [60.08 79.95 50.98 79.845 71.85 40.31 56.08 30.99 47.1];
 Nw = Ox.data./Mw;
 sNw = sum(Nw);
 
 Nw = Nw / sNw;

 % y are cation frac: Si, Fe, Mg, Ca, Na+K 
y(1,:) = Nw(1);
y(2,:) = Nw(4) + Nw(5);
y(3,:) = Nw(6);
y(4,:) = Nw(7);
y(5,:) = Nw(8) + Nw(9);

 
 %name of the variable, sapce, time and k*component
 name_variables = {'Space' ; 'Time' ; 'SiO_2' ; 'Fe^2 + Fe^3' ; 'MgO' ; 'CaO'...
      ; 'Na + K' ; 'Al + Ti'}; 
 dt=dt0;
 charge(:)=1;
 
f = figure('Units','normalized','Visible','off','Position',[.10,.10,.8,.8]);
hsurf = uicontrol('Units','normalized','Style','pushbutton','String','Stop',...
    'Position',[0.925,0.9,0.05,0.025],...
    'Callback',{@surfbutton_Callback});
a0 = uicontrol(f,'Units','normalized','Style','text',...
    'Position',[0.925,0.85,0.05,0.025],...
    'String','NGraph');
a1 = uicontrol(f,'Units','normalized', 'Style','edit',...
    'Position',[0.925,0.8,0.05,.025],...
    'CallBack',@callb1,...
    'String',Nfile);
a0 = uicontrol(f,'Units','normalized', 'Style','text',...
    'Position',[0.925,0.75,0.05,0.025],...
    'String','dt,t/Tfin');
a2 = uicontrol(f,'Units','normalized', 'Style','edit',...
    'Position',[0.925,0.7,0.05,0.025],...
    'CallBack',@callb2,...
    'String',dt);
a3 = uicontrol(f, 'Units','normalized','Style','edit',...
    'Position',[0.925,0.65,0.05,0.025],...
    'CallBack',@callb2,...
    'String',t);

ha = axes('Units','normalized','Position',[0.05,.38,.40,0.57]);
xlim([0 1])
ha1 = axes('Units','normalized','Position',[0.05,0.05,0.40,0.28]);
xlim([0 1])
ha2 = axes('Units','normalized','Position',[0.5,0.65,0.40,0.30]);
ha3 = axes('Units','normalized','Position',[0.5,0.45,0.40,0.15]);
ha4 = axes('Units','normalized','Position',[0.5,0.25,0.40,0.15]);
ha5 = axes('Units','normalized','Position',[0.5,0.05,0.40,0.15]);

xlabel(ha,'Normalised Distance')
ylabel(ha,'Major Elements')
xlabel(ha1,'Normalised Distance')
ylabel(ha1,'Trace Elements')
xlabel(ha2,'Time (hours)')
ylabel(ha2,'Crystal Composition')
xlabel(ha3,'Time (hours)')
ylabel(ha3,'Growth Rate x10^8, m.s^{-1}')
xlabel(ha4,'Time (hours)')
ylabel(ha4,'Liquidus T^o, ^oC')
xlabel(ha5,'Time (hours)')
ylabel(ha5,'Undercooling, ^oC')


set(f,'CurrentAxes',ha)
set(f,'Name','DIFFUSOR')
movegui(f,'center')
set(f,'Visible','on');

hold on
