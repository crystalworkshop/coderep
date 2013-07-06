 global n; if isempty(n), n=512; end;
 global t;  t=0;    % start time
 global T;  T=1178;    % start temperature
 global Nfile;  Nfile=1; % frequency of graphs
 global Kd; Kd=1d-1; % for trace element
 global Xs, Xs=0.1d0; %initial crystal position
 global V; V=0d0;  %initial growth rate
 
 global Clen; %dimentional length of the domain (m)
 global D0; D0=1d-12;  %dimentional diffusion coefficient.
 global Tscale;  %time scale
 
 global Tfin;   % final time;
 global dt0;   % time step; 

 
 global Ds; Ds=1d-2;  %itrace diffusion coefficient in the crystal
 global Dm; Dm=1d-1;  %itrace diffusion coefficient in the melt
global FileName name ext
 global Ox 
 global VT Tmin

 FileName = uigetfile('*.dat','Input Composition from:','input.dat');
[pathstr, name, ext] = fileparts(FileName);

Ox=importdata(FileName);

INparam=importdata(sprintf('%s_param%s',name,ext));
Tfin_d=INparam.data(1)*3600;
dt0=Tfin*1d-4;
VT=INparam.data(2);
Tmin=INparam.data(3);
dt0=INparam.data(4);
n=(INparam.data(5));
Clen=INparam.data(6);

Tscale=Clen^2/D0;

VT=VT/3600*Tscale;
Tfin=Tfin_d/Tscale;

%FileName = 'res1.dat' ; 
%FileName = uiputfile('*.dat','Save results to:','res1.dat');




 fid = fopen(strcat(name,'_Majore.dat'), 'w+');
fid1 = fopen(strcat(name,'_t.dat'), 'w+');
fid2 = fopen(strcat(name,'_trace.dat'), 'w+');
if fid == -1
    error('File is not opened');
end
if fid1 == -1
    error('File is not opened');
end

