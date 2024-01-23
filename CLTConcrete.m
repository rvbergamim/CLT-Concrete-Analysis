% CLT-CONCRETE DESIGN
%
% MSc. Eng. Ramon Vilela
% Version 1.0
% Updated December 27, 2023
% All rigth reserved

clc; clearvars; close all; tic

%% INPUT DATA
b   = 1.00;   % width, m
L   = 7.00;   % lenght, m
g   = 10;     % gravity, m/s^2
bt  = 0.14;   % lamella width, m  

% Tickness (m)
hc  = 0.08;     % concrete layer
t   = 0.00;     % spacing between materials
h1  = 0.04;     % CLT layer 1
h12 = 0.04;     % CLt layer 12
h2  = 0.04;     % CLT layer 2
h23 = 0.04;     % CLT layer 23
h3  = 0.04;     % CLT layer 3

% Concrete properties
fck = 40e6;       % characteristic strength (Pa)
alphaE = 1.2;     % factor according to NBR 6118:2023 (-)  

% Longitudinal elastic modulus of wood (Pa)
E(1) = 11000e6;    % layer 1
E(2) = 11000e6;    % layer 2
E(3) = 11000e6;    % layer 3

% Transversal elastic modulus (Pa)
G12     = 68e6;     % layer 12
G23     = 68e6;     % layer 23

% Characteristic Strength of wood (Pa)
fbk     = 24e6;     % in bending
ft0k    = 14e6;     % tension parallel
ft90k   = 0.4e6;    % tension perpendicular
fc0k    = 21e6;     % compression parallel
fc90k   = 2.5e6;    % compression perpendicular
fvk     = 4e6;      % shear
fvrk    = 1.6e6;    % rolling shear

% Mean density (kg/m^3)
rhoc    = 2400;     % concrete
rhow    = 550;      % wood

% Coefficient of creep and redution factor
phic    = 3;       % concrete
phiw    = 0.8;     % wood
psi2    = 0.4;     % redution factor according to NBR 8681:2004

% Connectors
Kser    = 49.5e6;   % service slip modulus (N/m)
Ku      = 2/3*Kser; % ultimate slip modulus (N/m)
FRks    = 25.23e3;  % characteristic strenght (N)
smin    = 0.50;     % minimum spacing (m)
smax    = 1.00;     % maximum spacing (m)
nx      = 3;        % number of connectors in x-direction of mid plate
ny      = 5;        % nunber of connectors in y-direction in width b

% Loads (N/m^2)
gk      = 1e3;      % permanent
qk      = 3e3;      % live

%% CALCULATIONS
fid = fopen('Results.txt', 'w');    % create a results file
% Sizing
hclt  = h1+h12+h2+h23+h3;           % CLT (m)
h     = hclt + hc + t;              % CLT-Concrete (m)
nclt  = length(E);                  % quantity of layers
nlam  = floor(b/bt);                % number of lamellae

% Loads (N/m)
g1kclt = b*hclt*rhow*g;             % dead load of CLT
g1kc   = b*hc*rhoc*g;               % dead load of concrete
g2k   = gk*b;                       % permanent
q1k   = qk*b;                       % live
qSd   = 1.25*g1kclt + 1.30*g1kc + 1.35*g2k + 1.5*q1k;

% Design bending moment (N.m)
MSd = qSd*L^2/8;  

% Design shear (N)
VSd = qSd*(L/2-2*h);  

% CLT areas (m^2)
A(1)      = b*h1;                   % Area da camada 1
A(2)      = b*h2;                   % Area da camada 2
A(3)      = b*h3;                   % Area da camada 3

% CLT moment of inertia (m^4)
I(1)      = b*h1^3/12;              % Area da camada 1
I(2)      = b*h2^3/12;              % Area da camada 2
I(3)      = b*h3^3/12;              % Area da camada 3

% Composition coefficient (-)
yclt(1)   = 1/(1+pi^2*E(1)*A(1)*h12/(L^2*b*G12));
yclt(2)   = 1;
yclt(3)   = 1/(1+pi^2*E(3)*A(3)*h23/(L^2*b*G23));

% Position CG's layers in relation of base(m)
zclt(1) = h1/2;
zclt(2) = h1+h12+h2/2;
zclt(3) = h1+h12+h2+h23+h3/2;

yEAz = 0;
yEA = 0;
for i=1:nclt
    yEAz = yclt(i)*E(i)/E(1)*A(i)*zclt(i) + yEAz;
    yEA = yclt(i)*E(i)/E(1)*A(i) + yEA;
end
zCG_CLT = yEAz/yEA;

aclt = zeros(1,nclt);
for i=1:nclt
    aclt(i) = zclt(i)-zCG_CLT;
end

EI(2) = 0;
for i=1:nclt
    EI(2) = E(i)*(I(i)+yclt(i)*A(1)*aclt(i)^2)+ EI(2);
end

% Connectors
sef   = 0.33*smin + 0.67*smax;      % spacing (m)
s = zeros(1,nx);
for i = 1:nx
    s(i) = (smin-smax)/(1-nx)*(i-nx)+smax;
end
x = zeros(1,nx);                    % position of connectors (m)
x(1) = s(1);
for i=2:nx
    x(i) = x(i-1)+s(i);
end
stot = sum(s);
V = qSd*(L/2-x);                    % shear in x-position (N)

Titulo(fid,'GENERAL DATA',70)
fprintf(fid,'- DIMENSIONS\n');
fprintf(fid,'Span (L): \t\t\t\t%4.0f mm\n',L*1000);
fprintf(fid,'Width (b): \t\t\t\t%4.0f mm\n',b*1000);
fprintf(fid,'Concrete tickness (hc): \t\t%4.0f mm\n',hc*1000);
fprintf(fid,'Spacing Concrete/CLT (t): \t\t%4.0f mm\n',t*1000);
fprintf(fid,'CLT tickness (hCLT): \t\t\t%4.0f mm\n',hclt*1000);
fprintf(fid,'Total tickness (h): \t\t\t%4.0f mm\n',h*1000);

fprintf(fid,'\n- LOADS\n');
fprintf(fid,'Dead load CLT (g1kCLT):\t\t\t%2.2f kN/m\n',g1kclt/1000);
fprintf(fid,'Dead load Conc. (g1kc):\t\t\t%2.2f kN/m\n',g1kc/1000);
fprintf(fid,'Permanent (g2k):\t\t\t%2.2f kN/m\n',g2k/1000);
fprintf(fid,'Live (qk):\t\t\t\t%2.2f kN/m\n',qk/1000);

Titulo(fid,'CLT STIFFNESS',70)
fprintf(fid,'Coef. of Composition (gamma1):\t\t%1.3f\n',yclt(1));
fprintf(fid,'Coef. of Composition (gamma2):\t\t%1.3f\n',yclt(2));
fprintf(fid,'Coef. of Composition (gamma3):\t\t%1.3f\n',yclt(3));
fprintf(fid,'Centroide (zCG):\t\t\t%3.0f mm\n',zCG_CLT*1000);
fprintf(fid,'Bending Stiffness (EIclt):\t\t%1.3fe6 Nm^2\n',EI(2)*10^-6);

Titulo(fid,'SPACING OF CONNECTORS',70)
fprintf(fid,'Minimum spacing (smin): \t\t%2.0f mm\n',smin*1000);
fprintf(fid,'Maximum spacing (smin): \t\t%2.0f mm\n',smax*1000);
fprintf(fid,'smax <= 4smin: \t\t\t\t\t');
if smax <= 4*smin
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'smax <= 1000 mm: \t\t\t');
if smax <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'Effective spacing (sef): \t\t%2.0f mm\n',sef*1000);

%% ULTIATE LIMIT STATE (ULS)
% Conector
K = ny*Ku;                          % Stiffness of conector (N/m)
FRds = FRks/1.4;                    % Design strengh of conector (N)

% Concrete
if fck <= 50e6                      % mean tension strength (Pa)
    fctm = 0.3*(fck*10^-6)^(2/3)*10^6;
else
    fctm = 2.12*log(1 + 0.11*fck*10^-6)*10^6;
end
fctkinf = 0.7*fctm;         % inferior characteristic strength (Pa)
fctd = fctkinf/1.4;         % design tension strenght (Pa)
fcd = fck/1.4;              % desing compression strength (Pa)
Eci = alphaE*5600*sqrt(fck*10^-6)*10^6; 
                            % initial Young's modulus (Pa)
alphai = 0.8 + 0.2*fck/80e6;  
                            % calculation factor
Ecs = alphai*Eci;           % secant Young's modulus (Pa)
rho1 = 0;                   % reinforcement rate 
tauRd = 0.25*fctd;          % design shear strength (Pa)
sigmacp = 0;                % normal stress (Pa)
VRdc1 = (tauRd*(1.2+40*rho1)+0.15*sigmacp)*b*hc;
                            % shear force strength (N)
Titulo(fid,'MECHANICAL PROPERTIES',70)
fprintf(fid,'- Tensile strength of concrete:\n');
fprintf(fid,'Medium (fctm): \t\t\t\t%2.2f MPa\n',fctm/10^6);
fprintf(fid,'Inferior (fctk,inf): \t\t\t%2.2f MPa\n',fctkinf/10^6);
fprintf(fid,'Design (fctd): \t\t\t\t%2.2f MPa\n',fctd/10^6);
fprintf(fid,'- Shear strength of concrete:\n');
fprintf(fid,'Design (tauRd): \t\t\t\t%2.2f MPa\n',tauRd/10^6);
fprintf(fid,'\n- Elasticit modulus of concrete:\n');
fprintf(fid,'Initial (Eci): \t\t\t\t%2.2f GPa\n',Eci/10^9);
fprintf(fid,'Secant (Ecs): \t\t\t\t%2.2f GPa\n',Ecs/10^9);

% Madeira
kmod = 0.7*1*0.95;          % modification factor
ksys  = (nlam + 34)/35;     % modification system factor
ft0d = ksys*kmod*ft0k/1.4;  % design parallel tension strengh (Pa)
fbd = ksys*kmod*fbk/1.4;    % design bending strength (Pa)
A0liq = sum(A);             % Net area (m^2)
fvrd = kmod*fvrk/1.8;       % rolling shear strength (Pa)

% CLT-Concrete
hcef=hc;                    % effective thickness (m)
sigma_cinf = fctd;        % stress inderside of concrete (Pa)
it=0;                      % iteration counter
EA(2) = 0;                  % axial stiffness of CLT (N)
    for i=1:nclt
        EA(2) = E(i)*A(i)+EA(2);
    end

% Loop to reduce effective concrete thickness
while sigma_cinf >= fctd
    if sigma_cinf <= fctd || it==0
        hcef=hc;
    else
        hcef=hcef*(1-.005);
    end  
    
    % Axial stiffness (N)
    EA(1) = Ecs*b*hcef;         % Concrete
    EA(2) = 0;                  % CLT
    for i=1:nclt
        EA(2) = E(i)*A(i)+EA(2);
    end
    
    % Composition coefficient (-)
    y(1) = 1/(1+pi^2*EA(1)*sef/(K*L^2)); % Concreto
    y(2) = 1;                            % CLT
    
    % Distance between CG (m)
    at = h-(hcef+hclt)/2;
    a(1) = y(2)*EA(2)/(y(1)*EA(1)+y(2)*EA(2))*at;
    a(2) = y(1)*EA(1)/(y(1)*EA(1)+y(2)*EA(2))*at;
    
    % Effective bending stiffness (N.m^2)
    EI(1) = Ecs*b*hcef^3/12;
    EIef = 0;
    for i=1:2
        EIef = EI(i)+ y(i)*EA(i)*a(i)^2 + EIef;
    end
       
    % Stress (Pa)
    sigma_csup = -MSd*Ecs/EIef*(y(1)*a(1)+0.5*hcef);	% top concrete
    sigma_cinf = -MSd*Ecs/EIef*(y(1)*a(1)-0.5*hcef);	% underside conc.
    sigma_cltsup = MSd*E(3)/EIef*(y(2)*a(2)-0.5*hclt);	% top CLT
    sigma_cltinf = MSd*E(1)/EIef*(y(2)*a(2)+0.5*hclt);	% underside CLT
    figure(1)
    hold on
    plot(it+1,sigma_cinf*10^-6,'o','Color',[0.4 0.4 0.4])
    if it<=8
        fprintf(fid,'In concrete (It. %.0f): \t\t\t%.2f MPa\n'...
            ,it+1,sigma_cinf/10^6);
    else
        fprintf(fid,'In concrete (It. %.0f): \t\t\t%.2f MPa\n'...
            ,it+1,sigma_cinf/10^6);
    end
    it=it+1;      
end
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Iterations')
ylabel('Normal stress (MPa)')
xlim([0 it+2])
title('TENSILE STRESS ON CONCRETE - ULS')
saveas(gcf,'StressConcrete.png')
tef = t+ hc-hcef;

% Discretization along the thickness
z(1) = 0;
z(2) = 0;
z(3) = h1;
z(4) = h1;
z(5) = h1+h12;
z(6) = h1+h12;
z(7) = h1+h12+h2;
z(8) = h1+h12+h2;
z(9) = h1+h12+h2+h23;
z(10) = h1+h12+h2+h23;
z(11) = hclt;
z(12) = hclt;
z(13) = h-hcef;
z(14) = h-hcef;
z(15) = h;
z(16) = h;

% Intermediate stress (Pa)
sigma = zeros(1,16);
for i=2:11
    sigma(i) = (sigma_cltsup-sigma_cltinf)*z(i)/hclt + sigma_cltinf;
end
sigma(1) = 0;
sigma(4) = 0;
sigma(5) = 0;
sigma(8) = 0;
sigma(9) = 0;
sigma(12) = 0;
sigma(13) = 0;
sigma(14) = sigma_cinf;
sigma(15) = sigma_csup;
sigma(16) = 0;

zLNclt = -sigma_cltinf*hclt/(sigma_cltsup-sigma_cltinf);
zLNc = h-hcef-sigma_cinf*hcef/(sigma_csup-sigma_cinf);

% Strain (m/m)
strain = sigma./E(1);
strain(14) = sigma_cinf/Ecs;
strain(15) = sigma_csup/Ecs;

% Figure of stress along the section
figure(2)
hold on
plot(sigma*10^-6,z*10^3)
plot([min(sigma) max(sigma)]*10^-6,[zLNclt zLNclt]*10^3,'--',...
    'Color',[.6,.3,0])
if zLNc <= hclt
else
    plot([min(sigma) max(sigma)]*10^-6,[zLNc zLNc]*10^3,'--',...
        'Color',[.5,.5,.5])
end
plot([0 0],[0 h*10^3],'k')
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
ylim([0 h*10^3])
xlabel('Normal stress (MPa)')
ylabel('CLT-Concrete tickness (mm)')
legend('Stress', 'NA CLT', 'NA Concrete','EdgeColor','none')
title('NORMAL STRESS AT MIDSPAN - ULS')
saveas(gcf,'StressResults.png')

%% RESISTANCE BENDING MOMENT
MRdc = fcd*EIef/(Ecs*(0.5*hc+y(1)*a(1)));   
                                     % bending moment resistence of the
                                     % concrete
MRdclt = EIef/(0.5*E(1)*hclt/fbd + a(2)*EA(2)/(ft0d*A0liq));
                                     % bending moment resistence of the
                                     % CLT
MRd = min([MRdc,MRdclt]);            % bending moment resistence
etab = MSd/MRd;                      % design factor
etabc = MSd/MRdc;                    % design factor - concrete
etabclt = MSd/MRdclt;                % design factor - CLT

fprintf(fid,'Concrete tickness (hcef): \t\t%.2f mm\n',hcef*1000);
fprintf(fid,'Eff. bending stiffness (EI)ef:');
fprintf(fid,'\t%2.3fe6 N.m^2\n',EIef*10^-6);
fprintf(fid,'Coef. of Composition (gammac)\t%.3f\n',y(1));
fprintf(fid,'\n- BENDING MOMENT\n');
fprintf(fid,'Design bending moment (MSd): \t%.2f kN.m\n',MSd/1000);
fprintf(fid,'Resistant bending CLT (MRdCLT): %.2f kN.m\n',MRdclt/1000);
fprintf(fid,'Resistant bending Conc. (MRdc): %.2f kN.m\n',MRdc/1000);
fprintf(fid,'Design factor (etab): \t\t\t%.2f %%\n',etab*100);

%% RESISTANCE SHEAR EFFORT
VRds = EIef/(EI(2)+0.5*EA(2)*a(2)*hclt)*ny*FRds;   
                                     % in function of the connector
VRclt = EI(2)*b*fvrd/(E(1)*A(1)*(h1/2+h12+h2/2));
VRdclt = EIef/(EI(2)+0.5*EA(2)*a(2)*hclt)*VRclt;
                                     % in function of the CLT
VRdc = EIef/(EI(1)+0.5*y(1)*EA(1)*a(1)*(2*hc-hcef+t))*VRdc1;
                                     % in function of the concrete                                 
VRd = min([VRds,VRdclt,VRdc]);       % shear resistence of concrete
etav = VSd/VRd;                      % design factor
etavs = VSd/VRds;                    % design factor - connector
etavclt = VSd/VRdclt;                % design factor - CLT
etavc = VSd/VRdc;                    % design factor - Concrete

fprintf(fid,'\n- SHEAR EFFORT\n');
fprintf(fid,'Design shear effort (VSd): \t\t%.2f kN\n',VSd/1000);
fprintf(fid,'Resistant shear connec. (VRds): %.2f kN\n',VRds/1000);
fprintf(fid,'Resistant shear CLT (VRdCLT):\t%.2f kN\n',VRdclt/1000);
fprintf(fid,'Resistant shear Conc. (VRdc):\t%.2f kN\n',VRdc/1000);
fprintf(fid,'Design factor (etav): \t\t\t%.2f %%\n',etav*100);

%% SERVICE LIMIT STATE - INSTANTANEOUS DISPLACEMENT
K = ny*Kser;                        % slip modulus of connector (N/m)
qSk   = g1kclt + g1kc + g2k + q1k;  % characteristic load (N/m)                                    
MSk = qSk*L^2/8;                    % design bending moment (N.m)

% CLT-Concrete
hcef=hc;                    % effective thickness (m)
sigma_cinf = 2*fctkinf;     % stress inderside of concrete (Pa)
it=0;                      % iteration counter
EA(2) = 0;                  % axial stiffness of CLT (N)
    for i=1:nclt
        EA(2) = E(i)*A(i)+EA(2);
    end

% Loop to reduce effective concrete thickness
while sigma_cinf >= fctkinf
    if sigma_cinf <= fctd || it==0
        hcef=hc;
    else
        hcef=hcef*(1-.005);
    end  
    
    % Axial stiffness (N)
    EA(1) = Ecs*b*hcef;         % Concrete
    EA(2) = 0;                  % CLT
    for i=1:nclt
        EA(2) = E(i)*A(i)+EA(2);
    end
    
    % Composition coefficient (-)
    y(1) = 1/(1+pi^2*EA(1)*sef/(K*L^2)); % Concreto
    y(2) = 1;                            % CLT
    
    % Distance between CG (m)
    at = h-(hcef+hclt)/2;
    a(1) = y(2)*EA(2)/(y(1)*EA(1)+y(2)*EA(2))*at;
    a(2) = y(1)*EA(1)/(y(1)*EA(1)+y(2)*EA(2))*at;
    
    % Effective bending stiffness (N.m^2)
    EI(1) = Ecs*b*hcef^3/12;
    EIef = 0;
    for i=1:2
        EIef = EI(i)+ y(i)*EA(i)*a(i)^2 + EIef;
    end
       
    % Stress (Pa)
    sigma_csup = -MSk*Ecs/EIef*(y(1)*a(1)+0.5*hcef);	% top concrete
    sigma_cinf = -MSk*Ecs/EIef*(y(1)*a(1)-0.5*hcef);	% underside conc.
    sigma_cltsup = MSk*E(3)/EIef*(y(2)*a(2)-0.5*hclt);	% top CLT
    sigma_cltinf = MSk*E(1)/EIef*(y(2)*a(2)+0.5*hclt);	% underside CLT
  
    it=it+1;      
end
t = t+ hc-hcef;             % space between concrete and CLT (m)

% Load vector (N/m)
p = [g1kclt    
    g1kc
    g2k
    q1k];

winst = sum(5.*p.*L^4./(384*EIef)); % instantaneous displacement (m)
winstlim = L/500;                  	% limit of instantaneous  
                                  	% displamcement (m)
etawinst = winst/winstlim;        	% design factor

% Print numerical results
fprintf(fid,'\n:::::::::::::SERVICE LIMIT STATE::::::::::::::\n');
fprintf(fid,'Number of iterations: \t\t\t%.0f\n',it);
fprintf(fid,'Concrete tickness (hcef): \t\t%.2f mm\n',hcef*1000);
fprintf(fid,'Eff. bending stiffness (EI)ef:');
fprintf(fid,'\t%2.3fe6 N.m^2\n',EIef*10^-6);
fprintf(fid,'Coef. of Composition (gammac):\t%.3f\n',y(1));
fprintf(fid,'\n- INSTANTANEOUS DISPLACEMENT\n');
fprintf(fid,'Displacement (winst): \t\t\t%.1f mm\n',winst*1000);
fprintf(fid,'Limit (winstlim): \t\t\t\t%.1f mm\n',winstlim*1000);
fprintf(fid,'Design factor (etawinst): \t\t%.2f %%\n',etawinst*100);

%% Service Limit State - Vibration
ml = rhow*b*hclt+rhoc*b*hc;         % linear mass (kg/m)
f1 = pi/(2*L^2)*sqrt(EIef/ml);      % natural frequency (Hz)
w1kN = 1000*L^3/(48*EIef);          % displacement duo to 1 kN load
w1kNlim = f1^1.43/39000;            % limit displacement
etaf = w1kN/w1kNlim;                % design factor

fprintf(fid,'\n- EXCESSIVE VIBRATION\n');
fprintf(fid,'Natural frequency (f1): \t\t%.2f Hz\n',f1);
fprintf(fid,'Displacement 1kN (w1kN): \t\t%.2f mm\n',w1kN*1000);
fprintf(fid,'Disp. limit (w1kNlim): \t\t\t%.2f mm\n',w1kNlim*1000);
fprintf(fid,'Design factor (etaf): \t\t\t%.2f %%\n',etaf*100);

%% Service Limit State - Final displacement
K = ny*Kser/(2*(1+phiw));           % creep connector stiffness
Ec = Ecs/(1+phic);                  % creep Young's modulus of concrete
E = E./(1+phiw);                    % creep elastic modulus in bendind 
                                    % of wood layers
G12 = G12./(1+phiw);                % creep transversal elastic modulus                                   
                                    % of 12 wood layer   
G23 = G23./(1+phiw);                % creep transversal elastic modulus                                   
                                    % of 23 wood layers                                      

% Composition coefficient (-)
yclt(1)   = 1/(1+pi^2*E(1)*A(1)*h12/(L^2*b*G12));
yclt(2)   = 1;
yclt(3)   = 1/(1+pi^2*E(3)*A(3)*h23/(L^2*b*G23));

% Distance of CGs in relation to base (m)
zclt(1) = h1/2;
zclt(2) = h1+h12+h2/2;
zclt(3) = h1+h12+h2+h23+h3/2;

yEAz = 0;
yEA = 0;
for i=1:nclt
    yEAz = yclt(i)*E(i)/E(1)*A(i)*zclt(i) + yEAz;
    yEA = yclt(i)*E(i)/E(1)*A(i) + yEA;
end
zCG_CLT = yEAz/yEA;

aclt = zeros(1,nclt);
for i=1:nclt
    aclt(i) = zclt(i)-zCG_CLT;
end

EI(2) = 0;
for i=1:nclt
    EI(2) = E(i)*(I(i)+yclt(i)*A(1)*aclt(i)^2)+ EI(2);
end

% Axial stiffness (N)
EA(1) = Ec*b*hc;            % Concrete
EA(2) = 0;                  % CLT
for i=1:nclt
    EA(2) = E(i)*A(i)+EA(2);           
end

% Composition coefficient (-)
y(1) = 1/(1+pi^2*EA(1)*sef/(K*L^2)); % Concrete
y(2) = 1;                            % CLT

% Distance among CGs (m)
a(1) = y(2)*EA(2)/(y(1)*EA(1)+y(2)*EA(2))*at;
a(2) = y(1)*EA(1)/(y(1)*EA(1)+y(2)*EA(2))*at;

% Effective bending stiffness (N.m^2)
EI(1) = Ec*b*hc^3/12;
EIef = 0;
for i=1:2
    EIef = EI(i)+ y(i)*EA(i)*a(i)^2 + EIef;
end
hcef = hc;
% Loop to reduce effective concrete thickness
sigma_cinf = 2*fctkinf;
while sigma_cinf >= fctkinf
    if sigma_cinf <= fctd || it==0
        hcef=hc;
    else
        hcef=hcef*(1-.005);
    end  
    
    % Axial stiffness (N)
    EA(1) = Ecs*b*hcef;         % Concrete
    EA(2) = 0;                  % CLT
    for i=1:nclt
        EA(2) = E(i)*A(i)+EA(2);
    end
    
    % Composition coefficient (-)
    y(1) = 1/(1+pi^2*EA(1)*sef/(K*L^2)); % Concreto
    y(2) = 1;                            % CLT
    
    % Distance between CG (m)
    at = h-(hcef+hclt)/2;
    a(1) = y(2)*EA(2)/(y(1)*EA(1)+y(2)*EA(2))*at;
    a(2) = y(1)*EA(1)/(y(1)*EA(1)+y(2)*EA(2))*at;
    
    % Effective bending stiffness (N.m^2)
    EI(1) = Ec*b*hcef^3/12;
    EIef = 0;
    for i=1:2
        EIef = EI(i)+ y(i)*EA(i)*a(i)^2 + EIef;
    end

    % Stress (Pa)
    sigma_csup = -MSk*Ec/EIef*(y(1)*a(1)+0.5*hcef);	% top concrete
    sigma_cinf = -MSk*Ec/EIef*(y(1)*a(1)-0.5*hcef);	% underside conc.
    sigma_cltsup = MSk*E(3)/EIef*(y(2)*a(2)-0.5*hclt);	% top CLT
    sigma_cltinf = MSk*E(1)/EIef*(y(2)*a(2)+0.5*hclt);	% underside CLT
  
    it=it+1;     
end
t = t+ hc-hcef;                     % space between concrete and CLT (m)

p(4) = p(4)*psi2;                   % reduced live load (N/m)
wfin = sum(5.*p.*L^4./(384*EIef));  % final displacement (m)
wfinlim = L/300;                    % limit displacement (m)
etawfin = wfin/wfinlim;             % design factor

fprintf(fid,'\n- FINAL DISPLACEMENT\n');
fprintf(fid,'Number of iterations: \t\t\t%.0f\n',it);
fprintf(fid,'Concrete tickness (hcef): \t\t%.2f mm\n',hcef*1000);
fprintf(fid,'Eff. bending stiffness (EI)ef:');
fprintf(fid,'\t%2.3fe6 N.m^2\n',EIef*10^-6);
fprintf(fid,'Coef. of Composition (gammac):\t%.3f\n',y(1));
fprintf(fid,'Displacement (wfin): \t\t\t%.1f mm\n',wfin*1000);
fprintf(fid,'Limit (wfinlim): \t\t\t\t%.1f mm\n',wfinlim*1000);
fprintf(fid,'Design factor (etawfin): \t\t%.2f %%\n',etawfin*100);

%% PLOT TABLE OF RESULTS
% ULS-Bending moment
fprintf(fid,'\nRESUME');
fprintf(fid,'\n---------------------------------');
fprintf(fid,'-------------------------------\n');
fprintf(fid,'CRITERION \t\t\t\t | DF \t\t | GOVERN \t | ACCEPTANCE\n');
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');
% Concrete
fprintf(fid,'ULS - Bending moment\t | %.2f %%\t | ',etabc*100);
fprintf(fid,'Concrete  | ');
if etabc <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
% CLT
fprintf(fid,'\t\t\t\t\t\t | %.2f %%\t | ',etabclt*100);
fprintf(fid,'CLT \t\t | ');
if etabclt <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% if MRd == MRdc
%     fprintf(fid,'Concrete  | ');
% else
%     fprintf(fid,'CLT \t\t | ');
% end
% if etab <= 1
%     fprintf(fid,'OK!\n');
% else
%     fprintf(fid,'NOT ACCEPTED!\n');
% end

% ULS-Shear Effort
fprintf(fid,'ULS - Shear effort\t\t | %.2f %%\t | ',etavs*100);
% Connector
fprintf(fid,'Connector | ');
if etavs <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
% Concrete
fprintf(fid,'\t\t\t\t\t\t | %.2f %%\t | ',etavc*100);
fprintf(fid,'Concrete  | ');
if etavc <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
% CLT
fprintf(fid,'\t\t\t\t\t\t | %.2f %%\t | ',etavclt*100);
fprintf(fid,'CLT \t\t | ');
if etavclt <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');
% if VRd == VRds
%     fprintf(fid,'Connector | ');
% elseif VRd == VRdclt
%     fprintf(fid,'CLT \t\t | ');
% else
%     fprintf(fid,'Concrete  | ');
% end
% if etav <= 1
%     fprintf(fid,'OK!\n');
% else
%     fprintf(fid,'NOT ACCEPTED!\n');
% end

% SLS-Instantanous displacement 
fprintf(fid,'SLS - Inst. Displacement | %.2f %%\t | ',etawinst*100);
fprintf(fid,'--- \t\t | ');
if etawinst <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% SLS-Final displacement
fprintf(fid,'SLS - Final displacement | %.2f %%\t | ',etawfin*100);
fprintf(fid,'--- \t\t | ');
if etawinst <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% SLS-Excessive vibration
fprintf(fid,'SLS - Vibration \t\t | %.2f %%\t | ',etaf*100);
fprintf(fid,'--- \t\t | ');
if etawinst <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

fclose(fid);                    % close fid
%clearvars                       % Clear variables

%% OPEN RESULTS
open Results.txt

fprintf('\nTotal time:\t%.2f s\n',toc)
