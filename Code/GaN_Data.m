
%% Data for GaAs
Na = 1e15;
Nd = 1e18;

ni = 1.9e-10;
npo = ni^2/Na;
pno = ni^2/Nd;


%Assuming Low level injection
% deln = 1e15;


%Assuming mid level injection
deln = 1e17;

n = Nd+deln;
p = pno+deln;


% cm2/Vs
mu_n = 1800;
mu_p = 30;

% SI
K_B = 1.38e-23;
T = 300;
q = 1.6e-19;


Eg = 3.4*q;


% cm2/s
Dn = (K_B*T*mu_n)/q;
Dp = (K_B*T*mu_p)/q;
Br = 1e-10;
tau_n = 1/(Br*(Nd+pno+deln));
tau_p = 1/(Br*(npo+Na+deln));

Ln = (Dn*tau_n)^0.5;
Lp = (Dp*tau_p)^0.5;

%cm to SI m
Ln_SI = Ln*1e-2;
Lp_SI = Lp*1e-2;


mo = 9.1e-28;
me = 0.2*mo;
mh = 0.8*mo;
mr =(me*mh)/(me+mh);

sr = 1e-14;
h = 6.626e-27;
h_cut = h/(2*pi);
c = 3e10;
Vol = 0.5e-3;
K_cm = 1.38e-16;

% Typical trap density of GaN
NT = 1e13;



% Refractive index
nr1 = 2.55; %GaN
nr2 = 1;   %Air