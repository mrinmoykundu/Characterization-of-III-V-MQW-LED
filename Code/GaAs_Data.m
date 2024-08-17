
%% Data for GaAs
eps0    = 8.854e-12;
eps     = 12.9*eps0;

% All units are in cm
Na = 1e15;
Nd = 5e17;

ni = 1.79e6;
npo = ni^2/Na;
pno = ni^2/Nd;


%Assuming Low level injection
deln = 1e15;


%Assuming mid level injection
% deln = 1e17;

n = Nd+deln;
p = pno+deln;

% cm2/Vs
mu_n = 8500;
mu_p = 400;

% SI
K_B = 1.38e-23;
q = 1.6e-19;

Eg = 1.43*q;

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

% Kg
mo = 9.1e-31;
me = 0.067*mo;
mh = 0.48*mo;
mr = (me*mh)/(me+mh);

sr = 1e-14;

% Typical trap density of GaAs
NT = 1e13;


% Refractive index
nr1 = 3.68; %GaAs
nr2 = 1;   %Air
