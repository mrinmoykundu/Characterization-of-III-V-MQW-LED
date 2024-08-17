function [phi, P] = Calculate_phi(d,N_mqw,J)
%% External environment

T = 300;

%% Choose data

i = 1;
switch (i)
    case 1 
        GaAs_Data()
    case 2 
        GaN_Data()
end


% SI unit
h = 6.626e-34;
h_cut = h/(2*pi);
c = 3e8;
k_B = 1.38e-23;

%% MQW Parameters

% d = 5e-9;       %Well width
% N_mqw = 50;     %Number of QW

N_state = 100;    %Number of states to consider

% tau_r is dependent of J
tau_r = 1/(sqrt((Br/(d*q)))*sqrt(J)); 


%% 3i_1
lambda = linspace(100,900,1000)*1e-9;
E = (h*c)./(lambda);
delE = (E-Eg);



% Calculating probability of radiative recombination/emission
Pem = 1/tau_r;



% Calculating Density of states
s = zeros(1,length(E));
for n = 1:N_state
    Ecn = h_cut^2/(2*me) * (n*pi/d)^2;
    Evn = h_cut^2/(2*mh) * (n*pi/d)^2;
    s = s + (E>(Eg+Ecn+Evn));
end
Nj_2D = N_mqw*mr/(pi*h_cut^2).*s/d;     % (1/eV .1/cm3)


% Calulating Distribution function (Weak Injection assumption)
Ef_n = k_B*T*log((Nd+deln)/ni);     % Ef_n-Ef_i
Ef_p = k_B*T*log((Na+deln)/ni);     % Ef_i-Ef_p

Efn_Efp = Ef_n + Ef_p;
f = exp((Efn_Efp-Eg)/(k_B*T)).*exp(-delE./(k_B*T));


% Total rsp(E)
R_sp = Pem .* Nj_2D .* f;



% SI unit
Vol = 20e-6*20e-6*d*N_mqw; %m3
phi = Vol*sum(R_sp)*(E(1)-E(2)); % 1/s

P = Vol*sum(R_sp .* E)*(E(1)-E(2)); % J/s




end