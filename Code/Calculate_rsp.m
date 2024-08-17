function [peak_lambda, phi, P, R_sp_cgs,lambda] = Calculate_rsp(d,N_mqw,J,T)
%% External environment

% T = 300;

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
tau_r = 1./(sqrt((Br/(d*q)))*sqrt(J)); 


%% 3i_1
lambda = linspace(100,900,1000)*1e-9;
E = (h*c)./(lambda);
delE = (E-Eg);



% Calculating probability of radiative recombination/emission
Pem = 1./tau_r;



% Calculating Density of states
s = zeros(1,length(E));
for n = 1:N_state
    Ecn = h_cut^2/(2*me) * (n*pi/d)^2;
    Evn = h_cut^2/(2*mh) * (n*pi/d)^2;
    s = s + (E>(Eg+Ecn+Evn));
end
Nj_2D = N_mqw*mr/(pi*h_cut^2).*s/d;     % (1/eV .1/cm3)


% Calulating Distribution function (Weak Injection assumption)
Nw = sqrt(1/(Br*d*q)).*sqrt(J);


% Ef_n = k_B*T*log((Nd+deln)/ni);     % Ef_n-Ef_i
% Ef_p = k_B*T*log((Na+deln)/ni);     % Ef_i-Ef_p

Ef_n = k_B*T*log(Nw/ni);     % Ef_n-Ef_i
Ef_p = k_B*T*log(Nw/ni);     % Ef_i-Ef_p


Efn_Efp = Ef_n + Ef_p;
f = exp((Efn_Efp-Eg)/(k_B*T)).*exp(-delE./(k_B*T));



% Total rsp(E)
R_sp = Pem .* Nj_2D .* f;
R_sp_cgs = R_sp*q/100^3;



% SI unit
% Per active area
% Vol = 20e-6*20e-6*d; %m3
Vol = 1e-6*1e-6*d; %m3


phi = Vol*sum(R_sp)*(E(1)-E(2)); % 1/s
P = Vol*sum(R_sp .* E)*(E(1)-E(2)); % J/s

peak_lambda = lambda(find(R_sp==max(R_sp),1));


% idx1 = find(R_sp < max(R_sp)/2,"last",1);
% idx2 = find(R_sp < max(R_sp)/2,"first",1);
% LW = lambda(idx2) - lambda(idx1);


end