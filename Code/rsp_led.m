clc;
clear all;

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
%% 3i_1
lambda = linspace(750,900,1000)*1e-9;
E = (h*c)./(lambda);
delE = (E-Eg);

%Recombination happens in p-GaAs region
tau_r = tau_p; 
Pem = 1/tau_r;

Ef_n = k_B*T*log((Nd+deln)/ni);     % Ef_n-Ef_i
Ef_p = k_B*T*log((Na+deln)/ni);     % Ef_i-Ef_p

Efn_Efp = Ef_n + Ef_p;

Nj_bulk = real((2*mr)^1.5/(2*pi^2*h_cut^3).*(delE.^0.5));
f = exp((Efn_Efp-Eg)/(k_B*T)).*exp(-delE./(k_B*T));

R_sp = Pem .* Nj_bulk .* f;

% R_sp = ((((2*mr)^1.5)/(2*pi^2*(h_cut^3)*tau_r)) ...
%     *exp((Efn_Efp-Eg)/(k_B*T)).*(delE.^0.5).*exp(-delE./(k_B*T)));

%SI to 1/s. 1/eV . 1/cm3 unit
R_sp_cgs = R_sp*q/100^3;


fig = figure();
subplot(211);
plot(lambda/1e-9,real(R_sp_cgs), 'LineWidth',2);
xlabel('\lambda (nm)');
ylabel('R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaAs');
grid on;

subplot(212);
plot(E/q,real(R_sp_cgs), 'LineWidth',2);
xlabel('E (eV)');
ylabel('R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaAs');
grid on;



% SI unit
Vol = 0.5/1000^3;
phi = Vol*sum(R_sp)*(E(1)-E(2))
% exportgraphics(fig,'i_1.png','Resolution',600);

%% 3i_2

lambda = linspace(300,400,1000)*1e-9;
E = (h*c)./(lambda);
delE = (E-Eg);

%Recombination happens in p-GaAs region
tau_r = tau_p; 

Ef_n = k_B*T*log((Nd+deln)/ni);     % Ef_n-Ef_i
Ef_p = k_B*T*log((Na+deln)/ni);     % Ef_i-Ef_p

Efn_Efp = Ef_n + Ef_p;


R_sp = ((((2*mr)^1.5)/(2*pi^2*(h_cut^3)*tau_r)) ...
    *exp((Efn_Efp-Eg)/(k_B*T)).*(delE.^0.5).*exp(-delE./(k_B*T)));

%SI to 1/s. 1/eV . 1/cm3 unit
R_sp_cgs = R_sp*q/100^3;


fig = figure();
subplot(211);
plot(lambda/1e-9,real(R_sp_cgs), 'LineWidth',2);
xlabel('\lambda (nm)');
ylabel('R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaN');
grid on;

subplot(212);
plot(E/q,real(R_sp_cgs), 'LineWidth',2);
xlabel('E (eV)');
ylabel('R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaN');
grid on;



% SI unit
Vol = 0.5/1000^3;
phi = (Vol/(sqrt(2)*h_cut^3*tau_r))*((mr*k_B*T/pi)^1.5)*exp((Efn_Efp-Eg)/(k_B*T))