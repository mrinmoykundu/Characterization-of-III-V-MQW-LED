clc;
clear all;
close all;

set(0,'DefaultAxesFontName', 'Latex');
set(0,'DefaultAxesFontSize', 15);

%% External environment

T = 300; %K
J = 100; %A/cm2 


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

d = 4e-9;         %Well width (3nm - 8nm)
N_mqw =4;       %Number of QW

N_state = 100;    %Number of states to consider

%Recombination happens in p-GaAs region
% tau_r = 0.5e-9;
tau_r = 1/(sqrt((Br/(d*q)))*sqrt(J));


%% 3i_1
E = linspace(Eg-0.5*q,Eg+2.5*q,5000);

lambda = (h*c)./E;
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
Nw = sqrt(1/(Br*d*q))*sqrt(J);


% Ef_n = k_B*T*log((Nd+deln)/ni);     % Ef_n-Ef_i
% Ef_p = k_B*T*log((Na+deln)/ni);     % Ef_i-Ef_p

Ef_n = k_B*T*log(Nw/ni);     % Ef_n-Ef_i
Ef_p = k_B*T*log(Nw/ni);     % Ef_i-Ef_p


Efn_Efp = Ef_n + Ef_p;
f = exp((Efn_Efp-Eg)/(k_B*T)).*exp(-delE./(k_B*T));


% Total rsp(E)
R_sp = Pem .* Nj_2D .* f;
R_sp_cgs = R_sp*q/100^3;


fig = figure();
subplot(211);
plot(lambda/1e-9,real(R_sp_cgs), 'LineWidth',2);
xlabel('\lambda (nm)');
ylabel('R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaAs MQW');
grid on;
hold on;
yline(max(real(R_sp_cgs))/2)

subplot(212);
plot(E/q,real(R_sp_cgs), 'LineWidth',2); hold on;

xlabel('E (eV)');
ylabel('R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaAs MQW');
grid on;



% SI unit
Vol = 20e-6*20e-6*d*N_mqw; %m3, 20umx20um
phi = Vol*sum(R_sp)*abs((E(2)-E(1))) % 1/s

P = Vol*sum(R_sp .* E)*abs((E(2)-E(1))) % J/s
% exportgraphics(fig,'i_1.png','Resolution',600);


%%
peak_lambda = lambda(find(R_sp==max(R_sp),1))

% idx1 = find(R_sp > max(R_sp)/2, 1, "first");
% idx2 = find(R_sp < max(R_sp)/2, 2, "last");
% LW = (lambda(idx1) - lambda(idx2))/1e-9

%% Compare density of states

Nj_bulk = (2*mr)^1.5/(2*pi^2*h_cut^3).*(delE.^0.5);

s = zeros(1,length(E));
for n = 1:N_state
    Ecn = h_cut^2/(2*me) * (n*pi/d)^2;
    Evn = h_cut^2/(2*mh) * (n*pi/d)^2;
    s = s + (E>(Eg+Ecn+Evn));
end
Nj_2D_s = mr/(pi*h_cut^2).*s/d;

figure();
plot(E/q, real(Nj_bulk),'LineWidth',2); hold on;
plot(E/q, Nj_2D_s,'LineWidth',2); hold on;
xlabel("E(eV)");
ylabel("N_J(1/eV 1/cm^3)");
grid on;

title("Comparison of Joint density of states for d=5nm");

legend("Bulk N_J(E)","2D QW N_J(E)");


%%
figure();
plot(E/q,f/max(f),'-', 'LineWidth',1); hold on;
plot(E/q,Nj_2D./max(Nj_2D),'-', 'LineWidth',1);