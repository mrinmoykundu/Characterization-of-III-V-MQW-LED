clc;
clear all;
close all;

set(0,'DefaultAxesFontName', 'Latex');
set(0,'DefaultAxesFontSize', 15);

%% External environment

T = 300; %K
J_MQW = 100; %A/cm2 


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
N_mqw = 4;       %Number of QW

N_state = 100;    %Number of states to consider

%Recombination happens in p-GaAs region
% tau_r = 0.5e-9;
tau_r = 1/(sqrt((Br/(d*q)))*sqrt(J_MQW));


%% Single QW

Vbi = k_B*T*log(Nd*Na/ni^2)/q;


xn = (2*eps*(Vbi/q)*(Na/Nd)*(1/(Nd+Na)))^0.5*100;      %in cm, 1um

Ec = 0;
Ef = -(Eg/2 - K_B*T*log(Nd/ni));

n = @(x) real(sqrt(x-Ec))./exp((x-Ef)/(k_B*T));

Va = linspace(0,Vbi,100);
J_SQW = zeros(1,length(Va));

%MQW
for i = 1:length(Va)

    %Divide by 2 beacause of QW in the middle
    x_start = q*(Vbi-Va(i))/2;
    x_end = x_start + 10*k_B*T;
    
    % E = linspace(0,10*k_B*T,100);
    % plot(E/q,n(E))
    Q = integral(@(x) n(x), x_start, x_end);
    n_inj = me*sqrt(2*me)/(pi^2*h_cut^3)*Q; %1/m3
    
    n_inj_cgs = n_inj/(100^3);
   
    J_SQW(i) = q*Dn*n_inj_cgs/xn; %A/cm2

end


A = 20e-4*20e-4; %20umx20um, in cm


% But No current should flow when Va=0, beacause space charge barrier 
% prevents further electron to diffuse
figure(1);
plot(Va,J_SQW-J_SQW(1),'LineWidth',2);
xlabel('V (V)');
ylabel('J (A/cm^2)');
title('I-V curve, GaAs SQW');
grid on;

hold on;


%% MQW


Ec = 0;
Ef = -(Eg/2 - K_B*T*log(Nd/ni));

n = @(x) real(sqrt(x-Ec))./exp((x-Ef)/(k_B*T));

Va = linspace(0,Vbi,100);
J_MQW = zeros(1,length(Va));



%MQW
for i = 1:length(Va)
    J_MQW(i) = 0;
    for j = 1:N_mqw
        %Divide by 2 beacause of QW in the middle
        x_start = q*(Vbi-Va(i))*(j/(N_mqw+1));
        x_end = x_start + 10*k_B*T;
        
        % E = linspace(0,10*k_B*T,100);
        % plot(E/q,n(E))
        Q = integral(@(x) n(x), x_start, x_end);
        n_inj = me*sqrt(2*me)/(pi^2*h_cut^3)*Q; %1/m3
        
        n_inj_cgs = n_inj/(100^3);
        J_MQW(i) = J_MQW(i)+ q*Dn*n_inj_cgs/xn; %A/cm2
    end
end


A = 20e-4*20e-4; %20umx20um, in cm


% But No current should flow when Va=0, beacause space charge barrier 
% prevents further electron to diffuse
figure(1);
plot(Va,J_MQW-J_MQW(1),'LineWidth',2);
xlabel('V (V)');
ylabel('J (A/cm^2)');
title('I-V curve, GaAs MQW');
grid on;
legend("SQW","MQW, N=4");

xlim([0, Vbi]);
ylim([0, 5e2]);