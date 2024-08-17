clear;


%% Plotting p vs J



Br = 1e-10; %1/cm3
q = 1.6e-19; %C

n_r = 0.92;
n_ex = 0.0015;

%% MQW Parameters
d = 4e-9;       %Well width
n_mqw = 4;     %Number of QW

T = 300;

J = linspace(1,1000,100); %A/cm2

phi = zeros(1,length(J));
P = zeros(1,length(J));

for i = 1:length(J)
      [~, phi(i), P_temp,~,~] = Calculate_rsp(d,n_mqw,J(i),T);
      P(i) = n_r*P_temp;
end
semilogy(J,P/1e-6,'LineWidth',2,'DisplayName',sprintf("Internally Generated")); 
hold on;

semilogy(J,P*n_ex/1e-6,'LineWidth',2,'DisplayName',sprintf("Externally Output")); 
hold on;



xlabel('J (A/cm^2)');
ylabel('P (\muW) per \mum^2');
title('Comparison Between Internally Generated and Extracted Power');
subtitle('n_{ex} = 0.0015, GaAs MQW');
grid on;
legend();
legend box off;