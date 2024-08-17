clear;


%% Plotting p vs J



Br = 1e-10; %1/cm3
q = 1.6e-19; %C

%% MQW Parameters
d = 4e-9;       %Well width
N_mqw = [1,5,10];     %Number of QW

T = 300;

for n_mqw = N_mqw
    J = linspace(1,1000,100); %A/cm2
    
    phi = zeros(1,length(J));
    P = zeros(1,length(J));
    
    for i = 1:length(J)
          [~, phi(i), P(i),~,~] = Calculate_rsp(d,n_mqw,J(i),T);
    end
    
    semilogy(J,P/1e-6,'LineWidth',2,'DisplayName',sprintf("N_{QW} = %d",n_mqw)); 
    hold on;
end

xlabel('J (A/cm^2)');
ylabel('P (\muW) per \mum^2');
title('Light output power (without loss) of GaAs MQW');
grid on;
legend();
legend box off;