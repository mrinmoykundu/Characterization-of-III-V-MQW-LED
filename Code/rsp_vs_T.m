clear;


%% Plotting rsp vs T


Br = 1e-10; %1/cm3
q = 1.6e-19; %C

%% MQW Parameters
d = 15e-9;       %Well width
N_mqw = 50;     %Number of QW
J = 100; %A/cm2

T = [300,400,500];
figure(1);
for i = 1:length(T)
      [~, ~,~,R_sp_cgs,lambda] = Calculate_rsp(d,N_mqw,J,T(i));
      plot(lambda/1e-9,R_sp_cgs/max(R_sp_cgs),'LineWidth',1.5,'DisplayName',sprintf("T = %dK",T(i))); 
hold on;

end

xlabel('\lambda (nm)');
ylabel('Norm R_{sp} (1/s. 1/eV . 1/cm^3)');
title('Emission Spectra of GaAs MQW');
xlim([600,1000]);
grid on;
legend();
legend box off;

%%
d = 15e-9;       %Well width
N_mqw = 50;     %Number of QW
J = 100; %A/cm2

T = linspace(300,500,100);

P = zeros(1,length(T));

figure(2);
for i = 1:length(T)
    [~, ~,P(i),~,~] = Calculate_rsp(d,N_mqw,J,T(i));
end
plot(T,P/1e-6,'LineWidth',2); 


xlabel('T(K)');
ylabel('P (\muW)');
title('Light output power (without loss) of GaAs MQW');
grid on;