clear;


%% Plotting peak lambda vs d



Br = 1e-10; %1/cm3
q = 1.6e-19; %C

%% MQW Parameters

N_mqw = [55];       %Number of QW
J = 100;            %A/cm2

T = 300;

for n_mqw = N_mqw
    d = linspace(3e-9,8e-9,100);       %Well width
    
    peak_lambda = zeros(1,length(J));
    P = zeros(1,length(J));
    
    for i = 1:length(d)
        [peak_lambda(i), phi(i), P(i),~,~] = Calculate_rsp(d(i),n_mqw,J,T);
    end
    
    figure(1);
    plot(d/1e-9,peak_lambda/1e-9,'LineWidth',2,'DisplayName',sprintf("N_{QW} = %d",n_mqw)); 
    hold on;

    figure(2);
    plot(d/1e-9,P,'LineWidth',2,'DisplayName',sprintf("N_{QW} = %d",n_mqw)); 
    hold on;
end


figure(1);
xlabel('d (nm)');
ylabel('\lambda_{peak} (nm)');
title('Peak \lambda vs QW width of GaAs MQW');
grid on;
legend();

figure(2);
xlabel('d (nm)');
ylabel('P (\muW)');
title('Light output power vs QW width GaAs MQW');
grid on;
legend();