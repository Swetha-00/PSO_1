function obj = Analyse_Grid(wf_var_hor_array)
global wind_farm wind_farm_x velocity_farm farm_power Efficiency cost_per_kW;
u01 = 8; u02 = 12; u03 = 17; z0 = 0.3; h0 = 60; r0 = 20;
total_power = zeros(); total_power_cap = zeros();sum_total_power = zeros();sum_total_power_cap = zeros();
velocity1_farm = zeros();velocity2_farm = zeros();velocity3_farm = zeros();
position = 1;
%Partha's best layout
%wf_var_hor_array = [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 1 0 1 0 1 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 0 1 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1 1 1];
%Grady's layout
%wf_var_hor_array = [1 0 1 1 1 0 1 1 0 1 1 0 0 0 0 0 0 0 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 0 0 0 0 1 0 0 1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 0 0 0 1 0 1 0 1 0 0 0 1 0 0 1 1 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1];
%Evolutive algorithm
%wf_var_hor_array = [1 0 1 1 1 0 1 0 1 1 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 1 1 0 0 0 1 0 1 0 1 0 0 1 0 1 0 0 0 0 1 1 0 0 0 0 0 1 0 1 0 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 0 0 1 0 1 1 0 0 0 1 0 0 0 0 1 1 0 1 0 1 0 1 1 0 1];
%Pookpunt layout - Binary PSO
%wf_var_hor_array = [1 1 1 1 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 1 1 0 0 1 0 1 0 0 0 1 1 0 1 0 1 0 0 0 0 1 1 0 0 1 0 0 1 0 0 1 1 0 0 0 0 0 0 1 0 1 1 0 1 0 1 0 1 0 0 1 0 1 0 0 0 0 0 1 0 1 0 1 0 1 0 1 0 0 0 1 1 1 1 0 1 1 1 1 1 1];
%wf_var_hor_array = [0 1 0 1 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 0 0 0 1 0 1 0 0 1 1 0 0 0 0 1 0 1 0 1 1 0 0 1 0 0 1 0 0 1 1 0 1 0 0 0 0 0 0 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 0 0 0 1 0 1 0 0 0 1 0 1 0 1 0 1 1 1 1 1 1 0 1 1 1];
wf_var_hor_array = (round(wf_var_hor_array))';
x_loc = repmat((100:200:2000)',10,1);
y_loc = repelem((100:200:2000)',10,1);
wf_t = horzcat(wf_var_hor_array,x_loc,y_loc);
wf_temp = wf_t(all(wf_t,2),:);
wf_var = ...
    [20 60];
wf_var_sel = wf_var(wf_temp(:,1),:);
wf_temp(:,1)=[];
wind_farm = horzcat(wf_temp,wf_var_sel);
wind_farm = sortrows(wind_farm,[2 1]);
N = size(wind_farm,1);
velocity_farm = zeros(N,1);
wind_farm_x = [];
wind_farm_x(:,3)= wind_farm(:,3);
wind_farm_x(:,4)= wind_farm(:,4);
   
for theta = 0:10:350
    x = wind_farm(:,1);
    y = wind_farm(:,2);        
    wind_farm_x(:,1)= sqrt(x.^2+y.^2).*cosd(atand(y./x)+theta);
    wind_farm_x(:,2)= sqrt(x.^2+y.^2).*sind(atand(y./x)+theta);
    wind_farm_x = round(wind_farm_x);
    wind_farm_x = sortrows(wind_farm_x,[2 1]);
    u011 = u01*log(wind_farm_x(:,4)./z0)/log(h0/z0);
    u022 = u02*log(wind_farm_x(:,4)./z0)/log(h0/z0);
    u033 = u03*log(wind_farm_x(:,4)./z0)/log(h0/z0);
    for i = 1:1:N
    x1 = wind_farm_x(i,1);
    y1 = wind_farm_x(i,2);
    velocity = check_wake(x1,y1,i);
    if velocity == 1
        velocity1_farm(i) = u011(i);
        velocity2_farm(i) = u022(i);
        velocity3_farm(i) = u033(i);
    else
        velocity1_farm(i) = u011(i)*velocity;
        velocity2_farm(i) = u022(i)*velocity;
        velocity3_farm(i) = u033(i)*velocity;
    end
    if theta <= 260
        %power = 0.3*(0.004*velocity1_farm(i)^3+0.008*velocity2_farm(i)^3+0.012*min(12.8058,velocity3_farm(i))^3);
        %power_cap = 0.3*(0.004*u011(i)^3+0.008*u022(i)^3+0.012*min(12.8058,u033(i))^3);
        power = 0.3*(0.004*velocity1_farm(i)^3+0.008*velocity2_farm(i)^3+0.012*velocity3_farm(i)^3);
        power_cap = 0.3*(0.004*u011(i)^3+0.008*u022(i)^3+0.012*u033(i)^3);
    elseif (theta==270) || (theta==350)
        %power = 0.3*(0.004*velocity1_farm(i)^3+0.01*velocity2_farm(i)^3+0.012*min(12.8058,velocity3_farm(i))^3);
        %power_cap = 0.3*(0.004*u011(i)^3+0.01*u022(i)^3+0.012*min(12.8058,u033(i))^3);
        power = 0.3*(0.004*velocity1_farm(i)^3+0.01*velocity2_farm(i)^3+0.012*velocity3_farm(i)^3);
        power_cap = 0.3*(0.004*u011(i)^3+0.01*u022(i)^3+0.012*u033(i)^3);
    elseif (theta==280) || (theta==340)
        %power = 0.3*(0.004*velocity1_farm(i)^3+0.012*velocity2_farm(i)^3+0.016*min(12.8058,velocity3_farm(i))^3);
        %power_cap = 0.3*(0.004*u011(i)^3+0.012*u022(i)^3+0.016*min(12.8058,u033(i))^3);
        power = 0.3*(0.004*velocity1_farm(i)^3+0.012*velocity2_farm(i)^3+0.016*velocity3_farm(i)^3);
        power_cap = 0.3*(0.004*u011(i)^3+0.012*u022(i)^3+0.016*u033(i)^3);
    elseif (theta==290) || (theta==330)
        %power = 0.3*(0.004*velocity1_farm(i)^3+0.015*velocity2_farm(i)^3+0.019*min(12.8058,velocity3_farm(i))^3);
        %power_cap = 0.3*(0.004*u011(i)^3+0.015*u022(i)^3+0.019*min(12.8058,u033(i))^3);
        power = 0.3*(0.004*velocity1_farm(i)^3+0.015*velocity2_farm(i)^3+0.019*velocity3_farm(i)^3);
        power_cap = 0.3*(0.004*u011(i)^3+0.015*u022(i)^3+0.019*u033(i)^3);
    elseif (theta==300) || (theta==320)
        %power = 0.3*(0.004*velocity1_farm(i)^3+0.016*velocity2_farm(i)^3+0.03*min(12.8058,velocity3_farm(i))^3);
        %power_cap = 0.3*(0.004*u011(i)^3+0.016*u022(i)^3+0.03*min(12.8058,u033(i))^3);
        power = 0.3*(0.004*velocity1_farm(i)^3+0.016*velocity2_farm(i)^3+0.03*velocity3_farm(i)^3);
        power_cap = 0.3*(0.004*u011(i)^3+0.016*u022(i)^3+0.03*u033(i)^3);
    elseif (theta==310)
        %power = 0.3*(0.004*velocity1_farm(i)^3+0.02*velocity2_farm(i)^3+0.036*min(12.8058,velocity3_farm(i))^3);
        %power_cap = 0.3*(0.004*u011(i)^3+0.02*u022(i)^3+0.036*min(12.8058,u033(i))^3);
        power = 0.3*(0.004*velocity1_farm(i)^3+0.02*velocity2_farm(i)^3+0.036*velocity3_farm(i)^3);
        power_cap = 0.3*(0.004*u011(i)^3+0.02*u022(i)^3+0.036*u033(i)^3);
    end
    total_power(i) = power;
    total_power_cap(i) = power_cap;
    end
    sum_total_power(position)= sum(total_power);
    sum_total_power_cap(position) = sum(total_power_cap);
    position = position+1;
end

farm_power = sum(sum_total_power);
farm_power_cap = sum(sum_total_power_cap);
Efficiency = farm_power/farm_power_cap;

cost_base = N*(2/3+(1/3)*exp(-0.00174*N^2));
cost = cost_base*(1+(1/N)*sum(0.39825*(wind_farm(:,3)./r0-1))+(1/N)*sum(0.16425*(wind_farm(:,4)./h0-1)));
cost_per_kW = cost/farm_power;

obj = cost_per_kW;
% fprintf('obj %0.8f',obj);
% fprintf('\n total_power %.6f',farm_power);
% fprintf('\n Efficiency %.6f',Efficiency);
% 
%     figure('DefaultAxesFontSize',9);
%     x1 = wind_farm(:,1);
%     y1 = wind_farm(:,2);
%     plot(x1,y1,'ks');
%     set(gca, 'XTickLabel', []);
%     set(gca, 'YTickLabel', []);
%     hold on
%     plot(x1,y1,'ks','MarkerSize',14, 'MarkerFaceColor', [0 0 0])
%     set(gca,'DefaultTextFontSize',12);
%     hold on;
%     grid on;
