function vel = calculate_velocity(j,location)
global wind_farm_x;
z0 = 0.3; a = 1/3;
counter = length(location);
Temp = wind_farm_x(location,:);

for i=1:counter-1
    ydistance = abs(Temp(i+1,2)-Temp(i,2));
    alpha = 0.5/log(Temp(i,4)/z0);
    radius = Temp(i,3)*(sqrt((1-a)/(1-2*a)))+alpha*ydistance;
    xmin = Temp(i,1) - radius;
    xmax = Temp(i,1) + radius;
    if (xmin < Temp(i+1,1)) && (xmax > Temp(i+1,1))
        Temp(i,:)= 0;
    end
end
Temp(all(~Temp,2),:)=[];
ydistance1 = abs(wind_farm_x(j,2) - Temp(:,2));
alpha1 = 0.5./log(Temp(:,4)./z0);
radius1 = Temp(:,3).*(sqrt((1-a)/(1-2*a)));
denominator = ((alpha1.* ydistance1)./ radius1 + 1).^ 2;
velr = (1-(1-(2*a./denominator))).^2;

vel = 1 - (sum(velr)^0.5);

