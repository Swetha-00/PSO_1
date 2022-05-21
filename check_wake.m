function ff = check_wake(x,y,j)
global wind_farm_x;
a = 1/3; z0 = 0.3;
max_wake_rad = 200;
x_low = x-max_wake_rad; x_high = x+max_wake_rad;

temp_farm = horzcat(wind_farm_x(1:j-1,:),(1:j-1)');
temp_farm(temp_farm(:,1)<x_low & temp_farm(:,1)>x_high)=[];
count = size(temp_farm,1);

counter = 0; location = zeros();

for i = 1:1:count
ydistance = abs(y - temp_farm(i,2));    
alpha = 0.5/log(temp_farm(i,4)/z0);
radius = temp_farm(i,3)*sqrt((1-a)/(1-2*a))+(alpha * ydistance);

%wake influence radius
xmin = temp_farm(i,1) - radius;
xmax = temp_farm(i,1) + radius;

if (xmin < x) && (xmax > x) % Checking for wake by radius
% Turbine in wake
counter = counter + 1;
location(counter) = temp_farm(i,5);
end
end

if counter == 0
ff = 1;
else
% Call calculate velocity
velocity = calculate_velocity(j,location);
ff = velocity;
end
