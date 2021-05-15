clear all;
close all;
clc

t_max=100;
t_interval=0:0.01:t_max;

init_cond=[1.46,0.8661,0.7964];

[t,y]=ode45(@(t,y)equations(t,y),t_interval,init_cond);
%[t,y]=ode45(@(t,y)equations(t,y),t_interval,init_cond);
results.t=t;
results.y=y;

while t(end) < t_max
  init_cond=results.y(end,:); 
  %init_cond(2)=init_cond(2)*0.5;
  %init_cond(3)=init_cond(3)*0.5;
  [t,y,]=ode45(@(t,y)wntmodel(t,y),[results.t(end),t_max],init_cond);  
  results.t=[results.t;t];
  results.y=[results.y;y];
  results.t(end)
end
M=(sin(2*pi/24*t))/2+0.5;
plot(t,y,'-');
hold on;
%plot(t,M,'-');

title('Concentration change of Wnt, destruction complex,and beta-catenin');
xlabel('Time (hour)');
xlim([20 100]);
%ylim([0 5]);
ylabel('Concentration [arbitrary unit]');
legend('wnt','Dest','Bcat');