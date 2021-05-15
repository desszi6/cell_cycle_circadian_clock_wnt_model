clear all;
close all;
clc

%[t,y]=ode23(@(t,y)equations(t,y,rhythm),tspan,Bcat)
t_max=600;
t_interval=[0 t_max];

%init_cond=[0,0]; % first one is dest complex the second element is beta catenin
%init_cond=zeros(1,28);
%init_cond=[0.0025,0.0028,2.3,0.077,0.016,0.012181,0.90053,0.012876,0.4343,0.016934,0.48284,0.0014623,0.040383,0.84912,0.99,1,1.7e-05,1e-05,9e-05,4.9572,0.9381,1.076,1.4,0.037,0.046,0.130];%1.4,0.037,0.046,0.13
init_cond=[0.0025,0.0028,2.3,0.077,0.016,0.012181,0.90053,0.012876,0.4343,0.016934,0.48284,0.0014623,0.040383,0.84912,0.99,1,1.7e-05,1e-05,9e-05,4.9572,0.9381,1.076,0,0,0];%23.=1.4,24=0.037,25=0.046,26.=0.13
%init_cond=[2.3,0.077,0.016,0.0025,0.0028,0.012181,0.90053,0.012876,0.4343,0.016934,0.48284,0.0014623,0.040383,0.84912,0.99,1,1.7e-05,1e-05,9e-05,4.9572,0.9381,1.076,1.4,0.037,0.046,0.13,0,10,10];
%solution
options=odeset('Events', @myevent);
[t,y,~,~,~]=ode45(@(t,y)wntmodel(t,y),t_interval,init_cond,options);
%[t,y,~,~,~]=ode45(@(t,y)wnt_cc_cr_model(y),t_interval,init_cond,options);
results.t=t;
results.y=y;

while t(end) < t_max
  init_cond=results.y(end,:); 
  init_cond(2)=init_cond(2)*0.5;
  init_cond(3)=init_cond(3)*0.5;
  [t,y,~,~,~]=ode45(@(t,y)wntmodel(t,y),[results.t(end),t_max],init_cond,options);  
  results.t=[results.t;t];
  results.y=[results.y;y];
  results.t(end)
end

%[t,y]=ode45(@(t,y)equations(t,y),t_interval,init_cond);
%[t,y,~,~,~]=ode45(@(t,y)wntmodel(y),t_interval,init_cond,options);
%event detection needed-> while-be adom meg amig el nem érem a szimulációs idõt
%wnt=(sin(2*pi/24*t))/2+0.5;%24 hours sinusiod
%y=zeros(1,28);
%y=[2.3,0.077,0.016,0.0025,0.0028,0.012181,0.90053,0.012876,0.4343,0.016934,0.48284,0.0014623,0.040383,0.84912,0.99,1,1.7e-05,1e-05,9e-05,4.9572,0.9381,1.076,1,8,1.4,0.037,0.046,0.13,0,10,10];
%for i=1:100
 %  dydt(i+1)=wntmodel(y(i));
%end

M=(sin(2*pi/24*t))/2+0.5;
%plot(results.t,results.y,'-')

plot(results.t,results.y(:,1),'-') %CycB
hold on
%plot(results.t,results.y(:,2),'-') %CycBp
%hold on
%plot(results.t,results.y(:,3),'-') %Wee1
%hold on
%plot(results.t,results.y(:,4),'-') %Wee1P
%hold on
%plot(results.t,results.y(:,5),'-') %Cdc25a
%hold on
%plot(results.t,results.y(:,6),'-') %ERG
%hold on
%plot(results.t,results.y(:,7),'-') %DRG
%hold on
plot(results.t,results.y(:,8),'-')%CycD
hold on
%plot(results.t,results.y(:,9),'-') %CD
%hold on
plot(results.t,results.y(:,10),'-')%CycE
hold on
%plot(results.t,results.y(:,11),'-') %CE
%hold on
plot(results.t,results.y(:,12),'-')%CycA
hold on
%plot(results.t,results.y(:,13),'-') %Ca
%hold on
%plot(results.t,results.y(:,14),'-') %P27
%hold on
%plot(results.t,results.y(:,15),'-') %Cdh1 !?!?!
%hold on
%plot(results.t,results.y(:,16),'-') %PPX
%hold on
%plot(results.t,results.y(:,17),'-') %IEP
%hold on
%plot(results.t,results.y(:,18),'-') %Cdc20
%hold on
%plot(results.t,results.y(:,19),'-') %Cdc20i
%hold on
%plot(results.t,results.y(:,20),'-') %E2F
%hold on
%plot(results.t,results.y(:,21),'-') %GM
%hold on
%plot(results.t,results.y(:,22),'-') %Mass
%hold on
%plot(results.t,results.y(:,23),'-') %M
%hold on
%plot(results.t,results.y(:,24),'-') %CP
%hold on
%plot(results.t,results.y(:,25),'-') %CP2
%hold on
%plot(results.t,results.y(:,26),'-') %TF
%hold on
plot(results.t,results.y(:,23),'-') %wnt
hold on
%plot(results.t,results.y(:,24),'-') %destruction_complex
%hold on
plot(results.t,results.y(:,25),'-') %Bcat
hold on

%plot (M,t,'-')
%hold on;

title('Change of CycB, CycD, CycE, CycA, Bcat when Wnt is periodically expressed');
xlabel('Time (hour)');
xlim([100 600]);
%ylim([0 5]);
ylabel('Concentration [arbitrary unit]');
legend('CycB', 'CycD','CycE','CycA','Wnt','Bcat');

