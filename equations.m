function dydt=equations(t,y)

%parameters from:
%https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2981-4/tables/1


M=(sin(2*pi/24*t))/2+0.5;
ksynWnt=0.5;
kdegWnt=0.3;
Wnt=y(1);
%wnt=(sin(2*pi/24*t))/2+0.5;%24 hours sinusiod

ksynDest=0.15; %destruction complex synthesis rate
kdegD1=0.05; %destruction complex degradation rate
kdegD2=0.3; %destruction complex degradation rate when wnt is present 
Dest=y(2); %destruction complex cc

ksynB=0.11; % synthesis rate
kdegB1=0.110257; % basic degradation rate
kdegB2=0.12; % degradtion rate when wnt is present
Bcat=y(3); %beta_catenin cc

%dWnt/dt
dydt(1)=ksynWnt*M-kdegWnt*Wnt;

%dDest/dt
dydt(2)=ksynDest-(kdegD1+kdegD2*Wnt)*Dest;

%dBcat/dt
dydt(3)=ksynB*Wnt-(kdegB1+kdegB2*Dest)*Bcat;
dydt=dydt';


end
