function dydt=wntmodel(t,y)

%parameters from:
%https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2981-4/tables/1

ksynDest=0.2909; %destruction complex synthesis rate
kdegD1=0.05; %destruction complex degradation rate
kdegD2=0.3; %destruction complex degradation rate when wnt is present 

ksynB=0.423; % synthesis rate
kdegB1=0.1195; % basic degradation rate
kdegB2=0.1383; % degradtion rate when dest is present

ksynWnt=0.4;
kdegWnt=0.3;

%Chris's model:
%parameters for clock-cell cycle coupling via Wnt

% initial conditions
CYCB=y(1);%0.0025;
 CYCBP=y(2);%0.0028;
 Wee1=y(3);% 2.3;
 Wee1P=y(4);%0.077;
 Cdc25a=y(5);%0.016;
 ERG=y(6);%0.012181;
 DRG=y(7);%0.90053;
 CYCD=y(8);%0.012876;
 CD=y(9);%0.4343;
 CYCE=y(10);%0.016934;
 CE=y(11);%0.48284;
 CYCA=y(12);%0.0014623;
 CA=y(13);%0.040383;
 P27=y(14);%0.84912;
 Cdh1=y(15);%0.99;
 PPX=y(16);%1;
 IEP=y(17);%1.7e-05;
 Cdc20=y(18);%1e-05;
 Cdc20i=y(19);%9e-05;
 E2F=y(20);%4.9572;
 GM=y(21);%0.9381;
 MASS=y(22);%1.076;
 
% eps=y(23);%1;
 %oszt=y(24);%8;
 
   %M=y(23);%1.4;
   %CP=y(24);%0.037;
   %CP2=y(25);%0.046;
   %TF=y(26);%0.13;
   
  Wnt=y(23);
  Dest=y(24); %destruction complex cc
  Bcat=y(25); %beta_catenin cc
M=(sin(2*pi/24*t))/2+0.5;
%no coupling 
%kw4(1)=2;
%kw4(2)=0;

%weak coupling 
%kw4(1)=0;
%kw4(2)=0.25;

%strong coupling 
kw4_1=0;
kw4_2=3;

%parameters for clock-cell cycle coupling via Wee1

%no coupling 
kw5_1=1;
kw5_2=0;
kw6=1;

%weak coupling 
%kw5_1=1.,kw5_2=0.25,kw6=1.

%strong coupling 
%kw5(1)=0.25;
%kw5(2)=2;
%kw6=1;

%parameters for Wee1P
kw2_1=0.2;
kw2_2=2;
Jw2=0.2;
kw1=0.4;
Jw1=0.2;
kwd=1;

%parameters for the Cdc25a
kc3_1=0.1;
kc3_2=1;
Jc3=0.05;
Jc4=0.05;
kc4=0.4;

%parameters for CYCB
K1_1=0.1;
K1=0.6;
J1=0.1;
K2_1=0.05;
K2=20;
K2_2=1;

%parameters for the CYCBP
kwee1_1=0.08;
kwee1_2=10;
kcdc25_1=0.05;
kcdc25_2=10;

%parameters for ERG
k15=0.25;
k16=0.25;
J15=0.1;

%parameters for DRG
k17_1=0.35;
k17=10;
J17=0.3;
k18=10;

%parameters for CycD
%K9=2.5;
K10=5;
k24=1000;
k24r=10;

%parameters for CycE
k7_1=0;
K7=0.6;
k8_1=0.1;
K8=2;
K25=1000;
K25R=10;
J8=0.1;
YE=1;
YB=0.05;

%parameters for CycA
K29=0.05;
K30=20;

% parameters for p27
K5=20;
k6_1=10;
K6=100;
HE=0.5;
HB=1;
HA=0.5;

% parameters for Rb
RBT=10;
LD=3.3;
LE=5;
LB=5;
LA=3;
K20=10;
K19_1=0;
K19=20;

% parameters for PP1A
K21=1;
PP1T=1;
FE=25;
FB=2;

% parameters for Cdh1
k3_1=7.5;
K3=140;
J3=0.01;
J4=0.01;
K4=40;
GE=0;
GB=1;
GA=0.3;

% parameters for PPX
K33=0.05;
K34=0.05;

% parameters for IEP activation/inactivation
K31=0.7;
K32=1.8;
J31=0.01;
J32=0.01;

% parameters for Cdc20 activation/inactivation
K13=5;
K14=2.5;
J13=0.005;
J14=0.005;

% parameters for Cdc20T synthesis/degradation
K11_1=0;
K11=1.5;
K12=1.5;

% parameters for E2F
E2FT=5;
K22=1;
K23_1=0.005;
K23=1;

% parameters for E2F:Rb
K26=10000;
K26R=200;

% parameters for GM
k27=0.2;
K28=0.2;

% specific growth rate
%MU=0.074;
%MU=0.0192; %36 hours link time
%MU= 0.0288;%24 hours link time
%MU=0.0433;%16 hours link time
MU=0.0577;%12 hours link time

% parameter for the time of growth factor deprivation
%DEPRIV=-10;

%CYCDT=CYCD+CD; %--> cseréljük ki inkább cycdt-t az összegre
%CYCET=CYCE+CE;



% max rate for Cdh1 inactivation
V4=K4*(GE*CYCE+GA*CYCA+GB*CYCB);

% rate of p27 degradation
V6=k6_1+K6*(HE*CYCE+HA*CYCA+HB*CYCB);

% rate of CycE degradation
V8=k8_1+K8*(YE*(CYCE+CYCA)+YB*CYCB)/(J8+(CYCE+CE));

% calculation of active form of PP1
PP1A=PP1T/(K21*(FE*(CYCE+CYCA)+FB*CYCB)+1);

% calculation of hypo-phosphorylated form of Rb
RBH=RBT/(K20*(LD*(CYCD+CD)+LE*CYCE+LA*CYCA+LB*CYCB)/(K19_1*(PP1T-PP1A)+K19*PP1A)+1);

% calculation of E2F and Rb complex
L = (K26R+K20*(LD*(CYCD+CD)+LE*CYCE+LA*CYCA+LB*CYCB))/K26;
E2RBC  = 2*E2FT*RBH/(E2FT+RBH+L+sqrt((E2FT+RBH+L)^2-4*E2FT*RBH));

% calculation of active form of E2F
E2FA=(E2FT-E2RBC)*E2F/E2FT;

% rate of CYCB degradation
V2=K2_1*(1-Cdh1) + K2*Cdh1 + K2_2*Cdc20;

%aux 
%CYCET=CYCE+CE;
%aux 
%(CYCD+CD)=CYCD+CD;
%aux 
%CYCAT=CYCA+CA;
%aux
%P27T=P27+CD+CE+CA;
%aux RBH=RBH



% simplified circadian rhythm model

% Scaling factor
kt=1;

% parameters for the messenger RNA of the clock protein
n=2;
J=0.3; 
kms=1; 
kmd=0.1;

% parameters for the clock protein
kcps=0.5; 
kcpd=0.525;
ka=100;
kd=0.01;
kp1=10;

% parameters for the dimer clock protein
kcp2d=0.0525; 
kicd=0.01; 
kica=20;

% parameters for the transcription factor of M
kp2=0.1;
Jp=0.05;

%TFtot=0.5;
%IC= TFtot - TF;

eps=1;
%inactive complex of CP2 and TF
%IC= TFtot - TF;

%aux IC=IC


%CYCB'
dydt(1)= eps*(K1_1+K1*(CYCB/J1)^2/(1+(CYCB/J1)^2))*MASS - V2*CYCB - (kwee1_1+kwee1_2*Wee1)*CYCB + (kcdc25_1+kcdc25_2*Cdc25a)*CYCBP; 


%CYCBP'
dydt(2)= (kwee1_1+kwee1_2*Wee1)*CYCB - (kcdc25_1+kcdc25_2*Cdc25a)*CYCBP - V2*CYCBP; 


%Wee1'
dydt(3)= (kw5_1+kw5_2*M)-kw6*Wee1 - ((kw2_1+kw2_2*CYCB)*Wee1)/(Jw2+Wee1) + (kw1*Wee1P)/(Jw1+Wee1P) ;


%Wee1P'
dydt(4)= ((kw2_1+kw2_2*CYCB)*Wee1)/(Jw2+Wee1) - (kw1*Wee1P)/(Jw1+Wee1P) - kwd*Wee1P ;


%Cdc25a' 
dydt(5)= ((kc3_1+kc3_2*CYCB)*(1-Cdc25a))/(Jc3+(1-Cdc25a)) - (kc4*Cdc25a)/(Jc4+Cdc25a);


%ERG'
dydt(6)= eps*k15/(1+(DRG/J15)^2) - k16*ERG; 


%DRG' 
dydt(7)= eps*(k17_1*ERG+k17*(DRG/J17)^2/(1+(DRG/J17)^2)) - k18*DRG ;

%CYCD'
dydt(8)= (kw4_1+kw4_2*(Bcat*0.5)) + V6*CD + k24r*CD - k24*CYCD*P27 - K10*CYCD ;
%atskalazas
%M should be replaced with beta_catenin

%CD'
dydt(9)= k24*CYCD*P27 - k24r*CD - V6*CD - K10*CD ;


%CYCE' 
dydt(10)= eps*(k7_1+K7*E2FA)*MASS - V8*CYCE - K25*CYCE*P27 + K25R*CE + V6*CE ;


%CE'
dydt(11)= K25*CYCE*P27 - K25R*CE - V8*CE - V6*CE ;


%CYCA' 
dydt(12)= eps*K29*E2FA - K30*Cdc20*CYCA - K25*CYCA*P27 + K25R*CA + V6*CA ;


%CA'
dydt(13)= K25*CYCA*P27 - K25R*CA - (V6+K30*Cdc20)*CA ;


%P27'
dydt(14)= eps*K5 + K10*CD + k24r*CD - V6*P27 - K25*P27*(CYCE+CYCA) - k24*CYCD*P27 + K25R*(CE+CA) + V8*CE + K30*Cdc20*CA ;


%Cdh1' 
dydt(15)= (k3_1+K3*Cdc20)*(1-Cdh1)/(J3+1-Cdh1) - V4*Cdh1/(J4+Cdh1);


%PPX'
dydt(16)= eps*K33 - K34*PPX ;

%IEP' 
dydt(17)= K31*CYCB*(1-IEP)/(J31+1-IEP) - K32*PPX*IEP/(J32+IEP) ;


%Cdc20'
dydt(18)= K13*IEP*(Cdc20i)/(J13+Cdc20i) - K14*Cdc20/(J14+Cdc20) - K12*Cdc20 ;


%Cdc20i' 
dydt(19)= eps*(K11_1+K11*CYCB) - K13*IEP*(Cdc20i)/(J13+Cdc20i) + K14*Cdc20/(J14+Cdc20) - K12*Cdc20i ;

Cdc20T = Cdc20 + Cdc20i;

%E2F'
dydt(20)= K22*(E2FT-E2F) - (K23_1+K23*(CYCA+CYCB))*E2F ;

%GM'
%dydt(23)= k27*MASS*HEAV(0.8-RBH/RBT) - K28*GM ;
n2=8;
dydt(21)= k27*MASS*((RBH/RBT)^n2/(0.8^n2 +  (RBH/RBT)^n2)) - K28*GM ;

%MASS'
dydt(22)= eps*MU*GM;

%messenger RNA of the clock proteins
%M'
%dydt(23)=kt*(kms*TF^n/(J^n+TF^n)-kmd*M);
%M=(sin(2*pi/24*time))/2+0.5;

%monomer clock proteins
%CP'
%dydt(24)=kt*(kcps*M-kcpd*CP-2*ka*CP^2+2*kd*CP2-kp1*CP/(Jp+CP+2*CP2+2*IC));

%dimer clock protein
%CP2'
%dydt(25)=kt*(ka*CP^2-kd*CP2-kcp2d*CP2+kicd*IC-kica*CP2*TF-kp2*CP2/(Jp+CP+2*CP2+2*IC));

%transcription factor of M
%TF'
%dydt(26)=kt*(kcp2d*IC+kicd*IC-kica*TF*CP2+kp2*IC/(Jp+CP+2*CP2+2*IC));

%eps'= 0;
%Cdc20T = Cdc20 + Cdc20i;

%dWnt/dt
dydt(23)=ksynWnt*M-kdegWnt*Wnt;
%Wnt=0.2; % when Wnt is consitutive with 0.2 (minimum)
%Wnt=0.8; % Wnt constitutive with 0.8 (medium)
%Wnt=1.4; % consitutive with 1.4 (high)
%dDest/dt
dydt(24)=ksynDest-(kdegD1+kdegD2*Wnt)*Dest;

%dBcat/dt
dydt(25)=ksynB*Wnt-(kdegB1+kdegB2*Dest)*Bcat;


dydt=dydt';


end
