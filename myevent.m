function [value,isterminal,direction] = myevent(t,x)
value=x(1)-0.2;
isterminal=1;
direction=-1;
end

%GLOBAL -1 {CYCB-oszt} {mass=0.5*mass; oszt=-8; GM=0.5*GM}
%GLOBAL -1 {CYCB-(kez+0.1)} {oszt=kez}
%oszt' = 0

%position:0.2
%direction: csökken

%The following two statements implement growth factor deprivation for DEPRIV < t < DEPRIV+10
%when depriv falls below t, set epsilon 0.5
%GLOBAL +1 {t-DEPRIV} {eps=0.5}

%when depriv-10 falls below t, set epsilon 1
%GLOBAL +1 {t-DEPRIV-10} {eps=1}