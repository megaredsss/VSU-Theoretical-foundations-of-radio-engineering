S0=1; TOL=1E-4;tau_s=[3,8];t=[-20:1E-2:20];
s=@(t,tau_s)S0.*exp(-(t./tau_s).^4);
for i=1:length(t)
s1(i)=s(t(i),tau_s(1));s2(i)=s(t(i),tau_s(2));
end
figure
plot(t,s1,'k-','LineWidth',2.5);ylim([0,1]);xlim([-20,20]);
hold on;
plot(t,s2,'k--','LineWidth',2.5);ylim([0,1]);xlim([-20,20]);
hleg=legend('First','Second');
hold off;

Es=@(t1,t2,tau_s)(integral(@(t)2.*s(t,tau_s).^2,t1,t2,'ArrayValued',true));
DEs1=@(t)(Es(t/2,t,tau_s(1))/Es(0,t,tau_s(1))-TOL);
DEs2=@(t)(Es(t/2,t,tau_s(2))/Es(0,t,tau_s(2))-TOL);
T1=fzero(DEs1,tau_s(1));T2=fzero(DEs2,tau_s(2));
T=[T1,T2]

S=(@(w,tau_s)(integral(@(t)2.*s(t,tau_s).*cos(w.*t),0,inf,'ArrayValued',true)));
W=[-2.5:0.01:2.5];
for i=1:length(W)
    S1(i)=abs(S(W(i),tau_s(1)));S2(i)=abs(S(W(i),tau_s(2)));
end
mS1=max(S1);mS2=max(S2);
figure
plot(W,S1./mS1,'-k','LineWidth',2.5);
hold all;
plot(W,S2./mS2,'--k','LineWidth',2.5);
hleg=legend('First','Second');
hold off;

Ef=(@(W,tau_s)integral(@(w)abs(S(w,tau_s)).^2/pi,0,W,'ArrayValued',true));
DEf1=@(WW)(Ef(WW,tau_s(1))/Es(0,T1,tau_s(1))-0.95);
DEf2=@(WW)(Ef(WW,tau_s(2))/Es(0,T1,tau_s(2))-0.95);
w1=1/T1;
WW1=fzero(DEf1,w1)
w2=1/T2;
WW2=fzero(DEf2,w2)

s_tau=@(t,tau_s,tau_t)s(t-tau_t,tau_s);
S_tau=(@(w,tau_s,tau_t,T)(integral(@(t)s_tau(t,tau_s,tau_t).*exp(-1j.*w.*t),tau_t-T,tau_t+T)));
tau_t=10;
for i=1:length(W)
    ABS_S_tau1(i)=abs(S_tau(W(i),tau_s(1),tau_t,T(1)));
    ABS_S_tau2(i)=abs(S_tau(W(i),tau_s(2),tau_t,T(2)));
end;
for i=1:length(W)
    
    ARG_S_tau1(i)=angle(S_tau(W(i),tau_s(1),tau_t,T(1)));
    ARG_S_tau2(i)=angle(S_tau(W(i),tau_s(2),tau_t,T(2)));
end;
figure
plot(W,ABS_S_tau1,'-k',W,ABS_S_tau2,'--k','LineWidth',2.5);
hleg=legend('First','Second');
figure
plot(W,ARG_S_tau1,'-k',W,ARG_S_tau2,':ko');xlim([-1,1]);
hleg=legend('First','Second');

syms tt;
ss1=@(tt)S0.*exp(-(tt./tau_s(1)).^2);
    
    
    
    
    
    
    
    






