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
ss1=@(tt)S0.*exp(-(tt./tau_s(1)).^4);
ss2=@(tt)S0.*exp(-(tt./tau_s(2)).^4);
s_d1=matlabFunction(tau_s(1).*diff(ss1(tt)));
s_d2=matlabFunction(tau_s(2).*diff(ss2(tt)));
for i=1:length(t)
    sd1(i)=s_d1(t(i));sd2(i)=s_d2(t(i));
end
figure
plot(t,sd1,'k-','LineWidth',2.5);
hold all;
plot(t,s1,'k--','LineWidth',2.5);
plot(t,sd2,'k-.','LineWidth',2.5);
plot(t,s2,'k:','LineWidth',2.5);
hleg=legend('First','Second','Third','Fourth');
hold off;

S_d1=(@(w,T)(integral(@(t)s_d1(t).*exp(-1j.*w.*t),-T,T)));
S_d2=(@(w,T)(integral(@(t)s_d2(t).*exp(-1j.*w.*t),-T,T)));
for i=1:length(W)
    ABS_S_d1(i)=abs(S_d1(W(i),T(1)));
    ABS_S_d2(i)=abs(S_d2(W(i),T(2)));
    ARG_S_d1(i)=angle(S_d1(W(i),T(1)));
    ARG_S_d2(i)=angle(S_d2(W(i),T(2)));
    ARG_S1(i)=angle(S(W(i),tau_s(1)));
    ARG_S2(i)=angle(S(W(i),tau_s(2)));
end
figure
plot(W,ABS_S_d1,'k-',W,ABS_S_d2,'k--',W,S1,'k:',W,S2,'k-.','LineWidth',2.5);
hleg=legend('First','Second','Third','Fourth');
figure
plot(W,ARG_S_d1,'k-',W,ARG_S1,'k--','LineWidth',2.5);
hleg=legend('First','Second');

W0=1.5;
s_W=@(t,tau_s)s(t,tau_s).*cos(W0.*t);
for i=1:length(t)
    s_W1(i)=s_W(t(i),tau_s(1));
    s_W2(i)=s_W(t(i),tau_s(2));
end;
figure
subplot(2,1,1);
plot(t,s_W1,'k-',t,s1,'k--','LineWidth',2.5);
hleg=legend('First','Second');
subplot(2,1,2);
plot(t,s_W2,'k-',t,s2,'k--','LineWidth',2.5);
hleg=legend('First','Second');

S_W=(@(w,tau_s,T)(integral(@(t)s_W(t,tau_s).*exp(-1j.*w.*t),-T,T)));
for i=1:length(W)
    ABS_S_W1(i)=abs(S_W(W(i),tau_s(1),T(1)));
    ABS_S_W2(i)=abs(S_W(W(i),tau_s(2),T(2)));
    ABS_S1(i)=abs(S(W(i),tau_s(1)));
    ABS_S2(i)=abs(S(W(i),tau_s(2)));
end
figure
subplot(2,1,1);
plot(W, ABS_S_W1,'k-',W,ABS_S1,'k--','LineWidth',2.5);
hleg=legend('First','Second');
figure
subplot(2,1,1);
plot(W, ABS_S_W2,'k-',W,ABS_S2,'k--','LineWidth',2.5);
hleg=legend('First','Second');
