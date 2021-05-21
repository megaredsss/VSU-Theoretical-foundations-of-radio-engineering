clear all; close all; clc;
%Задание 3.2
tau_K=[2,5];
K=(@(w,tau_K)1/(1+1j*w*tau_K));
W=[-5:1E-2:5];
for i=1:length(W)
    K1(i)=K(W(i),tau_K(1));K2(i)=K(W(i),tau_K(2));
end
figure
plot(W,abs(K1),'k-',W,abs(K2),'k--','LineWidth',2.5);
hleg=legend('First','Second');
figure
plot(W,angle(K1),'k-',W,angle(K2),'k--','LineWidth',2.5);
hleg=legend('First','Second');
h1=@(t)dirac(t);
h2=@(t,tau_K)(1./tau_K).*exp(-t./tau_K).*heaviside(t);
h=@(t,tau_K)h1(t)+h2(t,tau_K);
t=[-5:1E-2:20];
for i=1:length(t)
    hh1(i)=h(t(i),tau_K(1));hh2(i)=h(t(i),tau_K(2));
end
figure
plot(t,hh1,'k-',t,hh2,'k--','LineWidth',2.5);
hleg=legend('First','Second');
KK=(@(w,tau_K)(1+integral(@(t)h2(t,tau_K).*exp(-1j.*w*t),0,20)));
for i=1:length(W)
    KK1(i)=KK(W(i),tau_K(1));KK2(i)=KK(W(i),tau_K(2));
end
figure
plot(W,abs(KK1),'k-',W,abs(KK2),'k--','LineWidth',2.5);
hleg=legend('First','Second');
S0=1;
s_input=@(t,tau_s)S0.*exp(-(t./tau_s).^2).*heaviside(t);
tau_s=[3,8];
for i=1:length(t)
    s1(i)=s_input(t(i),tau_s(1)),s2(i)=s_input(t(i),tau_s(2));
end
figure
plot(t,s1,'k-',t,s2,'k--','LineWidth',2.5);
hleg=legend('First','Second');
S_input=(@(w,tau_s)(integral(@(t)s_input(t,tau_s).*exp(-1j.*w.*t),0,20)));
for i=1:length(W)
    S1_input(i)=abs(S_input(W(i),tau_s(1)));
    S2_input(i)=abs(S_input(W(i),tau_s(2)));
end
figure
plot(W,S1_input,'k-',W,S2_input,'k--','LineWidth',2.5);
hleg=legend('First','Second');
S_output=@(w,tau_s,tau_K)S_input(w,tau_s).*K(w,tau_K);
for i=1:length(W)
    S1_output(i)=abs(S_output(W(i),tau_s(1),tau_K(1)));
    S2_output(i)=abs(S_output(W(i),tau_s(1),tau_K(2)));
end
figure
plot(W,S1_output,'k-',W,S2_output,'k--','LineWidth',2.5);
hleg=legend('First','Second');
s1_output=@(t,tau_s)s_input(t,tau_s);
s2_output=@(t,tau_s,tau_K)integral(@(tt)s_input(tt,tau_s).*h2(t-tt,tau_K),0,t);
s_output=@(t,tau_s,tau_K)s1_output(t,tau_s)+s2_output(t,tau_s,tau_K);
for i=1:length(t)
    ss1_input(i)=s_input(t(i),tau_s(1));
    ss1_output(i)=s_output(t(i),tau_s(1),tau_K(1));
    ss2_output(i)=s_output(t(i),tau_s(1),tau_K(2));
end
figure
plot(t,ss1_input,'k-',t,ss1_output,'k--',t,ss2_output,'k:','LineWidth',2.5);
hleg=legend('First','Second','Third');


