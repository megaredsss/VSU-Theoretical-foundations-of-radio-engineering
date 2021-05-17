S0=1; TOL=1E-4;tau_s=[3,8];t=[-20:1E-2:20];
s=@(t,tau_s)S0.*exp(-(t./tau_s).^2);
for i=1:length(t)
s1(i)=s(t(i),tau_s(1));s2(i)=s(t(i),tau_s(2));
end
figure
plot(t,s1,'k-','LineWidth',2.5);ylim([0,1]);xlim([-20,20]);
hold on;
plot(t,s2,'k--','LineWidth',2.5);ylim([0,1]);xlim([-20,20]);
hleg=legend('First','Second');
hold off;