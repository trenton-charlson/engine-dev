function controlscrap(data)
function dydt = biteme(t,y)
u = -(1*y(1) + .5*y(1)*t + 2.1596*y(2));
c1 = 5;
c2=1000;
c3=7;
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2)=c1*u^2-c2*cosd(y(1))-c3*(y(2));
end
A = [0,1,0;0,-1000,0;1,0,0];
B = [0;7;0];
q = [1,0,0;0,1,0;0,0,1];
R = 1;
gains = lqr(A,B,q,R);
time=.001:.001:.001*length(data{:,1});
[t,y]=ode45(@biteme,[0 1.4],[-57,0]);
% [t,x]=ode45(@(t,x)biteme(t,x,1),[0 .45],[0,-50]);
% [t,a]=ode45(@(t,a)biteme(t,a,.5),[0 .45],[0,-50]);
% [t,b]=ode45(@(t,b)biteme(t,b,.75),[0 .45],[0,-50]);
% plot(t,y(:,1),'-',data{:,1},data{:,2},'o',t,x(:,1),t,a(:,1),t,b(:,1))
plot(t,y(:,1),'-',time,data{:,1})
legend('model','experimental')
% plot(time,data{:,1})
% legend('experimental')
title('Drop Test')
xlabel('Time(s)')
ylabel('angle(theta)')
end