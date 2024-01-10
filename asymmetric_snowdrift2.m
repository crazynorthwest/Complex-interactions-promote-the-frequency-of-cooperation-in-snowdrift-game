clear
clc



% scale=100;
% x1 = linspace(0, 1, 10);
% x2 = linspace(0, 1, 10);
% 
%     [X1, X2] = meshgrid(x1, x2);
%     u = zeros(size(X1));
%     v = zeros(size(X2));
% 
%     t = 0; 
%     for i = 1:numel(X1)
%         X_DOT = rigid1(t,[X1(i); X2(i)]);
%         Vmod = sqrt(X_DOT(1)^2 + X_DOT(2)^2);
%         u(i) = X_DOT(1)/Vmod;
%         v(i) = X_DOT(2)/Vmod;
%     end
% 
%     % Drawing
% %     h = figure;
%     hq = quiver(X1, X2, u, v, 'r'); 
%     hq.AutoScaleFactor = scale;
%     hold on;
%     xlabel('$x$', 'interpreter', 'latex')
%     ylabel('$y$', 'interpreter', 'latex')
%     axis tight equal;
%     xlim([0 1]);
%     ylim([0 1]);
%     hold on
% 
% s=5;
% 
% for ss=1:s
% for sss=1:s
% x0=[0.2*ss 0.2*sss]
% [t,x]=ode45(@rigid1, [0 10000], x0);
% plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
% hold on
% end
% end



% x0=[0.01 1]
% [t,x]=ode45(@rigid1, [0 10000], x0);
% plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
% hold on


b=1.2;
c=1;

figure(1)
p=linspace(0,1,100);
k=linspace(0,1,100);
x=(b./k-c)./(b./k-c/2);
plot((b-c)/(b-c/2),k,'b-')
hold on
plot(x,k,'r')
xlim([0 1]);
ylim([0 1]);



figure(2)
x0=[0.01 0.2]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
hold on
x0=[0.01 0.4]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.2 0.99]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.6 0.99]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.99]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.9 0.99]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.95]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.3]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.1]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.04]
[t,x]=ode45(@rigid1, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
b=1.2;
c=1;
k=0.5;
p=0.2 ;
y=((b-c)-p*(b-c/2))/((b-c/2)*(1-p));
plot(1, y, 'om');
xlim([0 1]);
ylim([0 1]);





figure(3)
x0=[0.01 0.1]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
hold on
x0=[0.01 0.3]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.5]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.9]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.5 0.99]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.9 0.99]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.99]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.93]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.03]
[t,x]=ode45(@rigid2, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
plot(1, 0, 'om');
xlim([0 1]);
ylim([0 1]);




figure(4)
x0=[0.01 0.1]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
hold on
x0=[0.01 0.03]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.3]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.7]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.01 0.92]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.4 0.99]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.8 0.99]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.99]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.94 0.99]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.5]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
x0=[0.99 0.9]
[t,x]=ode45(@rigid3, [0 10000], x0);
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2);
b=1.2;
c=1;
k=0.5;
p=0.8 ;
y=(p*(b-c)+(1-p)*(b/k-c))/(p*(b-c/2));
plot(y, 0, 'om');








xlim([0 1]);
ylim([0 1]);