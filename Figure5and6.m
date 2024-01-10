load('C:\Users\jxhan\Desktop\论文\实力不对称\实力不对称3\作图（final2)\程序\figure11.mat')
figure(1)
plot(cc,XX,'r')
hold on
plot(cc,YY(:,2),'k')

load('C:\Users\jxhan\Desktop\论文\实力不对称\实力不对称3\作图（final2)\程序\figure21.mat')
figure(2)
subplot(1,2,1)
plot(kkk,XX1(1,:),'rpentagram')
hold on
plot(kkk,XX2(1,:),'ko')
plot(kkk,ZZ(1,:),'g^')


subplot(1,2,2)
plot(kkk,XX1(2,:),'rpentagram')
hold on
plot(kkk,XX2(2,:),'ko')
plot(kkk,ZZ(2,:),'g^')

