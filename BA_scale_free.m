clear
clc

s=50;
sss=3;
N=100;%最终节点数
k=0.1;%选择强度
T=1000;%代数
TT=50;%统计代数
TTTT=5;%一张图重复次数
TTTTT=5;%图数
c=1;
kkk=linspace(0.02,1,s);
% cc=linspace(0.02,1,s);
cc=[0.2,0.5,0.8];
bb=(c./cc+c)/2;

X=zeros(1,T);
X1=zeros(1,T);
X2=zeros(1,T);
Y=zeros(1,T);
Z=zeros(1,T);
XX=zeros(sss,s);
XX1=zeros(sss,s);
XX2=zeros(sss,s);
YY=zeros(sss,s);
ZZ=zeros(sss,s);
for ttttt=1:TTTTT

m0=4;%初始节点数
m=2;%每次引入新节点时生成的边数
%m应不大于m0
x=100*rand(1,m0);
y=100*rand(1,m0);
A=ones(m0);
A(1:m0+1:m0^2)=0;
for k=m0+1:N
    x(k)=100*rand;
    y(k)=100*rand;
    p=(sum(A)+1)/sum(sum(A)+1);
    pp=cumsum(p);%累加
    A(k,k)=0;
    ind=[];
    while length(ind)<m
        jj=find(pp>rand);
        jj=jj(1);
        ind=union(ind,jj);
    end
    A(k,ind)=1;
    A(ind,k)=1;
end
k=0.1;%选择强度
for tttt=1:TTTT

for ss=1:s
for ssss=1:sss
b=bb(ssss);
kk=kkk(ss);

C=ones(1,N);%实力集，0到1。
AA=sum(A);
[AAA,AAAA]=sort(AA);
CC=zeros(1,N);
CC(AAAA((fix((5*N)/10)+1):N))=1;
C=CC*(1/kk-1)+C;

D=zeros(1,N);%收益。

B=ones(1,N);%%策略阵。1为合作，0为背叛。
BB=rand(1,N)-1/2;
BB=sign(sign(BB)+1);
B=B-BB;

for t=1:T
X(t)=sum(B)/N;
X1(t)=sum(B.*CC)/(sum(CC));%strong
X2(t)=sum(B.*(1-CC))/(sum(1-CC));%weak
for tt=1:N   
i=fix(rand*N)+1;
iii=find(A(i,:));
jj=fix(rand*length(iii))+1;
j=iii(jj);
jjj=find(A(j,:));
DD=B(i)*((b*C(i)-c)*length(iii)+sum(B(iii))*c/2)+(1-B(i))*b*C(i)*sum(B(iii));
D=B(j)*((b*C(j)-c)*length(jjj)+sum(B(jjj))*c/2)+(1-B(j))*b*C(j)*sum(B(jjj));
d=D-DD;
dd=sign(sign(d)+1);
p=dd/(1+exp((-1)*d/k));
if p>rand
    B(i)=B(j);
end
end
end
XX(ssss,ss)=XX(ssss,ss)+sum(X((T-TT+1):T))/TT;
XX1(ssss,ss)=XX1(ssss,ss)+sum(X1((T-TT+1):T))/TT;
XX2(ssss,ss)=XX2(ssss,ss)+sum(X2((T-TT+1):T))/TT;
end
ss
end

for ss=1:s
for ssss=1:sss
b=bb(ssss);
kk=kkk(ss);
b=b/kk;
D=zeros(1,N);%收益。

B=ones(1,N);%%策略阵。1为合作，0为背叛。
BB=rand(1,N)-1/2;
BB=sign(sign(BB)+1);
B=B-BB;

for t=1:T
Y(t)=sum(B)/N;
for tt=1:N   
i=fix(rand*N)+1;
iii=find(A(i,:));
jj=fix(rand*length(iii))+1;
j=iii(jj);
jjj=find(A(j,:));
DD=B(i)*((b-c)*length(iii)+sum(B(iii))*c/2)+(1-B(i))*b*sum(B(iii));
D=B(j)*((b-c)*length(jjj)+sum(B(jjj))*c/2)+(1-B(j))*b*sum(B(jjj));
d=D-DD;
dd=sign(sign(d)+1);
p=dd/(1+exp((-1)*d/k));
if p>rand
    B(i)=B(j);
end
end
end
YY(ssss,ss)=YY(ssss,ss)+sum(Y((T-TT+1):T))/TT;
end
end

for ss=1:s
for ssss=1:sss
b=bb(ssss);
kk=kkk(ss);
D=zeros(1,N);%收益。

B=ones(1,N);%%策略阵。1为合作，0为背叛。
BB=rand(1,N)-1/2;
BB=sign(sign(BB)+1);
B=B-BB;

for t=1:T
Z(t)=sum(B)/N;
for tt=1:N   
i=fix(rand*N)+1;
iii=find(A(i,:));
jj=fix(rand*length(iii))+1;
j=iii(jj);
jjj=find(A(j,:));
DD=B(i)*((b-c)*length(iii)+sum(B(iii))*c/2)+(1-B(i))*b*sum(B(iii));
D=B(j)*((b-c)*length(jjj)+sum(B(jjj))*c/2)+(1-B(j))*b*sum(B(jjj));
d=D-DD;
dd=sign(sign(d)+1);
p=dd/(1+exp((-1)*d/k));
if p>rand
    B(i)=B(j);
end
end
end
ZZ(ssss,ss)=ZZ(ssss,ss)+sum(Z((T-TT+1):T))/TT;
end
end
tttt
end
ttttt
end
XXX=XX/(TTTT*TTTTT);%横k纵c
XXX1=XX1/(TTTT*TTTTT);
XXX2=XX2/(TTTT*TTTTT);
YYY=YY/(TTTT*TTTTT);
ZZZ=ZZ/(TTTT*TTTTT);

subplot(2,2,1)
plot(kkk,XXX1(2,:),'rpentagram')
hold on
plot(kkk,XXX2(2,:),'ksquare')
plot(kkk,YYY(2,:),'g^')
plot(kkk,ZZZ(2,:),'co')
subplot(2,2,2)
plot(kkk,XXX1(3,:),'rpentagram')
hold on
plot(kkk,XXX2(3,:),'ksquare')
plot(kkk,YYY(3,:),'g^')
plot(kkk,ZZZ(3,:),'co')


% BA
