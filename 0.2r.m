clear
clc

c=1;
ss=2;
% kkk=linspace(0.02,1,ss);
kkk=[0.2,1];
s=50;
cc=linspace(0.02,1,s);
bb=(c./cc+c)/2;

XX=zeros(s,ss);%横k纵b
XX1=XX;
XX2=XX;

k=0.1;%选择强度
N=20;
T=N*N;
TT=2;%重复
TTTT=500;%代数
TTT=50;%统计代数

for sss=1:s
b=bb(sss);
for ssss=1:ss
kk=kkk(ssss);
X=zeros(TT,TTT);
X1=X;%weak
X2=X;%strong

for tt=1:TT
A=ones(N,N);%策略阵。1为合作，0为背叛。
AAA=rand(N,N)-1/2;
AAA=sign(sign(AAA)+1);
A=A-AAA;
AA=ones(N+2,N+2);

C=ones(N,N);%实力集，0为弱1为强。
CCC=rand(N,N)-1/5;
CCC=sign(sign(CCC)+1);
C=(1/kk)*C-((1/kk)-1)*CCC;

B=zeros(N,N);%收益阵。
BB=zeros(N+2,N+2);
for tttt=1:TTTT
for t=1:T
    

    
 AA(2:(N+1),2:(N+1))=A;
AA(1,2:(N+1))=A(N,:);
AA((N+2),2:(N+1))=A(1,:);
AA(2:(N+1),1)=A(:,N);
AA(2:(N+1),N+2)=A(:,1);   


At=A-AA(1:N,2:(N+1));
Ad=A-AA(3:(N+2),2:(N+1));
Al=A-AA(2:(N+1),1:N);
Ar=A-AA(2:(N+1),3:(N+2));

    


Bt=(b*C-At*c/2-c/2).*abs(At)+(1-abs(At)).*(A.*(C*b-c/2));
Bd=(b*C-Ad*c/2-c/2).*abs(Ad)+(1-abs(Ad)).*(A.*(C*b-c/2));
Bl=(b*C-Al*c/2-c/2).*abs(Al)+(1-abs(Al)).*(A.*(C*b-c/2));
Br=(b*C-Ar*c/2-c/2).*abs(Ar)+(1-abs(Ar)).*(A.*(C*b-c/2));
B=Bt+Bd+Br+Bl;
    

BB(2:(N+1),2:(N+1))=B;
BB(1,2:(N+1))=B(N,:);
BB((N+2),2:(N+1))=B(1,:);
BB(2:(N+1),1)=B(:,N);
BB(2:(N+1),N+2)=B(:,1);

    ii=fix(rand*N)+1;
    jj=fix(rand*N)+1;
    p=fix(rand*4)+1;
    if p==1
        d=BB(ii,jj+1)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj+1);
        end
    elseif p==2
        d=BB(ii+2,jj+1)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj+1);
        end
    elseif p==3
        d=BB(ii+1,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj);
        end
    else
        d=BB(ii+1,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj+2);
        end
    end  
end
ttt=tttt-TTTT+TTT;
if ttt>0
    X1(tt,ttt)=(sum(sum(A.*CCC)))/(sum(sum(CCC)));%weak
    X2(tt,ttt)=(sum(sum(A.*(1-CCC))))/(N*N-sum(sum(CCC)));%strong
    X(tt,ttt)=sum(sum(A))/(N*N);
end
%tt
end
end
XX(sss,ssss)=sum(sum(X))/(TT*TTT);
XX1(sss,ssss)=sum(sum(X1))/(TT*TTT);
XX2(sss,ssss)=sum(sum(X2))/(TT*TTT);
ssss
end
sss
end


for sss=1:s
for ssss=1:ss
b=bb(sss);
kk=kkk(ssss);
b=b/kk;
kk=1;
Y=zeros(TT,TTT);

for tt=1:TT
A=ones(N,N);%策略阵。1为合作，0为背叛。
AAA=rand(N,N)-1/2;
AAA=sign(sign(AAA)+1);
A=A-AAA;
AA=ones(N+2,N+2);

C=ones(N,N);%实力集，0为弱1为强。
CCC=rand(N,N)-1/2;
CCC=sign(sign(CCC)+1);
C=(1/kk)*C-((1/kk)-1)*CCC;

B=zeros(N,N);%收益阵。
BB=zeros(N+2,N+2);
for tttt=1:TTTT
for t=1:T
AA(2:(N+1),2:(N+1))=A;
AA(1,2:(N+1))=A(N,:);
AA((N+2),2:(N+1))=A(1,:);
AA(2:(N+1),1)=A(:,N);
AA(2:(N+1),N+2)=A(:,1);   

At=A-AA(1:N,2:(N+1));
Ad=A-AA(3:(N+2),2:(N+1));
Al=A-AA(2:(N+1),1:N);
Ar=A-AA(2:(N+1),3:(N+2));

Bt=(b*C-At*c/2-c/2).*abs(At)+(1-abs(At)).*(A.*(C*b-c/2));
Bd=(b*C-Ad*c/2-c/2).*abs(Ad)+(1-abs(Ad)).*(A.*(C*b-c/2));
Bl=(b*C-Al*c/2-c/2).*abs(Al)+(1-abs(Al)).*(A.*(C*b-c/2));
Br=(b*C-Ar*c/2-c/2).*abs(Ar)+(1-abs(Ar)).*(A.*(C*b-c/2));
B=Bt+Bd+Br+Bl;
    
BB(2:(N+1),2:(N+1))=B;
BB(1,2:(N+1))=B(N,:);
BB((N+2),2:(N+1))=B(1,:);
BB(2:(N+1),1)=B(:,N);
BB(2:(N+1),N+2)=B(:,1);

ii=fix(rand*N)+1;
jj=fix(rand*N)+1;
p=fix(rand*4)+1;
    if p==1
        d=BB(ii,jj+1)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj+1);
        end
    elseif p==2
        d=BB(ii+2,jj+1)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj+1);
        end
    elseif p==3
        d=BB(ii+1,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj);
        end
    else
        d=BB(ii+1,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj+2);
        end
    end
    
end
ttt=tttt-TTTT+TTT;
if ttt>0
    Y(tt,ttt)=(sum(sum(A)))/(N*N);
end
end
end
YY(sss,ssss)=sum(sum(Y))/(TT*TTT);
ssss
end
sss
end
subplot(2,3,1)
plot(cc,XX1(:,1),'rpentagram')
hold on
plot(cc,XX2(:,1),'ksquare')
plot(cc,YY(:,1),'g^')
plot(cc,YY(:,2),'co')
XXX=(YY(:,1)+4*YY(:,2))/5;
subplot(2,3,2)
plot(cc,XX(:,1),'r')
hold on
plot(cc,XXX,'y')