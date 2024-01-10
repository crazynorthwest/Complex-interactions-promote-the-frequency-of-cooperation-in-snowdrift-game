clear
clc

c=1;
ss=50;
kkk=linspace(0.02,1,ss);
s=3;
% cc=linspace(0.02,1,s);
cc=[0.2,0.5,0.8];
bb=(c./cc+c)/2;

XX=zeros(s,ss);%横k纵b
XX1=XX;
XX2=XX;

k=0.1;%选择强度
N=20;
T=N*N;
TT=10;%重复
TTTT=500;%代数
TTT=20;%统计代数

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
AA(1,1)=A(N,N);
AA(1,N+2)=A(N,1);
AA(N+2,1)=A(1,N);
AA(N+2,N+2)=A(1,1);

At=A-AA(1:N,2:(N+1));
Ad=A-AA(3:(N+2),2:(N+1));
Al=A-AA(2:(N+1),1:N);
Ar=A-AA(2:(N+1),3:(N+2));
Atl=A-AA(1:N,1:N);
Atr=A-AA(1:N,3:(N+2));
Adl=A-AA(3:(N+2),1:N);
Adr=A-AA(3:(N+2),3:(N+2));

    


Bt=(b*C-At*c/2-c/2).*abs(At)+(1-abs(At)).*(A.*(C*b-c/2));
Bd=(b*C-Ad*c/2-c/2).*abs(Ad)+(1-abs(Ad)).*(A.*(C*b-c/2));
Bl=(b*C-Al*c/2-c/2).*abs(Al)+(1-abs(Al)).*(A.*(C*b-c/2));
Br=(b*C-Ar*c/2-c/2).*abs(Ar)+(1-abs(Ar)).*(A.*(C*b-c/2));
Btl=(b*C-Atl*c/2-c/2).*abs(Atl)+(1-abs(Atl)).*(A.*(C*b-c/2));
Btr=(b*C-Atr*c/2-c/2).*abs(Atr)+(1-abs(Atr)).*(A.*(C*b-c/2));
Bdl=(b*C-Adl*c/2-c/2).*abs(Adl)+(1-abs(Adl)).*(A.*(C*b-c/2));
Bdr=(b*C-Adr*c/2-c/2).*abs(Adr)+(1-abs(Adr)).*(A.*(C*b-c/2));
B=Bt+Bd+Br+Bl+Btl+Btr+Bdl+Bdr;
    

BB(2:(N+1),2:(N+1))=B;
BB(1,2:(N+1))=B(N,:);
BB((N+2),2:(N+1))=B(1,:);
BB(2:(N+1),1)=B(:,N);
BB(2:(N+1),N+2)=B(:,1);
BB(1,1)=B(N,N);
BB(1,N+2)=B(N,1);
BB(N+2,1)=B(1,N);
BB(N+2,N+2)=B(1,1);

    ii=fix(rand*N)+1;
    jj=fix(rand*N)+1;
    p=fix(rand*8)+1;
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
    elseif p==4
        d=BB(ii+1,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj+2);
        end
        elseif p==5
        d=BB(ii+2,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj+2);
        end
        elseif p==6
        d=BB(ii,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj+2);
        end
        elseif p==7
        d=BB(ii+2,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj);
        end
    else
        d=BB(ii,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj);
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
AA(1,1)=A(N,N);
AA(1,N+2)=A(N,1);
AA(N+2,1)=A(1,N);
AA(N+2,N+2)=A(1,1);

At=A-AA(1:N,2:(N+1));
Ad=A-AA(3:(N+2),2:(N+1));
Al=A-AA(2:(N+1),1:N);
Ar=A-AA(2:(N+1),3:(N+2));
Atl=A-AA(1:N,1:N);
Atr=A-AA(1:N,3:(N+2));
Adl=A-AA(3:(N+2),1:N);
Adr=A-AA(3:(N+2),3:(N+2));

    


Bt=(b*C-At*c/2-c/2).*abs(At)+(1-abs(At)).*(A.*(C*b-c/2));
Bd=(b*C-Ad*c/2-c/2).*abs(Ad)+(1-abs(Ad)).*(A.*(C*b-c/2));
Bl=(b*C-Al*c/2-c/2).*abs(Al)+(1-abs(Al)).*(A.*(C*b-c/2));
Br=(b*C-Ar*c/2-c/2).*abs(Ar)+(1-abs(Ar)).*(A.*(C*b-c/2));
Btl=(b*C-Atl*c/2-c/2).*abs(Atl)+(1-abs(Atl)).*(A.*(C*b-c/2));
Btr=(b*C-Atr*c/2-c/2).*abs(Atr)+(1-abs(Atr)).*(A.*(C*b-c/2));
Bdl=(b*C-Adl*c/2-c/2).*abs(Adl)+(1-abs(Adl)).*(A.*(C*b-c/2));
Bdr=(b*C-Adr*c/2-c/2).*abs(Adr)+(1-abs(Adr)).*(A.*(C*b-c/2));
B=Bt+Bd+Br+Bl+Btl+Btr+Bdl+Bdr;
    

BB(2:(N+1),2:(N+1))=B;
BB(1,2:(N+1))=B(N,:);
BB((N+2),2:(N+1))=B(1,:);
BB(2:(N+1),1)=B(:,N);
BB(2:(N+1),N+2)=B(:,1);
BB(1,1)=B(N,N);
BB(1,N+2)=B(N,1);
BB(N+2,1)=B(1,N);
BB(N+2,N+2)=B(1,1);

    ii=fix(rand*N)+1;
    jj=fix(rand*N)+1;
    p=fix(rand*8)+1;
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
    elseif p==4
        d=BB(ii+1,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj+2);
        end
        elseif p==5
        d=BB(ii+2,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj+2);
        end
        elseif p==6
        d=BB(ii,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj+2);
        end
        elseif p==7
        d=BB(ii+2,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj);
        end
    else
        d=BB(ii,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj);
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


subplot(2,2,3)
plot(kkk,XX1(1,:),'rpentagram')
hold on
plot(kkk,XX2(1,:),'ksquare')
plot(kkk,YY(1,:),'g^')
plot(kkk,ZZ(1,:),'co')
subplot(2,2,4)
plot(kkk,XX1(3,:),'rpentagram')
hold on
plot(kkk,XX2(3,:),'ksquare')
plot(kkk,YY(3,:),'g^')
plot(kkk,ZZ(3,:),'co')







c=1;
ss=50;
kkk=linspace(0.02,1,ss);
s=3;
% cc=linspace(0.02,1,s);
cc=[0.2,0.5,0.8];
bb=(c./cc+c)/2;

XX=zeros(s,ss);%横k纵b
XX1=XX;
XX2=XX;

k=0.1;%选择强度
N=20;
T=N*N;
TT=4;%重复
TTTT=400;%代数
TTT=40;%统计代数

for sss=1:s
b=bb(sss);
for ssss=1:ss
kk=kkk(ssss);
X=zeros(TT,TTT);
for tt=1:TT
A=ones(N,N);%策略阵。1为合作，0为背叛。
AAA=rand(N,N)-1/2;
AAA=sign(sign(AAA)+1);
A=A-AAA;
AA=ones(N+2,N+2);

B=zeros(N,N);%收益阵。
BB=zeros(N+2,N+2);
for tttt=1:TTTT
for t=1:T 
AA(2:(N+1),2:(N+1))=A;
AA(1,2:(N+1))=A(N,:);
AA((N+2),2:(N+1))=A(1,:);
AA(2:(N+1),1)=A(:,N);
AA(2:(N+1),N+2)=A(:,1);   
AA(1,1)=A(N,N);
AA(1,N+2)=A(N,1);
AA(N+2,1)=A(1,N);
AA(N+2,N+2)=A(1,1);

At=A-AA(1:N,2:(N+1));
Ad=A-AA(3:(N+2),2:(N+1));
Al=A-AA(2:(N+1),1:N);
Ar=A-AA(2:(N+1),3:(N+2));
Atl=A-AA(1:N,1:N);
Atr=A-AA(1:N,3:(N+2));
Adl=A-AA(3:(N+2),1:N);
Adr=A-AA(3:(N+2),3:(N+2));

    


Bt=(b-At*c/2-c/2).*abs(At)+(1-abs(At)).*(A.*(b-c/2));
Bd=(b-Ad*c/2-c/2).*abs(Ad)+(1-abs(Ad)).*(A.*(b-c/2));
Bl=(b-Al*c/2-c/2).*abs(Al)+(1-abs(Al)).*(A.*(b-c/2));
Br=(b-Ar*c/2-c/2).*abs(Ar)+(1-abs(Ar)).*(A.*(b-c/2));
Btl=(b-Atl*c/2-c/2).*abs(Atl)+(1-abs(Atl)).*(A.*(b-c/2));
Btr=(b-Atr*c/2-c/2).*abs(Atr)+(1-abs(Atr)).*(A.*(b-c/2));
Bdl=(b-Adl*c/2-c/2).*abs(Adl)+(1-abs(Adl)).*(A.*(b-c/2));
Bdr=(b-Adr*c/2-c/2).*abs(Adr)+(1-abs(Adr)).*(A.*(b-c/2));
B=Bt+Bd+Br+Bl+Btl+Btr+Bdl+Bdr;
    

BB(2:(N+1),2:(N+1))=B;
BB(1,2:(N+1))=B(N,:);
BB((N+2),2:(N+1))=B(1,:);
BB(2:(N+1),1)=B(:,N);
BB(2:(N+1),N+2)=B(:,1);
BB(1,1)=B(N,N);
BB(1,N+2)=B(N,1);
BB(N+2,1)=B(1,N);
BB(N+2,N+2)=B(1,1);

    ii=fix(rand*N)+1;
    jj=fix(rand*N)+1;
    p=fix(rand*8)+1;
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
    elseif p==4
        d=BB(ii+1,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+1,jj+2);
        end
        elseif p==5
        d=BB(ii+2,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj+2);
        end
        elseif p==6
        d=BB(ii,jj+2)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj+2);
        end
        elseif p==7
        d=BB(ii+2,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii+2,jj);
        end
    else
        d=BB(ii,jj)-B(ii,jj);
        dd=sign(sign(d)+1);%理性
        pp=dd/(1+exp((-1)*d/k));
        if rand<pp
           A(ii,jj)=AA (ii,jj);
        end
    end  
end
ttt=tttt-TTTT+TTT;
if ttt>0
    X(tt,ttt)=sum(sum(A))/(N*N);
end
%tt
end
end
XX(sss,ssss)=sum(sum(X))/(TT*TTT);
ssss
end
sss
end