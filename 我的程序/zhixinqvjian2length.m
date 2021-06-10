function  C=zhixinqvjian2length(Sigma,Theta,Theta0,n,alpha )%n是样本容量,alpha为置信水平,Theta0表示真正的精度矩阵
p=size(Sigma,1);
I=eye(p);
Gamma=(1/2)*(Sigma*Theta+Theta*Sigma);
Gamma_1=Gamma-I;
A=Gamma_1(:);
T=Theta(:)-A;
U=zeros(p);
for i=1:p
    for j=1:p
        U(i,j)=(1/4)*Sigma(i,i)*Theta(j,j)+(1/4)*Sigma(j,j)*Theta(i,i)+(1/4)*Sigma(i,j)*Theta(i,j)+(1/2)*I(i,j)+1/2;
    end
end
V=sqrt(U(:));
a=norminv(1-alpha/2,0,1);
b=length(V);
for i=1:b
    C1(i)=T(i)-a*V(i)/sqrt(n);
    C2(i)=T(i)+a*V(i)/sqrt(n);
end
G=[Theta0(:),C1',C2'];
C=zeros(b,1);
for i=1:b
    if G(i,1)>=G(i,2) & G(i,1)<=G(i,3)
        C(i)=G(i,3)-G(i,2);
    else
        C(i)=0;
    end
end

