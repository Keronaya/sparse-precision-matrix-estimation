A=zeros(1,10);
B=zeros(1,10);
c_new=0;
d_new=0;
e_new=0;
f_new=0;
for k=1:1:30
N=50;p=100;n=100;para=0.2;
C=zeros(p*p,N);
for i=1:N
    Theta0= tridiag(p,1,para,para);
    lambda=sqrt(log(p)/n);
    rho=1;
    positive=1;
    u=zeros(p,1)';
    R=inv(Theta0);
    Y=multivrandn(u,R,n);
    Sigma=yangbenxiefancha(Y);
    [X Y iter positivity_fail]=inverse_cov(Sigma,lambda,rho,positive);
    Theta=X;
    C(:,i)=zhixinqvjian2length(Sigma,Theta,Theta0,n,0.01);
end
F=Theta0(:);
S1=0;S2=0;
for i=1:p*p
    for j=1:N
        if F(i)~=0
            S1=S1+C(i,j);
        else 
            S2=S2+C(i,j);
        end
    end
end
l1=S1/((3*p-2)*N);
l2=S2/((p*p-3*p+2)*N);
A(k)=l1;
B(k)=l2;
end
for i=1:1:10
 c_new=A(i)+c_new;   
 d_new=B(i)+d_new;
end
for i=1:1:10
 e_new=(A(i)-c_new/10)^2+e_new;   
 f_new=(B(i)-d_new/10)^2+f_new;
end
c_new/10
d_new/10
a_new=e_new/10
b_new=f_new/10





        

