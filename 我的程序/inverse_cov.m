function [X Y iter positivity_fail]=inverse_cov(A,lambda,rho,positive)
% A: sample covariance
% lambda: the tuning parameter for the sparsity, i.e., the weight of the l1 norm 
% rho: tuning parameter in ADMM. Theoretically, the algorithm always converges tp the same solution for arbitrary chosen rho>0; however, empirically the algorithm converges faster with appropriately chosen rho
% positive: 1 if we requires the output to be positive definite; 0 if there is no such requirement
% output X and Y are the estimated precision matrix.  X should be the same as Y (a property of ADMM algorithm). 
% iter: number of iterations 
maxIter=1000;
D=size(A,1);
%if init==1
X=inv(diag(diag(A))+eye(D));    %将A对角元素组成向量a拉成对角矩阵+A的行数构成的单位阵为b  然后取%
%else
%X=diag((diag(A)).^-1);
%end

oldX=X+0.01*eye(D); %%
Y=X;
oldY=Y+0.01*eye(D); 
Z=zeros(D,D);
iter=0;
[U,S]=svd(A);%svd是奇异值分解函数  U特征向量(正交单位)  c特征值
temp=rho/2+diag(S)/2;
temp1=repmat(temp,1,D)+repmat(temp,1,D)';
%allX{1}=1;allobjX=zeros(1,maxIter);allobjY=zeros(1,maxIter);
while iter<maxIter & (mod(iter,20)~=0 || norm(X-oldX,'fro')/max([1,norm(X,'fro'),norm(oldX,'fro')])>10^-7 || norm(Y-oldY,'fro')/max([1,norm(Y,'fro'),norm(oldY,'fro')])>10^-7)
oldX=X;
oldY=Y;
oldZ=Z;
X=U*((U'*(eye(D)-Z+rho*Y)*U)./temp1)*U';
Y=sign(Z/rho+X).*max(abs(Z/rho+X)-(lambda/rho)+(lambda/rho)*eye(D),0);
Z=Z+rho*(X-Y);
iter=iter+1;
% if mod(iter,100)==0
% X(1:4,1:4)
% end
end
[S1]=eig(Y);
positivity_fail=0;
if positive==1 & min(S1)<0
'positivity_fail'
positivity_fail=1;
oldX=X-1;
iter=1;
while iter<maxIter & (mod(iter,20)~=0 || norm(X-Y,'fro')>10^-5)
oldX=X;
oldY=Y;
oldZ=Z;
X=solve_G(A+rho*eye(D),eye(D)+rho*Y-Z,U,diag(S)+rho,0);
Y=shrink_G(Z/rho+X,lambda/rho,1);
Z=Z+rho*(X-Y);
iter=iter+1;
end

end