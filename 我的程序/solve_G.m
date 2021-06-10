function [X iter]=solve_G(A,B,UA,SigmaA,positive)
D=size(A,1);
C=2./(repmat(SigmaA,1,D)+(repmat(SigmaA,1,D))');
X=UA*((UA'*B*UA).*C)*UA';
if positive==1 & min(eig(X))<0
oldX=X-1;
iter=1;
hatA=max(SigmaA)*eye(D)-A;
while norm(oldX-X,'fro')>10^-5 
iter=iter+1;
oldX=X;
X=max(SigmaA)^-1*(B-(X*hatA+hatA*X)/2);
[UX,SX]=eig(X);
X=UX*max(SX,0)*UX';
end
end