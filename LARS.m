function LARS(Y,X)
[m,n]=size(X);
mua=zeros(m,1);
c=X'*(Y-mua);
[cm,A]=max(abs(c));
nA=length(A);
while nA<n
  XA=X(:,A);
  GA=XA'*XA;
  AA=(ones(nA,1)'*GA^(-1)*ones(nA,1))^(-0.5);
  ua=XA*AA*GA^(-1)*ones(nA,1);
  a=X'*ua;
  s=setdiff(1:n,A);
  c1=c(s);
  a1=a(s);
  ga0=[(cm-c1)./(AA-a1);(cm+c1)./(AA+a1)];
  [ga,j]=min(ga0(find(ga0>0)));
  if j>length(c1)
      j=j-length(c1);
  end
  j=s(j);
  A=sort([A,j]);
  mua=mua+ga*ua;
  nA=length(A);
end
