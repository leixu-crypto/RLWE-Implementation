[A,S,W]=Basis(256,14,859);

function S=BasicS(n,k,q)
S=zeros(n*k,n*k);
S0=zeros(k,k);
for i=1:k-1
   S0(i+1,i)=-1;
   S0(i,i)=2;
end
t=q;
while t>=2
    l=floor(log2(t));
    S0(l+1,k)=1;
    t=t-2^l;
end
S0(1,k)=t;
for i=1:n
    S(k*(i-1)+1:k*i,k*(i-1)+1:k*i)=S0;
end
end

function G=PrimG(n,k,q)
G=zeros(n,n*k);
for i=1:k
    g(i)=mod(2^(i-1),q);
end
for i=1:n
    G(i,k*(i-1)+1:k*i)=g;
end
end

function [A,R]=GenTrap(n,k,q)
m1=n;
w=n*k;
A0=randi([0,q-1],n,m1);
r=randi([0,q-1],1,k);
for i=1:n
    R(i,k*(i-1)+1:k*i)=r;
end
H=eye(n,n);
G=PrimG(n,k,q);
A1=H*G-A0*R;
A=mod([A0,A1],q);
end

function [A,S,W]=Basis(n,k,q)
m1=n;
w=n*k;
m=w+m1;
S0=BasicS(n,k,q);
G=PrimG(n,k,q);
H=eye(n,n);
[A,R]=GenTrap(n,k,q);
I1=eye(m1,m1);
I2=eye(w,w);
Z1=zeros(w,m1);
I3=eye(m1,m1);
Z2=zeros(m1,w);
B=-H*A*[eye(m1,m1);zeros(w,m1)];
W=solequ(G,B,q);
S=mod([I1,R;Z1,I2]*[I3,Z2;W,S0],q);
end

function C=GSNorm(S)
[m,n]=size(S);
C=zeros(m,m);
[Q,R] = qr(S');
for i=1:m
    C(:,i)=R(i,i)*Q(i,:);
end
end
