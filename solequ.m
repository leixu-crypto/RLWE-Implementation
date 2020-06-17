function y=solequ(G,b,q)
%format rat
[n1,m1]=size(G);
[n2,m2]=size(b);
k=m1/n1;
%A=randi([0,17-1],3,4)
%d=det(A);
%b=[1;1;1];
y=zeros(m1,m2);
for j=1:m2
    for i=1:n1
        y((i-1)*k+1,j)=b(i,j);
    end 
end
end



function y=moni(x,q)
y=mod(x^(q-2),q);
end

% function y=moni(x,p)
%  if x<0
%     x=mod(x,p);
%      x=x+p;
%  end
%  for i=1:p-1
%     t=x*i;
%     if  mod(t,p)==1
%     y=i;
%     end    
%  end
%  if mod(x,p)==0
%      y=0;
%  end     
% end