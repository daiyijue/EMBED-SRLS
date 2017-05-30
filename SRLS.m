function result=SRLS(dis,all_nodes)

%     load coordinates.mat;
    warning('off');
   
    A=zeros(all_nodes.num,all_nodes.dim+1);b=zeros(all_nodes.num,1);
    I=eye(all_nodes.dim);
    D=[I zeros(all_nodes.dim,1);zeros(1,all_nodes.dim) 0];
    f=[zeros(all_nodes.dim,1);-0.5];
    A=[-2*all_nodes.loc,ones(all_nodes.num,1)];
    for bh=1:all_nodes.num 
        b(bh,:)=dis(bh)^2-norm(all_nodes.loc(bh,:))^2;
    end
    AT=A';
    save Am.mat A;save Dm.mat D;save bm.mat b;save fm.mat f;
    g=eig((AT*A)^(-0.5)*D*(AT*A)^(-0.5));%求所有的特征值
    a=-1/max(g);
%     b1=a+1000;
%     while fun(b1)>0
%         b1=b1+1000;
%     end
    [lamx,err,yc,k]=bisect1(@fun1,a,10000);   %用二分法求namuda              
    A1(1,:)=(inv(AT*A+lamx*D)*(AT*b-lamx*f))';%SRLS用子问题得到的解
    result=A1(1,1:all_nodes.dim);
    
    
function [c,err,yc,k]=bisect1(f,a,b,delta)
%f是要求解的函数?
%a和b分别是有根区间的左右限
%delta是允许的误差界?
%c为所求近似解?
%yc为函数f在c的误差估计?
if nargin<4 
    delta=1e-5;
end

ya=fun1(a);
yb=fun1(b);
if yb==0
    c=b;
    return;
end
if ya*yb>0
    disp('(a,b)不是有根区间');
    return;
end
max1=1+round((log(b-a)-log(delta))/log(2));
for k=1:max1
    c=(a+b)/2;
%     yc=feval(f,c);
    yc=fun1(c);

    if yc==0
        a=c;
        b=c;
    break;
    elseif yb*yc>0
        b=c;yb=yc;
    else
        a=c;ya=yc;
    end
    if (b-a)<delta
        return;
    end

c=(a+b)/2;
err=abs(b-a);
% yc=feval(f,c);
yc=fun1(c);
end
 
function fam=fun1(lam)
load Am.mat;
load Dm.mat;
load bm.mat;
load fm.mat;
AT=A';
%y=inv(AT*A+lam*D)*(AT*b-lam*f);
y=(AT*A+lam*D)\(AT*b-lam*f);
yt=y';
ft=f';
fam=yt*D*y+2*ft*y;
