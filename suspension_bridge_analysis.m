tic
clear all;clc;
format short; 
global L;syms x;
global EcAc EhAh ji Lh k p H0 H2 delt3;
L=1385;
n=30;
% EcAc=103730000; H0=563167;
EcAc=103730000; H0=561754;
ji=16;EhAh=3280000/2;Lh=143.376;
k=131.905;p=500; limir=50;
a1=(8*k^2/L^2+1-H0/EcAc-16*k^2*H0/(EcAc*L^2));
a2=(32*k^2*H0/(EcAc*L^3)-16*k^2/L^3);
a3=(32*k^2/(3*L^4)-64*k^2*H0/(3*EcAc*L^4));

 
% delt1=3.257;   %500kn
% delt1=2.638;   %400kn
% delt1=1.998;   %300kn
% delt1=1.3467;  %200kn
% delt1=0.6786;  %100kn
delt1=0.3402;    %50kn 

 %delt1 is the total bend of the two towers, and it should interatively
 %determined in section 1 and 2 of this code

a=8*k^2/(3*L)+L-H0*L/EcAc-16*k^2*H0/(3*EcAc*L);
L2=L-delt1;
m=8/(3*L2);
l=H0*k*L2^3/(EcAc*L^2);
s=16*H0*k*L2/(3*EcAc*L^2);
f2=0.1*L2;
for ii=1:10
err1=a-(m*f2^2-s*f2+L2-l/f2) 
df2=-2*m*f2+s-l/f2^2;
if abs(err1)<1e-5
break
end
f2=f2-err1/df2;
end
H2=H0*k*L2^2/(f2*L^2);
k=f2;
L=L2;
 %%
c=rand(3*n,1); delt=0.01;d=sus(c);
for i=1:limir
for j=1:3*n
    c(j)=c(j)+delt;
    dv=sus(c);
    A(:,j)=d-dv;
    c(j)=c(j)-delt;
end 
    c=c+A\d'*delt;
    d=sus(c);
    err=max(abs(d))
if  err<0.00001
    break
end
end 
t=0:L/2:L;
g=pic4(t,c);
% h=@(t) 4*k*(1/L-2*t/L^2);
% ys=@(t) sqrt(1+(4*k*(1/L-2*t/L^2)).^2);
% e=@(t) sqrt((1+b).^2+(h(t)+z).^2);
% T1=@(t) EcAc*(e(t)./ys(t)-1);
% plot(t,g);
K1=231381.77;K2=239544.28;
delt2=g(2)/K1+g(2)/K2;
err3=delt2-delt1;
a4=(8*k^2/L^2+1-H2/EcAc-16*k^2*H2/(EcAc*L^2));
a5=(32*k^2*H2/(EcAc*L^3)-16*k^2/L^3);
a6=(32*k^2/(3*L^4)-64*k^2*H2/(3*EcAc*L^4));
delt3=g(2)/K1;
%%
q=0;
for jj=20.5:16:1364.5
      q=q+1;
x1=jj;
Ls0=a1*x1+a2*x1.^2+a3*x1.^3;
x2=x1;
for ii=1:10
Ls02=a4*x2+a5*x2.^2+a6*x2.^3;
err4=Ls0-Ls02;
if abs(err4)<1e-4
    break
end
derr=-(a4+2*a5*x2+3*a6*x2.^2);
x2=x2-err4/derr;
end
res1(q)=x2;
res(q)=x2+delt3-x1;
end
 
% subplot(2,2,1);
% f=0:0.1:1385;
f=20.5:16:1364.5;
y=pic1(f,res1,c);
z=pic2(f,res1,c);
b=pic3(f,c)+res;
r=max(abs(z-y));
r1=y.';
r2=z.';
r3=b.';
r4=f.';
plot(f,y,'r');
% subplot(2,2,2)
hold on
plot(f,z,'g');
% subplot(2,2,3);
plot(f,b,'b');
legend('v(x)','w(x)','u(x)')
toc
disp(['time elaspe:',num2str(toc)])







 

