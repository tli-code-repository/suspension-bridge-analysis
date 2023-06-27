function d=sus(c)
 global L;syms x ;
 global EcAc EhAh ji Lh k p H2;
 n=length(c)/3;
 E1I1=1074234000;
 v=@(x) 0;w=@(x) 0;u=@(x) 0;
 diff4v=@(x) 0; diffw=@(x) 0;diffu=@(x) 0;
 diff2u=@(x) 0; diff2w=@(x) 0;
 for i=1:n
     v=@(x) v(x)+c(i)*sin(i*pi*x/L);
     w=@(x) w(x)+c(i+n)*sin(i*pi*x/L);
     u=@(x) u(x)+c(i+2*n)*sin(i*pi*x/L);
     diff4v=@(x) diff4v(x)+c(i)*(pi*i/L)^4*sin(i*pi*x/L);
     diffw=@(x)  diffw(x)+c(i+n)*(pi*i/L)*cos(i*pi*x/L);
     diffu=@(x) diffu(x)+c(i+2*n)*(pi*i/L)*cos(i*pi*x/L);
     diff2w=@(x) diff2w(x)-(c(i+n)*(i*pi/L)^2*sin(i*pi*x/L));
     diff2u=@(x) diff2u(x)-(c(i+2*n)*(i*pi/L)^2*sin(i*pi*x/L));
 end
 ys=@(x) sqrt(1+(4*k*(1/L-2*x/L^2)).^2);
 y1=@(x) 4*k*(1/L-2*x/L^2);
 y2=-8*k/L^2;
 diffys=@(x) y2*y1(x)./ys(x);
 sh=@(x) v(x)-w(x);
 kh=@(x) EhAh./(ji*(Lh-4*k*(x/L-x.^2/L^2)));
 es=@(x) kh(x).*sh(x); 
 e=@(x) sqrt((1+diffu(x)).^2+(4*k*(1/L-2*x./L^2)+diffw(x)).^2);
 T0=@(x) H2*ys(x);
 T1=@(x) EcAc*(e(x)./ys(x)-1);
 diffe=@(x) (diff2u(x)+diffu(x).*diff2u(x)+(y1(x)+diffw(x)).*(y2+diff2w(x)))./e(x);
 diffT0=@(x) H2.*diffys(x);
 diffT1=@(x) EcAc*(diffe(x).*ys(x)-diffys(x).*e(x))./(ys(x).^2);
 ab=@(x) (diffT0(x)+diffT1(x))./e(x)-diffe(x).*(T0(x)+T1(x))./(e(x).^2);
 de=@(x) (e(x).*(diffT0(x).*diffw(x)+diff2w(x).*T0(x))-diffe(x).*T0(x).*diffw(x))./(e(x).^2)+...
     (y1(x)+diffw(x)).*(diffT1(x).*e(x)-diffe(x).*T1(x))./(e(x).^2)+T1(x).*(y2+diff2w(x))./e(x);
 Q=@(x) de(x)+es(x);
 M=@(x) ab(x).*(1+diffu(x))+diff2u(x).*(T0(x)+T1(x))./e(x);
 p2=@(x)   50*(x<=L/2   );
 K=@(x) E1I1*diff4v(x)+es(x)-p2(x);  %p为荷载
 x1=[0:1:L];
 s1=K(x1);
 s2=Q(x1);
 s3=M(x1);
 for i=1:n
%       d(i)=integral(@(x) sin(i*pi*x./L).*K(x),0,L,'AbsTol',1e-5);
%       d(i+n)=integral(@(x) sin(i*pi*x./L).*Q(x),0,L,'AbsTol',1e-5);
%      d(i+2*n)=integral(@(x) sin(i*pi*x./L).*M(x),0,L,'AbsTol',1e-5);
  d(i)=int(0,L,sin(i*pi*x1./L).*s1);
  d(i+n)=int(0,L,sin(i*pi*x1./L).*s2); 
  d(i+2*n)=int(0,L,sin(i*pi*x1./L).*s3);
 end
%  d(3*n+1)=u(0)-delt1;d(3*n+2)=u(L)-delt2;
end
%%
function y=int(a,b,s)
% a:start pt; b: end pt
% s: series to be integrated
n=length(s)-1;  
h=(b-a)/n;
y=3*(s(1)+s(n+1))+6*sum(s(2:1:n));
y=y*h/6;
end