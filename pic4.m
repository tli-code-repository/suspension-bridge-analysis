function T1=pic4(x,c)
global L k EcAc;
n=length(c)/3;
h=4*k*(1/L-2*x/L^2);
ys=sqrt(1+(4*k*(1/L-2*x/L^2)).^2);
diffw=0;diffu=0;
for i=1:n
    diffw=diffw+c(i+n)*(pi*i/L)*cos(i*pi*x/L);
    diffu=diffu+c(i+2*n)*(pi*i/L)*cos(i*pi*x/L);
end
e=sqrt((1+diffu).^2+(h+diffw).^2);
T1=EcAc*(e./ys-1);
end