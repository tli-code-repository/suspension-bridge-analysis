function v=pic1(x,x2,c)
global L k;
y4=4*k*(x2/L-x2.^2/L^2);
y3=(26381*x)/69250-(26381*x.^2)/95911250;
n=length(c)/3;
% n=length(c)/2;
v=0;
for i=1:n
   v=v+c(i)*sin(i*pi*x/1385);
end
v=v+y4-y3;
end