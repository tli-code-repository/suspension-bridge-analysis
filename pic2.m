function w=pic2(x,x2,c)
global L k;
n=length(c)/3;
y4=4*k*(x2/L-x2.^2/L^2);
y3=(26381*x)/69250 - (26381*x.^2)/95911250;
w=0;
for i=1:n
   w=w+c(i+n)*sin(i*pi*x/1385);
end
w=w+y4-y3;
end
