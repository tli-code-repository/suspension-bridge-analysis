function u=pic3(x,c)
n=length(c)/3;
u=0;
for i=1:n
   u=u+c(i+2*n)*sin(i*pi*x/1385);
end
end