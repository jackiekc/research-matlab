function size_new = get_size_power2(sizeimg) 
%[m,n]-size_new
%[x,y]-sizeimg
x=sizeimg(1);y=sizeimg(2);
m=ceil(log2(x));
n=ceil(log2(y));
size_new=[2^m,2^n];