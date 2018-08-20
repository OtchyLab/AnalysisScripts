i=1;
a=0;

for m = 1:2529
    if ttx_j6(10,m)==1
        a(i) = ttx_j6(3,m);
        i=i+1;
    end
end
x = 0:20:2000;
h=hist(a,x);
figure; plot(h);