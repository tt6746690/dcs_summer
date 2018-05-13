% display lattice structure of lcg 


m = 2048;
a = 64;
c = 1;
seed = 0;
U = LCG(a,c,m,seed,2048);
subplot(2,1,1)
plot(U(1:m-1), U(2:m), '.')
subplot(2,1,2)
plot(U(1:511), U(2:512), '.');


a = 1365;
c=1;
U = LCG(a,c,m,seed,2048);
figure
plot(U(1:m-1), U(2:m), '.');