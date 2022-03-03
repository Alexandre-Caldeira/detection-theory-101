rng('default')  % For reproducibility
% y = exprnd(5,100,1);
y = 5+2*randn(1e6,1);
m = bootstrp(100,@mean,y');
[fi,xi] = ksdensity(m);
plot(xi,fi)