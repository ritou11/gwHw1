function histDiagPri(mpc)
N = size(mpc.bus,1);
[Y, ~, ~] = makeYbus(mpc);
dp = zeros(N,1);
for i=1:N
    l = abs(Y(:, i));
    d = full(l(i));
    l(i) = 0;
    m = full(max(l));
    dp(i) = d / m;
end
edges = 0.5:0.5:5;
histogram(dp, edges);
hc = histcounts(dp, edges);
xbins = (edges(1:end-1) + edges(2:end))/2;
hold on;
plot(xbins, hc, '.', 'MarkerSize', 20);
hold off;
end