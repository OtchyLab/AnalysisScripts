function h = plotCorrelationDistribution(Rvalues, binSize)

[N] = histc(Rvalues,-1:binSize:1);
N = N / sum(N*binSize);
sum(N*binSize)
h = plot(-1:2/(length(N)-1):1,N);