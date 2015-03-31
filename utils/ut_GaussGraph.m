% Unit test for Gaussian Graph

clear all;
close all;
p = 6;
d = 3;

G = GaussGraph(p, d);
G.printGraph();
G.invSigma, 

n = 5000;
X = G.sample(n);
S = X'*X /n; % sample covariance
S, inv(G.invSigma),

