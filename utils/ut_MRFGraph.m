% A unit test for MRF Graph
clc;
clear all;
close all;

p = 20;
d = 4;

G = MRFGraph(p, d);
G.printGraph();
G.printMRFStructure();
X = G.sample(1000);
G.testSampleMembership(X),
G.testSampleMembership(randn(10, p)),

