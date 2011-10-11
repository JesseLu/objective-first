% Script used to test the SETUP() function.

% Simple test, make sure we get plane waves.
eps = [ones(15,1); 12.25*ones(10,1); ones(15,1)];
spec = setup(0.15, eps, [], 2*ones(1,40), 1);

