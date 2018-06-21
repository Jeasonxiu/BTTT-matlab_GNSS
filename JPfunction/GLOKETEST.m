close all; clear all;

%% GLONASS Navigation to Keplian Elements
load('glomat.mat');
JD = EphGlo(1,18);
r = EphGlo(1,3:5);                                              % Position Vector
v = EphGlo(1,6:8);                                              % Velocity Vector

[KE] = rv2KE(JD,r,v);