clear all
close all

s = tf('s');
G = 1/(s^2+s);
bode(G)

