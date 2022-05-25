clc
clear all
close all
[D_GPS,P_GPS]= textread('data/GPS_D_P.txt','%f %f');
[X_GPS,Y_GPS]= textread('data/GPS_XY.txt','%f %f');
[X_InSAR,Y_InSAR,D_InSAR,P_InSAR,looke,lookn,looku]= textread('data/InSAR.txt','%f %f %f %f %f %f %f');
G_InSAR= load('data/G_InSAR.txt','-ascii');
r_InSAR= load('data/r_InSAR.txt','-ascii');
T=load('data/T.txt','-ascii');
G_GPS= load('data/G_GPS.txt','-ascii');
r_GPS= load('data/r_GPS.txt','-ascii');
P1=diag(P_GPS);P2=diag(P_InSAR);
D1=r_GPS.^2; D2=r_InSAR.^2;slip_sum= load('data/slip_sum.txt','-ascii');
[h1,h2]=size(T)
Tzeros=zeros(h1,1);
save data/data.mat D_GPS D_InSAR slip_sum G_GPS G_InSAR Tzeros T P_GPS P_InSAR D1 D2