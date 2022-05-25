%Simulate slip, observations, Green function,
%% 1 simulate TRUE fault slip
clc;
clear all;
close all;
[xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=textread('data\fault.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
y=-14:2:14;
x=-14:2:14;
[X,Y]=meshgrid(x,y);
[m,n]=size(X);
a=1.5;%hanshu zuidazhi
D=[30 0;0 30];%gaosi  hanshu  fangchajuzhen
for i=1:m
    for j=1:n
        slip_sim(i,j)=a.*exp(-0.5*[X(i,j) Y(i,j)]*inv(D)*[X(i,j) Y(i,j)]');
    end
end
[mm,nn]=size(slip_sim)
Slip=reshape(slip_sim,mm*nn,1);
slip1=Slip.*cos(30.*pi/180);
slip2=Slip.*sin(30.*pi/180);
save data/slip_str.txt -ascii slip1
save data/slip_dip.txt -ascii slip2
a=[slip1,slip2];
b=a';
[m1,n1]=size(a);
slip=reshape(b,m1*n1,1);
slip_sum = (slip1.^2+slip2.^2).^0.5; 
save data/slip_sum.txt -ascii slip_sum
dd=max(slip_sum)
ee=mean(slip_sum)
%% 2 calulate observations 
% 2.1 GPS observations points
x_min=-60;y_min=-60;x_max=60;y_max=60;dx=7;dy=7;
[Xe,Yn]=meshgrid(linspace(x_min,x_max,dx),linspace(y_min,y_max,dy));
[m,n]=size(Xe);
Xe=reshape(Xe,m*n,1);
Yn=reshape(Yn,m*n,1);
XY_GPS=[Xe,Yn];
for j = 1:size(XY_GPS(:,1))
       Xd(j,:) = XY_GPS(j,1);%x distance, in km
       Yd(j,:) = XY_GPS(j,2);%y distance, in km
end
for i = 1:length(xs)
    for j = 1:size(Xe)
        xrs(i,j) =Xd(j)-xs(i);%x distance, in km
        yrs(i,j) =Yd(j)-ys(i);%y distance, in km
        zrs(i,j) = -zs(i);%z distance, in km
    end
end
G_GPS= okada_green_GPS(xrs,yrs,zrs,strike,dip,width,len);
save('data/G_GPS.txt','G_GPS','-ascii');
d_GPS=G_GPS*slip;
num_gps=length(d_GPS);
media=0;%qiwang
deviazione1=0.*10^(-3);%zhongwucha
r_GPS1=normrnd(media,deviazione1,num_gps/3);
deviazione2=3*10^(-3);%zhongwucha
r_GPS2=normrnd(media,deviazione2,num_gps/3);
deviazione3=6*10^(-3);%zhongwucha
r_GPS3=normrnd(media,deviazione3,num_gps/3);
D_GPS=[];r_GPS=[];
for i=1:3:(length(d_GPS)-2)
    d1=d_GPS(i)+r_GPS1(i);
    d2=d_GPS(i+1)+r_GPS2(i);
    d3=d_GPS(i+2)+r_GPS3(i);
    D_GPS=[D_GPS;d1;d2;d3];
    r_GPS=[r_GPS;r_GPS1(i);r_GPS2(i);r_GPS3(i)]
end
save data/r_GPS.txt -ascii r_GPS
P_GPS=(r_GPS(1).^2)./(r_GPS.^2);
P1=diag(P_GPS);
GPS_D_P=[D_GPS,P_GPS];
GPS_XY=[XY_GPS(:,1),XY_GPS(:,2)];
save data/GPS_D_P.txt -ascii GPS_D_P
save data/GPS_XY.txt -ascii GPS_XY

% 2.2 InSAR observations points
x_min=-60;y_min=-60;x_max=60;y_max=60;dx=30;dy=30;
[Xe,Yn]=meshgrid(linspace(x_min,x_max,dx),linspace(y_min,y_max,dy));
[m,n]=size(Xe);
Xe=reshape(Xe,m*n,1);
Yn=reshape(Yn,m*n,1);
XY_InSAR=[Xe,Yn];
looke_n_u_asc=[-0.62,-0.11,0.72];
% looke_n_u_des=[0.62,-0.11,0.77];
looke=looke_n_u_asc(1).*ones(length(XY_InSAR(:,1)),1);
lookn=looke_n_u_asc(2).*ones(length(XY_InSAR(:,1)),1);
looku=looke_n_u_asc(3).*ones(length(XY_InSAR(:,1)),1);
for j = 1:length(XY_InSAR(:,1))
       Xd(j,:) = XY_InSAR(j,1);%x distance, in km
       Yd(j,:) = XY_InSAR(j,2);%y distance, in km
end

for i = 1:length(xs)
    for j = 1:length(XY_InSAR(:,1))
        xrs(i,j) =Xd(j)-xs(i);%x distance, in km
        yrs(i,j) =Yd(j)-ys(i);%y distance, in km
        zrs(i,j) = -zs(i);%z distance, in km
    end
end
G_InSAR= okada_green_InSAR(xrs,yrs,zrs,strike,dip,width,len,looke,lookn,looku);
save('data\G_InSAR.txt','G_InSAR','-ascii');
d_InSAR=G_InSAR*slip;
media=0;%qiwang
deviazione4=6*10^(-3);%zhongwucha
r_InSAR=normrnd(media,deviazione4,size(d_InSAR));
D_InSAR=d_InSAR+r_InSAR;
P_InSAR=(r_InSAR(1).^2)./(r_InSAR.^2);
P2=diag(P_InSAR);
InSAR=[XY_InSAR(:,1),XY_InSAR(:,2),D_InSAR,P_InSAR,looke,lookn,looku];
save data/InSAR.txt -ascii InSAR
save data/r_InSAR.txt -ascii r_InSAR
%% 3 Laplace constraint
k=1;
T=[];
num_fault_L=16;%the points along the strike direction
num_fault_W=16;%the points along the dip direction
for j = 1:num_fault_L-1
    for i = 1:num_fault_W -1
        for m = 1:2
            index1 = (j-1)*(num_fault_L-1)+i;
            index2 = (j-1)*(num_fault_L-1)+i-1;
            index3 = (j-1)*(num_fault_L-1)+i+1;
            index4 = (j-2)*(num_fault_L-1)+i;
            index5 = (j)*(num_fault_L-1)+i;
            dx = len(1);
            dy = width(1);
            if (index1 >=1&& index1 <= length(xs))
                   T(k,2*(index1-1)+m) = -2*(dx.^-2+dy.^-2);
            end
            if (index2 >= 1 && index2 <= length(xs))
               
                   T(k,2*(index2-1)+m) = dx.^-2;
                
            end
            if (index3 >= 1 && index3 <= length(xs))
                
                   T(k,2*(index3-1)+m) = dx.^-2;
            end
            if (index4 >= 1 && index4 <= length(xs) )
                 T(k,2*(index4-1)+m) = dy.^-2;
            end
            if (index5 >= 1 && index5 <= length(xs))
                  T(k,2*(index5-1)+m) = dy.^-2;
            end
            k=k+1;
        end
    end
end
[h1,h2]=size(T)
save data/T.txt -ascii T
Tzeros=zeros(h1,1);
%% 4 figure
% gps points
[m1,n1]=size(D_GPS);
D_GPS1=reshape(D_GPS,3,m1/3)';P_GPS=reshape(P_GPS,3,m1/3)';r_GPS=reshape(r_GPS,3,m1/3)';%e,n,up
D_GPSE=D_GPS1(:,1);r_GPSE=r_GPS(:,1);P_GPSE=P_GPS(:,1);
D_GPSN=D_GPS1(:,2);r_GPSN=r_GPS(:,2);P_GPSN=P_GPS(:,2);
D_GPSU=D_GPS1(:,3);r_GPSU=r_GPS(:,3);P_GPSU=r_GPS(:,3);
x=[XY_GPS(:,1);25];y=[XY_GPS(:,2);25];D_GPSE=[D_GPSE;0];D_GPSN=[D_GPSN;0.05];D_GPSU=[D_GPSU;0.05];
figure;
subplot(1,2,1);
quiver(x,y,D_GPSE,D_GPSN,2,'b');

axis([-60,60,-60,60]);
hold on
subplot(1,2,2);
num_GPS=length(XY_GPS(:,1));
T_vertical=zeros(num_GPS+1,1);

quiver(x,y,T_vertical,D_GPSU,2,'r');
axis([-60,60,-60,60]);
% figure;
% subplot(1,3,1);
% scatter(XY_GPS(:,1),XY_GPS(:,2),15,D_GPSE);
% xlabel('E-W/km'), ylabel('N-S/km');
% hold on
% subplot(1,3,2);
% scatter(XY_GPS(:,1),XY_GPS(:,2),15,D_GPSN);
% xlabel('E-W/km'), ylabel('N-S/km');
% hold on
% subplot(1,3,3);
% scatter(XY_GPS(:,1),XY_GPS(:,2),15,D_GPSU);
% xlabel('E-W/km'), ylabel('N-S/km');
% h=colorbar('eastoutside','ytick',-0.3:0.05:0.1);
% set(get(h,'Title'),'string','m');
% colormap(jet);
% InSAR points
figure;
scatter(XY_InSAR(:,1),XY_InSAR(:,2),40,D_InSAR);
colormap(jet)
xlabel('E-W/km'), ylabel('N-S/km');
h=colorbar('eastoutside','ytick',-0.3:0.05:0.1);
set(get(h,'Title'),'string','m');
colormap(jet);
save data2/data2.mat D_GPS D_InSAR slip_sum G_GPS G_InSAR Tzeros T P1 P2



