% ABIC for slip distribution inversion
% Edited DEC 11, 2021, Wuhan University, by Xiong Zhao
%code from funning GJI && Yabuki 1992
clc ;
clear all;
close all;
tic
format long
load mycolormap;
%% 1 read data
num_fault_L=16;%the points along the strike direction
num_fault_W=16;%the points along the dip direction
%fault parameters
[xc,yc,zc,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=...
textread('data/fault.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
[D_GPS,P_GPS]= textread('data/GPS_D_P.txt','%f %f');
[X_GPS,Y_GPS]= textread('data/GPS_XY.txt','%f %f');
G_GPS= load('data/G_GPS.txt','-ascii');
r_GPS= load('data/r_GPS.txt','-ascii');
num_gps=length(D_GPS)
[X_InSAR,Y_InSAR,D_InSAR,P_InSAR,looke,lookn,looku]= textread('data/InSAR.txt','%f %f %f %f %f %f %f');
G_InSAR= load('data/G_InSAR.txt','-ascii');
r_InSAR= load('data/r_InSAR.txt','-ascii');
num_InSAR=length(X_InSAR)
T=load('data/T.txt','-ascii');
num_T=(num_fault_L-1)*(num_fault_W -1)
%% 2 Define the range of weight values(w1=1,w2,w3……） and regularization parameters（alpha)
  % notes: "[LL,UU]=ilu(**) & logdet=real()" are used to cal the log
  % value of determinant. For different matlab environments formats 
  % for the use of "ilu" are different. These two columns of codes
  % will generate same output as "[LL,UU]=lu(**) & logdet=real()".
w1=1;w2=linspace(0.05,0.5,10);alpha=logspace(log10(0.000001),log10(0.5),500);
Para_loop = [];
for i_w2 = 1:length(w2)
   for i_alpha = 1:length(alpha)
   D=diag([r_GPS.^2;1/w2(i_w2).*r_InSAR.^2]); P=diag([P_GPS;w2(i_w2).*P_InSAR]);d_obs=[D_GPS;D_InSAR];
   G=[G_GPS;G_InSAR]; 
   m=inv(G'*P*G+alpha(i_alpha)*T'*T)*(G'*P*d_obs);
   mod1=G_GPS*m;mod2=G_InSAR*m;d_mod=[mod1;mod2];
   s=(d_obs-d_mod)'*P*(d_obs-d_mod)+alpha(i_alpha)*((T*m)'*(T*m)); 
   LL1 =G'*P*G+alpha(i_alpha).*(T'*T);
   LL2 =alpha(i_alpha).*(T'*T);
   ABIC= (num_gps+num_InSAR)*log (s)+logdet(LL1)-logdet(LL2)+logdet(D);
   RMSE_GPS = norm((D_GPS-mod1).*diag(sqrt(P_GPS))) / sqrt(num_gps);
   RMSE_InSAR = norm((D_InSAR-mod2).*diag(sqrt(P_InSAR))) / sqrt(num_InSAR);
   Roughness=sum(abs(T*m)) /num_T;
   Para_loop=[Para_loop;w2(i_w2),alpha(i_alpha),ABIC,RMSE_GPS,RMSE_InSAR,Roughness];
   end
 end
OUTFILE =['All_para_for_ABIC_loops.out'];
header = 'w2,alpha,ABIC,RMSE_GPS,RMSE_InSAR,Roughness';
dlmwrite(OUTFILE,header,'delimiter',''); 
dlmwrite(OUTFILE,Para_loop,'-append','delimiter','\t','precision','%.6f');
%% 3. find and save solutions with ABIC minimum
a=reshape(Para_loop(:,3),length(alpha),length(w2));
figure;
plot(alpha,a(:,1),'c','LineWidth',2);
hold on
plot(alpha,a(:,2),'m','LineWidth',2);
hold on
plot(alpha,a(:,3),'y','LineWidth',2)
hold on
plot(alpha,a(:,4),'-.','LineWidth',2);
hold on
plot(alpha,a(:,5),'k','LineWidth',2);
hold on
plot(alpha,a(:,6),':','LineWidth',2);
hold on
plot(alpha,a(:,7),'--','LineWidth',2);
hold on
plot(alpha,a(:,8),'b','LineWidth',2);
hold on
plot(alpha,a(:,9),'g','LineWidth',2);
hold on
plot(alpha,a(:,10),'k','LineWidth',2);
legend('w2=0.05','w2=0.1','w2=0.15','w2=0.2','w2=0.25','w2=0.3','w2=0.35','w2=0.4','w2=0.45','w2=0.5');

% plot(Para_loop(:,2),Para_loop(:,3),'LineWidth',2);
[optimal_kappa_column, optimal_kappa_row] = find( Para_loop(:,3)==min(Para_loop(:,3)) );
w2_o = Para_loop(optimal_kappa_column,1);
alpha_o = Para_loop(optimal_kappa_column,2);
RMSE_GPS_o = Para_loop(optimal_kappa_column,4);
RMSE_InSAR_o = Para_loop(optimal_kappa_column,5);
Roughness_o = Para_loop(optimal_kappa_column,6);
% w2_o =0.1;alpha_o = 1.66e-06;
Slip=inv(w1*G_GPS'*diag(P_GPS)*G_GPS+w2_o*G_InSAR'*diag(P_InSAR)*G_InSAR+alpha_o.*T'*T)*(G_GPS'*diag(P_GPS)*D_GPS+w2_o*G_InSAR'*diag(P_InSAR)*D_InSAR);
 for i = 1:length(xc)
    S1(i,1)=Slip((i-1)*2+1);
    S2(i,1)=Slip((i-1)*2+2);
    end
slip = (S1.^2+S2.^2).^0.5;
theta=rad2deg(atan(S2./S1));
for i=1:length(theta)
    if (S1(i)<0&&S2(i)>0)
        theta(i)=theta(i)+180;
    end
end
maxslip=max(slip)
meanslip=mean(slip)
  for j = 1:length(slip)
    if slip(j)==maxslip
        theta0=rad2deg(atan(S2(j,1)./S1(j,1)));
    end
  end
theta0
theta1=mean( theta)
Mo = sum(30e9.*slip.*area.*1000.*1000)/1e-7;
Mw=(2/3)*log10(Mo)-10.7
Mo=Mo/1e7
RMSE=sqrt(((D_GPS-G_GPS*Slip)'*diag(P_GPS)*(D_GPS-G_GPS*Slip)+(D_InSAR-G_InSAR*Slip)'*w2_o*diag(P_InSAR)*(D_InSAR-G_InSAR*Slip))/(num_gps+num_InSAR))
%% 4 figure
slipt= load('data/slip_sum.txt','-ascii');
S1t= load('data/slip_str.txt','-ascii');
S2t= load('data/slip_dip.txt','-ascii');
figure; 
 subplot(2,3,1)

 L=30;
 W=30;
num_fault_L=16;%the points along the strike direction
num_fault_W=16;%the points along the dip direction
    pat_l=linspace(0,L,num_fault_L);
    pat_w=linspace(-W,0,num_fault_W);
    F=[];
   for j=1:num_fault_W-1
       for i=1:num_fault_L-1
        F0=[pat_l(i) pat_w(j) pat_l(i+1) pat_w(j)  pat_l(i+1) pat_w(j+1) pat_l(i) pat_w(j+1)];
        F=[F;F0];
       end
   end
lowl=F(:,1:2);
lowr=F(:,3:4);
upr=F(:,5:6);
upl=F(:,7:8);
xc=(upr(:,1)+lowl(:,1))/2;
yc=(upr(:,2)+lowr(:,2))/2;
xfault=[lowl(:,1),lowr(:,1),upr(:,1),upl(:,1)];
yfault=[lowl(:,2),lowr(:,2),upr(:,2),upl(:,2)];
zfault=[zf1,zf2,zf3,zf4];
for i=1:length(xc)
    patch(xfault(i,:),zfault(i,:),slip(i),'EdgeColor','white','LineWidth',0.001);
end
hold on
axis([0,30,-21,0]);colormap(flipud(ans));colorbar;caxis([min(slipt) max(slipt)])


hold on
subplot (2,3,2)
 for i=1:length(xc)
    patch(xfault(i,:),zfault(i,:),S1(i),'EdgeColor','white','LineWidth',0.001);
end
hold on
axis([0,30,-21,0]);colormap(flipud(ans));colorbar;caxis([min(S1) max(S1)])

subplot (2,3,3)
 for i=1:length(xc)
    patch(xfault(i,:),zfault(i,:),S2(i),'EdgeColor','white','LineWidth',0.001);
end
hold on
axis([0,30,-21,0]);colormap(flipud(ans));colorbar;caxis([min(S2) max(S2)])

 

subplot (2,3,4)
slipe=abs(slipt-slip);
  for i=1:length(xc)
    patch(xfault(i,:),yfault(i,:),zfault(i,:),slipe(i),'EdgeColor','k','LineWidth',0.0001);
end
hold on
view(3), grid on, axis equal,rotate3d;colormap(flipud(ans));colorbar;caxis([min(slipt) max(slipt)])


hold on 
subplot (2,3,5)
S1E=abs(S1t-S1);
S2E=abs(S2t-S2);
 for i=1:length(xc)
    patch(xfault(i,:),yfault(i,:),zfault(i,:),S1E(i),'EdgeColor','k','LineWidth',0.0001);
end
hold on
view(3), grid on, axis equal,rotate3d;colormap(flipud(ans));colorbar;caxis([min(S1) max(S1)])




hold on 
subplot (2,3,6)
for i=1:length(xc)
    patch(xfault(i,:),yfault(i,:),zfault(i,:),S2E(i),'EdgeColor','k','LineWidth',0.0001);
end
hold on
view(3), grid on, axis equal,rotate3d;colormap(flipud(ans));colorbar;caxis([min(S2) max(S2)])

toc
