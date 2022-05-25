clear all
clc
close all
format short
load mycolormap;
tic;

G1=load('data/G_GPS.txt','-ascii');
[D_GPS,P_GPS]= textread('data/GPS_D_P.txt','%f %f');
def1=D_GPS;
P1=diag(P_GPS);
[n1,n]=size(def1);

[X_InSAR,Y_InSAR,D_InSAR,P_InSAR,looke,lookn,looku]= textread('data/InSAR.txt','%f %f %f %f %f %f %f');
G2=load('data/G_InSAR.txt','-ascii');
def2=D_InSAR;
[n2,n]=size(def2);
P2=diag(P_InSAR);

T=load('data/T.txt','-ascii');
[h1,h2]=size(T)
Tzeros=zeros(h1,1);
p1=1;
p2=1;
p3=1;
m=1;
while (1)  
slip=inv(G1'*P1*p1*G1+G2'*P2*p2*G2+T'*p3*T)*(G1'*p1*P1*def1+G2'*p2*P2*def2);
v1=G1*slip-def1;
v2=G2*slip-def2;
v3=T*slip-Tzeros ;
w=[v1'*p1*P1*v1;v2'*p2*P2*v2;v3'*p3*v3];
N1=G1'*p1*P1*G1;
N2=G2'*p2*P2*G2;
N3=T'*p3*T;
N=N1+N2+N3;
s1=2*trace(inv(N)*N1)-trace(inv(N)*N1*inv(N)*N1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
s2=trace(inv(N)*N1*inv(N)*N2);
s3=trace(inv(N)*N1*inv(N)*N3);
s4=trace(inv(N)*N2*inv(N)*N1);
s5=2*trace(inv(N)*N2)-trace(inv(N)*N2*inv(N)*N2);
s6=trace(inv(N)*N2*inv(N)*N3);
s7=trace(inv(N)*N3*inv(N)*N1);
s8=trace(inv(N)*N3*inv(N)*N2);
s9=2*trace(inv(N)*N3)-trace(inv(N)*N3*inv(N)*N3);
s=[n1-s1,s2,s3;s4,n2-s5,s6;s7,s8,h1-s9];
% c=inv(s'*s)*s'*w;%Variance components
c=inv(s)*w;
if (abs(c(1)-c(2))<10^(-6)) && (abs(c(2)-c(3))<10^(-6))&& (abs(c(1)-c(3))<10^(-6))
    break
else
    c1=c(1);
p1=c1/c(1)*p1
p2=c1/c(2)*p2
p3=c1/c(3)*p3
m=m+1; 
end
end
 save data/p2.txt -ascii p2
 save data/p3.txt -ascii p3
 
 %% inversion slip
 [xc,yc,zc,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=...
textread('data/fault.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
slipt= load('data/slip_sum.txt','-ascii');
S1t= load('data/slip_str.txt','-ascii');
S2t= load('data/slip_dip.txt','-ascii');
Slip=inv(G1'*P1*p1*G1+G2'*P2*p2*G2+T'*p3*T)*(G1'*p1*P1*def1+G2'*p2*P2*def2);
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
RMSE=sqrt(((def1-G1*Slip)'*p1*P1*(def1-G1*Slip)+(def2-G2*Slip)'*p2*P2*(def2-G2*Slip))/(n1+n2))
Mo = sum(30e9.*slip.*area.*1000.*1000)/1e-7;
Mw=(2/3)*log10(Mo)-10.7
Mo=Mo/1e7

%% figure
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
    patch(xfault(i,:),yfault(i,:),zfault(i,:),slipe(i),'EdgeColor','k','LineWidth',0.001);
end
hold on
view(3), grid on, axis equal,rotate3d;colormap(flipud(ans));colorbar;caxis([min(slipt) max(slipt)])


hold on 
subplot (2,3,5)
S1E=abs(S1t-S1);
S2E=abs(S2t-S2);
 for i=1:length(xc)
    patch(xfault(i,:),yfault(i,:),zfault(i,:),S1E(i),'EdgeColor','k','LineWidth',0.001);
end
hold on
view(3), grid on, axis equal,rotate3d;colormap(flipud(ans));colorbar;caxis([min(S1) max(S1)])




hold on 
subplot (2,3,6)
for i=1:length(xc)
    patch(xfault(i,:),yfault(i,:),zfault(i,:),S2E(i),'EdgeColor','k','LineWidth',0.001);
end
hold on
view(3), grid on, axis equal,rotate3d;colormap(flipud(ans));colorbar;caxis([min(S2) max(S2)])



