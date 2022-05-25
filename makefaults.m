% function makefaults
clear all
clc
close all
format short
%�ϲ㼸�β���
L=30;
W=30;
strike=60;
dip0=45;
dip=dip0*pi/180;
arfa=(strike)*pi/180;
depth=0;
%%%%%%%���ת��?���ֱ������µ�patchģ��
f2=[16,16];
ast=f2(1);
adi=f2(2);
w1=depth/sin(dip);
w2=W-w1;
zeta0=L/2;%%�����ϵת��Ϊ�����ϵ�������ϵԭ��λ�������ģ������ϵԭ��λ�����½ǣ�
eta0=w2;
pat_l=linspace(0,L,ast);
pat_w=linspace(0,W,adi);
F=[];
%￵�ÿһ��С�����½ǵ������?
for j=1:adi-1
    for i=1:ast-1
        F0=[pat_l(i),pat_w(j),pat_l(i+1),pat_w(j),pat_l(i+1),pat_w(j+1),pat_l(i),pat_w(j+1) strike dip0];
        F=[F;F0];
    end
end
%read corners
%�õ�һ��С���ĸ��ǵ����?
lowl=F(:,1:2);
lowr=F(:,3:4);
upr=F(:,5:6);
upl=F(:,7:8);
strike=F(:,9);
dip0=F(:,10);
%Get along strike and dip lengths
strike_length=L/(ast-1)*ones(size(strike));
dip_length=W/(adi-1)*ones(size(dip0));
area=(strike_length.*dip_length);

%ת�����󣬽��ϲ�������ת��Ϊ�ֲ��ѿ���ֱ�����ϵ�µ���ꣻ
R1=[cos(arfa) sin(arfa)*cos(dip);sin(arfa) -cos(arfa)*cos(dip);0 sin(dip)];
lowlXYZ=[R1*[lowl(:,1)'-zeta0;lowl(:,2)'-eta0]]';
lowrXYZ=[R1*[lowr(:,1)'-zeta0;lowr(:,2)'-eta0]]';
uprXYZ=[R1*[upr(:,1)'-zeta0;upr(:,2)'-eta0]]';
uplXYZ=[R1*[upl(:,1)'-zeta0;upl(:,2)'-eta0]]';
xc=(lowlXYZ(:,1)+uprXYZ(:,1))/2;
yc=(lowlXYZ(:,2)+uprXYZ(:,2))/2;
mzd=(lowlXYZ(:,3)+uprXYZ(:,3))/2-depth;
Fault=[xc yc mzd lowlXYZ(:,1),lowrXYZ(:,1),uprXYZ(:,1),uplXYZ(:,1),lowlXYZ(:,2),lowrXYZ(:,2),uprXYZ(:,2),uplXYZ(:,2),lowlXYZ(:,3)-depth,lowrXYZ(:,3)-depth, uprXYZ(:,3)-depth, uplXYZ(:,3)-depth, strike dip0 strike_length dip_length area];
%��Ϊ�����[ÿС������x y z ÿС������ ���� ���� ���� x y z(12��) ��λ�� ���?�ϲ��泤 �� ���]
xfault=[lowlXYZ(:,1),lowrXYZ(:,1),uprXYZ(:,1),uplXYZ(:,1)];
yfault=[lowlXYZ(:,2),lowrXYZ(:,2),uprXYZ(:,2),uplXYZ(:,2)];
zfault=[lowlXYZ(:,3)-depth,lowrXYZ(:,3)-depth, uprXYZ(:,3)-depth, uplXYZ(:,3)-depth];
figure(1)
for i=1:length(xc)
    patch(xfault(i,:),yfault(i,:),zfault(i,:),[1,0,1])
%     patch(xfault(i,:),yfault(i,:),zfault(i,:),slip(i)) %[1,0,1]Ϊ������ɫRGBֵ
end
hold on
view(3), grid on, axis equal,rotate3d
colormap(jet);
%xlabel('E to E/km');
%ylabel('N to N/km');
box on
save data/fault.txt Fault -ascii