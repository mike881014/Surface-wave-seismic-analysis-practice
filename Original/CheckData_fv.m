% This is a program for on-site quality control (Checking data)
clear; close all; clc;
cd('/Users/gaoqianhui/Desktop/模擬/上次');
first_fnu=1; % first file number
last_fnu=12;  % Last file number

% =========================Read the seismic data=========================
NF=last_fnu-first_fnu+1; % Number of files
   fn1='X';
   fn3='.DB.FXZ.semd';
    for i=1:1:NF
        fn2=num2str(first_fnu+i-1); % The following lines putting the pieces of file names together
        FileName=[fn1 fn2 fn3];
        RawData=load(FileName);
        y(:,i)=RawData(:,2);
    end

dx=0.1; % Geophone spacing
x1=0.5; % distance between source & 1st receiver
dt=RawData(2,1)-RawData(1,1); % time step

[N, Ch]=size(y);
x=[x1:dx:x1+(Ch-1)*dx]';
t=0:dt:(N-1)*dt;

%---------------------------------------------------------------------------------------
% Peform f-v transform
fmin=0; fmax=14000;df=2000; % Frequency Range
vmin=0; vmax=2000; dv=500; % Velocity Range
[f,v,A,U]=fv(y,dt,x,fmin,fmax,df,vmin,vmax,dv);

% Find the dispersion curve from the overone image
Vphe=zeros(length(f),1);

for i=1:length(f);
   [Amax Ai]=max(A(i,:));
   Vphe(i)=v(Ai);
end

%---------------------------------------------------------------------------------------
% Perform MSASW regression of phi-x
p_unwrap = 1.0;
UangleU=unwrap(angle(U),[p_unwrap*pi],2); % Unwrap phase angle

Vphe2=zeros(length(f),1);
R2=zeros(length(f),1);
for i=1:length(f)
   phi=UangleU(i,:);
   a=1; b=Ch; % Effective channel from a to b
   %b=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis
   [b bint r rint stats]=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis with R2 statistic
   Vphe2(i)=abs(b(1))*2*pi*f(i);
   R2(i)=stats(1); % R2 statistic
end

%---------------------------------------------------------------------------------------
% Plot seismic data in t-x domain
figure
% subplot(141)
wigb(y,1.,x,t);
set(gca,'FontSize',16);
title('Seismic data')
ylabel('Time (s)')
xlabel('Offset (m)')
set(gca,'Ydir','reverse');

% %plot in f-x domain amaplitude
% figure
% subplot(142);
% wigb(abs(U),1.,x,f);
% set(gca,'FontSize',14);
% ylabel('frequency (Hz)')
% xlabel('Offset (m)')
% set(gca,'Ydir','reverse');
% 
% % plot in f-x domain R2 of phi-x regression
% subplot(143);
% plot(R2,f,'b'); hold on;
% axis([0 1.0 0 100]);
% set(gca,'Ydir','reverse','FontSize',14);
% ylabel('Frequency (Hz)','FontSize',14);
% xlabel('R^2','FontSize',14);
% 
% % Plotting the f-v spectrum and dispersion curve
% subplot(144);
% imagesc(v,f,A);
% %colormap(meshgrid(1:-1/63:0,1:3)');
% c1=0; c2=0; c3=0; cs=3;
% cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
% cmap=cmap.^cs;
% colormap(cmap);
% colorbar;
% %set(gca,'xticklabel',{v(20),v(40),v(60),v(80),v(100),v(120)});
% %set(gca,'yticklabel',{f(10),f(20),f(30),f(40),f(50),f(60),f(70),f(80),f(90)});
% %title('Overtone Diagram');
% set(gca,'FontSize',14);
% xlabel('Phase Velocity, v_p_h (m/s)');
% ylabel('Frequency, f (Hz)');
% hold on;
% plot(Vphe,f,'r.');hold on;
% plot(Vphe2,f,'bo');
% legend('MWTSW','MSASW');

%th=load('2lay.txt');

figure;
imagesc(v,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=3;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
hold on;
plot(Vphe,f,'bo');hold on;
set(gca,'FontSize',16);
title('Dispersion curve')
legend('1st mode','2nd mode','3rd mode','MWTSW')
xlabel('Phase Velocity, v_p_h (m/s)');
ylabel('Frequency, f (Hz)');
