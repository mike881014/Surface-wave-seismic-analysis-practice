clear all;clc%close all;clc
cd D:\study\Hank\onelayercavitysize5mxyz20022525ht2.5mh2mmeshsize0.7m0.8sec0point5\
first_fnu = 5;
last_fnu =  46;
NF=last_fnu-first_fnu+1; % Number of files
fn1='DB.X';Direction1='Z';Offset = 50;
fn3_Z=['.FX',Direction1,'.semd'];

Direction2='X';
fn3_X=['.FX',Direction2,'.semd'];

Direction3='Y';
fn3_Y=['.FX',Direction3,'.semd'];

for i=1:1:NF
    fn2=num2str(first_fnu+i-1); % The following lines putting the pieces of file names together
    FileName1=[fn1 fn2 fn3_Z];
    RawData=load(FileName1);
    yZ(:,i)=RawData(:,2);
    
    FileName2=[fn1 fn2 fn3_X];
    RawData=load(FileName2);
    yX(:,i)=-RawData(:,2);
    
    FileName3=[fn1 fn2 fn3_Y];
    RawData=load(FileName3);
    yY(:,i)=-RawData(:,2);

    %y(:,i) = smooth(y(:,i),10);
    yZ_norm(:,i) = yZ(:,i)/max(abs(yZ(:,i)));
    yX_norm(:,i) = yX(:,i)/max(abs(yX(:,i)));
    yY_norm(:,i) = yY(:,i)/max(abs(yY(:,i)));
end



[N, Ch]=size(yZ);
dx=1; % Geophone spacing
x1=10; % distance between source & 1st receiver
y1=5;
x=[x1:dx:x1+(Ch-1)*dx]';
y=[y1:dx:y1+(Ch-1)*dx]';
% t=0:dt:(N-1)*dt;
t = RawData(:,1); dt=t(2)-t(1); % time step
% Dt=50;
[m,t0_ind] = min(abs(t)); % Time correction, the simulation time start from negitive time 
t = t(t0_ind:end);
yZ = yZ(t0_ind:end,:);
yX = yX(t0_ind:end,:);
yY = yY(t0_ind:end,:);
N = length(t);
N=round(1/dt);t = [0:dt:(N-1)*dt];
data_length = length(yX);
yX=[yX;zeros(N-data_length,length(x))];
yY=[yY;zeros(N-data_length,length(x))];
yZ=[yZ;zeros(N-data_length,length(x))];
%% backscatter
% -------------------------
% matrix = [1 2 3; 0 0 0; 4 5 6; 0 0 0];
% zero_rows = all(matrix == 0, 2);
% zero_index = find(zero_rows, 1);
% matrix(zero_index, :) = [];
% zero_rows = all(yZ == 0, 2);
% zero_index = find(zero_rows, 1);
% yZ(zero_index,:)=[];
ybackscatterZ=extract_backscatter(yZ, dt); 
figure; % Examine the original waveforms
subplot(211); 
wigb(yZ,1.,x,t);
set(gca,'FontSize',18);
set(gcf,'color','w');
ylabel('Time (s)'); xlabel('Offset (m)');
title('Vertical Original Signal');
subplot(212); 
wigb(ybackscatterZ,1.,x,t);
set(gca,'FontSize',18);
set(gcf,'color','w');
ylabel('Time (s)'); xlabel('Offset (m)');
title('Vertical extracted backscatter'); 

% zero_rows = all(yX == 0, 2);
% zero_index = find(zero_rows, 1);
% yX(zero_index,:)=[];
% ybackscatterX=extract_backscatter(yX, dt); 
% figure; % Examine the original waveforms
% subplot(211); 
% wigb(yX,1.,x,t);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% ylabel('Time (s)'); xlabel('Offset (m)');
% title('Radial Original Signal');
% subplot(212); 
% wigb(ybackscatterX,1.,x,t);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% ylabel('Time (s)'); xlabel('Offset (m)');
% title('Radial extracted backscatter'); 
