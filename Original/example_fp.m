% This is an example of Multichannel Analysis of Surface Wave
% using f-p analysis

% Read seismic data
[Data, dt]=read_org('mc403_sx.org');

% Peform f-p Analysis
df=1.0;fmin=5; fmax=80; % Frequency Range
pmin=1/500; pmax=1/100; dp=(pmax-pmin)/500; % Wavenumber Range
dx=1; % Geophone spacing
%dx=4; % Geophone spacing
%Data=[Data(:,1) Data(:,5) Data(:,9) Data(:,13) Data(:,17) Data(:,21)];

[N, Ch]=size(Data);x1=10; x=[x1:dx:x1+(Ch-1)*dx]';
[f,p,A]=fp(Data,dt,x,fmin,fmax,df,pmin,pmax,dp);

% Find the dispersion curve from the f-p image
fe=f;
Vphe=zeros(length(fe),1);
pe=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A(i,:));
   pe(i)=p(Ai);
   Vphe(i)=1/pe(i);
end

% Plotting the f-p image and dispersion curve
figure;
imagesc(p,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=1;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('slowness, p (s/m)');
ylabel('Frequency, f (Hz)');
hold on;
plot(pe,fe,'r.');hold on;


% Plotting the f-v image and dispersion curve
[pp,ff]=meshgrid(p,f);
vv=1./pp;

figure
surfc(vv,ff,A);
view([0,-90]);
shading flat;

colormap(meshgrid(1:-1/63:0,1:3)');
%c1=0; c2=0; c3=0; cs=3;
%cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
%cmap=cmap.^cs;
%colormap(cmap);
colorbar;
%set(gca,'xticklabel',{v(20),v(40),v(60),v(80),v(100),v(120)});
%set(gca,'yticklabel',{f(10),f(20),f(30),f(40),f(50),f(60),f(70),f(80),f(90)});
title('Overtone Diagram');
xlabel('Phase Velocity, v_p_h (m/s)');
ylabel('Frequency, f (Hz)');
hold on;
plot(Vphe,fe,'r.');hold on;
axis([0 1000 0 80]);


% Perform f-v Analysis
vmin=150; vmax=500; dv=1; % Velocity Range
fmin=5; fmax=80; df=1; % Frequency Range
dx=1; % Geophone spacing
[N, Ch]=size(Data);x1=28; x=[x1:dx:x1+(Ch-1)*dx]';
[f,v,A1]=fv(Data,dt,x,fmin,fmax,df,vmin,vmax,dv);

% Find the dispersion curve from the f-v image
fe1=f;
Vphe1=zeros(length(fe1),1);
for i=1:length(fe1);
   [Amax Ai]=max(A1(i,:));
   Vphe1(i)=v(Ai);
end

% Perform mutichanel SASW analysis
[fe2, Vphe2]=MSASW(Data,dt,x,fmin,fmax,df);

% Compare and plot f-p, f-v, and MSASW
figure
plot(fe, Vphe,'r*');hold on;
plot(fe1, Vphe1,'bo');

plot(fe2, Vphe2,'b.');
%axis([0 80 0 500]);
ylabel('Phase Velocity, v_p_h (m/s)');
xlabel('Frequency, f (Hz)');
%legend('f-p','f-v','MSASW');

