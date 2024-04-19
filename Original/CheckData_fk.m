% This is a program for on-site quality control (Checking data)

% Read seismic data
[y, dt]=read_org('Mc202_SX.org');
dx=1; % Geophone spacing
x1=15; 

%y=fliplr(y);
%y=y(:,4:24);

[N, Ch]=size(y);
x=[x1:dx:x1+(Ch-1)*dx]';
t=0:dt:(N-1)*dt;

%---------------------------------------------------------------------------------------
% Peform f-k transform
fmin=1; fmax=100; df=1; % Frequency Range
kmin=2*pi*5/500; kmax=2*pi*100/50; dk=(kmax-kmin)/200; % Wavenumber Range
[f,k,A,U]=fk(y,dt,x,fmin,fmax,df,kmin,kmax,dk);

% Find the dispersion curve from the overone image
ke=zeros(length(f),1);

for i=1:length(f);
   [Amax Ai]=max(A(i,:));
   ke(i)=k(Ai);
end

% Find the dispersion curve from the f-v image
fe=f;
Vcphe=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A(i,:));
   Vcphe(i)=v(Ai);
end

% Plotting the f-v image and dispersion curve
figure;
imagesc(v,f,Ac);
cmap = load('colormap08.txt');
colormap(cmap);
colorbar;
%---------------------------------------------------------------------------------------
% Perform MSASW regression of phi-x
p_unwrap = 1.0;
UangleU=unwrap(angle(U),[p_unwrap*pi],2); % Unwrap phase angle

ke2=zeros(length(f),1);
R2=zeros(length(f),1);
for i=1:length(f)
   phi=UangleU(i,:);
   a=1; b=Ch; % Effective channel from a to b
   %b=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis
   [b bint r rint stats]=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis with R2 statistic
   %Vphe2(i)=abs(b(1))*2*pi*f(i);
   ke2(i)=1./abs(b(1));
   R2(i)=stats(1); % R2 statistic
end

%---------------------------------------------------------------------------------------
% Plot seismic data in t-x domain
figure
subplot(141)
wigb(y,1.,x,t);
set(gca,'FontSize',10);
ylabel('Time (s)')
xlabel('Offset (m)')
set(gca,'Ydir','reverse');

%plot in f-x domain amaplitude
subplot(142);
wigb(abs(U),1.,x,f);
set(gca,'FontSize',10);
ylabel('frequency (Hz)')
xlabel('Offset (m)')
set(gca,'Ydir','reverse');

% plot in f-x domain R2 of phi-x regression
subplot(143);
plot(R2,f,'b'); hold on;
axis([0 1.0 0 100]);
set(gca,'Ydir','reverse','FontSize',10);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('R^2','FontSize',10);

% Plotting the f-k image and dispersion curve
subplot(144)
imagesc(k,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=1;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('Wavenumber');
ylabel('Frequency, f (Hz)');
hold on;
plot(ke,f,'r.');hold on;
plot(ke2,f,'bo');
legend('MWTSW','MSASW');
