% This is an example of Multichannel Analysis of Surface Wave
% using overtone analysis

% Read seismic data
[Data, dt]=read_org('mc403_sx.org');
%Data=fliplr(Data);
%Data=Data(:,2:24);

% Peform f-v Analysis
vmin=100; vmax=400; dv=2; % Velocity Range
fmin=1; fmax=200; df=1; % Frequency Range
dx=1; % Geophone spacing
[N, Ch]=size(Data);
x1=10; x=[x1:dx:x1+(Ch-1)*dx]';
t=[0:dt:(N-1)*dt]';

figure
wigb(Data,1.,x,t);
ylabel('Time in secs')
xlabel('Source-receiver distance in metres')



[f,v,A]=fv(Data,dt,x,fmin,fmax,df,vmin,vmax,dv);


% Find the dispersion curve from the f-v image
fe=f;
Vphe=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A(i,:));
   Vphe(i)=v(Ai);
end

% Plotting the f-v image and dispersion curve
figure;
imagesc(v,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=3;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('Phase Velocity, v_p_h (m/s)');
ylabel('Frequency, f (Hz)');
hold on;
plot(Vphe,fe,'r.');hold on;

%nf=50;
%[maxA, IA]=max(A(nf,:));
%plot(v,A(nf,:)/maxA*f(nf));
%plot(v(IA),f(nf),'b^');


% perform SASW analysis
%[f1, Vph1]=sasw(Data(:,3), Data(:,5), dt, 10,1/100,100);
%figure(1)
%plot (Vph1, f1, 'go');

% Perform mutichanel SASW analysis
%[fe2, Vphe2]=MSASW(Data(:,1:6),dt,x(1:6),fmin,fmax,df);
%figure(1)
%plot(Vphe2,fe2,'bo');

% Compare and plot f-v vs. MSASW
%figure
%plot(fe,Vphe,'r');hold on;
%plot(fe2,Vphe2,'b');
%xlabel('Phase Velocity, v_p_h (m/s)');
%ylabel('Frequency, f (Hz)');
%legend('Overtone','MSASW');

