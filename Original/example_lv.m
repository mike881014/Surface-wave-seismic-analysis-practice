% This is an example of Multichannel Analysis of Surface Wave
% using l-f analysis

% Read seismic data
[Data, dt]=read_org('mc403_sx.org');

% Peform l-f Analysis
vmin=100; vmax=400; dv=2; % Velocity Range
lmin=1; lmax=50; dl=1; % Wavelength Range
dx=1; % Geophone spacing

[N, Ch]=size(Data);x1=10; x=[x1:dx:x1+(Ch-1)*dx]';
[v,l,A,U]=lv(Data,dt,x,vmin,vmax,dv,lmin,lmax,dl);

% Find the dispersion curve from the f-l image
le=l; 
ve=zeros(length(le),1);
for i=1:length(le);
   [Amax Ai]=max(A(:,i));
   ve(i)=v(Ai);
end
Vphe=ve;

% Plotting the v-l image and dispersion curve
figure;
imagesc(l,v,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=1;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('Wavelength (m)');
ylabel('Velocity, v (m/s)');
hold on;
plot(le,ve,'r.');hold on;


% Plotting the f-v image and dispersion curve
[ll,vv]=meshgrid(l,v);
ff=vv./ll;
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
fe=ve./le;
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
%[fe2, Vphe2]=MSASW(Data,dt,x,fmin,fmax,df);

% Compare and plot f-l, f-v, and MSASW
figure
plot(fe, Vphe,'r*');hold on;
plot(fe1, Vphe1,'bo');

%plot(fe2, Vphe2,'b.');
%axis([0 80 0 500]);
ylabel('Phase Velocity, v_p_h (m/s)');
xlabel('Frequency, f (Hz)');
%legend('f-k','f-v','MSASW');

