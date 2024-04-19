Data=load('plot_source_time_function.txt');
dt=Data(2,1)-Data(1,1);
[N, Ch]=size(Data);
df=1/(N*dt);
f=[0:df:(N-1)*df]';
A=abs(fft(Data(:,2)));
figure
set(gca,'FontSize',16);
plot(f,A)
ylabel('Amplitude','FontSize',16);
xlabel('Frequency (Hz)','FontSize',16);


% figure
% set(gca,'FontSize',16);
% plot(Data(:,1),Data(:,2))
% xlabel('Time (s)','FontSize',16);
% title('Source time function')