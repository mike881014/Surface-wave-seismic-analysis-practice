function y = seismic_data_muting(x,y,dt,N)
    beep; fprintf('Adjust image window, then Enter. \n'); pause; 
    beep; fprintf('Select upper and lower bound...\n');
    [xu tu] = ginput(2); plot(xu,tu);
    [xl tl] = ginput(2); plot(xl,tl);
    % muting with Tukey windows
    Tu = interp1(xu,tu,x);
    Tl = interp1(xl,tl,x);
    TuI = floor(Tu/dt);
    TlI = ceil(Tl/dt);
    for i=1:length(x)
        Nwin = TlI(i)-TuI(i);
        w = tukeywin(Nwin,0.25);
        temp1 = zeros(TuI(i)-1,1);
        size(temp1);
        temp2 = zeros(N-TlI(i)+1,1);
        W = [temp1;w;temp2];
        y(:,i) = y(:,i).*W;
    end