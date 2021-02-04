[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[da,ca]= butter(3,[38/250 40/250], 'stop');


load('session1_SHAM.mat')      %%SHAM DEFINITION
sham_tacs_t1 = 320000;
sham_tacs_t2 = sham_tacs_t1 + 149999;
Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG_mine(1, sham_tacs_t1:sham_tacs_t2)));
%Sham_SMA = filtfilt(da,ca,Sham_SMA);

load('session2_DATA.mat')     %%RAW DATA DEFINITION
tacs_t1 = 90000;        
tacs_t2 = tacs_t1 + 149999; %5 MINUTES (5*60*500)

X_init = EEG_mine_ses2(4,tacs_t1:tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 1000;    
N = length(X);
t = (1:N)';

RAW_FILTERED = filtfilt(b,a,filtfilt(d,c,X_init));

Fs = 500;                       %FFT OF RAW FILTERED                          
L = 150000;            
Y = fft(RAW_FILTERED/100000);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Y = fft(Sham_SMA);             %FFT OF SHAM
P22 = abs(Y/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
f = Fs*(0:(150000/2))/150000;
figure(22);
plot(f,P1, 'k') 
xlabel('f(Hz)')
ylabel('|P1(f)|')
ylim([0 9])

load('RC_1000_FROM90K.MAT') 


figure(29)
for m=1:10
  subplot(5,2,m);
  plot(t,RC(:,m)/10000000,'r-');
  ylabel(sprintf('RC %d',m));
  xlim([61800 62800])
end

figure(120);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 21;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end

figure(140);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 41;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end

figure(160);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 61;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end

figure(180);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 81;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end
%}

RC_TO_QUIT = zeros(150000,1000);
val = zeros(1000);

inner_col = 1;
for i=1:1:length(RC(1,:))
    if max(abs(RC(:,i))) > 40000
        RC_TO_QUIT(:,inner_col) =  RC(:,i);
        inner_col = inner_col + 1;
        val(inner_col) = i;
    end
end

summed = sum(RC_TO_QUIT(:,:),2);
%summed = sum(RC(:,1:4),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));


figure(6);
subplot(2,1,1)
plot(Sham_SMA, 'k');
subplot(2,1,2)
plot(free_filtered, 'r');

figure(7);
plot(free_filtered, 'r');
hold on
plot(Sham_SMA, 'k');

L= 150000;
Fs = 500;
Y = fft(Sham_SMA);                %%FFT OF SHAM
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(150000/2))/150000;

Y = fft(free_filtered);            %%FFT OF FREE FILTERED       
P22 = abs(Y/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
figure(10);
plot(f,P1, 'k') 
hold on
plot(f,P12, 'r') 
title('sham (black) and filtered free (red)')
xlabel('f(Hz)')
ylabel('|P1(f)|')


load('Test_Sham_ALPHA_PH.mat') 
sham_tacs_t1 = 15030;   %%CHANGE sham boundary
sham_tacs_t2 = sham_tacs_t1 + 29999;
Sham_SMA_ALPHA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));


load('Test_Sham_ERP_PH.mat') 
sham_tacs_t1 = 15032;   %%CHANGE sham boundary
sham_tacs_t2 = sham_tacs_t1 + 29999;
Sham_SMA_ERP=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));

L = 30000;
Y = fft(Sham_SMA_ALPHA);
P22 = abs(Y/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
fsham = Fs*(0:(L/2))/L;


Y = fft(Sham_SMA_ERP);
P222 = abs(Y/L);
P122 = P222(1:L/2+1);
P122(2:end-1) = 2*P122(2:end-1);


figure(3);
plot(f, P1, 'k')
hold on
plot(fsham,P12, 'r') 
hold on
plot(fsham,P122, 'b') 
title('SHAM 40 MIN (black )vs SHAM ALPHA (red) vs SHAM ERP ')
xlabel('f(Hz)')
ylabel('|P1(f)|')

figure(9)
plot(Sham_SMA(1,1:30000), 'k')
hold on
plot(Sham_SMA_ALPHA, 'r')
hold on
plot(Sham_SMA_ERP, 'b')


Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG_mine(1, 1:1200000)));
%Sham_SMA = filtfilt(da,ca,Sham_SMA);
L= 1200000;
Fs = 500;
Y = fft(Sham_SMA);                %%FFT OF SHAM
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(12);
plot(f,P1, 'k') 
title('sham (black) full')
xlabel('f(Hz)')
ylabel('|P1(f)|')


