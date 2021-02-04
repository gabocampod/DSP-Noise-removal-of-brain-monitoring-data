[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 2001;       %16000 long
upper = 18000; 

lower_sham = 7001;  %16000
upper_sham = 23000;

load('SMA_TEST_ERP_Sham.mat')

Sham_t1=128+15000;
Sham_t2=Sham_t1+29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data_tACS_AF(1,Sham_t1:Sham_t2))); %30000 long

load('CEEMDAN_FAST_20000_50_500_002_RESULTS')
load('CEEMDAN_of_mode_4_results.mat')

modes = CEEMDAN_FAST_20000_50_500_002;
modes_of_mode = CEEMDAN_of_mode_4;

figure(1);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end 

figure(5);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes_of_mode(i,:))
end 

Y = fft(modes(4,:));
Fs = 500;
L = 20000;

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(10);
plot(f,P1) 
xlim([2 8])
ylim([0 6000])
title('MODE 4 ORIGINAL')

%modes_of_mode(4,:) = zeros(1, 20000);   %set artifact to zero in multi level
%modes_of_mode(6,:) = zeros(1, 20000);

%for i= 1:1:11
%        modes(4,:) = modes(4,:) + modes_of_mode(i,:);
%end

[n,o] = butter(3,[4/250 6/250],'stop');
modes(4,:) = filtfilt(n,o,modes(4,:));


Y = fft(modes(4,:));
Fs = 500;
L = 20000;

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(11);
plot(f,P1) 
xlim([2 8])
ylim([0 6000])
title('MODE 4 filtered')

reconstructed = zeros(1,20000);

for i= 1:1:11
        reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end

reconstructed_artif_free = reconstructed(1,:);

FREE_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_artif_free));

Y = fft(FREE_after_filter);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(9);
plot(f,P1) 
title('Free After filter')

figure(3);
plot(ERP_Sham_SMA(1, 5001:25000))
hold on
plot(FREE_after_filter(1,:))
hold off
title('AFTER FILTER')


%%
%%COEFFICIENT
R1 = corrcoef(ERP_Sham_SMA(1, lower_sham: upper_sham), FREE_after_filter(1,lower:upper));
Coefficient = R1(1,2);

%%
%%RMS
RMS = 10^-3.*rms(ERP_Sham_SMA(1,lower_sham:upper_sham) - FREE_after_filter(1,lower:upper));

%%
%%SNR 
noise = (ERP_Sham_SMA(1,lower_sham:upper_sham) - FREE_after_filter(1,lower:upper));
SNR = (rms(FREE_after_filter(1,lower:upper) / rms(noise) ))^2;
SNR_DB = db(SNR, 'power');
%%
%%Coherence
[cxy,f2] = mscohere(ERP_Sham_SMA(1,lower_sham:upper_sham), FREE_after_filter(1,lower:upper), [],[],[],500);

figure(2)
plot(f2, cxy);





