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
%%
%%original mode 4
load('CEEMDAN_FAST_20000_50_500_002_RESULTS')
load('CEEMDAN_of_mode_4_results');

modes = CEEMDAN_FAST_20000_50_500_002;
modes_of_mode = CEEMDAN_of_mode_4;

figure(1);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end 

figure(33);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes_of_mode(i,:))
end 

Y = fft(modes_of_mode(3,:));
Fs = 500;
L = 20000;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(55);
plot(f,P1) 
title('MODE ORIGINAL')

%%
%%filtering mode 4
[n,o] = butter(3,[4.8/250 5.2/250],'stop');

modes_of_mode(4,:) = filtfilt(n,o, modes_of_mode(4,:));
modes_of_mode(3,:) = filtfilt(n,o, modes_of_mode(3,:));
modes_of_mode(5,:) = filtfilt(n,o, modes_of_mode(5,:));

modes(4,:) = zeros(1,20000);

for i= 1:1:11
        modes(4,:) = modes(4,:) + modes_of_mode(i,:);
end

Y = fft(modes(4,:));
Fs = 500;
L = 20000;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(11);
plot(f,P1)
xlim([0 8])
title('MODE')

%% Reconstructed
reconstructed = zeros(1,20000);

for i= 1:1:11
        reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end

reconstructed_artif_free = reconstructed(1,:);

FREE_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_artif_free));

%% SHAM RECONSTRUCTED

load('CEEMDAN_SHAM_PHANTOM_50_500_002_results.mat')

modes_of_sham = CEEMDAN_SHAM_PHANTOM_50_500_002;

reconstructed_SHAM(1,:) = zeros(1,20000);

for i= 1:1:12
        reconstructed_SHAM(1,:) = reconstructed_SHAM(1,:) + modes_of_sham(i,:);
end

SHAM_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_SHAM));

Y2 = fft(ERP_Sham_SMA(1, 5001:25000));
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
f = Fs*(0:(L/2))/L;
%figure(9);
%plot(f,P12, 'k')
%hold on
%plot(f,P1, 'b' )
%title('SHAM VS RECONSTRUCTED after filter')



%% raw data filtered
load('ERP_sham_test1_test2_cc.mat')

RAW_DATA = ERPTEST1.E1(1, 5001:25000);
data = filtfilt(nb,ob,RAW_DATA);
data_filtered =  filtfilt(b,a,filtfilt(d,c,data));

Y2 = fft(data_filtered);
P222 = abs(Y2/L);
P122 = P222(1:L/2+1);
P122(2:end-1) = 2*P122(2:end-1);
f = Fs*(0:(L/2))/L;
figure(19);
plot(f,P12, 'b')
hold on
plot(f,P1, 'k' )
hold on
plot(f, P122, 'r')
xlim([0 8])
title('ALL TOGETHER')

%% COMPARISSON PLOT
figure(3);
plot(ERP_Sham_SMA(1, 5001: 25000))
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

%figure(2)
%plot(f2, cxy);







