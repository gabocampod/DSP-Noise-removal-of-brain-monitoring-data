[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5000;
upper = 24999;

load('SMA_TEST_ERP_Sham.mat')

Sham_t1=128+15000;
Sham_t2=Sham_t1+29999;

RAW_SHAM = EEG.Data_tACS_AF(1,Sham_t1:Sham_t2);

ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data_tACS_AF(1,Sham_t1:Sham_t2)));

load('ERP_sham_test1_test2_cc.mat');

RAW_DATA = ERPTEST1.E1(1,:);

Fs = 500;            % Sampling frequency                         
L = 30000;             % Length of signal
Y = fft(RAW_DATA);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
L=20000;
Y = fft(ERP_Sham_SMA(5000:25000));      %SHAM
P23 = abs(Y/L);
P13 = P23(1:L/2+1);
P13(2:end-1) = 2*P13(2:end-1);


[l,h] = butter(3,[4/250 6/250]);

data_filtered_to_send = filtfilt(l,h,RAW_DATA);



load('NEW_APPROACH_TEST1_RESULTS');

modes_freq = NEW_APPROACH_TEST1_RESULTS;

figure(3);
for i = 1:1:9
    subplot(3,3,i)
    plot(modes_freq(i,:))
end 

modes_freq(1,:) = zeros(1,30000);

reconstructed_freq = zeros(1,30000);

for i = 1:1:9
    reconstructed_freq(1,:) = modes_freq(i,:);
end 

data_free = RAW_DATA - data_filtered_to_send + reconstructed_freq;

filtered_data_free = filtfilt(b,a,filtfilt(d,c,data_free));
%%
R1 = corrcoef(ERP_Sham_SMA(1,lower:upper), filtered_data_free(1,lower:upper));
Coefficient_test1 = R1(1,2);

figure(4);
plot(ERP_Sham_SMA)
hold on
plot(filtered_data_free)

L=20000;
Y = fft(filtered_data_free(5000:25000));     %DATA FREE filtered
P24 = abs(Y/L);
P14 = P24(1:L/2+1);
P14(2:end-1) = 2*P14(2:end-1);
f = Fs*(0:(L/2))/L;
figure(64);
plot(f,P14, 'r')
hold on
plot(f, P13, 'k')
xlim([0 100])
ylim([0 20000])
title('FILTERED final and SHAM')



L=30000;
Y = fft(data_filtered_to_send);
P22 = abs(Y/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
f = Fs*(0:(L/2))/L;
figure(70);
plot(f,P12, 'r')
hold on
plot(f,P1) 
hold on
plot(f, P13, 'k')
xlim([3 7])
ylim([0 20000])
title('ORIGINAL and band pass')



