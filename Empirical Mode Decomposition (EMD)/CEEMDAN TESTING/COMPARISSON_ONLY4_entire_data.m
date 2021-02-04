[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5001;       %16000 long
upper = 25000; 

load('CEEMDAN_ONLY4_TEST1_ENTIRE.mat')

modes_only4 = CEEMDAN_ONLY4_TEST1_ENTIRE;

figure(1);
for i = 1:1:6
    subplot(3,2,i)
    plot(modes_only4(i,:))
end 


load('ERP_sham_test1_test2_cc.mat')

RAW_DATA = ERPTEST1.E1(1,:);

modes_only4_original = modes_only4(4,:);
modes_only4_original_5 = modes_only4(5,:);

[n,o] = butter(3,[4.8/250 5.2/250],'stop');
modes_only4(4,:) = filtfilt(n,o, modes_only4(4,:));
modes_only4(5,:) = filtfilt(n,o, modes_only4(5,:));

FREE_DATA = RAW_DATA - modes_only4_original(1,:) + modes_only4(4,:);
FREE_DATA = FREE_DATA - modes_only4_original_5(1,:) + modes_only4(5,:);

FILTERED_FREE = filtfilt(b,a,filtfilt(d,c,FREE_DATA));

load('SMA_TEST_ERP_Sham.mat')

Sham_t1=128+15000;
Sham_t2=Sham_t1+29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data_tACS_AF(1,Sham_t1:Sham_t2))); %30000 long

figure(2);
plot(ERP_Sham_SMA)
hold on
plot(FILTERED_FREE)


%%
RONLY4 = corrcoef(ERP_Sham_SMA(1, lower: upper), FILTERED_FREE(1,lower:upper));
Coefficient = RONLY4(1,2);
%%
RMS = 10^-3.*rms(ERP_Sham_SMA(1,lower:upper) - FILTERED_FREE(1,lower:upper));
%%
noise = (ERP_Sham_SMA(1,lower:upper) - FILTERED_FREE(1,lower:upper));
SNR_ONLY4 = (rms(FILTERED_FREE(1,lower:upper) / rms(noise) ))^2;
SNR_DB_ONLY4 = db(SNR_ONLY4, 'power');
%%
%%Coherence
[cxy,f2] = mscohere(ERP_Sham_SMA(1,lower:upper), FILTERED_FREE(1,lower:upper), [],[],[],500);

figure(3)
plot(f2, cxy);

