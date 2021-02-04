[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[n,o] = butter(3,[4.8/250 5.2/250],'stop');

lower = 5001;       %16000 long
upper = 25000; 

load('ERP_sham_test1_test2_cc.mat')

%RAW_DATA = SUBJECT1_RESULTS.DATA5HZ(1, :);
RAW_DATA = ERPTEST2.E2(1, :);
data = filtfilt(n,o,RAW_DATA);
FILTERED_FREE =  filtfilt(b,a,filtfilt(d,c,data));

load('SMA_TEST_ERP_Sham.mat')

Sham_t1=128+15000;
Sham_t2=Sham_t1+29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data_tACS_AF(1,Sham_t1:Sham_t2))); %30000 long


figure(1);
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

