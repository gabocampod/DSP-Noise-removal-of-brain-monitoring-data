[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 2001;       %16000 long
upper = 18000; 

lower_sham = 7001;  %16000
upper_sham = 23000;

load('CEEMDAN_ONLY4_TEST1_PHANTOM_results.mat')

modes_only4 = CEEMDAN_ONLY4_TEST1_PHANTOM;

load('ERP_sham_test1_test2_cc.mat')

RAW_DATA = ERPTEST1.E1(1,5001:25000);

modes_only4_original = modes_only4(4,:);

[n,o] = butter(3,[4.8/250 5.2/250],'stop');
modes_only4(4,:) = filtfilt(n,o, modes_only4(4,:));

FREE_DATA = RAW_DATA - modes_only4_original(1,:) + modes_only4(4,:);

FILTERED_FREE = filtfilt(b,a,filtfilt(d,c,FREE_DATA));

load('SMA_TEST_ERP_Sham.mat')

Sham_t1=128+15000;
Sham_t2=Sham_t1+29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data_tACS_AF(1,Sham_t1:Sham_t2))); %30000 long

%%
RONLY4 = corrcoef(ERP_Sham_SMA(1, lower_sham: upper_sham), FILTERED_FREE(1,lower:upper));
Coefficient = RONLY4(1,2);

%%
noise = (ERP_Sham_SMA(1,lower_sham:upper_sham) - FILTERED_FREE(1,lower:upper));
SNR_ONLY4 = (rms(FILTERED_FREE(1,lower:upper) / rms(noise) ))^2;
SNR_DB_ONLY4 = db(SNR_ONLY4, 'power');

