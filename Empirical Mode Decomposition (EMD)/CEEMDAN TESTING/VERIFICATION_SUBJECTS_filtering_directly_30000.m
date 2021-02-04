[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[n,o] = butter(3,[39.8/250 40.2/250],'stop');

lower = 5001;       %16000 long
upper = 25000; 

load('SUBJECT_RESULTS.mat')

RAW_DATA = SUBJECT2_40HZ_RESULTS.DATA2_40HZ(1, :);

Y = fft(RAW_DATA);
Fs = 500;
L = 20000;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(55);
plot(f,P1) 
xlim([0 80])
title('MODE ORIGINAL')

data = filtfilt(n,o,RAW_DATA);
FILTERED_FREE =  filtfilt(b,a,filtfilt(d,c,data));

load('Test_Sham_ERP_SUBJ2.mat')  

ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15028:45027)));

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

