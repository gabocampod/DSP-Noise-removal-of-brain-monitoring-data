[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 2001;      %16000 long
upper = 18000; 

lower_sham = 7001; %16000
upper_sham = 23000;

load('Test_Sham_ERP_SUBJ1.mat')  

ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15015:45014)));

load('CEEMDAN_SUBJ1_50_500_002_results.mat')

modes = CEEMDAN_SUBJ1_50_500_002_results;

figure(1);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end 

[n,o] = butter(3,[4.8/250 5.2/250],'stop');
modes(4,:) = filtfilt(n,o,modes(4,:));

reconstructed = zeros(1,20000);

for i= 1:1:11
        reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end

reconstructed_artif_free = reconstructed(1,:);

FREE_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_artif_free));

figure(3);
plot(ERP_Sham_SMA(1,5001 : 25000 ))
hold on
plot(FREE_after_filter)
hold off
title('AFTER FILTER')


%%
%%COEFFICIENT
R1 = corrcoef(ERP_Sham_SMA(1, 7000:23000 ), FREE_after_filter(1,2000:18000));
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
[cxy,f] = mscohere(ERP_Sham_SMA(1,lower_sham:upper_sham), FREE_after_filter(1,lower:upper), [],[],[],500);

figure(2)
plot(f, cxy);





