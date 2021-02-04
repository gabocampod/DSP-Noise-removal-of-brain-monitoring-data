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

load('CEEMDAN_SHAM_PHANTOM_50_500_002_results.mat')

modes_of_sham = CEEMDAN_SHAM_PHANTOM_50_500_002;

figure(6);
for i = 1:1:12
    subplot(6,2,i)
    plot(modes_of_sham(i,:))
end 
title('modes sham')

load('CEEMDAN_FAST_20000_50_500_002_RESULTS')
modes = CEEMDAN_FAST_20000_50_500_002;

figure(1);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end 
title('modes data')

mean_mode = mean(modes(4,:));

modes(4,:) = zeros(1,20000);            %set artifact to zero in single level
          

reconstructed = zeros(1,20000);

for i= 1:1:11
        reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end

reconstructed_artif_free = reconstructed(1,:);

FREE_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_artif_free));

figure(3);
plot(ERP_Sham_SMA(1, lower_sham: upper_sham ))
hold on
plot(FREE_after_filter(1, lower: upper))
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
[cxy,f] = mscohere(ERP_Sham_SMA(1,lower_sham:upper_sham), FREE_after_filter(1,lower:upper), [],[],[],500);

figure(2)
plot(f, cxy);





