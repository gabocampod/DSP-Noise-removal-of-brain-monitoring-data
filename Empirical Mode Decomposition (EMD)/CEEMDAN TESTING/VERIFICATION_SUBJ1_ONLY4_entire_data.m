[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5001;       
upper = 25000; 

load('Test_Sham_ERP_SUBJ1.mat')  

ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15015:45014)));

load('CEEMDAN_ONLY4_SUBJ1_ENTIRE.mat')

modes_only4 = CEEMDAN_ONLY4_SUBJ1_ENTIRE;

figure(1);
for i = 1:1:6
    subplot(3,2,i)
    plot(modes_only4(i,:))
end 

load('SUBJECT_RESULTS.mat')

RAW_DATA = SUBJECT1_RESULTS.DATA5HZ(1,:);

modes_only4_original = modes_only4(4,:);
modes_only4_original_5 = modes_only4(5,:);

[n,o] = butter(3,[4.8/250 5.2/250],'stop');
modes_only4(4,:) = filtfilt(n,o, modes_only4(4,:));
modes_only4(5,:) = filtfilt(n,o, modes_only4(5,:));

FREE_DATA = RAW_DATA - modes_only4_original(1,:) + modes_only4(4,:);
FREE_DATA = FREE_DATA - modes_only4_original_5(1,:) + modes_only4(5,:);

FILTERED_FREE = filtfilt(b,a,filtfilt(d,c,FREE_DATA));

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

%%
%%RUN IN ARTIFACT
low_in = 1;
high_in = 500;
similarity_in = 0;
treshhold = 0.8812; %%CHANGE!!!
while similarity_in < treshhold
    
    R_in = corrcoef(ERP_Sham_SMA(1,low_in:high_in), FILTERED_FREE(1,low_in:high_in));
    similarity_in = R_in(1,2);
    
    low_in = low_in +1;
    high_in = high_in +1;
end

%%
%%RUN OUT ARTIFACT
low_out = 29501;
high_out = 30000;
similarity_out = 0;

while similarity_out < treshhold  
    
    R_out = corrcoef(ERP_Sham_SMA(1,low_out:high_out), FILTERED_FREE(1,low_out:high_out));
    similarity_out = R_out(1,2);
    
    low_out = low_out - 1;
    high_out = high_out - 1;
end

high_out = 30000 - high_out;

figure(2);
plot(ERP_Sham_SMA(1,:))
hold on
plot(FILTERED_FREE(1,:))

