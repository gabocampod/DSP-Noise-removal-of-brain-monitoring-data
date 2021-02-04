[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5001;       %16000 long
upper = 25000; 

load('Test_Sham_ERP_SUBJ2.mat')  

ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15028:45027)));

load('CEEMDAN_ONLY4_SUBJ2_40HZ_ENTIRE.mat')

modes_only = CEEMDAN_ONLY4_SUBJ2_40HZ_ENTIRE;

figure(1);
for i = 1:1:4
    subplot(2,2,i)
    plot(modes_only(i,:)/100000)
end 

load('SUBJECT_RESULTS.mat')

RAW_DATA = SUBJECT2_40HZ_RESULTS.DATA2_40HZ(1,:);

modes_only4_original_2 = modes_only(2,:);
%modes_only4_original_3 = modes_only(3,:);
modes_only4_original_4 = modes_only(4,:);

Y = fft(modes_only4_original_2);
Fs = 500;            % Sampling frequency                         
L = 30000;             % Length of signal
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(60);
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


[n,o] = butter(3,[39.8/250 40.2/250],'stop');
modes_only(2,:) = filtfilt(n,o, modes_only(2,:));
%modes_only(3,:) = filtfilt(n,o, modes_only(3,:));
modes_only(4,:) = filtfilt(n,o, modes_only(4,:));

FREE_DATA = RAW_DATA - modes_only4_original_2(1,:) + modes_only(2,:);
%FREE_DATA = FREE_DATA - modes_only4_original_3(1,:) + modes_only(3,:);
FREE_DATA = FREE_DATA - modes_only4_original_4(1,:) + modes_only(4,:);

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
treshhold = 0.8492; %%CHANGE!!!
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



