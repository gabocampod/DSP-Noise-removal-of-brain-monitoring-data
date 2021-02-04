[d,c]=butter(3,50/250,'low');
[da,ca]=butter(3,[50/250 50.2/250],'stop');
[b,a]=butter(3,0.5/250,'high');

lower = 5001;       
upper = 25000; 

load('Test_Sham_ALPHA_PH.mat') 
sham_tacs_t1 = 15053;   %%CHANGE sham boundary
sham_tacs_t2 = sham_tacs_t1 + 29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));


load('Test_250uA_10Hz_ALPHA_PH.mat') %%CHANGE data !!!!!!!!!!
tacs_t1 = 15040+1;        %%CHANGE data boundary !!!!!!!!!!
tacs_t2 = tacs_t1 + 29999; 

RAW_DATA = EEG.Data(1,tacs_t1:tacs_t2);

load('PYTHON_ALPHA_250uA_10HZ_100I.mat')  %%CHANGE RESULTS!!

modes = PYTHON_ALPHA_250uA_10HZ_100I;

figure(1);
for i = 1:1:5
    subplot(3,2,i)
    plot(modes(i,:))
end 

modes_original_3 = modes(3,:);
modes_original_2 = modes(2,:);
modes_original_4 = modes(4,:);

[n,o] = butter(3,[9.8/250 10.2/250],'stop');
modes(3,:) = filtfilt(n,o, modes(3,:));
modes(2,:) = filtfilt(n,o, modes(2,:));
%modes(4,:) = filtfilt(n,o, modes(4,:));

FREE_DATA = RAW_DATA - modes_original_3(1,:) + modes(3,:);
FREE_DATA = FREE_DATA - modes_original_2(1,:) + modes(2,:);
FREE_DATA = FREE_DATA - modes_original_4(1,:) + modes(4,:);

FILTERED_FREE = filtfilt(b,a,filtfilt(d,c,FREE_DATA));

figure(4);
plot(ERP_Sham_SMA(1,:), 'b')
hold on
plot(FILTERED_FREE(1,:), 'r')

%%
RONLY4 = corrcoef(ERP_Sham_SMA(1, lower: upper), FILTERED_FREE(1,lower:upper));
Coefficient = RONLY4(1,2);
%%
RMS = 10^-3.*rms(ERP_Sham_SMA(1,lower:upper) - FILTERED_FREE(1,lower:upper));
%%
noise = (ERP_Sham_SMA(1,lower:upper) - FILTERED_FREE(1,lower:upper));
SNR_ONLY4 = (rms(FILTERED_FREE(1,lower:upper) / rms(noise) ))^2;
SNR_ONLY4_DB = db(SNR_ONLY4, 'power');
%%
%%Coherence
[cxy,f2] = mscohere(ERP_Sham_SMA(1,lower:upper), FILTERED_FREE(1,lower:upper), [],[],[],500);

%figure(3)
%plot(f2, cxy);

%%
%%COHERENCE VALUE

Fs = 500;            % Sampling frequency                         
L = 20000;             % Length of signal
FFT_SHAM_FILTERED = fft(ERP_Sham_SMA(1,lower:upper));
FFT_FREE_FILTERED = fft(FILTERED_FREE(1,lower:upper));

P2 = abs(FFT_SHAM_FILTERED/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P2_2 = abs(FFT_FREE_FILTERED/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

R_FFT = corrcoef(P2, P2_2);
Coherence_value = R_FFT(1,2);

R_FFT = corrcoef(FFT_SHAM_FILTERED, FFT_FREE_FILTERED);
Coherence_value_2ndoption = R_FFT(1,2);

%%
%%RUN IN ARTIFACT
low_in = 1;
high_in = 500;
similarity_in = 0;
treshhold = 5; %db
while similarity_in < treshhold
    
    noise_in = (ERP_Sham_SMA(1,low_in:high_in) - FILTERED_FREE(1,low_in:high_in));
    SNR_in= (rms(FILTERED_FREE(1,low_in:high_in) / rms(noise_in) ))^2;
    SNR_in_db = db(SNR_in, 'power');
    similarity_in = SNR_in_db;
    
    low_in = low_in +1;
    high_in = high_in +1;
end

RUN_IN_ARTIFACT = (low_in - 1)/500;


%%
%%RUN OUT ARTIFACT
low_out = 29501;
high_out = 30000;
similarity_out = 0;

while similarity_out < treshhold
    
    noise_out = (ERP_Sham_SMA(1,low_out:high_out) - FILTERED_FREE(1,low_out:high_out));
    SNR_out= (rms(FILTERED_FREE(1,low_out:high_out) / rms(noise_out) ))^2;
    SNR_out_db = db(SNR_out, 'power');
    similarity_out = SNR_out_db;
    
    low_out = low_out - 1;
    high_out = high_out - 1;
end
high_out_tt = high_out + 1;
high_out = 30000 - high_out + 1;

RUN_OUT_ARTIFACT = (high_out - 1)/500;

    noise_in = (ERP_Sham_SMA(1,2500:3500) - FILTERED_FREE(1,2500:3500));
    SNR_in= (rms(FILTERED_FREE(1,2500:3500) / rms(noise_in) ))^2;
    SNR_in_db = db(SNR_in, 'power');







