[d,c]=butter(3,30/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5000;
upper = 24999;

load('Test_Sham_ALPHA_PH.mat') 
sham_tacs_t1 = 15048;   
sham_tacs_t2 = sham_tacs_t1 + 29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));

load('EEG_SINGULAR_ERP_1MA_5HZ.mat')

SINGULAR_RESULT_SUBJECT2 = Data_to_save;

%%
%%COEFFICIENT
R1 = corrcoef(ERP_Sham_SMA(1,lower:upper), SINGULAR_RESULT_SUBJECT2(1,lower:upper));
Coefficient_SUBJECT2 = R1(1,2);

%%
%%RMS
RMS_SUBJECT2 = 10^-3.*rms(ERP_Sham_SMA(1,lower:upper) - SINGULAR_RESULT_SUBJECT2(1,lower:upper));

%%
%%SNR 
noise_SUBJECT2 = (ERP_Sham_SMA(1,lower:upper) - SINGULAR_RESULT_SUBJECT2(1,lower:upper));
SNR_SUBJECT2= (rms(SINGULAR_RESULT_SUBJECT2(1,lower:upper) / rms(noise_SUBJECT2) ))^2;
SNR_SUBJECT2_db = db(SNR_SUBJECT2, 'power');
%%
%%Coherence
[cxy_SUBJECT1,f] = mscohere(ERP_Sham_SMA(1,lower:upper), SINGULAR_RESULT_SUBJECT2(1,lower:upper), [],[],[],500);

figure(1)
plot(f, cxy_SUBJECT1);

%%
%%COHERENCE VALUE

Fs = 500;            % Sampling frequency                         
L = 20000;             % Length of signal
FFT_SHAM_FILTERED = fft(ERP_Sham_SMA(1,lower:upper));
FFT_FREE_FILTERED = fft(SINGULAR_RESULT_SUBJECT2(1,lower:upper));

P2 = abs(FFT_SHAM_FILTERED/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P2_2 = abs(FFT_FREE_FILTERED/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

f = Fs*(0:(L/2))/L;

figure(99)

plot(f,P1, 'r') 
hold on 
plot(f, P1_2, 'b')
xlabel('f(Hz)')
ylabel('|P1(f)|')

R_FFT = corrcoef(P2, P2_2);
Coherence_value = R_FFT(1,2);

%%
%%RUN IN ARTIFACT
low_in = 1;
high_in = 500;
similarity_in = 0;
treshhold = 5; %db
while similarity_in < treshhold
    
    noise_in = (ERP_Sham_SMA(1,low_in:high_in) - SINGULAR_RESULT_SUBJECT2(1,low_in:high_in));
    SNR_in= (rms(SINGULAR_RESULT_SUBJECT2(1,low_in:high_in) / rms(noise_in) ))^2;
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
    
    noise_out = (ERP_Sham_SMA(1,low_out:high_out) - SINGULAR_RESULT_SUBJECT2(1,low_out:high_out));
    SNR_out= (rms(SINGULAR_RESULT_SUBJECT2(1,low_out:high_out) / rms(noise_out) ))^2;
    SNR_out_db = db(SNR_out, 'power');
    similarity_out = SNR_out_db;
    
    low_out = low_out - 1;
    high_out = high_out - 1;
end
high_out_tt = high_out + 1;
high_out = 30000 - high_out + 1;

RUN_OUT_ARTIFACT = (high_out - 1)/500;



figure(2);
plot(ERP_Sham_SMA(1,:))
hold on
plot(SINGULAR_RESULT_SUBJECT2(1,:))




