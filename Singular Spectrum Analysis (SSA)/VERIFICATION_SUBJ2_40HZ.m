[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5000;
upper = 24999;

load('Test_Sham_ERP_SUBJ2.mat')

ERP_Sham_SMA = filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15028:45027)));

load('singular_result_subj2_40HZ.mat')

SINGULAR_RESULT_SUBJECT2 = free_single_subj2_40hz;

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
%%RUN IN ARTIFACT
low_in = 1;
high_in = 500;
similarity_in = 0;
treshhold = 0.8982;
while similarity_in < treshhold
    
    R_in = corrcoef(ERP_Sham_SMA(1,low_in:high_in), SINGULAR_RESULT_SUBJECT2(1,low_in:high_in));
    similarity_in = R_in(1,2);
    
    low_in = low_in +1;
    high_in = high_in +1;
end

%%
%%RUN OUT ARTIFACT
low_out = 29501;
high_out = 30000;
similarity_out = 0;

while similarity_out < treshhold  %%CHANGE!!!
    
    R_out = corrcoef(ERP_Sham_SMA(1,low_out:high_out), SINGULAR_RESULT_SUBJECT2(1,low_out:high_out));
    similarity_out = R_out(1,2);
    
    low_out = low_out - 1;
    high_out = high_out - 1;
end

high_out = 30000 - high_out;

figure(2);
plot(ERP_Sham_SMA(1,:))
hold on
plot(SINGULAR_RESULT_SUBJECT2(1,:))





