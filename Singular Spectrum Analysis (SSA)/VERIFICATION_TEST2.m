[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5000;
upper = 24999;

load('ERP_sham_test1_test2_cc')

SHAM = ERPSHAM.ES(1,:);
ERP_Sham_SMA = filtfilt(b,a,filtfilt(d,c,SHAM));

load('singular_result_test2.mat')

SINGULAR_RESULT_SUBJECT2 = free_single_test2;

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




