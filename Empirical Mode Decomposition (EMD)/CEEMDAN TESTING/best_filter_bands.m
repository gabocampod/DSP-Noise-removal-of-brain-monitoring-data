[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 2001;       %16000 long
upper = 18000; 

lower_sham = 7001;  %16000
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

initial_mode4 = modes(4,:);
coefficient_loop_best = 0;
SNR_LOOP_BEST = 0;

for H = 0:0.1:0.9
    k_low = 4+H;
    k_high = 6-H;
    
    [n,o] = butter(3,[(k_low/250) (k_high/250)],'stop');
    modes(4,:) = filtfilt(n,o,initial_mode4);
    
    reconstructed_loop = zeros(1,20000);

    for i= 1:1:11
        reconstructed_loop(1,:) = reconstructed_loop(1,:) + modes(i,:);
    end
    
    FREE_after_filter_loop = filtfilt(b,a,filtfilt(d,c,reconstructed_loop));
    
    R1_loop = corrcoef(ERP_Sham_SMA(1, lower_sham: upper_sham), FREE_after_filter_loop(1,lower:upper));
    Coefficient_loop = R1_loop(1,2);
    
    noise = (ERP_Sham_SMA(1,lower_sham:upper_sham) - FREE_after_filter_loop(1,lower:upper));
    SNR = (rms(FREE_after_filter_loop(1,lower:upper) / rms(noise) ))^2;
    
    
    if SNR > SNR_LOOP_BEST
        SNR_LOOP_BEST = SNR;
        SNR_DB_best = db(SNR, 'power');
        SNR_low_best_loop = k_low;
        SNR_high_best_loop = k_high;
    end    
    
    if Coefficient_loop > coefficient_loop_best
        coefficient_loop_best = Coefficient_loop;
        coefficient_low_best_loop = k_low;
        coefficient_high_best_loop = k_high;
    end    
    
end