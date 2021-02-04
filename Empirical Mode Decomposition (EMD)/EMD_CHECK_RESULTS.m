load('ERP_sham_test1_test2_cc')
load('CEEMDAN_FAST_20000_50_500_002_RESULTS')
load('CEEMDAN_FAST_20000points_50_500_02_RESULTS')
load('CEEMDAN_of_mode_4_results')
load('CEEMDAN_FAST_TEST2_50_300_002_results.mat')

[d,c]=butter(3,25/250,'low');
[b,a]=butter(3,0.5/250,'high');

X = ERPSHAM.ES(1,5000:24999);
%% ERP DATA MODES PLOT

modes = CEEMDAN_FAST_20000_50_500_002;

figure(1);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end 

reconstructed= zeros(1,20000);

for i= 1:1:11
    reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end


%% Create reconstructed_free and filters
reconstructed_artif_free = reconstructed(1,:) - modes(4,:);

SHAM_after_filter = filtfilt(b,a,filtfilt(d,c,X));
FREE_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_artif_free));
%% Plots for filters

figure(3);
plot(SHAM_after_filter)
hold on
plot(FREE_after_filter)
hold off
title('AFTER FILTER')

similarity_best = 0;

for i= 0:1:15000
    low = 1 + i;
    high = 5000 + i;
    
    R = corrcoef(SHAM_after_filter(1, low: high), FREE_after_filter(1, low: high));
    similarity = R(1,2);
    
    if similarity > similarity_best
        similarity_best = similarity;
        low_best = low;
        high_best = high;
    end
end
 
R_most = corrcoef(SHAM_after_filter(1, 2000:18000 ), FREE_after_filter(1,2000:18000));
similarity_most = R_most(1,2);

figure(3);
plot(SHAM_after_filter(1, low_best: high_best))
hold on
plot(FREE_after_filter(1, low_best: high_best))
hold off
xlim([0 5000])
title('AFTER FILTER BEST 5 SEG')

figure(4);
plot(SHAM_after_filter(1, 2000:18000 ))
hold on
plot(FREE_after_filter(1, 2000:18000 ))
hold off
title('all')
