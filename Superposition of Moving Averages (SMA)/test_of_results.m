
load('SUBJECT_RESULTS.mat');
load('PYTHON_SUBJECT1_5HZ.mat');
load('PYTHON_SUBJECT2_10HZ.mat');

MATLAB_SHAM_SUBJECT1 = SUBJECT1_RESULTS.SHAM(1,:);
PYTHON_DATA_SUBJECT1 = PYTHON_SUBJECT1_5HZ(1,:);

figure(1);
plot(MATLAB_SHAM_SUBJECT1);
hold on
plot(PYTHON_DATA_SUBJECT1);

R1 = corrcoef(MATLAB_SHAM_SUBJECT1(1, 1720 : 6719) , PYTHON_DATA_SUBJECT1(1, 1720 : 6719));

MATLAB_SHAM_SUBJECT2_10HZ = SUBJECT2_10HZ_RESULTS.SHAM(1,:);
PYTHON_DATA_SUBJECT2_10HZ = PYTHON_SUBJECT2_10HZ(6,:);

figure(2);
plot(MATLAB_SHAM_SUBJECT2_10HZ);
hold on
plot(PYTHON_DATA_SUBJECT2_10HZ);

R2 = corrcoef(FILTERED_FREE_EEG_SUBJECT(1,:),PYTHON_DATA_SUBJECT2_10HZ(1,:));




