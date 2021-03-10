clc
clear all
clear all

datanames %importing all the data informatiton

%% output to text file
compiled_analysispath = '/Users/athula/Dropbox/Experiments/TesserScan/Behavioral Results';
cd(compiled_analysispath)
    
mean_name = sprintf('Induct_Results_Mean_Overall'); 
mean_fid=fopen([mean_name '.txt'], 'w'); % open the file
fprintf(mean_fid, 'SUBJECT\tPrim_Acc\tPrim_RT\tBound1_Acc\tBound1_RT\tBound2_Acc\tBound2_RT\tBound_Acc\tBound_RT\tOverall_Acc\tOverall_RT\n');

biasname = sprintf('Induct_Results_Bias_Overall'); 
bias_fid=fopen([biasname '.txt'], 'w'); % open the file
fprintf(bias_fid, 'SUBJECT\tPrim_Bias\tBound1_Bias\tBound2_Bias\tBound_Bias\tOverall_Bias\n');

rtname = sprintf('Induct_Results_RT_Overall'); 
rt_fid=fopen([rtname '.txt'], 'w'); % open the file
fprintf(rt_fid, 'SUBJECT\tPrim_Corr_RT\tPrim_Incorr_RT\tBound1_Corr_RT\tBound1_Incorr_RT\tBound2_Corr_RT\tBound2_Incorr_RT\tBound_Corr_RT\tBound_Incorr_RT\tOverall_Corr_RT\tOverall_Incorr_RT\n');

induct_means_all = []; 
induct_bias_all = []; 
induct_rt_all = [];

for c = 1:length(inductmatload_all)

%change directory
cd(datapaths_all{c});

all = load(inductmatload_all{c}); %change this only 
subj = all.data.subjNum; 
induct_mat = all.data.rundata;

%Col1: SubjNum
%Col2: TrialNum
%Col3: QuestType
    %1 = Prim
    %2 = Bound1
    %3 = Bound2
%Col4: Environment
    %1 = ocean
    %2 = desert
    %3 = forest
%Col5: CueObject#
%Col6: Option1Object#
%Col7: Option2Object#
%Col8: Resp
%Col9: Acc
%Col10: RT


comm1 = [1,2,3,18,19,20,21]; %OCEAN
comm1pri = [1,2,19,20,21];
comm1bou = [18,3];

comm2 = [4,5,6,7,8,9,10]; %DESERT
comm2pri = [5,6,7,8,9];
comm2bou = [4,10];

comm3 = [11,12,13,14,15,16,17]; %FOREST
comm3pri = [12,13,14,15,16];
comm3bou = [11,17];

hamforcomm = [18, 4, 11];
hambckcomm = [3, 17, 10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate means by each question type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subj = nanmean(induct_mat(:, 1));
prim = induct_mat(induct_mat(:, 3) == 1, :);
prim_acc = nanmean(prim(:, 9));
prim_rt = nanmean(prim(:, 10));
    
bound1 = induct_mat(induct_mat(:, 3) == 2, :);
bound1_acc = nanmean(bound1(:, 9));
bound1_rt = nanmean(bound1(:, 10));
    
bound2 = induct_mat(induct_mat(:, 3) == 3, :);
bound2_acc = nanmean(bound2(:, 9));
bound2_rt = nanmean(bound2(:, 10));
    
%accumulating all the boundary trials (12 total)
bound = [bound1; bound2];
bound_acc = nanmean(bound(:, 9)); 
bound_rt = nanmean(bound(:, 10)); 

%overall 
overall_acc = nanmean(induct_mat(:, 9)); 
overall_rt = nanmean(induct_mat(:, 10)); 
    
induct_means = [subj, prim_acc, prim_rt, bound1_acc, bound1_rt, bound2_acc, bound2_rt, bound_acc, bound_rt, overall_acc, overall_rt];
induct_means_all = [induct_means_all; induct_means];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate bias scores/RT + temporally congruent RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OVERALL BIAS + RT 
overall_corr = induct_mat(induct_mat(:, 9) == 1, :); 
overall_incorr = induct_mat(induct_mat(:, 9) == 0, :); 

overall_corr_count = size(overall_corr, 1); 
overall_incorr_count = size(overall_incorr, 1); 

overall_corr_prop = overall_corr_count/size(induct_mat, 1);
overall_incorr_prop = overall_incorr_count/size(induct_mat, 1);

overall_bias_total = overall_corr_count - overall_incorr_count; 
overall_bias = overall_bias_total/size(induct_mat, 1); 

overall_corr_rt = nanmean(overall_corr(:, 10)); 
overall_incorr_rt = nanmean(overall_incorr(:, 10));
overall_bias_rt = overall_corr_rt - overall_incorr_rt; 
% 

%PRIMARY BIAS + RT
prim_corr = prim(prim(:, 9) == 1, :) ; 
prim_incorr = prim(prim(:, 9) == 0, :); 

prim_corr_count = size(prim_corr, 1); 
prim_incorr_count = size(prim_incorr, 1); 

prim_bias_total = prim_corr_count - prim_incorr_count; 
prim_bias = prim_bias_total/size(prim, 1); 

prim_corr_rt = nanmean(prim_corr(:, 10)); 
prim_incorr_rt = nanmean(prim_incorr(:, 10));
prim_bias_rt = prim_corr_rt - prim_incorr_rt; 
%

%BOUND-1 BIAS + RT
bound1_corr = bound1(bound1(:, 9) == 1, :) ; 
bound1_incorr = bound1(bound1(:, 9) == 0, :); 

bound1_corr_count = size(bound1_corr, 1); 
bound1_incorr_count = size(bound1_incorr, 1); 

bound1_bias_total = bound1_corr_count - bound1_incorr_count; 
bound1_bias = bound1_bias_total/size(bound1, 1); 

bound1_corr_rt = nanmean(bound1_corr(:, 10)); 
bound1_incorr_rt = nanmean(bound1_incorr(:, 10));
bound1_bias_rt = bound1_corr_rt - bound1_incorr_rt; 
%

%BOUND-2 BIAS + RT
bound2_corr = bound2(bound2(:, 9) == 1, :) ; 
bound2_incorr = bound2(bound2(:, 9) == 0, :); 

bound2_corr_count = size(bound2_corr, 1); 
bound2_incorr_count = size(bound2_incorr, 1); 

bound2_bias_total = bound2_corr_count - bound2_incorr_count; 
bound2_bias = bound2_bias_total/size(bound2, 1); 

bound2_corr_rt = nanmean(bound2_corr(:, 10)); 
bound2_incorr_rt = nanmean(bound2_incorr(:, 10));
bound2_bias_rt = bound2_corr_rt - bound2_incorr_rt; 
%

%BOUND OVERALL BIAS + RT
bound_corr = bound(bound(:, 9) == 1, :) ; 
bound_incorr = bound(bound(:, 9) == 0, :); 

bound_corr_count = size(bound_corr, 1); 
bound_incorr_count = size(bound_incorr, 1); 

bound_bias_total = bound_corr_count - bound_incorr_count; 
bound_bias = bound_bias_total/size(bound, 1); 

bound_corr_rt = nanmean(bound_corr(:, 10)); 
bound_incorr_rt = nanmean(bound_incorr(:, 10));
bound_bias_rt = bound_corr_rt - bound_incorr_rt; 
%

induct_bias = [subj, prim_bias, bound1_bias, bound2_bias, bound_bias, overall_bias];
induct_bias_all = [induct_bias_all; induct_bias];

induct_rt = [subj, prim_corr_rt, prim_incorr_rt, bound1_corr_rt, bound1_incorr_rt, bound2_corr_rt, bound2_incorr_rt, bound_corr_rt, bound_incorr_rt, overall_corr_rt, overall_incorr_rt];
induct_rt_all = [induct_rt_all; induct_rt]; 

cd(compiled_analysispath)
  
fprintf(mean_fid, '%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n',...
    subj, prim_acc, prim_rt, bound1_acc, bound1_rt, bound2_acc, bound2_rt, bound_acc, bound_rt, overall_acc, overall_rt); 

fprintf(bias_fid, '%d\t %d\t %d\t %d\t %d\t %d\n',...
    subj, prim_bias, bound1_bias, bound2_bias, bound_bias, overall_bias); 

fprintf(rt_fid, '%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n',...
    subj, prim_corr_rt, prim_incorr_rt, bound1_corr_rt, bound1_incorr_rt, bound2_corr_rt, bound2_incorr_rt, bound_corr_rt, bound_incorr_rt, overall_corr_rt, overall_incorr_rt); 
  

end
close all

%% output to mat file
save(sprintf('Induct_Results_Mean_Overall'),'induct_means_all');
save(sprintf('Induct_Results_Bias_Overall'),'induct_bias_all');
save(sprintf('Induct_Results_RT_Overall'),'induct_rt_all');