clc
clear all
clear all

datanames %importing all the data informatiton
matloadnames = struct_pt2_matload_all;

outfname = sprintf('Struct_acc_rt_dprime_Part2_OVERALL');
fid=fopen([outfname '.txt'], 'w'); % open the file
fprintf(fid, 'SUBJCOUNT\t SUBJECT\t overall_acc\t overall_rt\t overall_dprime\n'); %Create header

struct_acc_rt_ALL = [];
subjcount = 0;
for c = 1:length(struct_pt2_matload_all)
    
    %change directory
    cd(datapaths_all{c});
    
    subjcount = subjcount + 1;
    struct_data_name = matloadnames{c};
    loaded = load(struct_data_name); %change this only
    struct = loaded.all_tog;
    subj = all(1, 1); 
    
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
    
    %Figuring out where the boundary objects are in the structure
    boundall = [];
    for i = 1:length(struct(:, 2))
        if ismember(struct(i, 4), comm1bou) || ismember(struct(i, 4), comm2bou) || ismember(struct(i, 4), comm3bou)
            bound = 1;
        else
            bound = 0;
        end
        
        boundall = [boundall; bound];
    end
    
    %append struct with bound category
    struct = [struct, boundall];
    %Col1: Subj#
    %Col2: Block#
    %Col3: Trial#
    %Col4: WalkType (1 = struct; 2 = scrambled)
    %Col5: Object
    %Col6: Rot/2back appearance (1/0)
    %0 = (ROT) 2back
    %1 = (NOTROT) none
    %Col7: RotResp(1)/NotRotResp(2)/NoResp(0)
    %1 = RespNot
    %2 = RespRot
    %0 = NoResp
    %Col8: Rot Accuracy (1/0)
    %1 = Correct
    %0 = Incorrect
    %Col9: Rot RT
    %Col10: Boundary Item or Not
    %1 = Boundary
    %0 = Non-Boundary
    
    %From Karuza et al., 2017
    %1. -eliminated all incorrect trials
    %2. -then eliminated trials with rotated objects
    %-removed implausible RT
    %-removed outlier data points from each subject constituting >3 SD from M
    
    %to calculate d-prime of Hits (correctly detecting rotated) and
    %FA (incorrectly saying a non-rotated item is rotated)
    
    hit_count = [];
    fa_count = [];
    for d = 1:size(struct, 1)
        rot_status = struct(d, 6); %1 = not rot; 0 = rot
        rot_resp = struct(d, 7);   %1 = not rot resp; 2 = rot resp; 0 = no resp
        
        if rot_status == 0 && rot_resp == 2 %rotated item, rotated resp
            %how many hits of rotated items?
            hit_count = [1; hit_count];
        elseif rot_status == 1 && rot_resp == 2 %non-rotated item, rotated resp
            %how many were false alarms to non-rotated items?
            fa_count = [1; fa_count];
        end
    end
    
    %first, total number of rotated items, hits, FA
    total_rot = sum(struct(:, 6) == 0);
    total_notrot = sum(struct(:, 6) == 1);
    total_hit = sum(hit_count);
    total_fa = sum(fa_count);
    
    %hit prop
    hit_prop = total_hit/total_rot;
    hit_z = norminv(hit_prop);
    
    %Adjust only the extreme values by replacing rates of 0 with 0.5/n
    %and rates of 1 with (n?0.5)/n
    %where n is the number of signal or noise trials (Macmillan & Kaplan, 1985)
    if hit_prop == 0
        hit_z = 0.5/total_rot;
    elseif hit_prop == 1
        hit_z = (total_rot - 0.5)/total_rot;
    end
    
    %fa prop
    fa_prop = total_fa/total_notrot;
    fa_z = norminv(fa_prop);
    
    if fa_prop == 0
        fa_z = 0.5/total_notrot;
    elseif fa_prop == 1
        fa_z = (total_rot - 0.5)/total_notrot;
    end
    
    %overall accuracy
    overall_acc_avg = nanmean(struct(:, 8));
    
    %overall RT
    overall_rt_avg = nanmean(struct(:, 9));
    
    %ovrall d-prime
    d_prime = hit_z - fa_z;
    
    struct_acc_rt = [subjcount, subj, overall_acc_avg, overall_rt_avg, d_prime];
    struct_acc_rt_ALL = [struct_acc_rt_ALL; struct_acc_rt];
    
    cd('/Users/athula/Desktop/Tesser/8_Tesser_Morph/version2/results/results_2019/compiled_results')
    
    %% Overall output to text file
    fprintf(fid, '%d\t %d\t %d\t %d\t %d\n',...
        subjcount, subj, overall_acc_avg, overall_rt_avg, d_prime);
    
end

close all

%also save out a .mat file of the same .txt file
save(sprintf('Struct_acc_rt_dprime_Part2_OVERALL'),'struct_acc_rt_ALL');