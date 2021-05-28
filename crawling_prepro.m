clear;

% run('../session_list.m');
% 
% Ndir = length(sessions);

in_out_folder = 'D:\Andres';

sessions{1} = '0403E\TUNL\STAGE 1\S3\2 SEC\AUGUST_29_19_H9_M33_S18';

Ndir = 1;

% loop on folders
for d = 1:Ndir
    
    % status used to make sure we have created the files necessary to continue.
    % Note that each of the functions defined and used below only read files,
    % process their content, and then write their output in other files. There
    % are no variables that are used in common with these different steps, and
    % no need for global variables then.
    status.timestamp = NaN;
    status.schedules = NaN;
    status.synchronization = NaN;
    status.trajectory = NaN;
    status.deconvolved = NaN;
    
    my_folder = [in_out_folder '\' sessions{d}];
    
    disp([num2str(d) '\' num2str(Ndir) ' Processing files in folder ' my_folder]);
    
    % process timestamp file
    status.timestamp = process_timestamp(my_folder);
    disp(' ');
    
    % process schedules file
    status.schedules = process_schedules(my_folder);
    disp(' ');
    
    % synchronize the data
    status.synchronization = process_synchro(my_folder,status);
    disp(' ');
    
    % process the trajectory files
    status.trajectory = process_trajectory(my_folder);
    disp(' ');
    
    % process the deconvolved data
    % create msDeconvolved.mat if it does not exist
    if exist([my_folder '\ms.mat'],'file')==1 && exist([my_folder '\msDeconvolved.mat'],'file')==0
        status.deconvolved = process_deconvolved(my_folder);
    else
        status.deconvolved = 1;
    end
    
    % merge the data
    dataset = process_merging(my_folder,status);

    % shifts the trajectory/dcs by any lag between cage and behav video
    % output, corrects vertical and horizontal orientation of the
    % trajectory and normalizes traj to its real proportions.
    if status.timestamp == 1 && status.schedules == 1 && ...
            status.synchronization == 1 && status.trajectory == 1 && ...
            status.deconvolved == 1
        
        disp('Measuring the gaps between the behavioral videos and the cage output.');
        
        % extract the lag
        lags_videos_cage = investigate_houselight_timing_fcn_2(my_folder,dataset);
        
        % shift by trajectory and dcs by the lag computed above and save
        % the new dataset.
        dataset = correct_dataset(dataset,lags_videos_cage);
        
        % display the trajectory and contact points after adjustment.
        check_cage_vs_behav(dataset);
        img = getframe(gcf);
        imwrite(img.cdata,[my_folder '\qualitycheck\summary_cage_traj_synch_with_lag_shift.png']);
        close all
        
        % make a pdf report with the figures that were produced
        disp('making pdf report of all the summary pictures.');
        func_figs2pdf([my_folder '\qualitycheck']);
        
        % Add the synchronization matrix
        if exist([my_folder '\msTouchSync_new.mat'],'file')
            load([my_folder '\msTouchSync_new.mat']);
        else
            load([my_folder '\msTouchSync.mat']);
        end
        dataset.synchronization = synchronization.miniscopeMaster;
        
        % save the final dataset
        disp(['Saving the finalized dataset in ' my_folder '\dataset.mat']);
        save([my_folder '\dataset.mat'],'dataset');
        disp('All done. :)')
        disp(' ');
        
    end
    
    disp(' ')
    disp('---------------------------------------------------------')
    disp(' ')
    
    %pause
    
end

% =============================== functions ===============================

% read the timestamp.dat file, clean it up and convert in into .csv format.
function r = process_timestamp(fol)

datfiles = dir([fol '\timestamp.dat']);

if length(datfiles)~=1
    disp('Not just one dat file. Cannot convert timestamp dat file into csv format ...');
    datfiles
    disp('Moving on but dataset will not be synchronized');
    r = 0;
else
    timestamp_file = datfiles(1).name;
    disp(['Converting dat timestamp file ' timestamp_file ' to csv format']);
    
    % dat file: convert to csv file
    fid = fopen([fol '\' timestamp_file],'rt') ;
    X = fread(fid) ;
    fclose(fid) ;
    X = char(X.') ;
    % replace string S1 with string S2
    Y = strrep(X,char(9), ',') ;
    fid2 = fopen([fol '\' 'timestamp.csv'],'wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    
    % make the folder where the diagnostic figures will be saved.
    mkdir([fol '\qualitycheck']);
    
    disp('ok')
    
    r = 1;
end

end

% read schedules.csv, remove header if necessary, and displays some stats.
function r = process_schedules(fol)

% schedules file is .csv but not timestamp.csv.
files = dir([fol '\*.csv']);

% look for 'schedules.csv' or name = session name minus time of recording

% remove timestamp.csv
for i1=1:length(files)
    if strcmp(files(i1).name,'timestamp.csv')==1
        files(i1) = [];
    end
end

% if there is "schedules.csv", keep it and move on
% look for "schedules.csv"
found_file = 0;
file_name = '';
for i1=1:length(files)
    if strcmpi(files(i1).name,'schedules.csv')==1
        found_file = 1;
        file_name = files(i1).name;
    end
end

if found_file==1
    % we have found "schedules.csv" or some version with lower/uppercase
    schedules_name = file_name;
else
    % we use the session name minus time
    ind = strfind(fol,'\');
    name = fol(ind(end)+1:end);
    under = strfind(name,'_');
    schedules_name = [name(1:(under(5)-1)) '.csv'];
end

% csv ---------------------------------------------------------------------
if strcmp(schedules_name,[])
    disp('Cannot find schedule file. Cannot fix schedule csv file...');
    files
    disp('Moving on but dataset will not be synchronized');
    r = 0;
else
    schedule_file = schedules_name;
    disp(['Fixing schedules file ' schedule_file]);
    
    % csv file: remove all '"'
    fid = fopen([fol '\' schedule_file],'rt');
    X = fread(fid);
    fclose(fid);
    X = char(X.');
    % replace string S1 with string S2
    Y = strrep(X, '"', '');
    
    % remove the header we do not care about (up to Evnt_Time) if there is
    % one.
%     ind = strfind(Y,'Evnt_Time');
%     Y = Y(ind(1):end);
    
    fid2 = fopen([fol '\' 'schedules.csv'],'wt');
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    
    % measure the length of the delay (d1 + d2) just in case...
    ind = find(Y==char(10));
    Nind = length(ind);
    my_delays = [];
    interval = 0;
    for num=2:Nind-1
        myline = Y((ind(num)+1):(ind(num+1)-1));
        comma = strfind(myline,',');
        if (interval==0) && not(isempty(strfind(myline,'Start Delay')))
            interval = str2num(myline(1:(comma(1)-1)));
        end
        if (interval>0) && not(isempty(strfind(myline,'Delay End')))
            newtime = str2num(myline(1:(comma(1)-1)));
            my_delays = [my_delays newtime-interval];
            interval = 0;
        end
    end
    disp(['Mean length of delay d1+d2: ' num2str(mean(my_delays)) 's']);
    
    % also make sure that there are there is a light pulse (at the
    % beginning of the delay)
    nb_light_pulses = 0;
    for num=2:Nind-1
        myline = Y((ind(num)+1):(ind(num+1)-1));
        comma = strfind(myline,',');
        % look for "Pulse Output Event" in column 3, and "TrayLight #1" in
        % column 4
        item3 = myline((comma(2)+1):(comma(3)-1));
        item4 = myline((comma(3)+1):(comma(4)-1));
        
        if strcmp(item3,'Pulse Output Event') && strcmp(item4,'TrayLight #1')
            nb_light_pulses = nb_light_pulses + 1;
        end
    end
    if nb_light_pulses == 0
        disp(['>>>>>>>>>>>>> There were NO light pulses during the recording.'])
    else
        disp(['There were ' num2str(nb_light_pulses) ' light pulses during the recording.'])
    end
    
    disp('ok')
    
    r = 1;
end

end

% reads ms and bh videos, computes mean image luminosity, and computes the 
% synchronization matrix.
function r = process_synchro(fol,st)

if st.timestamp==1 && st.schedules==1
    
    disp('Synchronizing the data ...');
    single_folder_synchronization(fol);
    disp('done');
    
    r = 1;
    
    close all;
    
else
    disp('Missing files:')
    disp(['Timestamp.csv: ' num2str(st.timestamp)]);
    disp(['schedules.csv: ' num2str(st.schedules)]);
    
    r = 0;
end

end

% reads trajectory and bad_data variable (if it exists) and produces single
% generic file.
function r = process_trajectory(fol)

% look for "mouse_traj_corrected.mat". If it is found, nothing to do,
% otherwise produce this file from DLC files.

traj = '';

% look for trajectory in mouse_traj_corrected.mat 
trajfiles = dir([fol '\mouse_traj_corrected.mat']);
if not(isempty(trajfiles))
    
    disp('Found mouse_traj_corrected.mat...')
    
    % we have found mouse_traj_corrected.mat
    r = 1;
    
    % add bad_data to the trajectory file (including the p)
    load([fol '\mouse_traj_corrected.mat']);
    
    % read bad_data if it exists or instead use the p from DLC
    if exist([fol '\bad_data.mat'],'file')
        load([fol '\bad_data.mat']);
    else
        if isfield(position.bodypart1,'bad_data')
            bad_data = position.bodypart1.bad_data;
        end
        if isfield(position.bodypart1,'p')
            bad_data = 1-position.bodypart1.p;
        end
    end
    
    position.bodypart1.bad_data = bad_data;
    save([fol '\mouse_traj_corrected.mat'],'position');
    
    disp('Data saved...')
    disp(' ');
    
    traj = 'ok';
    
end

% if mouse_traj_corrected.mat not found, try DeepLabCut_Trajectories.mat
trajfiles = dir([fol '\DeepLabCut_Trajectories.mat']);
if strcmp(traj,'') && not(isempty(trajfiles))
    
    disp('Found DeepLabCut_Trajectories.mat...')
    
    load([fol '\DeepLabCut_Trajectories.mat']);
    position = [];
    position.bodypart1.x = bodypart1.x_raw;
    position.bodypart1.y = bodypart1.y_raw;
    position.bodypart1.bad_data = 1-bodypart1.pval;
    save([fol '\mouse_traj_corrected.mat'],'position');
    
    disp('Data saved...')
    disp(' ');
    
    traj = 'ok';
    
    r = 1;
    
end

% otherwise try DLC csv files
trajfiles = dir([fol '\behavCam*DeepCut_resnet50_TUNL_TaskJan20shuffle1_1030000.csv']);
if strcmp(traj,'') && length(trajfiles)>0
    
    disp('Found behavCam*DeepCut_resnet50_TUNL_TaskJan20shuffle1_1030000.csv...')
    
    Nfiles = length(trajfiles);
    
    % extract file numbers
    nums = [];
    for i1=1:Nfiles
        temp = trajfiles(i1).name;
        temp = temp(9:end);
        temp = temp(1:(end-51));
        nums = [nums str2num(temp)];
    end
    
    % sort the file numbers
    nums = unique(nums);
    Nfiles = max(nums);
    
    traj = [];
    prob = [];
    
    for f=1:Nfiles
        
        my_file = ['behavCam' num2str(nums(f)) 'DeepCut_resnet50_TUNL_TaskJan20shuffle1_1030000.csv'];
        
        disp(['Processing ' my_file])
        
        % read the trajectory and the detection probability
        fid = fopen([fol '\' my_file]);
        tline = fgetl(fid);
        trial = 0;
        
        % csv file line number
        line_number = 0;
        
        while ischar(tline)
            
            tline = fgetl(fid);
            
            % count the lines
            line_number = line_number + 1;
            
            if not(isequal(tline,-1))
                
                % analyse the content of each line. "," is the separator
                separator_pos = find(tline==',');
                
                if line_number>3
                    
                    x = tline(separator_pos(1)+1:separator_pos(2)-1);
                    y = tline(separator_pos(2)+1:separator_pos(3)-1);
                    p = tline(separator_pos(3)+1:separator_pos(4)-1);
                    traj = [traj; str2num(x) str2num(y)];
                    prob = [prob; str2num(p)];
                    
                end
                
            end
            
        end
        
    end
    
    % save the results :)
    position = [];
    position.bodypart1.x = max(traj(:,1)) - traj(:,1);
    position.bodypart1.y = max(traj(:,2)) - traj(:,2);
    % including the p
    bad_data = (1-prob)>0;
    position.bodypart1.bad_data = bad_data;
    save([fol '\mouse_traj_corrected.mat'],'position');
    
    disp('Data saved...')
    disp(' ');
    
    traj = 'ok';
    
    r = 1;
    
end
   
% we have not found anything
if isempty(traj)
    disp('Cannot find any trajectory files.')
    disp('Moving on but dataset will not be synchronized');
    r = 0;
end

end

% read the deconvolved calcium signal, and saves it in generic file.
function r = process_deconvolved(fol)

disp('Adjusting conventions in deconvolved calcium signal file.');

% read the ms file and store its content in msDeconvolved.mat
ms = [];
calcium = [];
environment = [];
pos = [];
processed = [];
properties = [];
load([fol '\' 'ms.mat']);

if not(isempty(ms))
    temp = ms;
    clear('ms');
    
    if isfield(temp,'deconvolution')==1
        ms = [];
        ms.deconvolution = [];
        ms.deconvolution.deconvolvedSig = temp.deconvolution.deconvolvedSig;
    else
        ms = [];
        ms.deconvolution = [];
        ms.deconvolution.deconvolvedSig = temp.deconvolvedSig;
    end
end

if not(isempty(calcium))
    ms.deconvolution.deconvolvedSig = calcium.trace;
end

save([fol '\' 'msDeconvolved.mat'],'ms');

r = 1;

disp('Done')

end

% merges the trajectory, calcium signal and cage output into a single file
% in a single synchronized frame.
function r = process_merging(fol,stat)

if stat.timestamp == 1 && stat.schedules == 1 && ...
        stat.synchronization == 1 && stat.trajectory == 1 && ...
        stat.deconvolved == 1
    
    disp('Merging the data streams into a single variable...')
    r = data_convergence(fol);
    disp('done')
else
    disp(['Some files are missing. Not performing the merging of data...']);
end

end






















































