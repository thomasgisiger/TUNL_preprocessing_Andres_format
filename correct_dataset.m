function r = correct_dataset(ds,lags)

% cage dimensions
cage_dims = [240 332];
cage_prop = cage_dims(1)/cage_dims(2);

% shift to be performed
dt = 1/30;
delay = round(mean(lags(:)/dt));

disp(['Shifting by ' num2str(delay) ' frames.']);

% shift trajectory and dcs
Nf = size(ds.trialNumber,1);
% traj
traj = [ds.headPosition.x ds.headPosition.y ds.headPosition.p];
if delay>0
    temp = [traj(1,1:3).*ones(delay,3); traj(1:(Nf-delay),:)];
elseif delay<0
    temp = [traj((abs(delay)+1):Nf,:); traj(Nf,1:3).*ones(abs(delay),3)];
else
    temp = traj;
end
ds.headPosition.x = temp(:,1);
ds.headPosition.y = temp(:,2);
ds.headPosition.p = temp(:,3);

% dcs
dcs = ds.dcs;
Ncells = size(dcs,2);
if delay>0
    temp = [dcs(1,1:Ncells).*ones(delay,Ncells); dcs(1:(Nf-delay),1:Ncells)];
elseif delay<0
    temp = [dcs((abs(delay)+1):Nf,1:Ncells); dcs(Nf,1:Ncells).*ones(abs(delay),Ncells)];
else
    temp = dcs;
end

% points that will be used for the registration
temp = extract_cage_contact_points_4classes(ds,'p');
ds.registration.pos1 = temp.pos1;
ds.registration.pos5 = temp.pos5;
ds.registration.posR1 = temp.posR1;
ds.registration.posR5 = temp.posR5;

% -----------------------------------------------------------------

% normalize the trajectory with the cage's proportions: y in [0 1],
% and x in [-240.332 240.332]
dataset_corrected = ds;

% rescale the trajectory
ax = max(ds.headPosition.x);
bx = min(ds.headPosition.x);
dataset_corrected.headPosition.x = ...
    ((ds.headPosition.x-bx)/(ax-bx)-1/2)*2*cage_prop;
ay = max(ds.headPosition.y);
by = min(ds.headPosition.y);
dataset_corrected.headPosition.y = ...
    (ds.headPosition.y-by)/(ay-by);
% and the contact points
dataset_corrected.registration.pos1 = [((ds.registration.pos1(:,1)-bx)/(ax-bx)-1/2)*2*cage_prop (ds.registration.pos1(:,2)-by)/(ay-by)];
dataset_corrected.registration.pos5 = [((ds.registration.pos5(:,1)-bx)/(ax-bx)-1/2)*2*cage_prop (ds.registration.pos5(:,2)-by)/(ay-by)];
dataset_corrected.registration.posR1 = [((ds.registration.posR1(:,1)-bx)/(ax-bx)-1/2)*2*cage_prop (ds.registration.posR1(:,2)-by)/(ay-by)];
dataset_corrected.registration.posR5 = [((ds.registration.posR5(:,1)-bx)/(ax-bx)-1/2)*2*cage_prop (ds.registration.posR5(:,2)-by)/(ay-by)];

% -----------------------------------------------------------------

% check that the pad 1 points are on the left (x<0)
mean_pos1 = mean(dataset_corrected.registration.pos1,1);
if mean_pos1(1)>0
    disp('Mirror transforming (x -> -x) the trajectory and contact points')
    
    % we need to mirror transform along x the trajectory and the
    % contact points
    % traj
    dataset_corrected.headPosition.x = - dataset_corrected.headPosition.x;
    % contact points
    temp5 = [-dataset_corrected.registration.pos5(:,1) dataset_corrected.registration.pos5(:,2)];
    temp1 = [-dataset_corrected.registration.pos1(:,1) dataset_corrected.registration.pos1(:,2)];
    tempR5 = [-dataset_corrected.registration.posR5(:,1) dataset_corrected.registration.posR5(:,2)];
    tempR1 = [-dataset_corrected.registration.posR1(:,1) dataset_corrected.registration.posR1(:,2)];
    % 1
    
    dataset_corrected.registration.pos1 = temp1;
    % 5
    dataset_corrected.registration.pos5 = temp5;
    % R1
    dataset_corrected.registration.posR1 = tempR1;
    % R5
    dataset_corrected.registration.posR5 = tempR5;
end

% -----------------------------------------------------------------

% check that the reward port is at the bottom (y < 0.5)
mean_posR = mean([dataset_corrected.registration.posR1; dataset_corrected.registration.posR5],1);
if mean_posR(2)>0.5
    disp('Mirror transforming (y -> 1-y) the trajectory and contact points')
    
    % we need to mirror transform along x the trajectory and the
    % contact points
    % traj
    dataset_corrected.headPosition.y = 1 - dataset_corrected.headPosition.y;
    % contact points
    temp5 = [dataset_corrected.registration.pos5(:,1) 1-dataset_corrected.registration.pos5(:,2)];
    temp1 = [dataset_corrected.registration.pos1(:,1) 1-dataset_corrected.registration.pos1(:,2)];
    tempR5 = [dataset_corrected.registration.posR5(:,1) 1-dataset_corrected.registration.posR5(:,2)];
    tempR1 = [dataset_corrected.registration.posR1(:,1) 1-dataset_corrected.registration.posR1(:,2)];
    % 1
    dataset_corrected.registration.pos1 = temp1;
    % 5
    dataset_corrected.registration.pos5 = temp5;
    % R1
    dataset_corrected.registration.posR1 = tempR1;
    % R5
    dataset_corrected.registration.posR5 = tempR5;
end

r = dataset_corrected;

end