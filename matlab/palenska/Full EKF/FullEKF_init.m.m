%       Initialization file


clear all

trajDir = 'C:\Users\Markéta\Documents\MATLAB\trajectories\';

  load ([trajDir '_dtDyn6_noLA.mat']);
  load ([trajDir '_trDyn6_noLA.mat']);


dtgen(:,1) = dtgen(:,1) - 0.01;  
trmodel(:,1) = trmodel(:,1) - 0.01;
%--------------------------------------------------------------

dT = 0.01;
startT = 0;
endT = 600;
startindex = min(find(dtgen(:,1)>=startT));
endindex = max(find(dtgen(:,1)<=endT));
dtgen = dtgen(startindex:endindex,:);
trmodel = trmodel(startindex:endindex,:);
initP=diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6 1 1 1 5e-6 5e-6 5e-6 1e-8 1e-8 1e-8]); 
    
% For succesfull running of model modify path 'trajDir' to directory when
% trajectories are located.


 
