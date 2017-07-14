
% -------------------------------------------------------------------------
% Honey X-Format decoder
% -------------------------------------------------------------------------

% Instructions
% 1. Choose a file to decode
% 2. Choose an amount of samples
% 3. Run the script
% 4. Decoded file is named /data/ in workspace

% ------------------------------------------------------------------------
clear all

len = 30000;%31000
trajStr = input('trajectory characteristic string ... ','s');

% Read and save dtgen
fileName = 'dtgen';
fid = fopen(fileName,'r');      % choose a file
dtgen = zeros(len,7);

lun=fread(fid,1,'int32'); 
nr=fread(fid,1,'int32');
nc=fread(fid,1,'int32'); 
nm=fread(fid,1,'int32');
nt=fread(fid,1,'int32'); 

for i=1:len                % choose samples
      temp(i,1)=fread(fid,1,'int32'); 
      temp(i,2)=fread(fid,1,'int32'); 
      for j=1:nc+1
            dtgen(i,j)=fread(fid,1,'double');
      end;
      
end;
fclose(fid);

save(['_dt' trajStr],'dtgen');

% Read and save trmodel
fileName = 'trmodel';
fid = fopen(fileName,'r');      % choose a file

trmodel = zeros(len,18);

lun=fread(fid,1,'int32'); 
nr=fread(fid,1,'int32');
nc=fread(fid,1,'int32'); 
nm=fread(fid,1,'int32');
nt=fread(fid,1,'int32'); 

for i=1:len                % choose samples
      temp(i,1)=fread(fid,1,'int32'); 
      temp(i,2)=fread(fid,1,'int32'); 
      for j=1:nc+1
            trmodel(i,j)=fread(fid,1,'double');
      end;
      
end;
fclose(fid);

save(['_tr' trajStr],'trmodel');

