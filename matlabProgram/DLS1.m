%First DLS stage including demodulation at 2w_R and other frequencies over multiple data sets
%filename of delimited file containing data set names 
%folder to store results in
%number of rotations to demodulate over
function out = DLS1(filename,delimiter,folder,rotations)
runtime1 = clock; %start time
fprintf('starttime %4.1f \n',clock);
[pathstr, name, ext, versn] = fileparts(filename);
%read in data file names
[datafiles dummy] = textread(filename,'%q %n','delimiter',',');
%no. of data files
[dataNo, dummy] = size(datafiles);
fprintf('dataNo %9.1f \n',dataNo);
%define output variables to store results in
out.params = zeros(1,5);
out.sigs = zeros(1,5);
out.centreTimes = zeros(1,1);
out.setSizes = zeros(1,1);
%loop through each data set

for k=1:dataNo
    %read in data [index,timetag,angle,freq]
    tempData = dlmread(char(datafiles(k)),delimiter);
    [temprows,dummy] = size(tempData);
    i=1;
    %define temp variables for results in loop
    params = zeros(1,5);
    sigs = zeros(1,5);
    centreTimes = zeros(1,1);
    setSizes = zeros(1,1);
    %for each data set demodulate in blocks
    %numberPtsRot = 12;
    while i<temprows;%-1.1*12*rotations)
        %collect a subset of the data set to demodulate
        %test for gaps in data, no more than 1 missed point
        subset = tempData(i,:);
        j=i;
        rotationSpan=0;
        while (rotationSpan<rotations) & (i<temprows)
               i=i+1;
               if (tempData(i,3) - tempData(i-1,3)) < -4
                   rotationSpan = rotationSpan + 1;
               end
               subset = cat(1,subset,tempData(i,:));
               %angspan=tempData(i,3)-tempData(j,3);
        end
        [ssrows,sscols]=size(subset);
        %calculate the mean time of the subset
        centreTime=(subset(ssrows,2)+subset(1,2))/2;
        %make zero matrix in which to store A, the model matrix
        A = zeros(ssrows,5);
        A(:,1) = 1;
        A(:,2) = sin(subset(:,3));
        A(:,3) = cos(subset(:,3));
        A(:,4) = sin(2*subset(:,3));
        A(:,5) = cos(2*subset(:,3));
        %perform OLS
        varcov = inv((A')*A);
        %X=(A'*A)^-1.A'.Y
        param = varcov*((A')*subset(:,4));
        %residuals
        res = subset(:,4) - A*param;
        dres = std(res);
        sigraw = sqrt(diag(varcov));
        %param estimate sigmas
        sig = dres*sigraw;
        %store the results
        params = cat(1,params,param');
        sigs = cat(1,sigs,sig');
        centreTimes = cat(1,centreTimes,centreTime);
        setSizes = cat(1,setSizes,ssrows);
        i=i+1;
    end
    %trim of the first empty entry
    params(1,:)=[];
    sigs(1,:)=[];
    centreTimes(1,:)=[];
    setSizes(1,:)=[];
    [rows,dummy]=size(params);
    %just output from here
    params(rows,:)=[];
    sigs(rows,:)=[];
    centreTimes(rows,:)=[];
    setSizes(rows,:)=[];
    
%     out.params = cat(1,out.params,params);
%     out.sigs = cat(1,out.sigs,sigs);
%     out.centreTimes = cat(1,out.centreTimes,centreTimes);
%     out.setSizes = cat(1,out.setSizes,setSizes);
    
    [pathstr2, name2, ext2, versn2] = fileparts(char(datafiles(k)));
    workfilename = fullfile(pathstr,folder,[name2,'.tmp']);
    dlmwrite(strrep(workfilename,'tmp','param'),params,delimiter);
    dlmwrite(strrep(workfilename,'tmp','sig'),sigs,delimiter);
    dlmwrite(strrep(workfilename,'tmp','centTim'),centreTimes,delimiter);
    dlmwrite(strrep(workfilename,'tmp','setSiz'),setSizes,delimiter);
    fprintf('file %2.1f \n',k);
end

%%%%%%%%%%%%%%%%%%%%%%
FinalFileString = 'copy/b ';

[pathstr3, name3, ext3, versn3] = fileparts(char(datafiles(1)));
workfilename = fullfile(pathstr,folder,[name3,'.tmp']);
FinalFileString = [FinalFileString,strrep(workfilename,'tmp','param')];

if dataNo > 1
    
    for z = 2:dataNo
        [pathstr3, name3, ext3, versn3] = fileparts(char(datafiles(z)));
        workfilename = fullfile(pathstr,folder,[name3,'.tmp']);
        FinalFileString = [FinalFileString,'+',strrep(workfilename,'tmp','param')];
    end
end

FinalFileString = [FinalFileString, ' ', strrep(workfilename,'tmp','paramAll')]
system(FinalFileString);
%%%%%%%%%%%%%%%%%%%%%%%%%%
FinalFileString2 = 'copy/b ';

[pathstr3, name3, ext3, versn3] = fileparts(char(datafiles(1)));
workfilename = fullfile(pathstr,folder,[name3,'.tmp']);
FinalFileString2 = [FinalFileString2,strrep(workfilename,'tmp','centtim')];

if dataNo > 1
    
    for z = 2:dataNo
        [pathstr3, name3, ext3, versn3] = fileparts(char(datafiles(z)));
        workfilename = fullfile(pathstr,folder,[name3,'.tmp']);
        FinalFileString2 = [FinalFileString2,'+',strrep(workfilename,'tmp','centtim')];
    end
end

FinalFileString2 = [FinalFileString2, ' ', strrep(workfilename,'tmp','centtimAll')]
system(FinalFileString2);

runtime = clock - runtime1;%total runtime
out.runtime = runtime;
