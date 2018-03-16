% Usage:
% A function that searches recursively through all subdirectories 
% of a given directory, collecting a list of all file names it finds.
% 
% OUTPUT ARGUMENTS
%       flnmList                The list of full names incl. path
%       nmList                  The list of names without path
% 
% Liyan modifies codes that is originally downloaded from:
% o http://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flnmList = getAllFiles(dirName)

  dirData = dir(dirName);                       % Get the data for the current directory
  dirIndex = [dirData.isdir];                   % Find the index for directories
  flnmList = {dirData(~dirIndex).name}';        % Get a list of the files
  
  if ~isempty(flnmList)
    flnmList = cellfun(@(x) fullfile(dirName,x),...      % Prepend path to files
                       flnmList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};           % Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});   % Find index of subdirectories
                                                % that are not '.' or '..'
  for iDir = find(validIndex)                   % Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});  % Get the subdirectory path
    flnmList = [flnmList; getAllFiles(nextDir)];% Recursively call getAllFiles
  end

end
%% test
%{
clear,clc

dirName = 'D:\WORKSPACE\data_tmp\Caltech256\256_ObjectCategories\';
% 
flnmList = getAllFiles(dirName);
%}
