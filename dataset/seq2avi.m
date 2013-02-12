addpath(genpath('C:\Users\Anirudh\Documents\CS676\project\toolbox'));
savepath;                 %% include the folder in matlab's path

files = dir('*.seq');     %% get all .seq files in the folder

for i=1:length(files)
    [path fName ext] = fileparts(files(i).('name'));
    images = seqIo([fName ext], 'toImgs');        %% convert to images
    writerObj = VideoWriter(fName);
    open(writerObj);
    writeVideo(writerObj,images);                 %% write to avi file
    close(writerObj);
end