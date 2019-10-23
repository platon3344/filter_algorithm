% inputData是一个cell单元，每个cell单元中包含了矩阵
%数据格式：n，x,y,x,y.....
%其中n代表目标个数，x和y组成一对目标坐标
function outputData = writeFilterInputData(filename,inputData)
fid = fopen(filename,"w");
for i = 1 : length(inputData)
    writeData = inputData{i};
    fprintf(fid,"%d , " ,length(writeData));
    for jj = 1 : size(writeData,2)
        for kk = 1 : size(writeData,1)
            fprintf(fid,"%f ,", writeData(kk,jj));
        end
    end
    fprintf(fid,"\n");
end
fclose(fid);
end