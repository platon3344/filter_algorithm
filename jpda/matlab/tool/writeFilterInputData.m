% inputData��һ��cell��Ԫ��ÿ��cell��Ԫ�а����˾���
%���ݸ�ʽ��n��x,y,x,y.....
%����n����Ŀ�������x��y���һ��Ŀ������
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