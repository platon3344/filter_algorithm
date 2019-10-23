function outputData = readFilterInputData(filename)
data = textread(filename,"%s","delimiter","\n");
for i = 1 : length(data)
    lineData = split(data{i},',');
    measureData = [];
    for j = 1: length(lineData)
        measureData = [measureData str2num(lineData{j})];
    end
    targetNum = measureData(1);
    useData = measureData(2:end);
    outputData{i} = reshape(useData,2,[]);
end