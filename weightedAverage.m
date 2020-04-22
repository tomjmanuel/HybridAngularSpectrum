function wAverage = weightedAverage(matrix1, matrix2)

wAverage = sum(sum(matrix1.*matrix2))./sum(sum(matrix1));
