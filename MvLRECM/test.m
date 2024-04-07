load('3D.mat')
[y] = MvLRECM(data, y0, options); 
res=EvidentialEvaluationMetrics_old(y0,y);


