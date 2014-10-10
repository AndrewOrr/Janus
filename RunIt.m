close all
clear all
clc


X_small = RGER(10,0.2);
Y_small = RGER(10,0.2);

save('X_small','X_small')
save('Y_small','Y_small')



tic;
[Network_Distance, Bijection, B] = NetworkCompare(X_small,Y_small)
Time = toc;

save('Network_Distance_small','Network_Distance')
save('Bijection_small','Bijection')
save('B_small','B')
save('Time_small','Time')




load X
load Y

tic;
[Network_Distance, Bijection, B] = NetworkCompare(X(1:50,1:50),Y(1:50,1:50));
Time = toc;

save('Network_Distance','Network_Distance')
save('Bijection','Bijection')
save('B','B')
save('Time','Time')




if matlabpool('size') == 0
%     matlabpool open
    parpool(2);
else
    matlabpool close
    parpool(2);
end

tic;
[Network_Distance, Bijection, B] = NetworkCompareParallel(X(1:50,1:50),Y(1:50,1:50));
Time = toc;

save('Network_Distance_parallel','Network_Distance')
save('Bijection_parallel','Bijection')
save('B_parallel','B')
save('Time_parallel','Time')



matlabpool close