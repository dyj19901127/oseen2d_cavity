clear;
close all;

n = 10;
[U,T,t,x] = importFoamData(n,'data');

vidName = 'Temperature.avi';

writerObj = VideoWriter(vidName);
open(writerObj);

X = reshape(x(:,1),n,n);
Y = reshape(x(:,2),n,n);

numTS = size(T,2);
Z = zeros(n);

for k = 1:numTS
  Z = reshape(T(:,k),n,n);
  surf(X,Y,Z);
  frame = getframe;
  writeVideo(writerObj,frame);
end

close(writerObj);

