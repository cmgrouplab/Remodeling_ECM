clear;

nImages = 31;
fig = figure;

for idx = 1:1000:30001
    T1=readtable(['/Users/yuzheng/Desktop/Fiber_model_Code/fiber/fiber_huan2/position',num2str(idx-1),'.csv']);
    x=T1.Var1;
    y=T1.Var2;
    x = x - 0.5;
    y = y - 0.5;
    distance = sqrt(x.^2 + y.^2);  
    [row,col] = size(distance);

    delta_r = 0.01;
    res = zeros(25,2);
    j = 1;
    for r = 0.02:0.02:0.5
        count = 0;
        for i = 1:row
            if r-delta_r <=distance(i,1) && distance(i,1) <=r + delta_r
                count = count + 1;
            end
        end
        res(j,1) = r;
        res(j,2) = count / (2*pi*r*2*delta_r);
        j = j + 1;
    end

    bar(res(:,1),res(:,2),1)
    title([num2str(idx-1)]);
    axis([0.01 0.51 0 5000])
    axis square;
    drawnow
    frame = getframe(fig); 
    im{idx}=frame2im(frame);
end
close;
filename = 'fiber/histogram_particle.gif'; % Specify the output file name
for idx = 1:1000:30001
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
    end
end
