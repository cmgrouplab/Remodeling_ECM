clear;

nImages = 31;
fig = figure;

for idx = 1:300:10001
    T1=readtable(['/Users/yuzheng/Desktop/Fiber_model_Code/fiber_data/rotation_16/angle20/density_0.3_0.01/fiberPosition',num2str(idx-1),'.csv']);
    x=T1.Var1;
    y=T1.Var2;
    density = T1.Var3;
    scatter(x,y,30,density,'filled');
    title(num2str(idx-1));
    colormap(jet)
    colorbar;
    caxis([0 1])
    axis square;
    drawnow
    frame = getframe(fig); 
    im{idx}=frame2im(frame);
end
close;
filename = 'fiber/density.gif'; % Specify the output file name
for idx = 1:300:10001
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
