rng(7)
load('data.mat')

aux_t = P.t_initial:P.delta_t:P.t_end;
ini = 10*60/P.delta_t;
t_end = length(aux_t);


for l=ini:(10/P.delta_t):t_end
    f1 = figure;
    set(gcf, 'Position',  [200, 200, 200, 200])
    axes('Position', [0 0 1 1]);
%     background box
    aux_box = [-1.6 -1.6;1.6 -1.6;1.6 1.6;-1.6 1.6];
    hold on
    fill(aux_box(:,1),aux_box(:,2),[0.7 0.7 0.7]);
    lll = aux_S{l};
    fill(lll(:,1),lll(:,2),[0.6 0.6 0.6]);
    plot([lll(:,1); lll(1,1)],[lll(:,2); lll(1,2)],'Color',[0.6 0.6 0.6])    
%     random points inside the spine
    px = -1.6 + (2*1.6).*rand(2250,1);
    py = -1.6 + (2*1.6).*rand(2250,1);
    in = inpolygon(px,py,lll(:,1),lll(:,2));

    aux_color = [0.5 0.5 0.5];
    plot(px(in),py(in),'o','Markersize',5,'MarkerEdgeColor',aux_color,'MarkerFaceColor',aux_color)
    axis([-1.6 1.6 -1.6 1.6])
    axis off
%     apply gaussian filter and save image 
    saveas(f1,sprintf('Spine%d.tif',(l-ini)*P.delta_t/10+1));
    I = imread(sprintf('Spine%d.tif',(l-ini)*P.delta_t/10+1));
    close(f1)
    Iblur = imgaussfilt(I,1.99);
    f1 = figure;
    set(gcf, 'Position',  [200, 200, 200, 200])
    axes('Position', [0 0 1 1]);
    imshow(Iblur);
    saveas(f1,sprintf('Spine%d.tif',(l-ini)*P.delta_t/10+1));
    close(f1)
end 


