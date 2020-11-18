% Data.mat contains aux_area,aux_index2,aux_index3,aux_S,aux_foci,aux_B, 
% ff_length,P
% load('data.mat')
% aux_t = P.t_initial:P.delta_t:P.t_end;

figure;
plot(aux_t/60,aux_area)
xlabel('time')
ylabel('Area')

for l=1:(10/P.delta_t):length(aux_t)
    figure(100);clf
    set(gca,'fontsize',20)
    hold on 
    lll = aux_S{l,1}; 
    plot([lll(:,1); lll(1,1)],[lll(:,2); lll(1,2)],'k','linewidth',2)
    plot(lll(aux_index2(l,1):aux_index2(l,2),1),lll(aux_index2(l,1):aux_index2(l,2),2),'r','linewidth',2)
    plot(lll(aux_index3(l,1):aux_index3(l,2),1),lll(aux_index3(l,1):aux_index3(l,2),2),'b','linewidth',2)
    title(['t = ' num2str(floor(aux_t(l)/60),'%.0f') ' min  ' num2str(rem(aux_t(l),60),'%.0f') ' s'])
    axis([-1.5 1.5 -1.5 1.5])
    xlabel('{\mu}m')
    ylabel('{\mu}m')
    set(gcf,'color','w');   
end

%obtain shape descriptors
Li_descriptors

figure
hold on 
plot(0:10:290,100*std(roi_area')'./mean(roi_area,2),'-.','color','k','linewidth',2) 
plot(0:10:290,100*std(S')'./mean(S,2),'color','k','linewidth',2)
plot(0:10:290,100*std(D')'./mean(D,2),':','color','k','linewidth',2)
plot(0:10:290,100*std(O')'./mean(O,2),'--','color','k','linewidth',2) 
set(gca,'fontsize',14)
ylabel('CV (%)')
xlabel('time (s)')
axis([-.25 290.25 0 50])
legend('Area','S','D','O')

