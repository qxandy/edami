
% =========================================================================
% -------------------------------------------------------------------------
%
%  SCRIPT TO PLOT TRAJECTORIES OF THE DIFFERENT METHODS FOR VISUALISATION
% -------------------------------------------------------------------------
% =========================================================================
%    Version : 1.0 (7 September 2013)
%    Since   : 2 March 2014
%    Authors : Momodou Lamin Sanyang
%              School of Computer Science, UoB
%    Contact : mls248@cs.bham.ac.uk
% =========================================================================
% You need to number the .mat file like BFC_N and BFG_N.
% N is the population size the programs use to get the results.
% For each Functions, you choose a dimension and vary the population size
% to include [300 1000 5000]
% For example: Function 1:
% Choose a dimension and run [300 1000 5000] population sizes. save the
% results with the population size as above. You should have a Folder with
% the following name as example, F1_D20, where F1 is the function 1 you
% tried and D20 means dimension 20 was used.Inside this folder, you should
% have BFC_300, BFC_1000, BFC_2000, BFG_300, BFG_1000 and BFG_2000
% =========================================================================

function Trajectory(pD,pNp,pFun)
%N = [20];
nX = floor((1*10^4*pD)/pNp);
h=figure;
lw=3;
clear B;
clear Bx;
load (['BFE_' num2str(pD) '_' num2str(pNp) '_F' num2str(pFun)]);
for i=1:25
    Bx(i,:) = B(i).fv(1:nX);
end
semilogy(mean(Bx), '-r','LineWidth',lw)
%meanc = mean(Bx);
% stdc  = std(Bx);

hold on

clear B;
clear Bx;
load (['BFG_' num2str(pD) '_' num2str(pNp) '_F' num2str(pFun)]);
for i=1:25
    Bx(i,:) = B(i).fv(1:nX);
end
semilogy(mean(Bx), '--b','LineWidth',lw)
meang = mean(Bx);
stdg  = std(Bx);

hold on 

clear B;
clear Bx;
load (['BFM_' num2str(pD) '_' num2str(pNp) '_F' num2str(pFun)]);
for i=1:25
    Bx(i,:) = B(i,:).fv(1:nX);
end
semilogy(mean(Bx), '-.g','LineWidth',lw)
means = mean(Bx);
stds  = std(Bx);
set(gca,'FontSize',16)
xlabel('Generations','fontsize',18)
ylabel('Fitness - Optimal fitness','fontsize',18)
%title(['Population size is:' num2str(N)],'fontsize',22)
legend('EDA-MCC','EDA-MCC-MI','EDA-MCC-GMI')
axis('tight')

%  if N ==20
%  print -depsc -tiff -r300 F20
%  elseif N== 2
%      print -depsc -tiff -r300 F2
%   elseif N== 4
%       print -depsc -tiff -r300 F4
%   elseif N== 6
%       print -depsc -tiff -r300 F6
%   elseif N== 9
%       print -depsc -tiff -r300 F9
%  else
%      print -depsc -tiff -r300 F13
% end
hold off
saveas(h,['result_F' num2str(pFun) '_D' num2str(pD) '_P' num2str(pNp)],'png');
clear h;
disp(['File "result_D' num2str(pD) 'P_' num2str(pNp) '_F' num2str(pFun) '.png" saved.']);
clear;
close;
end