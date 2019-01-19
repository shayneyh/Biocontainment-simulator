close all

%plot extinction time and probability vs alpha/beta
extinction_time_WT =  zeros(length(final),10);
extinction_time_Mut = zeros(length(final),10);

average_extinction_time_WT = zeros(1, length(final));
average_extinction_time_Mut = zeros(1, length(final));

for i = 1:length(final)
    for j = 1:10
        if final{1,i}.final_WT_pop(j) > 0 && final{1,i}.final_Mut_pop(j) > 0
            extinction_time_WT(i,j) = 0;
            extinction_time_Mut(i,j) = 0;
        elseif final{1,i}.final_Mut_pop(j) > 0
            extinction_time_WT(i,j) = final{1,i}.extinction_time(j);
            extinction_time_Mut(i,j) = 0;
        else
            extinction_time_WT(i,j) = 0;
            extinction_time_Mut(i,j) = final{1,i}.extinction_time(j);
        end
    end
%         average_extinction_prob_WT(i) = final{1,i}.extinction_time(1);
%         average_extinction_prob_Mut(i) = final{1,i}.extinction_time(2);
        average_extinction_time_WT(i) = final{1,i}.extinction_time(1);
        average_extinction_time_Mut(i) = final{1,i}.extinction_time(2);
end

figure
hold on
for i = 1:10
    scatter(alpha_range, extinction_time_WT(:,i), 40, [77/256 77/256 216/256], 'Marker', '*')
    scatter(alpha_range, extinction_time_Mut(:,i), 40, [216/256 59/256 59/256])
end

plot([2/3 1.5], [log(2)/2.5076557*20 log(2)/2.5076557*20])
% plot(alpha_range, average_extinction_time_WT, 'Color', 'cyan')
% plot(alpha_range, average_extinction_time_Mut, 'Color', 'black')


hold off
ylim([-10 150])
xlim([0.66 1.5])
title('Extinction time vs Interspecific competing factor')
xlabel('alpha, or 1/beta')
ylabel('Time to extinction (h)')
legend('Wild Type', 'Mutant')

% figure
% hold on
% plot(alpha_range, average_extinction_prob_WT)
% plot(alpha_range, average_extinction_prob_Mut)
% hold off
% title('Extinction Probability vs Interspecific competing factor');







%plot extinction time and probability vs efficiency
extinction_time_WT =  zeros(size(final_dilution_eff_summary));
extinction_time_Mut = zeros(size(final_dilution_eff_summary));

average_extinction_prob_WT = zeros(size(final_dilution_eff_summary));
average_extinction_prob_Mut = zeros(size(final_dilution_eff_summary));

for i = 1:size(final_dilution_eff_summary, 1)
    for j = 1:size(final_dilution_eff_summary, 2)
        extinction_time_WT(i,j) = final_dilution_eff_summary{i,j}.avg_extinction_time(1);
        extinction_time_Mut(i,j) = final_dilution_eff_summary{i,j}.avg_extinction_time(2);
        average_extinction_prob_WT(i,j) = final_dilution_eff_summary{i,j}.extinction_prob(1);
        average_extinction_prob_Mut(i,j) = final_dilution_eff_summary{i,j}.extinction_prob(2);
    end
end

color = {'black', 'red', [0.7500   0    0.7500], 'blue', 'green'};
marker_string = 'o+*s^';
% figure
% hold on
% for i = 1:size(final_dilution_eff_summary, 1)
% %     lines(10) use symbols for each species and colours for dilution
%     plot(efficiency_range, average_extinction_time_WT(i,:), 'Color', color(i,:), 'Marker', marker_string(i))
% end
% hold off
% ylim([-10 150])
% title('Wild type extinction time vs CRISPR efficiency')
% legend('D = 0.5', 'D = 0.135', 'D = 0.1', 'D = 0.01', 'D = 0.001')%WT never dies

figure
hold on
for i = 1:size(final_dilution_eff_summary, 1)
%     lines(10) use symbols for each species and colours for dilution
%     plot(efficiency_range, extinction_time_Mut(i,:), 'Color', [77/256 77/256 40*i/256], 'Marker', marker_string(i))
    plot(efficiency_range, extinction_time_Mut(i,:), 'Color', color{i}, 'Marker', marker_string(i))

end
hold off
ylim([-10 150])
title('Mutant extinction time vs CRISPR efficiency')
xlabel('CRISPR efficiency')
ylabel('Time to extinction (h)')
legend('D = 0.5', 'D = 0.135', 'D = 0.1', 'D = 0.01', 'D = 0.001')


% figure
% hold on
% for i = 1:size(final_dilution_eff_summary, 1)
%     plot(efficiency_range, average_extinction_prob_WT(i,:), 'Color', color(i,:), 'Marker', marker_string(i))
% end
% hold off
% title('Wild type extinction probability vs CRISPR efficiency')
% legend('D = 0.5', 'D = 0.135', 'D = 0.1', 'D = 0.01', 'D = 0.001')
% 
% figure
% hold on
% for i = 1:size(final_dilution_eff_summary, 1)
%     plot(efficiency_range, average_extinction_prob_Mut(i,:), 'Color', color(i,:), 'Marker', marker_string(i))
% end
% hold off
% title('Mutant extinction probability vs CRISPR efficiency')
% legend('D = 0.5', 'D = 0.135', 'D = 0.1', 'D = 0.01', 'D = 0.001')
% 
