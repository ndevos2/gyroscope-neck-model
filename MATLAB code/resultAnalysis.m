%%
%{



%}


%% Average

errA = [results.exp11A_RMSEn, results.exp12A_RMSEn, results.exp13A_RMSEn, results.exp21A_RMSEn, results.exp22A_RMSEn, results.exp51A_RMSEn, results.exp52A_RMSEn, results.exp53A_RMSEn, results.exp61A_RMSEn, results.exp62A_RMSEn, results.exp63A_RMSEn];
errAmean = mean(errA);
errAmedian = median(errA);

maxA = [results.exp11A_EMax, results.exp12A_EMax, results.exp13A_EMax, results.exp21A_EMax, results.exp22A_EMax, results.exp51A_EMax, results.exp52A_EMax, results.exp53A_EMax, results.exp61A_EMax, results.exp62A_EMax, results.exp63A_EMax];
maxAmean = mean(maxA);
maxAmedian = median(maxA);

errB = [results.exp11B_RMSEn, results.exp12B_RMSEn, results.exp13B_RMSEn, results.exp21B_RMSEn, results.exp22B_RMSEn, results.exp51B_RMSEn, results.exp52B_RMSEn, results.exp53B_RMSEn, results.exp61B_RMSEn, results.exp62B_RMSEn, results.exp63B_RMSEn];
errBmean = mean(errB);

maxB = [results.exp11B_EMax, results.exp12B_EMax, results.exp13B_EMax, results.exp21B_EMax, results.exp22B_EMax, results.exp51B_EMax, results.exp52B_EMax, results.exp53B_EMax, results.exp61B_EMax, results.exp62B_EMax, results.exp63B_EMax];
maxBmean = mean(maxB);

criteraA = [results.exp11A_RMSEn, results.exp13A_RMSEn, results.exp21A_RMSEn, results.exp22A_RMSEn, results.exp51A_RMSEn, results.exp52A_RMSEn, results.exp53A_RMSEn, results.exp61A_RMSEn, results.exp62A_RMSEn, results.exp63A_RMSEn];
criteriaAmean = mean(criteraA);

criteriaMaxA = [results.exp11A_EMax, results.exp12A_EMax, results.exp13A_EMax, results.exp22A_EMax, results.exp51A_EMax, results.exp52A_EMax, results.exp53A_EMax, results.exp61A_EMax, results.exp62A_EMax, results.exp63A_EMax];
criteriaMaxAmean = mean(criteriaMaxA);

criteriaB = [results.exp52B_RMSEn, results.exp53B_RMSEn];
criteriaBmean = mean(criteriaB);

criteriaMaxB = [results.exp52B_EMax, results.exp53B_EMax];
criteriaMaxBmean = mean(criteriaMaxB)