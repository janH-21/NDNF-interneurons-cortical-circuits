function [mean_dats, std_dats, sem_dats, median_dats] = statistics(dats, dim)
    mean_dats = mean(dats,dim,"omitnan");
    std_dats = std(dats, [], dim,"omitnan");
    sem_dats = std(dats, [], dim,"omitnan")./sqrt(size(dats,dim));
    median_dats = median(dats,dim,"omitnan");
end