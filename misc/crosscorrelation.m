function [h, lags_collect] = crosscorrelation(x,y,BinWidth)
    lags_collect = [];
    for xter = 1:length(x)
        lags = x(xter) - y;
        lags_collect = [lags_collect lags];
    end
    h = histcounts(lags_collect, "BinWidth", BinWidth);
end