function magnetogramInfo = save_Info(fname, Lp)
    magnetogramInfo.imname = fname;
    magnetogramInfo.label = {Lp.label};
    magnetogramInfo.Area = cellfun(@double, {Lp.area});
    magnetogramInfo.xc = arrayfun(@(x) double(x.center_x), Lp);
    magnetogramInfo.yc = arrayfun(@(x) double(x.center_y), Lp);
    %magnetogramInfo.edge = edge;
end
