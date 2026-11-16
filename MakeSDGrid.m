function X = MakeSDGrid(d,s_d)
    if d ~= 2
        error('Cette version est prévue pour d=2.');
    end
    x = linspace(0,1,s_d);
    [X1,X2] = meshgrid(x,x);
    X = [X1(:), X2(:)];
end
