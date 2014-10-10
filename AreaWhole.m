function [Area] = AreaWhole(X1,Curve1,X2,Curve2,Length_Xaxis)

%     Length_Xaxis = 1000;
    
    Xaxis = 0:1/(Length_Xaxis-1):1;
    Xaxis = Xaxis';
    
    Interpolation1 = interp1(X1,Curve1,Xaxis);
    Interpolation2 = interp1(X2,Curve2,Xaxis);

    Area = 0;
    for a = 1:length(Xaxis)-1
        Area = Area + (0.5 * 1/(Length_Xaxis-1) * (abs(Interpolation1(a)-Interpolation2(a)) + abs(Interpolation1(a+1)-Interpolation2(a+1))));
    end
    
end