
size2 = []; 
size2_10 = [];
size10_20 = [];
size20 = [];
diameter = 10:0.1:20;

for ii = 1:length(diameter)
    disp(diameter(ii))
    if diameter(ii) <= 2
        size2 = [size2;diameter(ii)];
    elseif diameter(ii)>= 2 & diameter(ii) <= 10
        size2_10 = [size2_10;diameter(ii)];
    elseif  diameter(ii)>= 10 & diameter(ii) <= 20
        size10_20 =[size10_20; diameter(ii)];
    elseif diameter(ii) >= 20 
        size20 =[ size20;diameter(ii)];
    end
 end