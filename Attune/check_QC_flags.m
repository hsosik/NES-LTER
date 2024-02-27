basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\';
cruises = dir([basepath '20*']);
cruises = {cruises.name};
cruises(ismember(cruises, {'20180810_SR2018' '20190303_AL' '20190725_HB1907' '20190921_AR38' '20210512_SG2105'})) = []; %no underway data
cruises(ismember(cruises, {'20180131_EN608' '20190414_AR34B' '20190705_TN368' '20200201_EN649' '20211108_AR61B' '20211117_AR62' '20211126_AR63' '20220216_AT46'})) = []; %QC_flowrates missing from table
cruises(ismember(cruises, {'20221117_AR70B' '20230111_EN695' '20231011_AR77' '20231023_AR78'})) = []; %not processed yet
for ii = 1:length(cruises)
    if exist([basepath cruises{ii} '\bead_calibrated\AttuneTable.mat'], 'file')
        load([basepath cruises{ii} '\bead_calibrated\AttuneTable.mat']);
    else
        load([basepath cruises{ii} '\ifcb_calibrated\AttuneTable.mat']);
    end
    iii = find(AttuneTable.QC_flag==0);
    figure
    plot(AttuneTable.QC_flowrates(:,1), '.'), hold on
    plot(iii,AttuneTable.QC_flowrates(iii,1), 'bo')
    line(xlim, [1.5 1.5], 'color', 'b')
    yyaxis right, plot(AttuneTable.QC_flowrates(:,2), '.')
    plot(iii,AttuneTable.QC_flowrates(iii,2), 'ro')
    line(xlim, [2 2], 'color', 'r')
    title(cruises{ii}, 'Interpreter','none')
end


%%
load('\\sosiknas1\Lab_data\Attune\cruise_data\20220806_EN688\bead_calibrated\AttuneTable.mat')
ii = find(AttuneTable.Syn_carbon>10000)
%check movie at these frames; strange noise? pile up at high FSC-W
%   {'OTZ_EN688_08Aug2022(4)_Group_day0_Syn(23).fcs' }
%     {'OTZ_EN688_09Aug2022(4)_Group_day0_Syn(19).fcs' }
%     {'OTZ_EN688_10Aug2022(14)_Group_day0_Syn(1).fcs' }
%     {'OTZ_EN688_11Aug2022(3)_Group_day0_Syn(4).fcs'  }
%     {'OTZ_EN688_13Aug2022(2)_Group_day0_Syn(22).fcs' }
%     {'OTZ_EN688_13Aug2022(3)_Group_day0_Syn(5).fcs'  }
%     {'OTZ_EN688_16Aug2022(10)_Group_day0_Syn(21).fcs'}

ii2 = find(AttuneTable.Syn_count<2000)
%     {'OTZ_EN688_07Aug2022(9)_Group_day0_Syn(15).fcs' } %burst of noise
%     {'OTZ_EN688_08Aug2022(4)_Group_day0_Syn(23).fcs' }
%     {'OTZ_EN688_09Aug2022(7)_Group_day0_Syn(11).fcs' }
%     {'OTZ_EN688_09Aug2022(13)_Group_day0_Syn(7).fcs' }
%     {'OTZ_EN688_10Aug2022(13)_Group_day0_Syn(17).fcs'}
%     {'OTZ_EN688_10Aug2022(15)_Group_day0_Syn(4).fcs' }
%     {'OTZ_EN688_11Aug2022(1)_Group_day0_Syn(21).fcs' }
%     {'OTZ_EN688_11Aug2022(2)_Group_day0_Syn(5).fcs'  }
%     {'OTZ_EN688_11Aug2022(3)_Group_day0_Syn(4).fcs'  }
%     {'OTZ_EN688_11Aug2022(5)_Group_day0_Syn(3).fcs'  }
%     {'OTZ_EN688_11Aug2022(5)_Group_day0_Syn(8).fcs'  }
%     {'OTZ_EN688_11Aug2022(6)_Group_day0_Syn(6).fcs'  }
%     {'OTZ_EN688_11Aug2022(6)_Group_day0_Syn(8).fcs'  }
%     {'OTZ_EN688_11Aug2022(7)_Group_day0_Syn.fcs'     }
%     {'OTZ_EN688_11Aug2022(10)_Group_day0_Syn(7).fcs' }
%     {'OTZ_EN688_13Aug2022(7)_Group_day0_Syn(8).fcs'  }
%     {'OTZ_EN688_13Aug2022(10)_Group_day0_Syn(21).fcs'}
%     {'OTZ_EN688_15Aug2022(5)_Group_day0_Syn.fcs'     }
%     {'OTZ_EN688_16Aug2022(5)_Group_day0_Syn(12).fcs' }
%     {'OTZ_EN688_16Aug2022(10)_Group_day0_Syn(21).fcs'}
%     {'OTZ_EN688_16Aug2022(25)_Group_day0_Syn(9).fcs' }

% AttuneTable(setdiff(ii2,ii),{'Filename' 'QC_flag'})
%                          Filename                         QC_flag
%     __________________________________________________    _______
% 
%     {'OTZ_EN688_07Aug2022(9)_Group_day0_Syn(15).fcs' }       1   
%     {'OTZ_EN688_09Aug2022(7)_Group_day0_Syn(11).fcs' }       1   
%     {'OTZ_EN688_09Aug2022(13)_Group_day0_Syn(7).fcs' }       0   
%     {'OTZ_EN688_10Aug2022(13)_Group_day0_Syn(17).fcs'}       1   
%     {'OTZ_EN688_10Aug2022(15)_Group_day0_Syn(4).fcs' }       1   
%     {'OTZ_EN688_11Aug2022(1)_Group_day0_Syn(21).fcs' }       0   
%     {'OTZ_EN688_11Aug2022(2)_Group_day0_Syn(5).fcs'  }       1   
%     {'OTZ_EN688_11Aug2022(5)_Group_day0_Syn(3).fcs'  }       0   
%     {'OTZ_EN688_11Aug2022(5)_Group_day0_Syn(8).fcs'  }       1   
%     {'OTZ_EN688_11Aug2022(6)_Group_day0_Syn(6).fcs'  }       1   
%     {'OTZ_EN688_11Aug2022(6)_Group_day0_Syn(8).fcs'  }       1   
%     {'OTZ_EN688_11Aug2022(7)_Group_day0_Syn.fcs'     }       0   
%     {'OTZ_EN688_11Aug2022(10)_Group_day0_Syn(7).fcs' }       1   
%     {'OTZ_EN688_13Aug2022(7)_Group_day0_Syn(8).fcs'  }       0   
%     {'OTZ_EN688_13Aug2022(10)_Group_day0_Syn(21).fcs'}       0   
%     {'OTZ_EN688_15Aug2022(5)_Group_day0_Syn.fcs'     }       1   
%     {'OTZ_EN688_16Aug2022(5)_Group_day0_Syn(12).fcs' }       0   
%     {'OTZ_EN688_16Aug2022(25)_Group_day0_Syn(9).fcs' }       0  

%Compare these two
C1 = load('\\sosiknas1\Lab_data\Attune\cruise_data\20220806_EN688\bead_calibrated\class\OTZ_EN688_07Aug2022(9)_Group_day0_Syn(13).mat')
C = load('\\sosiknas1\Lab_data\Attune\cruise_data\20220806_EN688\bead_calibrated\class\OTZ_EN688_07Aug2022(9)_Group_day0_Syn(15).mat')
figure
%plot(C1.ssc_value(C1.class==0), '+')
hold on
plot(C1.ssc_value(C1.class==2), '+')
hold on
%plot(C.ssc_value(C.class==0), '.')
plot(C.ssc_value(C.class==2), '.')
%--> somehow the syn are not being measured properly on SSC
% "Num_particles" value seems to be all triggers and it's about the same for these adjacent files
% but Syn count is much lower in the bad one and noise is much higher by same differences
%Look at this plot!!
semilogy(C1.ssc_value, '.')
hold on
semilogy(C.ssc_value, '.')


%{'OTZ_EN688_08Aug2022(12)_Group_day0_Syn(24).fcs'}
%This is an different(?) case with a portion of the Syn apparently
%mis-measured on GL2 and BL2
semilogy(C3.ssc_value, '.') %looks weird at end but movie shows bigger problem appear in GL2, BL2
