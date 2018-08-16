function syringe_acqtime_QC2(finv,ordered_acq,oind,ss,gca)
%plot expected acquisition times following a normal distribution 
n=length(ordered_acq);
plot(gca,finv,ordered_acq,'o'), hold on

if oind(end)==1 || oind(end)==n
    plot(finv(end),ordered_acq(end),'ro','markerface','r')
end

line([ordered_acq(1) ordered_acq(end)],[ordered_acq(1) ordered_acq(end)]) %expectation
title(['sum of squares:' num2str(ss)])

end %syringe acqtime plotting function


