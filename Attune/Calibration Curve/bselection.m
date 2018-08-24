function bselection(bg,event,F)
display(['Previous: ', event.OldValue.Text]);
display(['Current: ', event.NewValue.Text]);
display('------------------')
switch event.NewValue.Text
    case 'Bad'
        exp(F.index) = 3;
    case 'Concerning'
        exp(F.index) = 2;
    case 'Good'
        exp(F.index) = 0 ;
end