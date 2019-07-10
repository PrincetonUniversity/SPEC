function Compare_with_master


fidL1 = fopen('input_L1');

i = 1;
label = cell(1,1);
while ~feof(fidL1)
   tline = fgetl(fidL1);
   filename = char(tline);
    
   label{i} = filename;
   
   out_L1 = [filename, '.h5'];
   out_L3 = [filename(1:6), '3', filename(8:end), '.h5'];
   
   
   figure
   
   idataL1 = read_spec_iota(out_L1);
   fdataL1 = read_spec_field(out_L1);
   pdataL1 = read_spec_poincare(out_L1);
   plot_spec_iota(idataL1, pdataL1, fdataL1, 'i', 'f', 0)
   
   idataL3 = read_spec_iota(out_L3);
   fdataL3 = read_spec_field(out_L3);
   pdataL3 = read_spec_poincare(out_L3);
   plot_spec_iota(idataL3, pdataL3, fdataL3, 'i', 'f', 0)
   
   lines = findobj(gca, 'Type', 'Line');
   lines(1).Marker = 'o';
   lines(1).MarkerSize = 8;
   title(filename(1:end-3))
   
   set(gca,'FontSize',14)                                                                                                                                                                    
   grid on                                                                                                                                                                                   
   legend('Lconstraint=1', 'Lconstraint=3') 
   
   
%    dir = pwd;
%    cd Figures;
%    saveas(gcf, [filename(1:end-3), '.fig'], 'fig'); 
%    saveas(gcf, [filename(1:end-3), '.eps'], 'epsc');
%    cd(dir)
   
   
   [Dabs(i), Drel(i), fest(i)] = compare_spec_outputs(...
            out_L1, out_L3);
   i = i+1;
   
   disp('');
    
end


figure
b = bar(fest);
b(1).FaceColor = [0 0.4470 0.7410];
set(gca, 'XTickLabel', label)
set(gca, 'FontSize', 14)
ylabel('\Delta\mid f\mid')



end