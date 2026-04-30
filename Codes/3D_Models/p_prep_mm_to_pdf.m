function [pdf_terp, mm_terp, dmm_di] = p_prep_mm_to_pdf(pdfs, nmm)

%%% Go from different mm (posible model values) at each station to a common
%%% mm basis. 
nsta = length(pdfs); 

pdf_terp = zeros(nmm, nsta); 

all_mm = [pdfs(:).mm]; 
min_mm = min(min(all_mm)); 
max_mm = max(max(all_mm)); 
mm_terp = linspace(min_mm, max_mm, nmm)'; % Corresponds to pdf_terp dimension 1. 

dmm_di = mm_terp(2) - mm_terp(1); % How much does mm change per index. 

for ista = 1:nsta; 
    pdf_terp_i = interp1(pdfs(ista).mm, pdfs(ista).pm, mm_terp,...
        'cubic', 0 ); 
    pdf_terp(:, ista) = pdf_terp_i; 
end



end