function  Diff(refc,refc_true)
%DIFF Check error and estimate error
%   Input:
%           refc: numerical result
%           refc_true: exact answer
fprintf('relative $L_2$ error: %s\n\n', num2str(norm(refc-refc_true)/norm(refc_true)));
fprintf('relative $L_{\\infty}$ error: %s\n\n',num2str(norm(refc-refc_true,inf)/norm(refc_true,inf)));

end

