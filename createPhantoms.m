function created_Image = createPhantoms(Model, TE, T2, S0, Nrow, Ncol, c_gt)
%%%% The function "Create_Phantoms" creates 256 * 256 phantom images without noise; Mean = 0, Sigma = 0

%%%% Inputs:-
%%%% Model:- "Exponential + Constant", "Exponential" or "Biexponential"
%%%% T2:- T2* value "One element"

%%%% Outputs:-
%%%% created_Image:- 3D image, the third dimension is for the intensity with Time to Echo "TE".

% ============================================================= %


if strcmp(Model, 'exp+c')
    SI = S0 * exp(-TE / T2) + c_gt;                                         % Exponential + Const
    mkdir('exp_plus_c');
    Folder_Name = 'exp_plus_c\';
elseif  strcmp(Model, 'exp')
    SI = S0 * exp(-TE / T2);                 % Exponential
    mkdir('exp');
    Folder_Name = 'exp\';
else strcmp(Model, 'exp+exp')
    SI = 0.9 * S0 * exp(-TE / T2) + 0.1 * S0 * exp(-TE / 200);    % Exponential + Exponential
    mkdir('exp_exp');
    Folder_Name = 'exp_exp\';
end

for i = 1:length(TE)
    created_Image(:,:,i) = SI(i) .* ones(Nrow, Ncol);
end

end