function stop = record_max_fval_params(optimValues, state)
    persistent max_fval_store worst_params_store;
    stop = false;

    if nargin == 0
        % max_fval = max_fval_store;
        stop = worst_params_store;
        return;
    end
    
    if strcmp(state, 'init')
        max_fval_store = -Inf;
        worst_params_store = [];
    elseif strcmp(state, 'iter')
        if optimValues.localsolution.Fval > max_fval_store
            max_fval_store = optimValues.localsolution.Fval;
            worst_params_store = optimValues.localsolution.X;
        end
    end
end