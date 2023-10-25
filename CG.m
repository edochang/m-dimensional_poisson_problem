function [ x, niters ] = CG( A, b, x0 )
    % Set iteration limit
    itr_limit = 1001;
    
    %{
    fprintf( "The iteration limit (itr_limit) is %d \n", itr_limit );
    fprintf( "Execute CG \n" );
    %}

    x = x0;
    r = b;
    k = 0;
      
    % Implement 8.3.6.1 Stopping Criteria
    while norm( r, 2 ) >= eps * norm( b, 2 ) && k < itr_limit
        if k == 0
            p = r;
        else
            %gamma = -1 * ( p_prev.' * A * r ) / ( p_prev.' * A * p_prev );
            gamma = ( r.' * r ) / ( r_prev.' * r_prev);
            p = r + ( gamma * p_prev );
        end
        % Store current iteration p values as the next iteration's
        % previous iteration values.
        p_prev = p;

        alpha = ( r.' * r ) / ( p.' * A * p );
        %alpha = ( p.' * r ) / ( p.' * A * p );
        x = x + alpha * p;
        
        % Store current iteration r values as the next iteration's
        % previous iteration values.
        r_prev = r;

        r = r - ( alpha * A * p );
        k = k + 1;
    end

    niters = k;
    %{
    if niters == itr_limit
        fprintf( "The stopping critiera of 8.3.6.1 was not met and " + ...
            "ran for %d iterations.", niters );
    else
        fprintf( "Requires %d iterations before meeting the " + ...
            "stopping criteria of 8.3.6.1.", niters );
    end
    %}
end