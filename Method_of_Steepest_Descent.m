%{
    Uses the algorithm as stated in the textbook in section 8.2.4 Method
    of Steepest Descent:
    
    https://www.cs.utexas.edu/users/flame/laff/alaff/
        chapter08-steepest-descent.html

    This also implements the stopping criteria defined in 8.3.6.1 Stopping
    Criteria:

    https://www.cs.utexas.edu/users/flame/laff/alaff/
        chapter08-CGS-final-details.html
%}
function [ x, niters ] = Method_of_Steepest_Descent( A, b, x0 )
    % Set iteration limit
    itr_limit = 1001;
    
    %{
    fprintf( "The iteration limit (itr_limit) is %d \n", itr_limit );
    fprintf( "Execute Method_of_Steepest_Descent \n" );
    %}

    x = x0;
    k = 0;
    r = b - A * x;
    
    %norm( r, 2 ) <= eps * norm( b, 2 )  % debug statement
    %k < itr_limit  % debug statement
    
    %while r ~= 0
    % Implement 8.3.6.1 Stopping Criteria
    while norm( r, 2 ) >= eps * norm( b, 2 ) && k < itr_limit
        %fprintf("Iteration %d", k);
        p = r;
        q = A * p;
        alpha = ( p.' * r ) / ( p.' * q );
        x = x + alpha * p;
        r = r - alpha * q;
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
