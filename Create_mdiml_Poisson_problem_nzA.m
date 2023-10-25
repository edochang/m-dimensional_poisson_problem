function [ nzA, ir, ic ] = Create_mdiml_Poisson_problem_nzA( N , m )
    % Note: m >= 1 with N^m interior points.
    %
    % Assume with Poisson equation you will have at most ( 2 * m ) + 1
    % elements per row.
    max_row_elements = ( 2 * m ) + 1;
    
    % Number of rows in A is N^m
    a_rows = N^m;

    nzA = zeros( ( a_rows ) * max_row_elements, 1);
    ir = zeros( a_rows + 1, 1);
    ic = zeros( ( a_rows ) * max_row_elements, 1 );
    
    % Calculate the \tilde{a}_i for A ( N^m x N^m )
    for i = 1 : a_rows
        % Set ir(1) = 1
        ir(1) = 1;

        % Calculate the column indices for u_i.
        %a_row_column_indices = [ i-N, i-1, i, i+1, i+N ];
        
        % Generate the index pattern for the finite differential
        % pattern / stencil.
        % Add 1-D index stencil
        a_row_column_indices = zeros( 1, max_row_elements );
        a_row_middle = ceil( max_row_elements / 2 );

        %i % debug

        % Set i
        a_row_column_indices( a_row_middle ) = i;

        %a_row_column_indices  % debug

        % Positive Direction N^( m - 1 ) * 2 + 1 and N^( m - 1 ) * 5 + 1
        % and loop N^( m - 1) times where i + N is set to 0.
        for j = 1 : a_row_middle - 1
            power = m-j;
            np = N^power;            
            
            if j > 1
                % i - np
                np2 = np * N;
                valid_mod_neg = 1 : np;
                valid_mod_pos = zeros( 1, np );
                if np2 - np + 1 <= np2 - 1
                    valid_mod_pos( 1 : np - 1 ) = np2 - np + 1 : np2 - 1;
                end

                edge = mod( i , np2 );

                if ( ismember( edge , valid_mod_neg ) )
                    a_row_column_indices( j ) = 0;
                else
                    a_row_column_indices( j ) = ...
                        i - np;
                end

                % i + np 
                if ( ismember( edge , valid_mod_pos ) )
                    a_row_column_indices( max_row_elements + 1 - j ) = 0;
                else
                    a_row_column_indices( max_row_elements + 1 - j ) = ...
                        i + np;
                end
            else
                a_row_column_indices( j ) = i - np;
                a_row_column_indices( max_row_elements + 1 - j ) = ...
                    i + np;
            end
        end

        %a_row_column_indices  % debug
        
        % Remove the i-1 or i+1 that conceptually sits on the side mesh
        % boundaries.

        %{
        % Calculate the mod / remainder division to determine whether the
        % row you're in sits on a side boundary or not.
        edge = mod( i, N );

        if edge == 0 && i ~= 1
            a_row_column_indices( a_row_middle - 1 ) ...
                = [ i - 1 ];
            a_row_column_indices( a_row_middle + 1 ) = [];
        elseif edge == 1 && i ~= N^m
            a_row_column_indices( a_row_middle + 1 ) ...
                = [ i + 1 ];
            a_row_column_indices( a_row_middle - 1 ) = [];
        else
            a_row_column_indices( a_row_middle - 1 ) ...
                = [ i - 1 ];
            a_row_column_indices( a_row_middle + 1 ) ...
                = [ i + 1 ];
        end
        %}

        %a_row_column_indices  % debug

        % Ignore column indices outside of 1 and N^m that occurs with i-N,
        % i-1, i+1, i+N
        nz_indices = find( a_row_column_indices > 0 & ...
            a_row_column_indices <= N^m );

        % Use the indices from find() to only use valid nonzero column 
        % indices.
        a_row_nz_column_indices = a_row_column_indices( nz_indices );
        
        % Construct the row values of -1 and 2 * m based on the current i
        % (row) iteration, where column_indices == i is 2 * m and
        % everything else is -1.
        a_row_values = a_row_nz_column_indices;
        a_row_values ( a_row_values ~= i ) = -1;
        a_row_values ( a_row_values == i ) = 2 * m;
        
        %i  % debug
        %a_row_nz_column_indices  % debug
        %a_row_values  % debug

        % Set ir( i + 1 ) based on the number of elements in the current
        % row + ir( i ).  Also set ic and nzA.
        ir( i + 1 ) = ir( i ) + length( a_row_values );
        start_indx = ir( i );
        end_indx = ir( i + 1 ) - 1;
        ic( start_indx : end_indx ) = a_row_nz_column_indices;
        nzA( start_indx : end_indx ) = a_row_values;

        % Compress the oversized ic and nzA arrays.
        if i == a_rows
            ic = ic( 1 : ir(end) - 1 );
            nzA = nzA( 1 : ir(end) - 1 );
        end
    end
end