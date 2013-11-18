function triangulation_rcm ( prefix )

%*****************************************************************************80
%
%% MAIN is the main program for TRIANGULATION_RCM.
%
%  Discussion:
%
%    TRIANGULATION_RCM applies the RCM reordering to a triangulation.
%
%    The user supplies a node file and a triangle file, containing
%    the coordinates of the nodes, and the indices of the nodes that
%    make up each triangle.  Either 3-node or 6-node triangles may
%    be used.
%
%    The program reads the data, computes the adjacency information,
%    carries out the RCM algorithm to get the permutation, applies
%    the permutation to the nodes and triangles, and writes out
%    new node and triangle files that correspond to the RCM permutation.
%
%    Note that node data is normally two dimensional, that is,
%    each node has an X and Y coordinate.  In some applications, it
%    may be desirable to specify more information.  This program
%    will accept node data that includes DIM_NUM entries on each line,
%    as long as DIM_NUM is the same for each entry.  
%
%  Usage:
%
%    triangulation_rcm ( 'prefix' )
%
%    where
%
%    * 'prefix'_nodes.txt contains the node coordinates;
%    * 'prefix'_elements.txt contains the element definitions.
%    * 'prefix'_rcm_nodes.txt contains the reordered node coordinates;
%    * 'prefix'_rcm_elements.txt contains the reordered element definitions.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 September 2009
%
%  Author:
%
%    John Burkardt
%
  timestamp ( );

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TRIANGULATION_RCM\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Read a node dataset of NODE_NUM points in 2 dimensions.\n' );
  fprintf ( 1, '  Read an associated triangulation dataset of TRIANGLE_NUM\n' );
  fprintf ( 1, '  triangles using 3 or 6 nodes.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Apply the RCM reordering (Reverse Cuthill-McKee).\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Reorder the data and write it out to files.\n' );
%
%  The command line argument is the common filename prefix.
%
  if ( nargin < 1 )

    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRIANGULATION_RCM:\n' );

    prefix = input ( ...
      'Please enter the filename prefix:' );

  end
%
%  Create the filenames.
%
  node_filename = strcat ( prefix, '_nodes.txt' );
  element_filename = strcat ( prefix, '_elements.txt' );
  node_rcm_filename = strcat ( prefix, '_rcm_nodes.txt' );
  element_rcm_filename = strcat ( prefix, '_rcm_elements.txt' );
%
%  Read the node data.
%
  [ dim_num, node_num ] = r8mat_header_read ( node_filename );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read the header of ""%s".\n', node_filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Spatial dimension DIM_NUM = %d\n', dim_num );
  fprintf ( 1, '  Number of nodes NODE_NUM  = %d\n', node_num );

  node_xy(1:dim_num,1:node_num) = r8mat_data_read ( ...
    node_filename, dim_num, node_num );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read the data in ""%s".\n', node_filename );

  r8mat_transpose_print_some ( dim_num, node_num, node_xy, ...
    1, 1, dim_num, 5, '  Coordinates of first 5 nodes:' );
%
%  Read the element data.
%
  [ element_order, element_num ] = i4mat_header_read ( element_filename );

  if ( element_order ~= 3 && element_order ~= 6 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRIANGULATION_RCM - Fatal error!\n' );
    fprintf ( 1, '  Data is not for a 3-node or 6-node triangulation.\n' );
    error ( 'TRIANGULATION_RCM - Fatal error!' );
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read the header of ""%s".\n', element_filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Element order = %d\n', element_order );
  fprintf ( 1, '  Number of elements ELEMENT_NUM  = %d\n', element_num );

  element_node(1:element_order,1:element_num) = i4mat_data_read ( ...
    element_filename, element_order, element_num );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read the data in ""%s".\n', element_filename );

  i4mat_transpose_print_some ( element_order, element_num, element_node, ...
    1, 1, element_order, 5, '  First five elements:' );
%
%  Detect and correct 0-based indexing.
%
  element_node = mesh_base_one ( node_num, element_order, element_num, ...
    element_node );
%
%  Create the element neighbor array.
%
  element_neighbor = triangulation_neighbor_triangles ( element_order, ...
    element_num, element_node );
%
%  Following code depends on element order.
%
  if ( element_order == 3 )
%
%  Count the number of adjacencies, and set up the ADJ_ROW 
%  adjacency pointer array.
%
    [ adj_num, adj_row ] = triangulation_order3_adj_count ( node_num, ...
      element_num, element_node, element_neighbor );
%
%  Set up the ADJ adjacency array.
%
    adj = triangulation_order3_adj_set ( node_num, element_num, element_node, ...
      element_neighbor, adj_num, adj_row );

  elseif ( element_order == 6 )
%
%  Count the number of adjacencies, and set up the ADJ_ROW 
%  adjacency pointer array.
%
    [ adj_num, adj_row ] = triangulation_order6_adj_count ( node_num, ...
      element_num, element_node, element_neighbor );
%
%  Set up the ADJ adjacency array.
%
    adj = triangulation_order6_adj_set ( node_num, element_num, element_node, ...
      element_neighbor, adj_num, adj_row );

  end

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  ADJ bandwidth = %d\n', bandwidth );
%
%  Compute the RCM permutation.
%
  perm = genrcm ( node_num, adj_num, adj_row, adj );

  perm_inv = perm_inverse3 ( node_num, perm );

  bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, ...
    perm, perm_inv );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Permuted ADJ bandwidth = %d\n', bandwidth );
%
%  Permute the nodes according to the permutation vector.
%
  node_xy = r8col_permute ( dim_num, node_num, node_xy, perm );
%
%  Permute the node indices in the element array.
%
  for i = 1 : element_order
    for j = 1 : element_num
      node = element_node(i,j);
      element_node(i,j) = perm_inv ( node );
    end
  end
%
%  Write out the node and element data.
%
  r8mat_write ( node_rcm_filename, dim_num, node_num, node_xy );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Created the node file "%s".\n', node_rcm_filename );

  i4mat_write ( element_rcm_filename, element_order, element_num, ...
    element_node );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Created the element file "%s".\n', ...
    element_rcm_filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TRIANGULATION_RCM\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
function value = adj_bandwidth ( node_num, adj_num, adj_row, adj )

%*****************************************************************************80
%
%% ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%  Parameters:
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Output, integer VALUE, the bandwidth of the adjacency
%    matrix.
%
  band_lo = 0;
  band_hi = 0;

  for i = 1 : node_num

    for j = adj_row(i) : adj_row(i+1)-1
      col = adj(j);
      band_lo = max ( band_lo, i - col );
      band_hi = max ( band_hi, col - i );
    end

  end

  value = band_lo + 1 + band_hi;

  return
end
function value = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, ...
  perm, perm_inv )

%*****************************************************************************80
%
%% ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
%
%  Discussion:
%
%    The matrix is defined by the adjacency information and a permutation.
%
%    The routine also computes the bandwidth and the size of the envelope.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%  Parameters:
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Input, integer PERM(NODE_NUM), integer PERM_INV(NODE_NUM), the permutation
%    and inverse permutation.
%
%    Output, integer ADJ_PERM_BANDWIDTH, the bandwidth of the permuted
%    adjacency matrix.
%
  band_lo = 0;
  band_hi = 0;

  for i = 1 : node_num

    for j = adj_row(perm(i)) : adj_row(perm(i)+1)-1
      col = perm_inv(adj(j));
      band_lo = max ( band_lo, i - col );
      band_hi = max ( band_hi, col - i );
    end

  end

  value = band_lo + 1 + band_hi;

  return
end
function [ deg, iccsze, ls ] = degree ( root, adj_num, adj_row, adj, mask, ...
  node_num )

%*****************************************************************************80
%
%% DEGREE computes the degrees of the nodes in the connected component.
%
%  Discussion:
%
%    The connected component is specified by MASK and ROOT.
%    Nodes for which MASK is zero are ignored.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 January 2007
%
%  Author:
%
%    Original FORTRN77 version by Alan George, Joseph Liu.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%  Parameters:
%
%    Input, integer ROOT, the node that defines the connected component.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Input, integer MASK(NODE_NUM), is nonzero for those nodes which are
%    to be considered.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, integer DEG(NODE_NUM), contains, for each  node in the connected
%    component, its degree.
%
%    Output, integer ICCSZE, the number of nodes in the connected component.
%
%    Output, integer LS(NODE_NUM), stores in entries 1 through ICCSIZE the nodes
%    in the connected component, starting with ROOT, and proceeding
%    by levels.
%

%
%  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
%
  ls(1) = root;
  adj_row(root) = -adj_row(root);
  lvlend = 0;
  iccsze = 1;
%
%  LBEGIN is the pointer to the beginning of the current level, and
%  LVLEND points to the end of this level.
%
  while ( 1 )

    lbegin = lvlend + 1;
    lvlend = iccsze;
%
%  Find the degrees of nodes in the current level,
%  and at the same time, generate the next level.
%
    for i = lbegin : lvlend

      node = ls(i);
      jstrt = -adj_row(node);
      jstop = abs ( adj_row(node+1) ) - 1;
      ideg = 0;

      for j = jstrt : jstop

        nbr = adj(j);

        if ( mask(nbr) ~= 0 )

          ideg = ideg + 1;

          if ( 0 <= adj_row(nbr) )
            adj_row(nbr) = -adj_row(nbr);
            iccsze = iccsze + 1;
            ls(iccsze) = nbr;
          end

        end

      end

      deg(node) = ideg;

    end
%
%  Compute the current level width.
%
    lvsize = iccsze - lvlend;
%
%  If the current level width is nonzero, generate another level.
%
    if ( lvsize == 0 )
      break
    end

  end
%
%  Reset ADJ_ROW to its correct sign and return.
%
  for i = 1 : iccsze
    node = ls(i);
    adj_row(node) = -adj_row(node);
  end

  return
end
function column_num = file_column_count ( input_file_name )

%*****************************************************************************80
%
%% FILE_COLUMN_COUNT counts the columns in the first line of a file.
%
%  Discussion:
%
%    The file is assumed to be a simple text file.
%
%    Most lines of the file are presumed to consist of COLUMN_NUM words,
%    separated by spaces.  There may also be some blank lines, and some 
%    comment lines, which have a "#" in column 1.
%
%    The routine tries to find the first non-comment non-blank line and
%    counts the number of words in that line.
%
%    If all lines are blanks or comments, it goes back and tries to analyze
%    a comment line.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    21 February 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string INPUT_FILE_NAME, the name of the file.
%
%    Output, integer COLUMN_NUM, the number of columns in the file.
%
  FALSE = 0;
  TRUE = 1;
%
%  Open the file.
%
  input_unit = fopen ( input_file_name );

  if ( input_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'FILE_COLUMN_COUNT - Error!\n' );
    fprintf ( 1, '  Could not open the file "%s".\n', input_file_name );
    error ( 'FILE_COLUMN_COUNT - Error!' );
  end
%
%  Read one line, but skip blank lines and comment lines.
%  Use FGETL so we drop the newline character!
%
  got_one = FALSE;

  while ( 1 )

    line = fgetl ( input_unit );

    if ( line == -1 )
      break;
    end

    if ( s_len_trim ( line ) == 0 )

    elseif ( line(1) == '#' )

    else
      got_one = TRUE;
      break;
    end

  end

  fclose ( input_unit );

  if ( got_one == FALSE ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'FILE_COLUMN_COUNT - Warning!\n' );
    fprintf ( 1, '  The file does not seem to contain any data.\n' );
    column_num = -1;
    return;
  end

  column_num = s_word_count ( line );

  return
end
function row_num = file_row_count ( input_file_name )

%*****************************************************************************80
%
%% FILE_ROW_COUNT counts the number of row records in a file.
%
%  Discussion:
%
%    Each input line is a "RECORD".
%
%    The records are divided into three groups:
%    
%    * BLANK LINES (nothing but blanks)
%    * COMMENT LINES (begin with a '#')
%    * DATA RECORDS (anything else)
%
%    The value returned by the function is the number of data records.
%
%    By the way, if the MATLAB routine FGETS is used, instead of
%    FGETL, then the variable LINE will include line termination 
%    characters, which means that a blank line would not actually
%    have zero characters.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    31 December 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string INPUT_FILE_NAME, the name of the input file.
%
%    Output, integer ROW_NUM, the number of rows found. 
%
  input_unit = fopen ( input_file_name );

  if ( input_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'FILE_ROW_COUNT - Error!\n' );
    fprintf ( 1, '  Could not open the file "%s".\n', input_file_name );
    error ( 'FILE_ROW_COUNT - Error!' );
  end

  blank_num = 0;
  comment_num = 0;
  row_num = 0;
  
  record_num = 0;

  while ( 1 )

    line = fgetl ( input_unit );

    if ( line == -1 )
      break;
    end

    record_num = record_num + 1;
    record_length = s_len_trim ( line );
    
    if ( record_length <= 0 )
      blank_num = blank_num + 1;
    elseif ( line(1) == '#' )
      comment_num = comment_num + 1;
    else
      row_num = row_num + 1;
    end

  end

  fclose ( input_unit );

  return
end
function perm = genrcm ( node_num, adj_num, adj_row, adj )

%*****************************************************************************80
%
%% GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
%
%  Discussion:
%
%    For each connected component in the graph, the routine obtains
%    an ordering by calling RCM.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    03 January 2007
%
%  Author:
%
%    Original FORTRAN77 version by Alan George, Joseph Liu.
%    MATLAB version by John Burkardt
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%  Parameters:
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Output, integer PERM(NODE_NUM), the RCM ordering.
%
%  Local Parameters:
%
%    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
%    structure.  The level structure is stored in the currently unused
%    spaces in the permutation vector PERM.
%
%    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
%
  mask(1:node_num) = 1;

  num = 1;

  for i = 1 : node_num
%
%  For each masked connected component...
%
    if ( mask(i) ~= 0 )

      root = i;
%
%  Find a pseudo-peripheral node ROOT.  The level structure found by
%  ROOT_FIND is stored starting at PERM(NUM).
%
      [ root, level_num, level_row, level ] = root_find ( root, ...
        adj_num, adj_row, adj, mask, node_num );
%
%  RCM orders the component using ROOT as the starting node.
%
      [ mask, level, iccsze ] = rcm ( root, adj_num, adj_row, adj, mask, ...
        node_num );

      perm(num:num+iccsze-1) = level(1:iccsze);

      num = num + iccsze;
%
%  We can stop once every node is in one of the connected components.
%
      if ( node_num < num )
        return
      end

    end

  end

  return
end
function isgn = i4col_compare ( m, n, a, i, j )

%*****************************************************************************80
%
%% I4COL_COMPARE compares columns I and J of a integer array.
%
%  Example:
%
%    Input:
%
%      M = 3, N = 4, I = 2, J = 4
%
%      A = (
%        1  2  3  4
%        5  6  7  8
%        9 10 11 12 )
%
%    Output:
%
%      ISGN = -1
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 June 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns.
%
%    Input, integer A(M,N), an array of N columns of vectors of length M.
%
%    Input, integer I, J, the columns to be compared.
%    I and J must be between 1 and N.
%
%    Output, integer ISGN, the results of the comparison:
%    -1, column I < column J,
%     0, column I = column J,
%    +1, column J < column I.
%

%
%  Check.
%
  if ( i < 1)
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4COL_COMPARE - Fatal error!\n' );
    fprintf ( 1, '  Column index I = %d < 1.\n', i );
    error ( 'I4COL_COMPARE - Fatal error!' );
  end

  if ( n < i )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4COL_COMPARE - Fatal error!\n' );
    fprintf ( 1, '  N = %d < column index I = %d.\n', n, i );
    error ( 'I4COL_COMPARE - Fatal error!' );
  end

  if ( j < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4COL_COMPARE - Fatal error!\n' );
    fprintf ( 1, '  Column index J = %d < 1.\n', j );
    error ( 'I4COL_COMPARE - Fatal error!' );
  end

  if ( n < j )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4COL_COMPARE - Fatal error!\n' );
    fprintf ( 1, '  N = %d < column index J = %d.\n', n, j );
    error ( 'I4COL_COMPARE - Fatal error!' );
  end

  isgn = 0;

  if ( i == j )
    return
  end

  k = 1;

  while ( k <= m )

    if ( a(k,i) < a(k,j) )
      isgn = -1;
      return
    elseif ( a(k,j) < a(k,i) )
      isgn = +1;
      return
    end

    k = k + 1;

  end

  return
end
function a = i4col_sort_a ( m, n, a )

%*****************************************************************************80
%
%% I4COL_SORT_A ascending sorts an I4COL.
%
%  Discussion:
%
%    In lexicographic order, the statement "X < Y", applied to two real
%    vectors X and Y of length M, means that there is some index I, with
%    1 <= I <= M, with the property that
%
%      X(J) = Y(J) for J < I,
%    and
%      X(I) < Y(I).
%
%    In other words, the first time they differ, X is smaller.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the number of rows of A, and the length of
%    a vector of data.
%
%    Input, integer N, the number of columns of A.
%
%    Input, integer A(M,N), the array of N columns of M-vectors.
%
%    Output, integer A(M,N), the columns of A have been sorted in ascending
%    lexicographic order.
%
  if ( m <= 0 )
    return
  end

  if ( n <= 1 )
    return
  end
%
%  Initialize.
%
  indx = 0;
  isgn = 0;
%
%  Call the external heap sorter.
%
  while ( 1 )

    [ indx, i, j ] = sort_heap_external ( n, indx, isgn );
%
%  Interchange the I and J objects.
%
    if ( 0 < indx )

      a = i4col_swap ( m, n, a, i, j );
%
%  Compare the I and J objects.
%
    elseif ( indx < 0 )

      isgn = i4col_compare ( m, n, a, i, j );

    elseif ( indx == 0 )

      break

    end

  end

  return
end
function a = i4col_swap ( m, n, a, i, j )

%*****************************************************************************80
%
%% I4COL_SWAP swaps columns I and J of a integer array of column data.
%
%  Example:
%
%    Input:
%
%      M = 3, N = 4, I = 2, J = 4
%
%      A = (
%        1  2  3  4
%        5  6  7  8
%        9 10 11 12 )
%
%    Output:
%
%      A = (
%        1  4  3  2
%        5  8  7  6
%        9 12 11 10 )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns in the array.
%
%    Input, integer A(M,N), an array of N columns of length M.
%
%    Input, integer I, J, the columns to be swapped.
%
%    Output, integer A(M,N), the array, with columns I and J swapped.
%
  if ( i < 1 || n < i || j < 1 || n < j )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4COL_SWAP - Fatal error!\n' );
    fprintf ( 1, '  I or J is out of bounds.\n' );
    fprintf ( 1, '  I =    %d\n', i );
    fprintf ( 1, '  J =    %d\n', j );
    fprintf ( 1, '  N =    %d\n', n );
    error ( 'I4COL_SWAP - Fatal error!' );
  end

  if ( i == j )
    return
  end

  col(1:m) = a(1:m,i)';
  a(1:m,i) = a(1:m,j);
  a(1:m,j) = col(1:m)';

  return
end
function table = i4mat_data_read ( input_filename, m, n )

%*****************************************************************************80
%
%% I4MAT_DATA_READ reads data from an I4MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    27 January 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string INPUT_FILENAME, the name of the input file.
%
%    Input, integer M, N, the number of rows and columns in the data.
%
%    Output, integer TABLE(M,N), the point coordinates.
%
  table = zeros ( m, n );
%
%  Build up the format string for reading M real numbers.
%
  string = ' ';

  for i = 0 : m
    string = strcat ( string, ' %d' );
  end

  input_unit = fopen ( input_filename );

  if ( input_unit < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4MAT_DATA_READ - Error!\n' );
    fprintf ( 1, '  Could not open the input file.\n' );
    error ( 'I4MAT_DATA_READ - Error!' );
  end

  i = 0;

  while ( i < n )

    line = fgets ( input_unit );

    if ( line == -1 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'I4MAT_DATA_READ - Error!\n' );
      fprintf ( 1, '  End of input while reading data.\n' );
      error ( 'I4MAT_DATA_READ - Error!' );
    end

    if ( line(1) == '#' )

    elseif ( s_len_trim ( line ) == 0 )
      
    else

      [ x, count ] = sscanf ( line, string );

      if ( count == m )
        i = i + 1;
        table(1:m,i) = x(1:m);
      end

    end

  end

  fclose ( input_unit );

  return
end
function [ m, n ] = i4mat_header_read ( input_filename )

%*****************************************************************************80
%
%% I4MAT_HEADER_READ reads the header from an I4MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string INPUT_FILENAME, the name of the input file.
%
%    Output, integer M, the spatial dimension.
%
%    Output, integer N, the number of points.
%
  m = file_column_count ( input_filename );

  if ( m <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4MAT_HEADER_READ - Fatal error!\n' );
    fprintf ( 1, '  There was some kind of I/O problem while trying\n' );
    fprintf ( 1, '  to count the number of data columns in\n' );
    fprintf ( 1, '  the file %s.\n', input_filename );
  end

  n = file_row_count ( input_filename );

  if ( n <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4MAT_HEADER_READ - Fatal error!\n' );
    fprintf ( 1, '  There was some kind of I/O problem while trying\n' );
    fprintf ( 1, '  to count the number of data rows in\n' );
    fprintf ( 1, '  the file %s\n', input_filename );
  end

  return
end
function i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

%*****************************************************************************80
%
%% I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 June 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns.
%
%    Input, integer A(M,N), an M by N matrix to be printed.
%
%    Input, integer ILO, JLO, the first row and column to print.
%
%    Input, integer IHI, JHI, the last row and column to print.
%
%    Input, string TITLE, an optional title.
%
  incx = 10;

  if ( 0 < s_len_trim ( title ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, '%s\n', title );
  end

  for i2lo = max ( ilo, 1 ) : incx : min ( ihi, m )

    i2hi = i2lo + incx - 1;
    i2hi = min ( i2hi, m );
    i2hi = min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    fprintf ( 1, '\n' );
    fprintf ( 1, '  Row: ' );
    for i = i2lo : i2hi
      fprintf ( 1, '%7d  ', i );
    end
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Col\n' );
    fprintf ( 1, '\n' );

    j2lo = max ( jlo, 1 );
    j2hi = min ( jhi, n );

    for j = j2lo : j2hi

      fprintf ( 1, '%5d  ', j );
      for i2 = 1 : inc
        i = i2lo - 1 + i2;
        fprintf ( 1, '%7d  ', a(i,j) );
      end
      fprintf ( 1, '\n' );

    end

  end

  return
end
function i4mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% I4MAT_WRITE writes an I4MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, integer TABLE(M,N), the points.
%
%    Input, logical HEADER, is TRUE if the header is to be included.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'I4MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %12d', round ( table(i,j) ) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function a_out = i4vec_heap_d ( n, a )

%*****************************************************************************80
%
%% I4VEC_HEAP_D reorders an I4VEC into an descending heap.
%
%  Definition:
%
%    A descending heap is an array A with the property that, for every index J,
%    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
%    2*J and 2*J+1 are legal).
%
%  Diagram:
%
%                  A(1)
%                /      \
%            A(2)         A(3)
%          /     \        /  \
%      A(4)       A(5)  A(6) A(7)
%      /  \       /   \
%    A(8) A(9) A(10) A(11)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    24 February 2004
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Albert Nijenhuis, Herbert Wilf,
%    Combinatorial Algorithms,
%    Academic Press, 1978, second edition,
%    ISBN 0-12-519260-6.
%
%  Parameters:
%
%    Input, integer N, the size of the input array.
%
%    Input, integer A(N), an unsorted array.
%
%    Output, integer A_OUT(N), the array has been reordered into a heap.
%
  a_out(1:n) = a(1:n);
%
%  Only nodes N/2 down to 1 can be "parent" nodes.
%
  for i = floor ( n/2 ) : -1 : 1
%
%  Copy the value out of the parent node.
%  Position IFREE is now "open".
%
    key = a_out(i);
    ifree = i;

    while ( 1 )
%
%  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
%  IFREE.  (One or both may not exist because they exceed N.)
%
      m = 2 * ifree;
%
%  Does the first position exist?
%
      if ( n < m )
        break;
      end
%
%  Does the second position exist?
%
      if ( m + 1 <= n )
%
%  If both positions exist, take the larger of the two values,
%  and update M if necessary.
%
        if ( a_out(m) < a_out(m+1) )
          m = m + 1;
        end

      end
%
%  If the large descendant is larger than KEY, move it up,
%  and update IFREE, the location of the free position, and
%  consider the descendants of THIS position.
%
      if ( a_out(m) <= key )
        break;
      end

      a_out(ifree) = a_out(m);
      ifree = m;

    end
%
%  Once there is no more shifting to do, KEY moves into the free spot IFREE.
%
    a_out(ifree) = key;

  end

  return
end
function a = i4vec_reverse ( n, a )

%*****************************************************************************80
%
%% I4VEC_REVERSE reverses the elements of an I4VEC.
%
%  Example:
%
%    Input:
%
%      N = 5,
%      A = ( 11, 12, 13, 14, 15 ).
%
%    Output:
%
%      A = ( 15, 14, 13, 12, 11 ).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    17 April 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of entries in the array.
%
%    Input, integer A(N), the array to be reversed.
%
%    Output, integer A(N), the reversed array.
%
  b(1:n) = a(n:-1:1);
  a(1:n) = b(1:n);

  return
end
function a_sorted = i4vec_sort_heap_a ( n, a )

%*****************************************************************************80
%
%% I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    01 May 2004
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    A Nijenhuis and H Wilf,
%    Combinatorial Algorithms,
%    Academic Press, 1978, second edition,
%    ISBN 0-12-519260-6.
%
%  Parameters:
%
%    Input, integer N, the number of entries in the array.
%
%    Input, integer A(N), the array to be sorted;
%
%    Output, integer A_SORTED(N), the sorted array.
%
  if ( n < 1 )
    a_sorted = [];
    return
  end

  if ( n == 1 )
    a_sorted(1) = a(1);
    return
  end
%
%  1: Put A into descending heap form.
%
  a_sorted = i4vec_heap_d ( n, a );
%
%  2: Sort A.
%
%  The largest object in the heap is in A(1).
%  Move it to position A(N).
%
  temp = a_sorted(1);
  a_sorted(1) = a_sorted(n);
  a_sorted(n) = temp;
%
%  Consider the diminished heap of size N1.
%
  for n1 = n-1 : -1 : 2
%
%  Restore the heap structure of A(1) through A(N1).
%
    a_sorted(1:n1) = i4vec_heap_d ( n1, a_sorted );
%
%  Take the largest object from A(1) and move it to A(N1).
%
    temp = a_sorted(1);
    a_sorted(1) = a_sorted(n1);
    a_sorted(n1) = temp;

  end

  return
end
function [ mask, level_num, level_row, level ] = level_set ( root, adj_num, ...
  adj_row, adj, mask, node_num )

%*****************************************************************************80
%
%% LEVEL_SET generates the connected level structure rooted at a given node.
%
%  Discussion:
%
%    Only nodes for which MASK is nonzero will be considered.
%
%    The root node chosen by the user is assigned level 1, and masked.
%    All (unmasked) nodes reachable from a node in level 1 are
%    assigned level 2 and masked.  The process continues until there
%    are no unmasked nodes adjacent to any node in the current level.
%    The number of levels may vary between 2 and NODE_NUM.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 January 2007
%
%  Author:
%
%    Original FORTRAN77 version by Alan George, Joseph Liu.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%  Parameters:
%
%    Input, integer ROOT, the node at which the level structure
%    is to be rooted.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Input, integer MASK(NODE_NUM).  On input, only nodes with nonzero
%    MASK are to be processed.  On output, those nodes which were included
%    in the level set have MASK set to 1.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, integer MASK(NODE_NUM).  On input, only nodes with nonzero
%    MASK are to be processed.  On output, those nodes which were included
%    in the level set have MASK set to 1.
%
%    Output, integer LEVEL_NUM, the number of levels in the level
%    structure.  ROOT is in level 1.  The neighbors of ROOT
%    are in level 2, and so on.
%
%    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the rooted
%    level structure.
%
  mask(root) = 0;
  level(1) = root;
  level_num = 0;
  lvlend = 0;
  iccsze = 1;
%
%  LBEGIN is the pointer to the beginning of the current level, and
%  LVLEND points to the end of this level.
%
  while ( 1 )

    lbegin = lvlend + 1;
    lvlend = iccsze;
    level_num = level_num + 1;
    level_row(level_num) = lbegin;
%
%  Generate the next level by finding all the masked neighbors of nodes
%  in the current level.
%
    for i = lbegin : lvlend

      node = level(i);
      jstrt = adj_row(node);
      jstop = adj_row(node+1)-1;

      for j = jstrt : jstop

        nbr = adj(j);

        if ( mask(nbr) ~= 0 )
          iccsze = iccsze + 1;
          level(iccsze) = nbr;
          mask(nbr) = 0;
        end

      end

    end
%
%  Compute the current level width (the number of nodes encountered.)
%  If it is positive, generate the next level.
%
    lvsize = iccsze - lvlend;

    if ( lvsize <= 0 )
      break
    end

  end

  level_row(level_num+1) = lvlend + 1;
%
%  Reset MASK to 1 for the nodes in the level structure.
%
  mask(level(1:iccsze)) = 1;

  return
end
function element_node = mesh_base_one ( node_num, element_order, ...
  element_num, element_node )

%*****************************************************************************80
%
%% MESH_BASE_ONE ensures that the element definition is one-based.
%
%  Discussion:
%
%    The ELEMENT_NODE array contains nodes indices that form elements.
%    The convention for node indexing might start at 0 or at 1.
%    Since a MATLAB program will naturally assume a 1-based indexing, it is
%    necessary to check a given element definition and, if it is actually
%    0-based, to convert it.
%
%    This function attempts to detect 0-based node indexing and correct it.
%
%    Thanks to Feifei Xu for pointing out that I was subtracting 1 when I
%    should have been adding 1!  29 November 2012.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    29 November 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input/output, integer ELEMENT_NODE(ELEMENT_ORDE,ELEMENT_NUM), the element
%    definitions.
%
  node_min = min ( min ( element_node(1:element_order,1:element_num) ) );
  node_max = max ( max ( element_node(1:element_order,1:element_num) ) );

  if ( node_min == 0 && node_max == node_num - 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MESH_BASE_ONE:\n' );
    fprintf ( 1, '  The element indexing appears to be 0-based!\n' );
    fprintf ( 1, '  This will be converted to 1-based.\n' );
    element_node(1:element_order,1:element_num) = ...
      element_node(1:element_order,1:element_num) + 1;
  elseif ( node_min == 1 && node_max == node_num )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MESH_BASE_ONE:\n' );
    fprintf ( 1, '  The element indexing appears to be 1-based!\n' );
    fprintf ( 1, '  No conversion is necessary.\n' );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MESH_BASE_ONE - Warning!\n' );
    fprintf ( 1, '  The element indexing is not of a recognized type.\n' );
    fprintf ( 1, '  NODE_MIN = %d\n', node_min );
    fprintf ( 1, '  NODE_MAX = %d\n', node_max );
    fprintf ( 1, '  NODE_NUM = %d\n', node_num );
  end

  return
end
function perm_inv = perm_inverse3 ( n, perm )

%*****************************************************************************80
%
%% PERM_INVERSE3 produces the inverse of a given permutation.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of items permuted.
%
%    Input, integer PERM(N), a permutation.
%
%    Output, integer PERM_INV(N), the inverse permutation.
%
  perm_inv(perm(1:n)) = ( 1 : n );

  return
end
function a = r8col_permute ( m, n, a, p )

%*****************************************************************************80
%
%% R8COL_PERMUTE permutes an R8COL in place.
%
%  Discussion:
%
%    An R8COL is an M by N array of R8's, regarded
%    as an array of N columns of length M.
%
%    This routine permutes an array of real "objects", but the same
%    logic can be used to permute an array of objects of any arithmetic
%    type, or an array of objects of any complexity.  The only temporary
%    storage required is enough to store a single object.  The number
%    of data movements made is N + the number of cycles of order 2 or more,
%    which is never more than N + N/2.
%
%  Example:
%
%    Input:
%
%      M = 2
%      N = 5
%      P = (   2,    4,    5,    1,    3 )
%      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
%          (11.0, 22.0, 33.0, 44.0, 55.0 )
%
%    Output:
%
%      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
%             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the size of the objects.
%
%    Input, integer N, the number of objects.
%
%    Input, real A(M,N), the array to be permuted.
%
%    Input, integer P(N), the permutation.  P(I) = J means
%    that the I-th element of the output array should be the J-th
%    element of the input array.  P must be a legal permutation
%    of the integers from 1 to N, otherwise the algorithm will
%    fail catastrophically.
%
%    Output, real A(M,N), the permuted array.
%

%
%  Search for the next element of the permutation that has not been used.
%
  for istart = 1 : n

    if ( p(istart) < 0 )

      continue

    elseif ( p(istart) == istart )

      p(istart) = - p(istart);
      continue

    else

      a_temp(1:m) = a(1:m,istart);
      iget = istart;
%
%  Copy the new value into the vacated entry.
%
      while ( 1 )

        iput = iget;
        iget = p(iget);

        p(iput) = - p(iput);

        if ( iget < 1 || n < iget )
          fprintf ( 1, '\n' );
          fprintf ( 1, 'R8COL_PERMUTE - Fatal error!\n' );
          fprintf ( 1, '  A permutation index is out of range.\n' );
          fprintf ( 1, '  P(%d) = %d\n', iput, iget );
          error ( 'R8COL_PERMUTE - Fatal error!' );
        end

        if ( iget == istart )
          a(1:m,iput) = a_temp(1:m)';
          break
        end

        a(1:m,iput) = a(1:m,iget);

      end

    end

  end
%
%  Restore the signs of the entries.
%
%  p(1:n) = - p(1:n);

  return
end
function table = r8mat_data_read ( input_filename, m, n )

%*****************************************************************************80
%
%% R8MAT_DATA_READ reads data from an R8MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    27 January 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string INPUT_FILENAME, the name of the input file.
%
%    Input, integer M, N, the number of rows and columns of data.
%
%    Output, real TABLE(M,N), the point coordinates.
%
  table = zeros ( m, n );
%
%  Build up the format string for reading M real numbers.
%
  string = ' ';

  for i = 0 : m
    string = strcat ( string, ' %f' );
  end

  input_unit = fopen ( input_filename );

  if ( input_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_DATA_READ - Error!\n' );
    fprintf ( 1, '  Could not open the file.\n' );
    error ( 'R8MAT_DATA_READ - Error!' );
  end

  i = 0;

  while ( i < n )

    line = fgets ( input_unit );

    if ( line == -1 )
      break;
    end

    if ( line(1) == '#' )

    elseif ( s_len_trim ( line ) == 0 )
      
    else

      [ x, count ] = sscanf ( line, string );

      if ( count == m )
        i = i + 1;
        table(1:m,i) = x(1:m);
      end

    end

  end

  fclose ( input_unit );

  return
end
function [ m, n ] = r8mat_header_read ( input_filename )

%*****************************************************************************80
%
%% R8MAT_HEADER_READ reads the header from an R8MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string INPUT_FILENAME, the name of the input file.
%
%    Output, integer M, the spatial dimension.
%
%    Output, integer N, the number of points.
%
  m = file_column_count ( input_filename );

  if ( m <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_HEADER_READ - Fatal error!\n' );
    fprintf ( 1, '  There was some kind of I/O problem while trying\n' );
    fprintf ( 1, '  to count the number of data columns in\n' );
    fprintf ( 1, '  the file %s.\n', input_filename );
  end

  n = file_row_count ( input_filename );

  if ( n <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_HEADER_READ - Fatal error!\n' );
    fprintf ( 1, '  There was some kind of I/O problem while trying\n' );
    fprintf ( 1, '  to count the number of data rows in\n' );
    fprintf ( 1, '  the file %s\n', input_filename );
  end

  return
end
function r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

%*****************************************************************************80
%
%% R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 May 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns.
%
%    Input, real A(M,N), an M by N matrix to be printed.
%
%    Input, integer ILO, JLO, the first row and column to print.
%
%    Input, integer IHI, JHI, the last row and column to print.
%
%    Input, string TITLE, an optional title.
%
  incx = 5;

  if ( 0 < s_len_trim ( title ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, '%s\n', title );
  end

  for i2lo = max ( ilo, 1 ) : incx : min ( ihi, m )

    i2hi = i2lo + incx - 1;
    i2hi = min ( i2hi, m );
    i2hi = min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;
    
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Row: ' );
    for i = i2lo : i2hi
      fprintf ( 1, '%7d       ', i );
    end
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Col\n' );

    j2lo = max ( jlo, 1 );
    j2hi = min ( jhi, n );

    for j = j2lo : j2hi

      fprintf ( 1, '%5d ', j );
      for i2 = 1 : inc
        i = i2lo - 1 + i2;
        fprintf ( 1, '%12f', a(i,j) );
      end
      fprintf ( 1, '\n' );

    end

  end

  return
end
function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% R8MAT_WRITE writes an R8MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  For greater precision, try:
%
%     fprintf ( output_unit, '  %24,16f', table(i,j) );
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %14f', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function [ mask, perm, iccsze ] = rcm ( root, adj_num, adj_row, adj, ...
  mask, node_num )

%*****************************************************************************80
%
%% RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
%
%  Discussion:
%
%    The connected component is specified by a node ROOT and a mask.
%    The numbering starts at the root node.
%
%    An outline of the algorithm is as follows:
%
%    X(1) = ROOT.
%
%    for ( I = 1 to N-1)
%      Find all unlabeled neighbors of X(I),
%      assign them the next available labels, in order of increasing degree.
%
%    When done, reverse the ordering.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 January 2007
%
%  Author:
%
%    Original FORTRAN77 version by Alan George, Joseph Liu.
%    MATLAB version by John Burkardt
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%  Parameters:
%
%    Input, integer ROOT, the node that defines the connected component.
%    It is used as the starting point for the RCM ordering.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Input, integer MASK(NODE_NUM), a mask for the nodes.  Only
%    those nodes with nonzero input mask values are considered by the
%    routine.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, integer MASK(NODE_NUM), the nodes numbered by RCM have 
%    their mask values set to zero.
%
%    Output, integer PERM(NODE_NUM), the RCM ordering.
%
%    Output, integer ICCSZE, the size of the connected component
%    that has been numbered.
%
%  Local Parameters:
%
%    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
%    the degree of the nodes in the section graph specified by mask and root.
%

%
%  Find the degrees of the nodes in the component specified by MASK and ROOT.
%
  [ deg, iccsze, perm ] = degree ( root, adj_num, adj_row, adj, mask, ...
    node_num );

  mask(root) = 0;

  if ( iccsze <= 1 )
    return
  end

  lvlend = 0;
  lnbr = 1;
%
%  LBEGIN and LVLEND point to the beginning and
%  the end of the current level respectively.
%
  while ( lvlend < lnbr )

    lbegin = lvlend + 1;
    lvlend = lnbr;

    for i = lbegin : lvlend
%
%  For each node in the current level...
%
      node = perm(i);
      jstrt = adj_row(node);
      jstop = adj_row(node+1) - 1;
%
%  Find the unnumbered neighbors of NODE.
%
%  FNBR and LNBR point to the first and last neighbors
%  of the current node in PERM.
%
      fnbr = lnbr + 1;

      for j = jstrt : jstop

        nbr = adj(j);

        if ( mask(nbr) ~= 0 )
          lnbr = lnbr + 1;
          mask(nbr) = 0;
          perm(lnbr) = nbr;
        end

      end
%
%  If no neighbors, skip to next node in this level.
%
      if ( lnbr <= fnbr )
        continue
      end
%
%  Sort the neighbors of NODE in increasing order by degree.
%  Linear insertion is used.
%
      k = fnbr;

      while ( k < lnbr )

        l = k;
        k = k + 1;
        nbr = perm(k);

        while ( fnbr < l )

          lperm = perm(l);

          if ( deg(lperm) <= deg(nbr) )
            break
          end

          perm(l+1) = lperm;
          l = l-1;

        end

        perm(l+1) = nbr;

      end

    end

  end
%
%  We now have the Cuthill-McKee ordering.  Reverse it.
%
  perm = i4vec_reverse ( iccsze, perm );

  return
end
function [ root, level_num, level_row, level ] = root_find ( root, ...
  adj_num, adj_row, adj, mask, node_num )

%*****************************************************************************80
%
%% ROOT_FIND finds a pseudo-peripheral node.
%
%  Discussion:
%
%    The diameter of a graph is the maximum distance (number of edges)
%    between any two nodes of the graph.
%
%    The eccentricity of a node is the maximum distance between that
%    node and any other node of the graph.
%
%    A peripheral node is a node whose eccentricity equals the
%    diameter of the graph.
%
%    A pseudo-peripheral node is an approximation to a peripheral node;
%    it may be a peripheral node, but all we know is that we tried our
%    best.
%
%    The routine is given a graph, and seeks pseudo-peripheral nodes,
%    using a modified version of the scheme of Gibbs, Poole and
%    Stockmeyer.  It determines such a node for the section subgraph
%    specified by MASK and ROOT.
%
%    The routine also determines the level structure associated with
%    the given pseudo-peripheral node; that is, how far each node
%    is from the pseudo-peripheral node.  The level structure is
%    returned as a list of nodes LS, and pointers to the beginning
%    of the list of nodes that are at a distance of 0, 1, 2, ...,
%    NODE_NUM-1 from the pseudo-peripheral node.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 January 2007
%
%  Author:
%
%    Original FORTRAN77 version by Alan George, Joseph Liu.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Alan George, Joseph Liu,
%    Computer Solution of Large Sparse Positive Definite Systems,
%    Prentice Hall, 1981.
%
%    Norman Gibbs, William Poole, Paul Stockmeyer,
%    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
%    SIAM Journal on Numerical Analysis,
%    Volume 13, pages 236-250, 1976.
%
%    Norman Gibbs,
%    Algorithm 509: A Hybrid Profile Reduction Algorithm,
%    ACM Transactions on Mathematical Software,
%    Volume 2, pages 378-387, 1976.
%
%  Parameters:
%
%    Input, integer ROOT a node in the the component of the graph for 
%    which a pseudo-peripheral node is sought.
%
%    Input, integer ADJ_NUM, the number of adjacency entries.
%
%    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
%    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
%
%    Input, integer ADJ(ADJ_NUM), the adjacency structure.
%    For each row, it contains the column indices of the nonzero entries.
%
%    Input, integer MASK(NODE_NUM), specifies a section subgraph.  Nodes
%    for which MASK is zero are ignored by FNROOT.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, integer ROOT, the pseudo-peripheral node obtained.
%
%    Output, integer LEVEL_NUM, is the number of levels in the level structure
%    rooted at the node ROOT.
%
%    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
%    level structure array pair containing the level structure found.
%

%
%  Determine the level structure rooted at ROOT.
%
  [ mask, level_num, level_row, level ] = level_set ( root, adj_num, ...
    adj_row, adj, mask, node_num );
%
%  Count the number of nodes in this level structure.
%
  iccsze = level_row(level_num+1) - 1;
%
%  Extreme case:
%    A complete graph has a level set of only a single level.
%    Every node is equally good (or bad).
%
  if ( level_num == 1 )
    return
  end
%
%  Extreme case:
%    A "line graph" 0--0--0--0--0 has every node in its only level.
%    By chance, we've stumbled on the ideal root.
%
  if ( level_num == iccsze )
    return
  end
%
%  Pick any node from the last level that has minimum degree
%  as the starting point to generate a new level set.
%
  while ( 1 )

    mindeg = iccsze;

    jstrt = level_row(level_num);
    root = level(jstrt);

    if ( jstrt < iccsze )

      for j = jstrt : iccsze

        node = level(j);
        ndeg = 0;
        kstrt = adj_row(node);
        kstop = adj_row(node+1)-1;

        for k = kstrt : kstop
          nabor = adj(k);
          if ( 0 < mask(nabor) )
            ndeg = ndeg+1;
          end
        end

        if ( ndeg < mindeg )
          root = node;
          mindeg = ndeg;
        end

      end

    end
%
%  Generate the rooted level structure associated with this node.
%
    [ mask, level_num2, level_row, level ] = level_set ( root, adj_num, ...
      adj_row, adj, mask, node_num );
%
%  If the number of levels did not increase, accept the new ROOT.
%
    if ( level_num2 <= level_num )
      break
    end

    level_num = level_num2;
%
%  In the unlikely case that ROOT is one endpoint of a line graph,
%  we can exit now.
%
    if ( iccsze <= level_num )
      break
    end

  end

  return
end
function len = s_len_trim ( s )

%*****************************************************************************80
%
%% S_LEN_TRIM returns the length of a character string to the last nonblank.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 June 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string S, the string to be measured.
%
%    Output, integer LEN, the length of the string up to the last nonblank.
%
  len = length ( s );

  while ( 0 < len )
    if ( s(len) ~= ' ' )
      return
    end
    len = len - 1;
  end

  return
end
function word_num = s_word_count ( s )

%*****************************************************************************80
%
%% S_WORD_COUNT counts the number of "words" in a string.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 January 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string S, the string to be examined.
%
%    Output, integer WORD_NUM, the number of "words" in the string.
%    Words are presumed to be separated by one or more blanks.
%
  FALSE = 0;
  TRUE = 1;

  word_num = 0;
  s_length = length ( s );

  if ( s_length <= 0 )
    return;
  end

  blank = TRUE;

  for i = 1 : s_length

    if ( s(i) == ' ' )
      blank = TRUE;
    elseif ( blank == TRUE )
      word_num = word_num + 1;
      blank = FALSE;
    end

  end

  return
end
function [ indx, i, j ] = sort_heap_external ( n, indx, isgn )

%*****************************************************************************80
%
%% SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
%
%  Discussion:
%
%    The actual list of data is not passed to the routine.  Hence this
%    routine may be used to sort integers, reals, numbers, names,
%    dates, shoe sizes, and so on.  After each call, the routine asks
%    the user to compare or interchange two items, until a special
%    return value signals that the sorting is completed.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 February 2004
%
%  Author:
%
%    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
%    MATLAB version by John Burkardt
%
%  Reference:
%
%    Albert Nijenhuis, Herbert Wilf.
%    Combinatorial Algorithms,
%    Academic Press, 1978, second edition,
%    ISBN 0-12-519260-6.
%
%  Parameters:
%
%    Input, integer N, the number of items to be sorted.
%
%    Input, integer INDX, the main communication signal.
%    The user must set INDX to 0 before the first call.
%    Thereafter, the user should set the input value of INDX
%    to the output value from the previous call.
%
%    Input, integer ISGN, results of comparison of elements I and J.
%    (Used only when the previous call returned INDX less than 0).
%    ISGN <= 0 means I is less than or equal to J;
%    0 <= ISGN means I is greater than or equal to J.
%
%    Output, integer INDX, the main communication signal.
%    If INDX is
%
%      greater than 0, the user should:
%      * interchange items I and J;
%      * call again.
%
%      less than 0, the user should:
%      * compare items I and J;
%      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
%      * call again.
%
%      equal to 0, the sorting is done.
%
%    Output, integer I, J, the indices of two items.
%    On return with INDX positive, elements I and J should be interchanged.
%    On return with INDX negative, elements I and J should be compared, and
%    the result reported in ISGN on the next call.
%
  persistent i_save;
  persistent j_save;
  persistent k;
  persistent k1;
  persistent n1;
%
%  INDX = 0: This is the first call.
%
  if ( indx == 0 )
      
    i_save = -1;
    j_save = -1;
    k = floor ( n / 2 );
    k1 = k;
    n1 = n;
%
%  INDX < 0: The user is returning the results of a comparison.
%
  elseif ( indx < 0 )

    if ( indx == -2 )

      if ( isgn < 0 )
        i_save = i_save + 1;
      end

      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    end

    if ( 0 < isgn )
      indx = 2;
      i = i_save;
      j = j_save;
      return;
    end

    if ( k <= 1 )

      if ( n1 == 1 )
        i_save = 0;
        j_save = 0;
        indx = 0;
      else
        i_save = n1;
        n1 = n1 - 1;
        j_save = 1;
        indx = 1;
      end

      i = i_save;
      j = j_save;
      return;

    end

    k = k - 1;
    k1 = k;
%
%  0 < INDX, the user was asked to make an interchange.
%
  elseif ( indx == 1 )

    k1 = k;

  end

  while ( 1 )

    i_save = 2 * k1;

    if ( i_save == n1 )
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    elseif ( i_save <= n1 )
      j_save = i_save + 1;
      indx = -2;
      i = i_save;
      j = j_save;
      return;
    end

    if ( k <= 1 )
      break;
    end

    k = k - 1;
    k1 = k;

  end

  if ( n1 == 1 )
    i_save = 0;
    j_save = 0;
    indx = 0;
    i = i_save;
    j = j_save;
  else
    i_save = n1;
    n1 = n1 - 1;
    j_save = 1;
    indx = 1;
    i = i_save;
    j = j_save;
  end

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
function triangle_neighbor = triangulation_neighbor_triangles ( ...
  triangle_order, triangle_num, triangle_node )

%*****************************************************************************80
%
%% TRIANGULATION_NEIGHBOR_TRIANGLES determines triangle neighbors.
%
%  Discussion:
%
%    A triangulation of a set of nodes can be completely described by
%    the coordinates of the nodes, and the list of nodes that make up
%    each triangle.  However, in some cases, it is necessary to know
%    triangle adjacency information, that is, which triangle, if any,
%    is adjacent to a given triangle on a particular side.
%
%    This routine creates a data structure recording this information.
%
%    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
%    data items.
%
%    This routine was modified to use columns instead of rows.
%
%  Example:
%
%    The input information from TRIANGLE_NODE:
%
%    Triangle   Nodes
%    --------   ---------------
%     1         3      4      1
%     2         3      1      2
%     3         3      2      8
%     4         2      1      5
%     5         8      2     13
%     6         8     13      9
%     7         3      8      9
%     8        13      2      5
%     9         9     13      7
%    10         7     13      5
%    11         6      7      5
%    12         9      7      6
%    13        10      9      6
%    14         6      5     12
%    15        11      6     12
%    16        10      6     11
%
%    The output information in TRIANGLE_NEIGHBOR:
%
%    Triangle  Neighboring Triangles
%    --------  ---------------------
%
%     1        -1     -1      2
%     2         1      4      3
%     3         2      5      7
%     4         2     -1      8
%     5         3      8      6
%     6         5      9      7
%     7         3      6     -1
%     8         5      4     10
%     9         6     10     12
%    10         9      8     11
%    11        12     10     14
%    12         9     11     13
%    13        -1     12     16
%    14        11     -1     15
%    15        16     14     -1
%    16        13     15     -1
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 September 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer TRIANGLE_ORDER, the order of the triangles.
%
%    Input, integer TRIANGLE_NUM, the number of triangles.
%
%    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), the nodes that 
%    make up each triangle.
%
%    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the three triangles that are direct
%    neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I) is the index of the triangle
%    which touches side 1, defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I)
%    is negative if there is no neighbor on that side.  In this case, that
%    side of the triangle lies on the boundary of the triangulation.
%

%
%  Step 1.
%  From the list of nodes for triangle T, of the form: (I,J,K)
%  construct the three neighbor relations:
%
%    (I,J,3,T) or (J,I,3,T),
%    (J,K,1,T) or (K,J,1,T),
%    (K,I,2,T) or (I,K,2,T)
%
%  where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
%
  col = zeros ( 4, triangle_order * triangle_num );
  
  for tri = 1 : triangle_num

    i = triangle_node(1,tri);
    j = triangle_node(2,tri);
    k = triangle_node(3,tri);

    if ( i < j )
      col(1:4,1+3*(tri-1)) = [ i, j, 3, tri ]';
    else
      col(1:4,1+3*(tri-1)) = [ j, i, 3, tri ]';
    end

    if ( j < k )
      col(1:4,2+3*(tri-1)) = [ j, k, 1, tri ]';
    else
      col(1:4,2+3*(tri-1)) = [ k, j, 1, tri ]';
    end

    if ( k < i )
      col(1:4,3+3*(tri-1)) = [ k, i, 2, tri ]';
    else
      col(1:4,3+3*(tri-1)) = [ i, k, 2, tri ]';
    end

  end
%
%  Step 2. Perform an ascending dictionary sort on the neighbor relations.
%  We only intend to sort on rows 1 and 2; the routine we call here
%  sorts on rows 1 through 4 but that won't hurt us.
%
%  What we need is to find cases where two triangles share an edge.
%  Say they share an edge defined by the nodes I and J.  Then there are
%  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
%  we make sure that these two columns occur consecutively.  That will
%  make it easy to notice that the triangles are neighbors.
%
  col = i4col_sort_a ( 4, 3*triangle_num, col );
%
%  Step 3. Neighboring triangles show up as consecutive columns with
%  identical first two entries.  Whenever you spot this happening,
%  make the appropriate entries in TRIANGLE_NEIGHBOR.
%
  triangle_neighbor(1:3,1:triangle_num) = -1;

  icol = 1;

  while ( 1 )

    if ( 3 * triangle_num <= icol )
      break
    end

    if ( col(1,icol) ~= col(1,icol+1) || col(2,icol) ~= col(2,icol+1) )
      icol = icol + 1;
      continue
    end

    side1 = col(3,icol);
    tri1 =  col(4,icol);
    side2 = col(3,icol+1);
    tri2 =  col(4,icol+1);

    triangle_neighbor(side1,tri1) = tri2;
    triangle_neighbor(side2,tri2) = tri1;

    icol = icol + 2;

  end

  return
end
function [ adj_num, adj_col ] = triangulation_order3_adj_count ( node_num, ...
  tri_num, triangle_node, triangle_neighbor )

%*****************************************************************************80
%
%% TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
%
%  Discussion:
%
%    This routine is called to count the adjacencies, so that the
%    appropriate amount of memory can be set aside for storage when
%    the adjacency structure is created.
%
%    The triangulation is assumed to involve 3-node triangles.
%
%    Two nodes are "adjacent" if they are both nodes in some triangle.
%    Also, a node is considered to be adjacent to itself.
%
%  Diagram:
%
%       3
%    s  |\
%    i  | \
%    d  |  \
%    e  |   \  side 2
%       |    \
%    3  |     \
%       |      \
%       1-------2
%
%         side 1
%
%    The local node numbering
%
%
%   21-22-23-24-25
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%   16-17-18-19-20
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%   11-12-13-14-15
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%    6--7--8--9-10
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%    1--2--3--4--5
%
%    A sample grid.
%
%
%    Below, we have a chart that summarizes the adjacency relationships
%    in the sample grid.  On the left, we list the node, and its neighbors,
%    with an asterisk to indicate the adjacency of the node to itself
%    (in some cases, you want to count this self adjacency and in some
%    you don't).  On the right, we list the number of adjancencies to
%    lower-indexed nodes, to the node itself, to higher-indexed nodes,
%    the total number of adjacencies for this node, and the location
%    of the first and last entries required to list this set of adjacencies
%    in a single list of all the adjacencies.
%
%    N   Adjacencies                Below  Self   Above   Total First  Last
%
%   --  -- -- -- -- -- -- --           --    --      --      --   ---     0   
%    1:  *  2  6                        0     1       2       3     1     3
%    2:  1  *  3  6  7                  1     1       3       5     4     8
%    3:  2  *  4  7  8                  1     1       3       5     9    13
%    4:  3  *  5  8  9                  1     1       3       5    14    18
%    5:  4  *  9 10                     1     1       2       4    19    22
%    6:  1  2  *  7 11                  2     1       2       5    23    27
%    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
%    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
%    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
%   10:  5  9  * 14 15                  2     1       2       5    49    53
%   11:  6  7  * 12 16                  2     1       2       5    54    58
%   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
%   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
%   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
%   15: 10 14  * 19 20                  2     1       2       5    80    84
%   16: 11 12  * 17 21                  2     1       2       5    85    89
%   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
%   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
%   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
%   20: 15 19  * 24 25                  2     1       2       5   111   115
%   21: 16 17  * 22                     2     1       1       4   116   119
%   22: 17 18 21  * 23                  3     1       1       5   120   124
%   23: 18 19 22  * 24                  3     1       1       5   125   129
%   24: 19 20 23  * 25                  3     1       1       5   130   134
%   25: 20 24  *                        2     1       0       3   135   137
%   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 August 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer TRI_NUM, the number of triangles.
%
%    Input, integer TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
%    make up each triangle, in counterclockwise order. 
%
%    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
%    a triangle, lists the neighboring triangle, or -1 if there is
%    no neighbor.
%
%    Output, integer ADJ_NUM, the number of adjacencies.
%
%    Output, integer ADJ_COL(NODE_NUM+1).  Information about column J is stored
%    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
%

%
%  Set every node to be adjacent to itself.
%
  adj_col(1:node_num) = 1;
%
%  Examine each triangle.
%
  for triangle = 1 : tri_num

    n1 = triangle_node(1,triangle);
    n2 = triangle_node(2,triangle);
    n3 = triangle_node(3,triangle);
%
%  Add edge (1,2) if this is the first occurrence,
%  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
%  or if this triangle is the first of the pair in which the edge
%  occurs (TRIANGLE < TRIANGLE2).
%
    triangle2 = triangle_neighbor(1,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj_col(n1) = adj_col(n1) + 1;
      adj_col(n2) = adj_col(n2) + 1;
    end
%
%  Add edge (2,3).
%
    triangle2 = triangle_neighbor(2,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj_col(n2) = adj_col(n2) + 1;
      adj_col(n3) = adj_col(n3) + 1;
    end
%
%  Add edge (3,1).
%
    triangle2 = triangle_neighbor(3,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj_col(n1) = adj_col(n1) + 1;
      adj_col(n3) = adj_col(n3) + 1;
    end
      
  end
%
%  We used ADJ_COL to count the number of entries in each column.
%  Convert it to pointers into the ADJ array.
%
  adj_col(2:node_num+1) = adj_col(1:node_num);

  adj_col(1) = 1;
  for i = 2 : node_num+1
    adj_col(i) = adj_col(i-1) + adj_col(i);
  end

  adj_num = adj_col(node_num+1) - 1;

  return
end
function adj = triangulation_order3_adj_set ( node_num, ...
  tri_num, triangle_node, triangle_neighbor, adj_num, adj_col )

%*****************************************************************************80
%
%% TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
%
%  Discussion:
%
%    This routine is called to count the adjacencies, so that the
%    appropriate amount of memory can be set aside for storage when
%    the adjacency structure is created.
%
%    The triangulation is assumed to involve 3-node triangles.
%
%    Two nodes are "adjacent" if they are both nodes in some triangle.
%    Also, a node is considered to be adjacent to itself.
%
%    This routine can be used to create the compressed column storage
%    for a linear triangle finite element discretization of 
%    Poisson's equation in two dimensions.
%
%  Diagram:
%
%       3
%    s  |\
%    i  | \
%    d  |  \
%    e  |   \  side 2
%       |    \
%    3  |     \
%       |      \
%       1-------2
%
%         side 1
%
%    The local node numbering
%
%
%   21-22-23-24-25
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%   16-17-18-19-20
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%   11-12-13-14-15
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%    6--7--8--9-10
%    |\ |\ |\ |\ |
%    | \| \| \| \|
%    1--2--3--4--5
%
%    A sample grid
%
%
%    Below, we have a chart that summarizes the adjacency relationships
%    in the sample grid.  On the left, we list the node, and its neighbors,
%    with an asterisk to indicate the adjacency of the node to itself
%    (in some cases, you want to count this self adjacency and in some
%    you don't).  On the right, we list the number of adjancencies to
%    lower-indexed nodes, to the node itself, to higher-indexed nodes,
%    the total number of adjacencies for this node, and the location
%    of the first and last entries required to list this set of adjacencies
%    in a single list of all the adjacencies.
%
%    N   Adjacencies                Below  Self    Above  Total First  Last
%
%   --  -- -- -- -- -- -- --           --    --      --      --   ---     0   
%    1:  *  2  6                        0     1       2       3     1     3
%    2:  1  *  3  6  7                  1     1       3       5     4     8
%    3:  2  *  4  7  8                  1     1       3       5     9    13
%    4:  3  *  5  8  9                  1     1       3       5    14    18
%    5:  4  *  9 10                     1     1       2       4    19    22
%    6:  1  2  *  7 11                  2     1       2       5    23    27
%    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
%    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
%    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
%   10:  5  9  * 14 15                  2     1       2       5    49    53
%   11:  6  7  * 12 16                  2     1       2       5    54    58
%   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
%   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
%   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
%   15: 10 14  * 19 20                  2     1       2       5    80    84
%   16: 11 12  * 17 21                  2     1       2       5    85    89
%   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
%   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
%   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
%   20: 15 19  * 24 25                  2     1       2       5   111   115
%   21: 16 17  * 22                     2     1       1       4   116   119
%   22: 17 18 21  * 23                  3     1       1       5   120   124
%   23: 18 19 22  * 24                  3     1       1       5   125   129
%   24: 19 20 23  * 25                  3     1       1       5   130   134
%   25: 20 24  *                        2     1       0       3   135   137
%   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 September 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer TRI_NUM, the number of triangles.
%
%    Input, integer TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
%    make up each triangle in counterclockwise order.
%
%    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
%    a triangle, lists the neighboring triangle, or -1 if there is
%    no neighbor.
%
%    Input, integer ADJ_NUM, the number of adjacencies.
%
%    Input, integer ADJ_COL(NODE_NUM+1).  Information about column J is stored
%    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
%
%    Output, integer ADJ(ADJ_NUM), the adjacency information.
%
  adj(1:adj_num) = -1;
  adj_copy(1:node_num) = adj_col(1:node_num);
%
%  Set every node to be adjacent to itself.
%
  for node = 1 : node_num
    adj(adj_copy(node)) = node;
    adj_copy(node) = adj_copy(node) + 1;
  end
%
%  Examine each triangle.
%
  for triangle = 1 : tri_num

    n1 = triangle_node(1,triangle);
    n2 = triangle_node(2,triangle);
    n3 = triangle_node(3,triangle);
%
%  Add edge (1,2) if this is the first occurrence,
%  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
%  or if this triangle is the first of the pair in which the edge
%  occurs (TRIANGLE < TRIANGLE2).
%
    triangle2 = triangle_neighbor(1,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj(adj_copy(n1)) = n2;
      adj_copy(n1) = adj_copy(n1) + 1;
      adj(adj_copy(n2)) = n1;
      adj_copy(n2) = adj_copy(n2) + 1;
    end
%
%  Add edge (2,3).
%
    triangle2 = triangle_neighbor(2,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj(adj_copy(n2)) = n3;
      adj_copy(n2) = adj_copy(n2) + 1;
      adj(adj_copy(n3)) = n2;
      adj_copy(n3) = adj_copy(n3) + 1;
    end
%
%  Add edge (3,1).
%
    triangle2 = triangle_neighbor(3,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj(adj_copy(n1)) = n3;
      adj_copy(n1) = adj_copy(n1) + 1;
      adj(adj_copy(n3)) = n1;
      adj_copy(n3) = adj_copy(n3) + 1;
    end
      
  end
%
%  Ascending sort the entries for each node.
%
  for node = 1 : node_num
    k1 = adj_col(node);
    k2 = adj_col(node+1)-1;
    adj(k1:k2) = i4vec_sort_heap_a ( k2+1-k1, adj(k1:k2) );
  end

  return
end
function [ adj_num, adj_col ] = triangulation_order6_adj_count ( node_num, ...
  tri_num, triangle_node, triangle_neighbor )

%*****************************************************************************80
%
%% TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
%
%  Discussion:
%
%    This routine is called to count the adjacencies, so that the
%    appropriate amount of memory can be set aside for storage when
%    the adjacency structure is created.
%
%    The triangulation is assumed to involve 6-node triangles.
%
%    Two nodes are "adjacent" if they are both nodes in some triangle.
%    Also, a node is considered to be adjacent to itself.
%
%  Diagram:
%
%       3
%    s  |\
%    i  | \
%    d  |  \
%    e  6   5  side 2
%       |    \
%    3  |     \
%       |      \
%       1---4---2
%
%         side 1
%
%    The local node numbering
%
%
%   21-22-23-24-25
%    |\    |\    |
%    | \   | \   |
%   16 17 18 19 20
%    |   \ |   \ |
%    |    \|    \|
%   11-12-13-14-15
%    |\    |\    |
%    | \   | \   |
%    6  7  8  9 10
%    |   \ |   \ |
%    |    \|    \|
%    1--2--3--4--5
%
%    A sample grid.
%
%
%    Below, we have a chart that lists the nodes adjacent to each node, with 
%    an asterisk to indicate the adjacency of the node to itself
%    (in some cases, you want to count this self adjacency and in some
%    you don't).
%
%    N   Adjacencies
%
%    1:  *  2  3  6  7 11
%    2:  1  *  3  6  7 11
%    3:  1  2  *  4  5  6  7  8  9 11 12 13
%    4:  3  *  5  8  9 13
%    5:  3  4  *  8  9 10 13 14 15
%    6:  1  2  3  *  7 11
%    7:  1  2  3  6  *  8 11 12 13
%    8:  3  4  5  7  *  9 11 12 13
%    9:  3  4  5  8  * 10 13 14 15
%   10:  5  9  * 13 14 15
%   11:  1  2  3  6  7  8  * 12 13 16 17 21
%   12:  3  7  8 11  * 13 16 17 21
%   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
%   14:  5  9 10 13  * 15 18 19 23
%   15:  5  9 10 13 14  * 18 19 20 23 24 25
%   16: 11 12 13  * 17 21
%   17: 11 12 13 16  * 18 21 22 23
%   18: 13 14 15 17  * 19 21 22 23
%   19: 13 14 15 18  * 20 23 24 25
%   20: 15 19  * 23 24 25
%   21: 11 12 13 16 17 18  * 22 23
%   22: 13 17 18 21  * 23
%   23: 13 14 15 17 18 19 20 21 22  * 24 25
%   24: 15 19 20 23  * 25
%   25: 15 19 20 23 24  *    
%
%    Below, we list the number of adjancencies to lower-indexed nodes, to 
%    the node itself, to higher-indexed nodes, the total number of 
%    adjacencies for this node, and the location of the first and last 
%    entries required to list this set of adjacencies in a single list 
%    of all the adjacencies.
%
%    N   Below  Self   Above   Total First  Last
%
%   --      --    --      --      --   ---     0
%    1:      0     1       5       6     1     6
%    2:      1     1       4       6     7    12
%    3:      2     1       9      12    13    24
%    4:      1     1       4       6    25    30
%    5:      2     1       6       9    31    39
%    6:      3     1       2       6    40    45
%    7:      4     1       4       9    46    54
%    8:      4     1       4       9    55    63
%    9:      4     1       4       9    62    72
%   10:      2     1       3       6    73    78
%   11:      6     1       5      12    79    90
%   12:      4     1       4       9    91    99
%   13:      9     1       9      19   100   118
%   14:      4     1       4       9   119   127
%   15:      5     1       6      12   128   139
%   16:      3     1       2       6   140   145
%   17:      4     1       4       9   146   154
%   18:      4     1       4       9   155   163
%   19:      4     1       4       9   164   172
%   20:      2     1       3       6   173   178
%   21:      6     1       2       9   179   187
%   22:      4     1       1       6   188   193
%   23:      9     1       2      12   194   205
%   24:      4     1       1       6   206   211
%   25:      5     1       0       6   212   217
%   --      --    --      --      --   218   ---
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 August 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer TRI_NUM, the number of triangles.
%
%    Input, integer TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
%    make up each triangle.  The first three nodes are the vertices,
%    in counterclockwise order.  The fourth value is the midside
%    node between nodes 1 and 2; the fifth and sixth values are
%    the other midside nodes in the logical order.
%
%    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
%    a triangle, lists the neighboring triangle, or -1 if there is
%    no neighbor.
%
%    Output, integer ADJ_NUM, the number of adjacencies.
%
%    Output, integer ADJ_COL(NODE_NUM+1).  Information about column J is stored
%    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
%

%
%  Set every node to be adjacent to itself.
%
  adj_col(1:node_num) = 1;
%
%  Examine each triangle.
%
  for triangle = 1 : tri_num

    n1 = triangle_node(1,triangle);
    n2 = triangle_node(2,triangle);
    n3 = triangle_node(3,triangle);
    n4 = triangle_node(4,triangle);
    n5 = triangle_node(5,triangle);
    n6 = triangle_node(6,triangle);
%
%  For sure, we add the adjacencies:
%    43 / (34)
%    51 / (15)
%    54 / (45)
%    62 / (26)
%    64 / (46)
%    65 / (56)
%
    adj_col(n3) = adj_col(n3) + 1;
    adj_col(n4) = adj_col(n4) + 1;
    adj_col(n1) = adj_col(n1) + 1;
    adj_col(n5) = adj_col(n5) + 1;
    adj_col(n4) = adj_col(n4) + 1;
    adj_col(n5) = adj_col(n5) + 1;
    adj_col(n2) = adj_col(n2) + 1;
    adj_col(n6) = adj_col(n6) + 1;
    adj_col(n4) = adj_col(n4) + 1;
    adj_col(n6) = adj_col(n6) + 1;
    adj_col(n5) = adj_col(n5) + 1;
    adj_col(n6) = adj_col(n6) + 1;
%
%  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
%  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
%  or if this triangle is the first of the pair in which the edge
%  occurs (TRIANGLE < TRIANGLE2).
%
%  Maybe add
%    21 / 12
%    41 / 14
%    42 / 24
%
    triangle2 = triangle_neighbor(1,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj_col(n1) = adj_col(n1) + 1;
      adj_col(n2) = adj_col(n2) + 1;
      adj_col(n1) = adj_col(n1) + 1;
      adj_col(n4) = adj_col(n4) + 1;
      adj_col(n2) = adj_col(n2) + 1;
      adj_col(n4) = adj_col(n4) + 1;
    end
%
%  Maybe add
%    32 / 23
%    52 / 25
%    53 / 35
%
    triangle2 = triangle_neighbor(2,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj_col(n2) = adj_col(n2) + 1;
      adj_col(n3) = adj_col(n3) + 1;
      adj_col(n2) = adj_col(n2) + 1;
      adj_col(n5) = adj_col(n5) + 1;
      adj_col(n3) = adj_col(n3) + 1;
      adj_col(n5) = adj_col(n5) + 1;
    end
%
%  Maybe add
%    31 / 13
%    61 / 16
%    63 / 36
%
    triangle2 = triangle_neighbor(3,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj_col(n1) = adj_col(n1) + 1;
      adj_col(n3) = adj_col(n3) + 1;
      adj_col(n1) = adj_col(n1) + 1;
      adj_col(n6) = adj_col(n6) + 1;
      adj_col(n3) = adj_col(n3) + 1;
      adj_col(n6) = adj_col(n6) + 1;
    end
      
  end
%
%  We used ADJ_COL to count the number of entries in each column.
%  Convert it to pointers into the ADJ array.
%
  adj_col(2:node_num+1) = adj_col(1:node_num);

  adj_col(1) = 1;
  for i = 2 : node_num+1
    adj_col(i) = adj_col(i-1) + adj_col(i);
  end

  adj_num = adj_col(node_num+1) - 1;

  return
end
function adj = triangulation_order6_adj_set ( node_num, tri_num, ...
  triangle_node, triangle_neighbor, adj_num, adj_col )

%*****************************************************************************80
%
%% TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
%
%  Discussion:
%
%    This routine is called to count the adjacencies, so that the
%    appropriate amount of memory can be set aside for storage when
%    the adjacency structure is created.
%
%    The triangulation is assumed to involve 6-node triangles.
%
%    Two nodes are "adjacent" if they are both nodes in some triangle.
%    Also, a node is considered to be adjacent to itself.
%
%    This routine can be used to create the compressed column storage
%    for a quadratic triangle finite element discretization of 
%    Poisson's equation in two dimensions.
%
%  Diagram:
%
%       3
%    s  |\
%    i  | \
%    d  |  \
%    e  6   5  side 2
%       |    \
%    3  |     \
%       |      \
%       1---4---2
%
%         side 1
%
%    The local node numbering
%
%
%   21-22-23-24-25
%    |\    |\    |
%    | \   | \   |
%   16 17 18 19 20
%    |   \ |   \ |
%    |    \|    \|
%   11-12-13-14-15
%    |\    |\    |
%    | \   | \   |
%    6  7  8  9 10
%    |   \ |   \ |
%    |    \|    \|
%    1--2--3--4--5
%
%    A sample grid.
%
%
%    Below, we have a chart that lists the nodes adjacent to each node, with 
%    an asterisk to indicate the adjacency of the node to itself
%    (in some cases, you want to count this self adjacency and in some
%    you don't).
%
%    N   Adjacencies
%
%    1:  *  2  3  6  7 11
%    2:  1  *  3  6  7 11
%    3:  1  2  *  4  5  6  7  8  9 11 12 13
%    4:  3  *  5  8  9 13
%    5:  3  4  *  8  9 10 13 14 15
%    6:  1  2  3  *  7 11
%    7:  1  2  3  6  *  8 11 12 13
%    8:  3  4  5  7  *  9 11 12 13
%    9:  3  4  5  8  * 10 13 14 15
%   10:  5  9  * 13 14 15
%   11:  1  2  3  6  7  8  * 12 13 16 17 21
%   12:  3  7  8 11  * 13 16 17 21
%   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
%   14:  5  9 10 13  * 15 18 19 23
%   15:  5  9 10 13 14  * 18 19 20 23 24 25
%   16: 11 12 13  * 17 21
%   17: 11 12 13 16  * 18 21 22 23
%   18: 13 14 15 17  * 19 21 22 23
%   19: 13 14 15 18  * 20 23 24 25
%   20: 15 19  * 23 24 25
%   21: 11 12 13 16 17 18  * 22 23
%   22: 13 17 18 21  * 23
%   23: 13 14 15 17 18 19 20 21 22  * 24 25
%   24: 15 19 20 23  * 25
%   25: 15 19 20 23 24  *    
%
%    Below, we list the number of adjancencies to lower-indexed nodes, to 
%    the node itself, to higher-indexed nodes, the total number of 
%    adjacencies for this node, and the location of the first and last 
%    entries required to list this set of adjacencies in a single list 
%    of all the adjacencies.
%
%    N   Below  Self   Above   Total First  Last
%
%   --      --    --      --      --   ---     0
%    1:      0     1       5       6     1     6
%    2:      1     1       4       6     7    12
%    3:      2     1       9      12    13    24
%    4:      1     1       4       6    25    30
%    5:      2     1       6       9    31    39
%    6:      3     1       2       6    40    45
%    7:      4     1       4       9    46    54
%    8:      4     1       4       9    55    63
%    9:      4     1       4       9    62    72
%   10:      2     1       3       6    73    78
%   11:      6     1       5      12    79    90
%   12:      4     1       4       9    91    99
%   13:      9     1       9      19   100   118
%   14:      4     1       4       9   119   127
%   15:      5     1       6      12   128   139
%   16:      3     1       2       6   140   145
%   17:      4     1       4       9   146   154
%   18:      4     1       4       9   155   163
%   19:      4     1       4       9   164   172
%   20:      2     1       3       6   173   178
%   21:      6     1       2       9   179   187
%   22:      4     1       1       6   188   193
%   23:      9     1       2      12   194   205
%   24:      4     1       1       6   206   211
%   25:      5     1       0       6   212   217
%   --      --    --      --      --   218   ---
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 August 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer TRI_NUM, the number of triangles.
%
%    Input, integer TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
%    make up each triangle.  The first three nodes are the vertices,
%    in counterclockwise order.  The fourth value is the midside
%    node between nodes 1 and 2; the fifth and sixth values are
%    the other midside nodes in the logical order.
%
%    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
%    a triangle, lists the neighboring triangle, or -1 if there is
%    no neighbor.
%
%    Input, integer ADJ_NUM, the number of adjacencies.
%
%    Input, integer ADJ_COL(NODE_NUM+1).  Information about column J is stored
%    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
%
%    Output, integer ADJ(ADJ_NUM), the adjacency information.
%
  adj(1:adj_num) = -1;
  adj_copy(1:node_num) = adj_col(1:node_num);
%
%  Set every node to be adjacent to itself.
%
  for node = 1 : node_num
    adj(adj_copy(node)) = node;
    adj_copy(node) = adj_copy(node) + 1;
  end
%
%  Examine each triangle.
%
  for triangle = 1 : tri_num

    n1 = triangle_node(1,triangle);
    n2 = triangle_node(2,triangle);
    n3 = triangle_node(3,triangle);
    n4 = triangle_node(4,triangle);
    n5 = triangle_node(5,triangle);
    n6 = triangle_node(6,triangle);
%
%  For sure, we add the adjacencies:
%    43 / (34)
%    51 / (15)
%    54 / (45)
%    62 / (26)
%    64 / (46)
%    65 / (56)
%
    adj(adj_copy(n3)) = n4;
    adj_copy(n3) = adj_copy(n3) + 1;
    adj(adj_copy(n4)) = n3;
    adj_copy(n4) = adj_copy(n4) + 1;

    adj(adj_copy(n1)) = n5;
    adj_copy(n1) = adj_copy(n1) + 1;
    adj(adj_copy(n5)) = n1;
    adj_copy(n5) = adj_copy(n5) + 1;

    adj(adj_copy(n4)) = n5;
    adj_copy(n4) = adj_copy(n4) + 1;
    adj(adj_copy(n5)) = n4;
    adj_copy(n5) = adj_copy(n5) + 1;

    adj(adj_copy(n2)) = n6;
    adj_copy(n2) = adj_copy(n2) + 1;
    adj(adj_copy(n6)) = n2;
    adj_copy(n6) = adj_copy(n6) + 1;

    adj(adj_copy(n4)) = n6;
    adj_copy(n4) = adj_copy(n4) + 1;
    adj(adj_copy(n6)) = n4;
    adj_copy(n6) = adj_copy(n6) + 1;

    adj(adj_copy(n5)) = n6;
    adj_copy(n5) = adj_copy(n5) + 1;
    adj(adj_copy(n6)) = n5;
    adj_copy(n6) = adj_copy(n6) + 1;
%
%  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
%  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
%  or if this triangle is the first of the pair in which the edge
%  occurs (TRIANGLE < TRIANGLE2).
%
%  Maybe add
%    21 / 12
%    41 / 14
%    42 / 24
%
    triangle2 = triangle_neighbor(1,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj(adj_copy(n1)) = n2;
      adj_copy(n1) = adj_copy(n1) + 1;
      adj(adj_copy(n2)) = n1;
      adj_copy(n2) = adj_copy(n2) + 1;
      adj(adj_copy(n1)) = n4;
      adj_copy(n1) = adj_copy(n1) + 1;
      adj(adj_copy(n4)) = n1;
      adj_copy(n4) = adj_copy(n4) + 1;
      adj(adj_copy(n2)) = n4;
      adj_copy(n2) = adj_copy(n2) + 1;
      adj(adj_copy(n4)) = n2;
      adj_copy(n4) = adj_copy(n4) + 1;
    end
%
%  Maybe add
%    32 / 23
%    52 / 25
%    53 / 35
%
    triangle2 = triangle_neighbor(2,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj(adj_copy(n2)) = n3;
      adj_copy(n2) = adj_copy(n2) + 1;
      adj(adj_copy(n3)) = n2;
      adj_copy(n3) = adj_copy(n3) + 1;
      adj(adj_copy(n2)) = n5;
      adj_copy(n2) = adj_copy(n2) + 1;
      adj(adj_copy(n5)) = n2;
      adj_copy(n5) = adj_copy(n5) + 1;
      adj(adj_copy(n3)) = n5;
      adj_copy(n3) = adj_copy(n3) + 1;
      adj(adj_copy(n5)) = n3;
      adj_copy(n5) = adj_copy(n5) + 1;
    end
%
%  Maybe add
%    31 / 13
%    61 / 16
%    63 / 36
%
    triangle2 = triangle_neighbor(3,triangle);

    if ( triangle2 < 0 || triangle < triangle2 )
      adj(adj_copy(n1)) = n3;
      adj_copy(n1) = adj_copy(n1) + 1;
      adj(adj_copy(n3)) = n1;
      adj_copy(n3) = adj_copy(n3) + 1;
      adj(adj_copy(n1)) = n6;
      adj_copy(n1) = adj_copy(n1) + 1;
      adj(adj_copy(n6)) = n1;
      adj_copy(n6) = adj_copy(n6) + 1;
      adj(adj_copy(n3)) = n6;
      adj_copy(n3) = adj_copy(n3) + 1;
      adj(adj_copy(n6)) = n3;
      adj_copy(n6) = adj_copy(n6) + 1;
    end
      
  end
%
%  Ascending sort the entries for each node.
%
  for node = 1 : node_num
    k1 = adj_col(node);
    k2 = adj_col(node+1)-1;
    adj(k1:k2) = i4vec_sort_heap_a ( k2+1-k1, adj(k1:k2) );
  end

  return
end
