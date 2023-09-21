% Grid class - 3D Cartesian (Anisotropic) Grids 
% G = Grid(NX,NY,NZ,dx,dy,dz) generates a NX x NY x NZ rectangular cuboid
%   finescale width dx, dy, and dz in the x, y, and z dimensions.
% 
% INPUT:
%       nx, ny, nz = number of cells(grid blocks) by dimension 
%       dx, dy, dz = size of each cell, for each dimension dx = x2-x1 
% OUTPUT:
%       G = Grid Class Object with the properties:
%
% PROPERTIES
%   .ncells:   total number of cells (including corner & ghost cells)
%   .nicells:  total number of interior cells
%   .nfaces:   total number of faces
%   .cells:    array of cell centers, dimensions, etc.
%   .icells: row indices corresponding to interior cells
%   .faces: nfaces-by-5-matrix
% 
% G.cells: ncells-by-16 matrix
%       1: x-coordinate of cell center
%       2: y-coordinate of cell center
%       3: z-coordinate of cell center
%       ------------------------------------------------ 
%       4: size of cell, x-dimension
%       5: size of cell, y-dimension
%       6: size of cell, z-dimension
%       ------------------------------------------------ 
%       7: refinement level in x
%       8: refinement level in y
%       9: refinement level in z
%       ------------------------------------------------ 
%       Lower-Left Coarse cell index (i,j,k)
%       10: coarse cell x-index
%       11: coarse cell y-index
%       12: coarse cell z-index
%       ------------------------------------------------ 
%       13: boundary cell bitmask (interior cells = 0)
%           bits: 
%               1=left boundary, 2=right boundary, = x 
%               3=bottom boundary, 4=top boundary, =  y
%               5=down boundary, 6=up boundary = z
%               e.g., 20 = 010100 = bottom/down boundary
%               bitmask used since each cell can be multiple
%               boundaries
%       14: index of interior cell neighbor, if bdry cell, but not edge or corner (0 otherwise, ~0 = legit bdry cell)
%       15: row index of interior cell, else = 0 (i.e., 0 if not interior)
%       ------------------------------------------------ 
%       16: pressure of cell i as flow driven in x-direction
%       17: pressure of cell i .............. in y-direction
%       18: pressure of cell i .............. in y-direction
%       ------------------------------------------------ 
%       19: ID (unique referential id)      
%
% G.faces: nfaces-by-5-matrix
%   columns:
%       1: index of cell on left/bottom/down side of face
%       2: index of cell on right/top/up side of face
%       3: orthogonal to which axis 
%          =1  x-axis, left/right side
%          =2  y-axis, bottom/top side 
%          =3  z-axis, down/up side
%       4: FaceID (row number)
%       5: Boundary Cell (>0).  Boundary Face Indicator 
%
% Example:  G = Grid(4,4,4,1,1,1)
%
% See also

%{
-----------------------------------------------------------------------
Copyright (C) 2021 J. Quinlan

This file is part of Variable Compact Multipoint Method (VCMP).

VCMP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VCMP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VCMP.  If not, see <http://www.gnu.org/licenses/>.

-----------------------------------------------------------------------
MATLAB Version:        MATLAB Version: 9.10.0.1602886 (R2021a)
Author:                James Quinlan                                  
Modified:              10/18/2021
 
Real-cell detector: sum(G.cells(:,1:3)==0,2)<2
----------------------------------------------------------------------- 
%}

classdef Grid < handle

    properties
        cells         % Cells (array)
        ncells        % Total Number of Cells
        icells        % Interior Cell indices
        nicells       % Number of Interior Cells
        fine=[1,1,1]  % Fine-scale Size
        domain        % Grid Domain  
        faces         % Faces (array)
        nfaces        % Number of Faces
        nbfaces       % Number of Bdry Faces
        ghosts        % Boundary cells ordered by cells (no corners or edge cells)
        bcells        % Boundary cells ordered by face
        
        coarse        % Coarse Grid Size

        % Prop5  double {mustBePositive} = 15
    end
    
    
    methods
        % Constructor   
        % ---------------------- 
        function G = Grid(nx,ny,nz,dx,dy,dz)     
            if nargin<4                          
                dx=1;dy=1;dz=1; % Default Grid Width
            end
            
            % Coarse grid size (if fine, then coarse = fine)
            G.coarse = [dx dy dz];
                          
            % Initialize Cells
            G.ncells = (nx+2)*(ny+2)*(nz+2);    % Total cells (w/ghosts)
            G.cells = zeros(G.ncells,19);       % Grid Cell array
            G.nicells = nx*ny*nz;               % Num. interior cells
            G.icells = zeros(G.nicells,1);      % Interior Cell array
            

            % Initialize Faces
            G.nfaces = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1);
            G.faces = zeros(G.nfaces,5);

            % Unit vectors,[1,...,1]'   
            ex = ones(nx,1);
            ey = ones(ny,1);
            ez = ones(nz,1);

            % Add boundary positions
            ex2 = [ 1; ex; 1 ];
            ey2 = [ 1; ey; 1 ];
            ez2 = [ 1; ez; 1 ];

            % Grid index by dimension
            i1 = (1:nx+1)';               
            j1 = (1:ny+1)';               
            k1 = (1:nz+1)';               

            i = i1(1:nx);  % nx*dx
            j = j1(1:ny);
            k = k1(1:nz);

            ix = [ 0; i-1/2; nx ];      % centers (by dimension) w/brdy coordinates
            iy = [ 0; j-1/2; ny ];
            iz = [ 0; k-1/2; nz ];
                        
            %% CELL defintions
            
            % ----------------------- 1:3 cell centers ------------------ %
            % Cell centroids (x,y,z), including ghost cells 
            G.cells(:,1) = kron(ez2,kron(ey2,dx*ix)); % = x
            G.cells(:,2) = kron(ez2,kron(dy*iy,ex2)); % = y
            G.cells(:,3) = kron(dz*iz,kron(ey2,ex2)); % = z

             
            % ------------------ 4:6 cell dims -------------------------- %
            % cell dims (ordered triple cells(:,4:6) = dimensions of ith cell/row)
            G.cells(:,4) = kron(ez2,kron(ey2,dx*[ 0; ex; 0 ]));   % x-width
            G.cells(:,5) = kron(ez2,kron(dy*[ 0; ey; 0 ],ex2));   % y-width
            G.cells(:,6) = kron(dz*[ 0; ez; 0 ],kron(ey2,ex2));   % z-width

            % ------------------ 7:9 refinement ------------------------- %
            % Cell refinement level(s) 
            % Number of refinements - backwards?
            % 0 = no refinement
            % 1 = one refinement 
            % 2 = two refinements
            G.cells(:,7) = log2(nx); % e.g., nx=256, log2(nx)=8
            G.cells(:,8) = log2(ny);
            G.cells(:,9) = log2(nz);

            % ----------------- 10:12 coarse cell indices --------------- %
            G.cells(:,10) = kron(ez2,kron(ey2,[ 0; i; nx+1 ]));
            G.cells(:,11) = kron(ez2,kron([ 0; j; ny+1 ],ex2));
            G.cells(:,12) = kron([ 0; k; nz+1 ],kron(ey2,ex2));

            % ----------------- 13 boundary indicators ------------------ %  
            % x      
            G.cells(:,13) = bitset(G.cells(:,13),1,log2(G.cells(:,10))<0);
            G.cells(:,13) = bitset(G.cells(:,13),2,log2(G.cells(:,10))>G.cells(:,7));

            % y
            G.cells(:,13) = bitset(G.cells(:,13),3,log2(G.cells(:,11))<0);
            G.cells(:,13) = bitset(G.cells(:,13),4,log2(G.cells(:,11))>G.cells(:,8));

            % z
            G.cells(:,13) = bitset(G.cells(:,13),5,log2(G.cells(:,12))<0);
            G.cells(:,13) = bitset(G.cells(:,13),6,log2(G.cells(:,12))>G.cells(:,9));

            % ----------------------------- 14:16-18 ----------------------------
            % 14: index of interior cell neighbor, if bdry cell (0 otherwise)
            % 15: row index of interior cell, 0 if boundary cell
            % 16-18: pressure of ith cell.  Used for Local-Global iteration
            % -------------------------------------------------------------------
            G.cells(:,16:18) = 0.5;  % set to 0, then set all fake cells NaN
            % 0.5 indicates interior cells
            % NaN = Fake cells (corners or edges)
            % -1 = No-flow boundary
            % 1 / 0 = Dirichlet boundary pressure value

            % -------------------- 19 cell number (row index) sequential --------  
            G.cells(:,19) = 1:G.ncells;
            
            
            %% FACES: IDENTIFY & INDEX

            % -------- x-faces ---------- %
            base=(nx+1)*ny*nz;

            I = kron(i1,ones(ny*nz,1));
            J = kron(ones(nx+1,1),kron(j,ones(nz,1))); % J=kron(ones(base/ny,1),kron(j,ones(ny,1)));
            k = (1:nz)'; K=kron(ones((nx+1)*ny,1),k);
            fid = (nx+1)*ny*(K-1)+(nx+1)*(J-1)+I;
            idx1 = (nx+2)*(ny+2)*K+(nx+2)*J+I;
            idx2 = K*(nx+2)*(ny+2)+J*(nx+2)+(I+1);

            xfaces(fid,:) = [ idx1 idx2 ones(base,1) fid];

            G.faces(1:base,1:4) = xfaces(:,1:4);

            % Indicate boundary faces in cells matrix
            X = xfaces(:,4);         % index for xfaces

            idx1 = X(G.cells(xfaces(:,1),13)==1);      % row index of xfaces that have col. 13 = 1
            G.cells(xfaces(idx1,1),14) = xfaces(idx1,2); % Update Col.14 = index of interior cell nbr.

            idx2 = X(G.cells(xfaces(:,2),13)==2);
            G.cells(xfaces(idx2,2),14) = xfaces(idx2,1); % Update Col.14 = index of interior cell nbr.


            % -------- y-faces ---------- %
            numYfaces = nx*(ny+1)*nz;

            I = kron(i,ones((ny+1)*nz,1));
            J = kron(ones(nx,1),kron(j1,ones(nz,1)));
            K = kron(ones((nx)*(ny+1),1),k);

            fid = (nx)*(ny+1)*(K-1)+(nx)*(J-1)+I;
            jdx1 = (nx+2)*(ny+2)*K+(nx+2)*(J-1)+(I+1);
            jdx2 = (nx+2)*(ny+2)*K+(nx+2)*J+(I+1);

            yfaces(fid,:) = [ jdx1 jdx2 2*ones(numYfaces,1) base+fid];

            G.faces(base+1:base+numYfaces,1:4) = yfaces(:,1:4);

            % Indicate boundary faces in cells matrix
            Y = yfaces(:,4)-base;                             % index for yfaces

            jdx1 = Y(G.cells(yfaces(:,1),13)==4);           % row index of xfaces that have col. 13 = 1
            G.cells(yfaces(jdx1,1),14) = yfaces(jdx1,2);    % Update Col.14 = index of interior cell nbr.

            jdx2 = Y(G.cells(yfaces(:,2),13)==8);
            G.cells(yfaces(jdx2,2),14) = yfaces(jdx2,1);    % Update Col.14 = index of interior cell nbr.

            base = base+numYfaces; % =base+nx*(ny+1)*nz;


            % -------- z-faces ---------- %
            numZfaces = nx*ny*(nz+1);

            I = kron(i,ones(ny*(nz+1),1));
            J = kron(ones(nx,1),kron(j,ones((nz+1),1)));
            K = kron(ones((nx)*(ny),1),k1);
            
            fid = nx*ny*(K-1)+nx*(J-1)+I;
            
            kdx1 = (nx+2)*(ny+2)*(K-1)+(nx+2)*J+(I+1);
            kdx2 = (nx+2)*(ny+2)*K+(nx+2)*J+(I+1);
            
            zfaces(fid,:) = [ kdx1 kdx2 3*ones(numZfaces,1) base+fid];
            
            G.faces(base+1:base+numZfaces,1:4) = zfaces(:,1:4);
            
            % Indicate boundary faces in cells matrix
            Z = zfaces(:,4)-base;   % index for zfaces
            
            kdx1 = Z(G.cells(zfaces(:,1),13)==16);        % row index of xfaces that have col. 13 = 1
            G.cells(zfaces(kdx1,1),14) = zfaces(kdx1,2);    % Update Col.14 = index of interior cell nbr.
            
            kdx2 = Z(G.cells(zfaces(:,2),13)==32);
            G.cells(zfaces(kdx2,2),14) = zfaces(kdx2,1);    % Update Col.14 = index of interior cell nbr.
            
            % Boundary Indicator (cell_id = bdry face, 0 = interior)
            G.faces(:,5) = G.bindicator();
            
            %% Grid Domain
            %G.domain=max(G.cells(:,4:6).*2.^(G.cells(:,7:9)));   % 12s/1000 calls
            G.domain = max(G.cells(:,1:3));

            %% Set Generic Dirichlet Boundary Conditions (-1 indicates no-flow) 
            % (i.e., constant-pressure /  no-flow)

            % IDENTIFIES left and right boundary CELLS
            % x-boundary cells
            x1 = G.cells(G.cells(:,1)==0 & G.cells(:,14)>0,19);
            x0 = G.cells(G.cells(:,1)==G.domain(1) & G.cells(:,14)>0,19);

            % y-boundary cells
            y1 = G.cells(G.cells(:,2)==0 & G.cells(:,14)>0,19);
            y0 = G.cells(G.cells(:,2)==G.domain(2) & G.cells(:,14)>0,19);

            % z-boundary cells
            z1 = G.cells(G.cells(:,3)==0 & G.cells(:,14)>0,19);
            z0 = G.cells(G.cells(:,3)==G.domain(3) & G.cells(:,14)>0,19);


            % Generic Dirichlet (Future release consider linear Dirichlet BC's)
            G.cells(x0,16) = 0;
            G.cells(x1,16) = 1;
            G.cells(y0,17) = 0;
            G.cells(y1,17) = 1;
            G.cells(z0,18) = 0;
            G.cells(z1,18) = 1;

            % Neumann BC (No-flow)
            G.cells([y0;y1;z0;z1],16) = -1;   
            G.cells([x0;x1;z0;z1],17) = -1;   
            G.cells([x0;x1;y0;y1],18) = -1;   

           % Assign Interior Pressure Values (Linearly)


            %% Assign Grid Properties

            G.nbfaces = sum(G.faces(:,5)>0);
            % G.cells = G.cells(sum(G.cells(:,4:6),2)>1,:);
            % G.ncells = (nx+2)*(ny+2)*(nz+2) - 4*(nx+ny+nz) - 8;
            G.icells = find(prod(G.cells(:,4:6),2)); % finds index of nonzero entries in product of elements.  0 in any dimension (x,y,y = 4:6) is a boundary
            G.cells(G.icells,15) = (1:G.nicells);  % G.cells(G.icells,15)=(1:G.nicells)'; 14s/1000
            % G.cells(:,19) = 1:G.ncells;
            G.ghosts = G.getGhosts;            % boundary cells ordered by cells
            G.bcells = G.faces(G.bfaces(),5);% boundary cells ordered by face % =G.faces(G.faces(:,5)>0,5);
            %G.setfake;
            
            
        end % Grid Constructor Method 
          % / ------------------------------------------------- 
          

            
        %% LOCAL FUNCTIONS
        %  ------------------------------------------------- 
        % Boundary functions
        % 
        % ---------------------- %
        % GET boundary pressure values AND type indicators
        function [bv, dn] = bvals(G,bc)
             % G.bvals(bc) returns bdry vals in cell 16:18 and dn
             % ORDER by FACE
             % bv = boundary (pressure) values driving flow in direction = bc
             % dn = logical (=1 if Dirichlet, 0 if Neumann)
             dn = G.cells(G.bcells,15+bc)>=0;  % size = # boundaries
             % dn = ~isnan(G.cells(G.bcells,15+bc));
             bv = G.cells(G.bcells,15+bc).*dn; %
        end


        %% Get (all) Dirichlet cells    
        function dcells = dirichlet(G)
            dcells = G.cells(logical(sum(~isnan(G.cells(:,16:18)),2)>0),19);
        end


        %% Ghost + Interior cell IDs
        % Get (ALL) non-Edge and non-Corner cells (i.e., ghost + interior)    
        function bc = realcells(G)
            bc = G.cells(sort([G.ghosts; G.icells]),19);
            % ordered by CELLS
            % G.cells(G.cells(:,14)>0 | G.cells(:,15)>0,19); %
            %    = bdry & interior = real cells
        end


        %% Get (all) Boundary Cells    
        function bc = getGhosts(G)
            bc = G.cells(G.cells(:,14)>0,19);
            % cells(:,15) = boundries, corners, & edges
        end

         %% Set Edge and Corner Cell to NaN    
        function setfake(G)
            % Fake cells are edges and corner cells 
            edges_and_corners = setdiff((1:G.ncells)',G.realcells);
            G.cells(edges_and_corners,16:18) = NaN;
            % cells(:,15) = boundries, corners, & edges
        end

        %% Fake = corner or edge
        function fake = enc(G)
            edges_and_corners=isnan(G.cells(:,16));
            fake = G.cells(edges_and_corners,19);
            % cells(:,15) = boundries, corners, & edges
        end

        %% Get (all) Boundary Faces (ordered by Face ID)   
        function bf = bfaces(G)
            bf = G.faces(G.faces(:,5)>0,4);
        end

        %% Get (all) Interior Faces (ordered by Face ID)   
        function intfaces = ifaces(G)
            intfaces = G.faces(G.faces(:,5) == 0,4);
        end

        %% Left and Right Boundary Faces (orderd by Face)
        % ---------------------- %
        % Get Left Boundary Faces (ordered by Face ID)   
        function bf = lfaces(G)
            bf = G.faces(G.faces(:,1)==G.faces(:,5),4);
        end


        %% Get Right Boundary Faces (ordered by Face ID)   
        function bf = rfaces(G)
            bf = G.faces(G.faces(:,2)==G.faces(:,5),4);
        end

        %% Get face_ids by orientation (ordered by Face ID)   
        % Need faces by orientation?

        % x-faces
        function xf = xfaces(G)
            xf = G.faces(G.faces(:,3)==1,4);
        end

        % y-faces     
        function yf = yfaces(G)
            yf=G.faces(G.faces(:,3)==2,4);
        end

        % z-faces     
        function zf = zfaces(G)
            zf=G.faces(G.faces(:,3)==3,4);
        end

        %% Boundary indicators
        function bi = bindicator(G)
            % RETURNS cell ID if boundary, else 0
            % Selects (by face) the cell associated with a boundary face
            % ORDERED BY FACE
            % 0 = interior cell
            bdry_cells = G.getGhosts;
            X = [ismember(G.faces(:,1),bdry_cells),ismember(G.faces(:,2),bdry_cells)]; % X is recycled (dummy) variable
            bi = sum(X.*G.faces(:,1:2),2);
        end


        %% Boundary Map (Faces --> Cells)
        function face2cell=bmap(G)
            % BMAP returns mapping between boundary faces and cells
            % Boundary Face ID --> Boundary Cell ID

            bcid = bindicator(G);

            bcid = bcid(bcid > 0);  % Cells Ordered by Face
            bfid = G.bfaces();
            face2cell = [bfid, bcid];
        end


        %% Left-Right Boundary Cell IDs
        function bcell_LR=rbc(G,side)                   
            % Select cells from Left or right side (ordered by Face)
            if side == 1
                bcell_LR = G.faces(G.lfaces,5);
            elseif side == 2
                bcell_LR = G.faces(G.rfaces,5);
            else
                bcell_LR = G.faces(G.lfaces,5);
            end
        end

        %% 
        function y = get_cell_neighbors(G,list)
            % GET_NEIGHBORS returns an array of cell ids
            % given a list of cells.  
            % For example, list=[44, 3, 9] as a list of cells.
            % CALL: G.get_neighbors(list) return an array of the 
            % original list (1st column) and rows 2:24 its neigboring
            % cells.  y = cells ordered by face id.

            % CURRENTLY list = singleton 

            [m,n] = size(list);
                if m > n
                    list = list';
                end

            y = zeros(24,length(list),'uint32');

            f1 = G.faces(G.faces(:,1)==list,4);         % faces bottom/left/down
            f2 = G.faces(G.faces(:,2)==list,4);         % top/right/up
            n = [G.faces(f1,2); G.faces(f2,1)];
            y(1:length(n)) = n;
            y = y(y > 0);
        end
          
    end % Methods
     
end % Class


%% NOTES:

%{

% Selects real cells.  Call prior to building all cells array
sum(G.cells(:,1:3)==0,2)<2



Question 3/20/2020 - Grid domain

 % Refinement/Update  (not implemented here, see refine_grid)    

% Refinement adds 1 to G.cells(:,7:9).  The dx, dy, dz play no role in this
% value.  Think of starting with 1 x 1 x 1 grid (1 block).  
%         Then columns 7:9 tell how many times grid has been refined in that direction.    
%   For example, Gc=Grid(2,2,2,4,4,4);  Gc.cells(Gc.icells,7:9) --> 1 1 1
%                Qc=Grid(4,4,4,2,2,2);  Qc.cells(Qc.icells,7:9) --> 2 2 2
%   In the first example, since there are 2 cells in each dimension, there
%   effectively is 1 level of refinement.  In example 2, to go from (1x1x1)
%   to (4x4x4), there would be 2 levels of refinement in each dimension.
% Another example, G=Grid(256,64,16,dx,dy,dz) --> 8 6 4 as refinement
% levels 7:9 for each cells.  Note that dx, dy, and dz do not matter.  

% Mapping: Cell number of interior cell to Matrix Row
          % G.cells(G.icells,[15,19])
%
% BoundryCellsByFace=G.faces(G.faces(:,5)>0,[1:2,4])

% NOTES: may want integer arithmetic to wrap as opposed to saturate (which is what MATLAB does) on overflow.
%
%   Vector order matters!  That is when cell indexing, 
%   X([1,2,3]) ~= X([2,3,1]) as an ordered array (vector), but equal as
%   sets.  So, when obtaining cell ID for computation, G.cells(X, 1), X
%   Makes difference.

Modification History:
01/08/2020
     1. The bug in function bvals where dn returns all true fixed.  
             Applied Hadamard product to boundary values.
             dn.*G.cells(:,15+bc), since -1 indicates Neumann
             boundary condition, and we ultimately want this value
             to be 0 in matrix setup and calculations.
             This was caused by assigning: 
               dn=~isnan(G.cells(G.bcells,15+bc));
             after code was changed not to assign NaN to 
             boundary values in Section "Set Generic Dirichlet Boundary
             Conditions"

12/05/2019 
     1. Implemented linear pressure values per direction when initiating Grid.
     2. Added 2 Methods LFACES and RFACES to return left/right boundary faceIDs
          - G.lfaces = G.faces(G.faces(:,1)==G.faces(:,5),4);
          - G.lfaces = Boundary faces on left side of grid
             
10/22/2019 - Grid needs to be rethought.  For example, the following are
potential inefficiencies that may need redesigned and recoded.  
    1. Corner points and edge points are stored in the structure.  Not
    necessary and adds large amount to storage.  Every call to retrieve
    boundary cells requires filtering out corners and edges.

    2. Cell structure stores large amount of data as well.  Perhaps it
    should be redesigned.  For instance, the number of times that cell
    bitmask is needed could use a function to call and compute it on the
    fly.

    3. Number of dummy variables seems excessive.

    4. kron is extremely slow.  Try redesigning using meshgrid instead.

    5. Consider either making a subclass of Grid for local grid compuations
    since much of the data is NOT used in local calculuations OR making
    local setup based on conditions.  For example, cells 7:9 refinement
    level is not needed/used in local calculuations.  Same with 10:12
    finescale indices.  Possibly same for 13 and 16:18.

    6. Do you store or compute?  For example, does finding and storing
    ghost cells, then recalling them more efficient than finding them 
    each time needed?  

        a. Store referential id vs. using the FIND command.
 
    7. Parallel computing issues with associativity?

    8.  Perhaps ID should be a column itself.  So instead of G.cells(:,19)
    used as self-referential id to avoid FIND command, there should be an
    ID property, i.e., G.id=(1:ncells)';  Then when looking for the cell
    number with a certain condition, use G.id(G.cells(:,4)>0).  G.cell_id

    9.  Similiarly, a new class for Pressure should be created to store
    P.id and P.vals (size numel(G.cells) x 3 (or by number of bc)).

10/31/2019

    1. Added mapping between boundary cells and boundary faces.  From the
    face structure column 5 contains the cell from column 1:2 which is a
    boundary cell (if either).  

    2. Gridology:  Grid needs to be redesigned (potential student project).
     No need to store corner points and edges.

    3. Look to add events and listeners
https://www.mathworks.com/help/matlab/matlab_oop/events-and-listeners-syntax-and-techniques.html

    4. pack - Consolidate workspace memory.






          
          % x-boundary cells
          % G.cells(G.cells(:,16)>=0,19);  % Dirichlet BC all bdry cells
          % G.cells(~isnan(G.cells(:,16)),19); % bdry cells of x-faces
          % G.cells(~isnan(G.cells(:,18)),19); % bdry cells of z-faces
          % G.cells(~isnan(G.cells(:,16)),19)
          % [dn, G.bcells G.cells(G.bcells,16)]; x-boundary info
          % dn =1 if Dirichlet, 0 if Neumann

          
          % Not Implemented Yet
          
          % ---------------------- %
          % Get Boundary Faces                      % Currently Not Class Method
          %function bf=bfaces(G)
          %      bf=G.faces(G.faces(:,5)==1,4);
          % end
          % BorderFaces=G.faces(G.faces(:,5)>0,[1:2,4])
          % BoundryCellsByFace=G.faces(G.faces(:,5)>0,[1:2,4])
          
          % bdfaces     % G.faces(G.faces(:,5)==1,4) % ordered by FaceID
          
          % Cells involved in x-faces
          % G.faces(G.faces(:,5)==1 & G.faces(:,3)==1,1:2)
          %
          % x-faces
          % G.faces(G.faces(:,5)==1 & G.faces(:,3)==1,4)
          
          % ---------------------- %
          % Get Number of Boundary Faces            % Currently Not Class Method
          
          % Mapping: Cell number of interior cell to Matrix Row
          % G.cells(G.icells,[15,19])
          
          % Viewgrid                                % Currently Not Class Method




%}

