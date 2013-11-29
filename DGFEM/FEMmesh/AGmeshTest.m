    %% Example 
    % by Manuel Diaz, NTU, 2013.11.28
    
    clear all; close all; cla;
    Lx=[0,1];   Ly=[0,1];
    N=10;        M=8;
    
    % Generate mesh
    %[VX VY EtoV nV nE] = AGmeshGenerator(Lx,Ly,N,M,'Horseshoe');
    [VX VY EtoV nV nE] = AGmeshGenerator(Lx,Ly,N,M,'PolarCoordinates');
    
    % Display data arrays
    fprintf('Number of elements = %1.0f\n', nE );
    fprintf('Number of nodes = %1.0f\n', nV );
    disp('nodes Coordinates'); disp([VX,VY]); 
    disp('Connectivity Table'); disp(EtoV);

    % mistake element one in purpose
    i = 13; EtoV(i,:) = EtoV(i,[1 4 3 2]); clear i; 
    
    % Plot elements
    subplot(1,2,1); patch('Faces',EtoV,'Vertices',[VX,VY],...
        'FaceColor','none','EdgeColor','k'); hold on;   
    
        % Plot nodes and elements info
        for i = 1:nV
            text(VX(i),VY(i),int2str(i),'fontsize',8,....
                'fontweight','bold','Color','r'); 
        end 
        for e = 1:nE
            pos = [mean(VX(EtoV(e,:))),mean(VY(EtoV(e,:)))];
            text(pos(1),pos(2),int2str(e),'fontsize',8,...
                'fontweight','bold','color','b');
        end

	% test elements (one by one)
    subplot(1,2,2);
    for e = 1:nE
    ax = VX(EtoV(e,1)); ay = VY(EtoV(e,1));
    bx = VX(EtoV(e,2)); by = VY(EtoV(e,2));
    cx = VX(EtoV(e,3)); cy = VY(EtoV(e,3));
    dx = VX(EtoV(e,4)); dy = VY(EtoV(e,4));
    Area = 0.5*(-ay.*bx+ax.*by-by.*cx+bx.*cy-cy.*dx+cx.*dy);
    %Area = 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
        if Area>0
            patch([ax bx cx dx],[ay by cy dy],'b'); hold on;
            pause(0.1)
        else
            patch([ax bx cx dx],[ay by cy dy],'r'); hold on;
            pause(0.1)
        end
    end
    
    % Correct elements
   	ax = VX(EtoV(:,1)); ay = VY(EtoV(:,1));
    bx = VX(EtoV(:,2)); by = VY(EtoV(:,2));
    cx = VX(EtoV(:,3)); cy = VY(EtoV(:,3));
    dx = VX(EtoV(:,4)); dy = VY(EtoV(:,4));
    Area = 0.5*(-ay.*bx+ax.*by-by.*cx+bx.*cy-cy.*dx+cx.*dy);
    %Area = 0.5*(ax-dx).*(cy-dy)-(bx-dx).*(by-dy)-(cx-dx).*(ay-dy);
    i = find(Area<0); fprintf('problem elements: %1.0f\n',i);
    EtoV(i,:) = EtoV(i,[1 4 3 2]);
    
    % test elements (one by one)
    subplot(1,2,2);
    for e = 1:nE
    ax = VX(EtoV(e,1)); ay = VY(EtoV(e,1));
    bx = VX(EtoV(e,2)); by = VY(EtoV(e,2));
    cx = VX(EtoV(e,3)); cy = VY(EtoV(e,3));
    dx = VX(EtoV(e,4)); dy = VY(EtoV(e,4));
    Area = 0.5*(-ay.*bx+ax.*by-by.*cx+bx.*cy-cy.*dx+cx.*dy);
    %Area = 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
        if Area>0
            patch([ax bx cx dx],[ay by cy dy],'g'); hold on;
            pause(0.1)
        else
            patch([ax bx cx dx],[ay by cy dy],'r'); hold on;
            pause(0.1)
        end
    end
    
    