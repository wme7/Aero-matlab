% Writing test
write_ans = 1;

%% Open a Files to store the Results
if write_ans == 1
    file = fopen(IDn,'w');
    % 'file' gets the handel for the file "case.plt".
    % 'w' specifies that it will be written.
    % similarly 'r' is for reading and 'a' for appending.
    fprintf(file, 'TITLE = "%s"\n',ID);
    fprintf(file, 'VARIABLES = "x" "y" "n" "E" "p" "t" "z"\n');
end

%% Write Results
if write_ans == 1 && (mod(tsteps,5*dt) == 0 || tsteps == time(end))
    fprintf(file, 'ZONE T = "time %0.4f"\n', tsteps);
    fprintf(file, 'I = %d, J = %d, K = 1, F = POINT\n\n',nx,ny);
    for j = 1:ny
        for i = 1:nx
            fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                x(j,i),y(j,i),rho(j,i),E(j,i),p(j,i),t(1,1,j,i),z(1,1,j,i));
        end
    end
end

%% Write Results
fprintf('Simulation has been completed succesfully!\n')
if write_ans == 1
    fclose(file);
    fprintf('All Results have been saved!\n')
end