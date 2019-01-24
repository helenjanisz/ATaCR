% function to read a sac pole zero file
% usage:
%   [zeros, poles, constant] = read_sac_pole_zero(pole_zero_file_name)
%   creates two output vectors (poles, zeros) and one scalar (constant)
% file is the format as output from rdseed if the option
% "output polezero file" is accepted
% format:
% ZEROS nzeros (only listed if not at 0+0i - this version only allows zeros at the origin!)
% zero_r zero_i
% zero_r zero_i
% ...    ...
% POLES npoles
% pole_r pole_i
% pole_r pole_i
% ...     ...
% CONSTANT constant
% example:
%    ZEROS 3
%    POLES 5
%   -0.0370  0.0370
%   -0.0370  -0.0370
%   -118.7520  423.4880
%   -118.7520  -423.4880
%   -251.3270  0.0000
%   CONSTANT 3.127953e+16
%
% RWP 5.30.2011
% added functionality to allow #, *, or % as comments in sac pole zero file



function [zz,pp,constant] = read_sac_pole_zero(pole_zero_file_name)
    % First open the file
    pz_fid = fopen(pole_zero_file_name);
    %initialize some flags and various variables
    buff = 0;
    nzeros = 0;
    npoles = 0;
    zero_scan_flag = -1;
    pole_scan_flag = -1;
    pole_count_flag = 1;
    zero_count_flag = 1;
    zz=0;
    pp=0;
    % loop over the entire file
    while (buff ~= -1)
        % read the next line in the file
        buff = fgets(pz_fid);
        % check to make sure it isn't the end of the file
        if (buff ~= -1)
            %assume first line is ZEROS nzeros
%           tmp = sscanf(buff, '%s" "%d');
            if (strcmp(buff(1),'*'))
            elseif (strcmp(buff(1),'#'))
            elseif (strcmp(buff(1),'%'))
            elseif strcmp(buff(1), 'Z')
                zero_scan_flag = 1;
                pole_scan_flag = 0;
                tmp = sscanf(buff, '%s %d');
                nzeros = tmp(6);
            elseif strcmp(buff(1), 'P')
                pole_scan_flag = 1;
                zero_scan_flag = 0;
                tmp = sscanf(buff,'%s %d');
                npoles = tmp(6);
            elseif strcmp(buff(1), 'C')
                pole_scan_flag = 0;
                tmp = sscanf(buff, '%s %e');
                constant = tmp(9);
            elseif (zero_scan_flag == 1)
                if (zero_count_flag < nzeros+1)
                    tmp = sscanf(buff,'%f %f');
                    zero_count_flag = zero_count_flag + 1;
                    if (zero_count_flag > 1)
                        n=zero_count_flag-1;
                        z1(n)=tmp(1);
                        z2(n)=tmp(2);
                    end
                end
            elseif (pole_scan_flag == 1) 
                if (pole_count_flag < npoles+1)
                    tmp = sscanf(buff, '%f %f');
                    pole_count_flag = pole_count_flag + 1;
                    if (pole_count_flag > 1) 
                        n=pole_count_flag-1;
                        p1(n)=tmp(1);
                        p2(n)=tmp(2);
                    end
                end
            end
        end
    end
    if nzeros == 0
        zz = [];
    else
    zz = complex(z1', z2');
    end
    if npoles == 0
        pp = [];
    else
    pp = complex(p1', p2');
    end
    % fill in missing poles and zeros
    if (length(zz) < nzeros)
        for n=length(zz)+1:nzeros
            zz(n) = 0+0i;
        end
    end
    if (length(pp) < npoles)
        for n=length(pp)+1:npoles
            pp(n) = 0+0i;
        end
    end
    
    
    fclose(pz_fid);

    return
