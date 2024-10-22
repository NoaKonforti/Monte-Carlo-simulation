classdef grids
    properties
    n_r % num of slots in r
    n_z % num of slots in z
    n_a % num of slots in a
    dr %dr- resolution for r ;where n_z*dz = dr (num slot*resolution = the
       %layer)
    dz %dz- resolution for z
    da %da- resolution for a
    ir %index r
    iz %index z
    ia %index a
    r %array r; ir starts from 0
    z %array z; iz starts from 0
    a %array a; ia starts from 0
    A_rz %absorbtion matrix - matrix r over z
    Az
    A
    A_not_scatter
    T_ra %transmittion matrix
    T
    Rd_ra %reflection matrix
    Rd
    end
    
    methods
        function obj = grids(n_r, n_z, n_a, dr, dz, da)
            obj.n_r = n_r;
            obj.n_z = n_z;
            obj.n_a = n_a;
            obj.dr = dr;
            obj.dz = dz;
            obj.da = da;
            obj.ir = 0:(n_r - 1);
            obj.iz = 0:(n_z - 1);
            obj.ia = 0:(n_a - 1);
            obj.r = (obj.ir+0.5+1./(12.*(obj.ir+0.5))).*dr;
            obj.z = (obj.iz+0.5).*dz;
            obj.a = (obj.ia+0.5)*da+(1-0.5*da*cot(da))*cot((obj.ia+0.5)*da);
            obj.A_rz = zeros(n_r, n_z);
            obj.T_ra = zeros(n_r, n_a);
            obj.Rd_ra = zeros(n_r, n_a);
            obj.A_not_scatter = zeros(n_r,1);
        end

        function [coor_r, num_coor] = get_coor_r(obj, ri)
            if ismember(ri, obj.r) %ri is the exact value in r array
                coor_r = find(obj.r == ri);
                num_coor = 1;
            elseif ri > obj.r(end)
                coor_r = obj.n_r;%should add -1
                num_coor = 1; 
            else
                differences = abs(obj.r - ri);
                [min_diff, ~] = min(differences);
                indexes = find(differences == min_diff);
                if length(indexes) == 1 %ri is close to a specific value in r also if passes the border
                    coor_r = indexes;
                    num_coor = 1;
                else %ri falls inbetween 2 coords and needs to split equally 
                    coor_r = indexes;
                    num_coor = 2;
                end
            end  
        end

        function [coor_z, num_coor] = get_coor_z(obj, zi)
            if ismember(zi, obj.z) %zi is the exact value in r array
                coor_z = find(obj.z == zi);
                num_coor = 1;
            else
                differences = abs(obj.z - zi);
                [min_diff, ~] = min(differences);
                indexes = find(differences == min_diff);
                if length(indexes) == 1 %zi is close to a specific value in z also if passes the border
                    coor_z = indexes;
                    num_coor = 1;
                else %zi falls inbetween 2 coords and needs to split equally 
                    coor_z = indexes;
                    num_coor = 2;
                end
            end  
        end
        
        function [coor_a, num_coor] = get_coor_a(obj, ai)
            if ismember(ai, obj.a) %ai is the exact value in r array
                coor_a = find(obj.a == ai);
                num_coor = 1;
            else
                differences = abs(obj.a - ai);
                [min_diff, ~] = min(differences);
                indexes = find(differences == min_diff);
                if length(indexes) == 1 %ai is close to a specific value in a also if passes the border
                    coor_a = indexes;
                    num_coor = 1;
                else %ai falls inbetween 2 coords and needs to split equally 
                    coor_a = indexes;
                    num_coor = 2;
                end
            end  
        end

        function obj = update_Arz(obj, ri, zi, dw)
            [coor_r, num_coor_r] = obj.get_coor_r(ri);
            [coor_z, num_coor_z] = obj.get_coor_z(zi);
            obj.A_rz(coor_r,coor_z) = obj.A_rz(coor_r,coor_z) + dw/(num_coor_r * num_coor_z);

        end
        
        function obj = update_Tra(obj, ri, ai, w)
            [coor_r, num_coor_r] = obj.get_coor_r(ri);
            [coor_a, num_coor_a] = obj.get_coor_a(ai);
            obj.T_ra(coor_r,coor_a) = obj.T_ra(coor_r,coor_a) + w/(num_coor_r * num_coor_a);

        end
        
        function obj = update_Rdra(obj, ri, ai, w)
            [coor_r, num_coor_r] = obj.get_coor_r(ri);
            [coor_a, num_coor_a] = obj.get_coor_a(ai);
            obj.Rd_ra(coor_r,coor_a) = obj.Rd_ra(coor_r,coor_a) + w/(num_coor_r * num_coor_a);
        end

        function obj = normalize_Arz(obj, dr, dz, n_packets, nr, nz)
            obj.A = sum(obj.A_rz(:))/n_packets;
            obj.Az = sum(obj.A_rz, 1)/(dz*n_packets);
            delta_a = 2*pi*dr^2*dz*n_packets*(repmat(0.5+(1:nr).',[1,nz]));
            obj.A_not_scatter = obj.A_not_scatter/(n_packets*dz);
            obj.A_rz(:,end) = transpose(obj.A_not_scatter);
            obj.A_rz = obj.A_rz./(delta_a);
        end
        
        function obj = update_A_not_scatter(obj, ri, dw)
            [coor_r, ~] = obj.get_coor_r(ri);
            obj.A_not_scatter(coor_r,1) = obj.A_rz(coor_r,1) + dw;

        end
        
        function obj = normalize_A_not_scatter(obj, ri, dw)
            
        end

        function obj = normalize_Rdra(obj, dr, da, n_packets, nr, na)
            obj.Rd = sum(obj.Rd_ra(:))/n_packets;
            delta_a = 2*pi*dr^2*(0.5+(1:nr));
            delta_sigma = 4*pi*sin((0.5+(1:na)).*da)*sin(0.5.*da);
            obj.Rd_ra = obj.Rd_ra./(delta_a.*cos(obj.a).*delta_sigma.*n_packets);
        end

        function obj = normalize_Tra(obj, dr, da, n_packets, nr, na)
            obj.T = sum(obj.T_ra(:))/n_packets;
            delta_a = 2*pi*dr^2*(0.5+(1:nr));
            delta_sigma = 4*pi*sin((0.5+(1:na)).*da)*sin(0.5.*da);
            obj.T_ra = obj.T_ra./(delta_a.*cos(obj.a).*delta_sigma.*n_packets);
        end
    end 
end 