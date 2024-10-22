classdef packet 
    
    properties
    x, y, z % coordinates of photon packet
    ux, uy, uz % direction cosines
    w % weight of photon packet
    dead % 0 if propagating, 1 if terminated 
    s % step size
    s_left
    scatters % # of scatter events
    layer % up to 2 layers in our case
    update %update tramsnit = 1, reflection =2 matrixes
    bio_layers % # of biological layers
    r_sp %specular reflection
    end
    
    methods

        function obj = packet(x, y, bio_layers, layers_lst)
            obj.x = x;
            obj.y = y;
            obj.z = 0;
            obj.ux = 0;
            obj.uy = 0;
            obj.uz = 1;
            obj.w = 1;
            obj.dead = 0;
            obj.scatters = 0;
            obj.layer = 2;
            obj.s = 0;
            obj.update = 0;% 0-no update; 1-update transmit matrix; 2-update reflection matrix
            obj.bio_layers = bio_layers;
            top_layer = layers_lst(1);
            current_layer = layers_lst(2);
            %bottom_layer = layers_lst(3); 
            n1 = top_layer.n;
            n2 = current_layer.n;
            obj.r_sp = 0; 
            if n1 ~= n2
                obj.r_sp = ((n1 - n2)/(n1 + n2))^2;
                obj.w = 1 - obj.r_sp;
            else
                obj.w = 1;
            end
        end

        function obj = step_sample(obj, layers_lst)
            % sampling the step size for a free path
            curr_layer = layers_lst(obj.layer);
            xi = rand();
            obj.s = (-log(xi)/ curr_layer.mu_t);
            %obj.s = (-log(xi)/ curr_layer.mu_s); %3.11
            obj.s_left = obj.s;
        end

        function obj = move_me(obj,layers_lst)
            % check if packet will cross boundry
            curr_layer = layers_lst(obj.layer);
            mu_t = curr_layer.mu_t;
            boundry = [curr_layer.z0, curr_layer.z1];
            if obj.uz > 0
                db = (boundry(2) - obj.z)/obj.uz;
            elseif obj.uz == 0
                db = inf;
            else 
                db = (boundry(1)- obj.z)/obj.uz;
            end
            if db*mu_t <= obj.s
                % cross boundry
                obj.x = obj.x + obj.ux*db;
                obj.y = obj.y + obj.uy*db;
                obj.z = obj.z + obj.uz*db;
                obj.s = obj.s - db*mu_t;
                obj.s_left = obj.s;
            else
                % didn't cross boundry
                obj.x = obj.x + obj.ux*obj.s; 
                obj.y = obj.y + obj.uy*obj.s;
                obj.z = obj.z + obj.uz*obj.s;
                obj.s = 0;
            end
        end

        function [r, at] = fresnel(obj, layers_lst)
            ai = acos(abs(obj.uz));
            current_layer = layers_lst(obj.layer);
            sg = sign(obj.uz);
            if sg < 0
                next_layer  = layers_lst(obj.layer - 1);
            else
                next_layer = layers_lst(obj.layer + 1);
            end
            ni = current_layer.n;
            nt = next_layer.n;
            at = asin((ni/nt) * sin(ai));
            if ni == nt %transfer between layers n0 == n1 
                r = 0;
            elseif ai == 0 %perpendicular Fresnel reflection
                r = ((ni - nt)/ (ni + nt))^2;
            elseif ni > nt && ai > asin(nt/ni) %n0 > n1: ai is critical angle
                r = 1; % total internal reflection
            else %general case of Fresnal reflection
                r = 0.5 *((sin(ai-at)/sin(ai+at))^2 + (tan(ai -at)/tan(ai + at))^2);
            end
             
        end 

        function obj = reflection_transmission(obj, layers_lst)
            [r, at] = fresnel(obj, layers_lst);
            xi = rand();
            sg = sign(obj.uz);
            ni = layers_lst(obj.layer).n;
            nt = layers_lst(obj.layer + sg).n;
            if xi > r %transmitted
                obj.ux = obj.ux * (ni/nt);
                obj.uy = obj.uy * (ni/nt);
                obj.uz = sg * cos(at);
                if obj.layer == 2 && sg == -1  
                    obj.dead = 1;
                    obj.update = 2; %update reflection matrix
                elseif obj.layer == (obj.bio_layers + 1) && sg == 1
                    obj.dead = 1;
                    obj.update = 1; %update transmit matrix
                else
                    obj.layer = obj.layer + sg;
                end
            else % reflected
                obj.uz = -obj.uz; 
            end 

        end

        function [obj, delta_w] = absorption(obj, layers_lst)
            % update weights according to absorbtion
            curr_layer = layers_lst(obj.layer);
            delta_w = (curr_layer.mu_a/curr_layer.mu_t)*obj.w;
            %delta_w = obj.w*(1-exp(-curr_layer.mu_a*obj.s_left));%3.11
            obj.w = obj.w - delta_w;            
        end
        
        function obj = scattering(obj, layers_lst)
            % update cosine direction according to scattering
            curr_layer = layers_lst(obj.layer);
            g = curr_layer.g;
            xi = rand();
            % sampling the polar scattering angle
            if g == 0
                cos_theta = 2*xi - 1;
            else
                cos_theta = (1/(2*g))*(1 + g^2 - ((1 - g^2)/(1 - g + 2*g*xi))^2);
            end
            % sampling the azimuthal scattering angle
            xi = rand();
            theta = acos(cos_theta);
            phi = 2* pi * xi;
            if abs(obj.uz) > 0.9999
                obj.ux = sin(theta) * cos(phi);
                obj.uy = sin(theta) * sin(phi);
                obj.uz = sign(obj.uz) * cos(theta);
            else
                obj.ux = (sin(theta) * obj.ux *obj.uz * cos(phi) - obj.uy *sin(phi))/ (sqrt(1 - obj.uz^2)) + obj.ux * cos(theta);
                obj.uy = (sin(theta) * obj.uy *obj.uz * cos(phi) + obj.ux *sin(phi))/ (sqrt(1 - obj.uz^2)) + obj.uy * cos(theta);
                obj.uz =  (-sqrt(1- obj.uz^2)) * sin(theta) * cos(phi) + obj.uz * cos(theta);
            end
            obj.scatters = obj.scatters +1;
        end

        function obj = termination(obj, m)
            % terminate photon packate
            xi = rand();
            if xi <= 1/m 
                obj.w = m * obj.w;
            else
                obj.w = 0;
                obj.dead = 1;
            end
        end
    end
end
