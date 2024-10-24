classdef layer
    properties
        n
        z0
        z1
        mu_a
        mu_s
        mu_t
        g
    end

    methods
        function obj = layer(n, z0, z1, mu_a, mu_s, g)
            obj.n = n;
            obj.z0 = z0;
            obj.z1 = z1;
            obj.mu_a = mu_a;
            obj.mu_s = mu_s;
            obj.mu_t = mu_a + mu_s;
            obj.g = g;
        end
    end
end
