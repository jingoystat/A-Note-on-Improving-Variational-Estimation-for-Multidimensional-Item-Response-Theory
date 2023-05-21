person_list = [200, 500];
item_list = [30, 55];
domain_list = [2, 5];
r_list = ["high", "low"];
within_list = [0, 1];


maxiter = 50; 
Nrep = 100;
lr = 0.01; 
S = 10;
M = 10; 

parfor p = 1:2
    for i = 1:2
        for rl = 1:2
            for w = 1:2
                person = person_list(p);
                item = item_list(i);
                domain = domain_list(i);
                r = r_list(rl);
                within = within_list(w);

                simulation_IS_parallel(person, item, domain, r, within, maxiter, Nrep, S, M);
            end
        end
    end
end

