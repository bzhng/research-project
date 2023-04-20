import xarray as xr
import numpy as np
import pandas
from random import sample,choice
from time import time 
import matplotlib as mpl
import matplotlib.pyplot as plt



ds_range = xr.open_dataset('data/processed/range-pd-sync3-1.nc')
ds_stream = xr.open_dataset('data/processed/stream-pd-int2.nc')

latitude = ds_range.lat.values
longitude = ds_range.lon.values
t = ds_range.time.values

r_pd = ds_range.pd.values
s_pd = ds_stream.pd.values

r_area = 3.8*3e6 # 3.8 +/- 0.4 km longitude
s_area = 150000

r_power = r_pd*r_area
s_power = s_pd*s_area




not_nan = pandas.read_csv('data/processed/not-nan.csv')
df_idx = list(not_nan.index)




### place a threshold on min power during initialisation
### or just try 2 pixels at first

def initialise(N,k):
    '''
    random initial population
    '''
    
    init_pop = np.zeros((N,k),dtype=int)
    
    for i in range(N):
        chrom = np.array(sample(df_idx,k))
        init_pop[i,:] = chrom
    
    return init_pop







def fitness(pop):
    '''
    '''
    
    f_pop = np.zeros((pop.shape[0]))
    
    for i in range(pop.shape[0]):
        df_chrom = not_nan.iloc[pop[i]]
        chrom_pow = np.zeros((pop.shape[1],r_power.shape[0]))
        
        for j in range(pop.shape[1]):
            if df_chrom.iloc[j,0] == 'r':
                chrom_pow[j,:] = r_power[:,df_chrom.iloc[j,1],df_chrom.iloc[j,2]]

            elif df_chrom.iloc[j,0] == 's':
                chrom_pow[j,:] = s_power[:,df_chrom.iloc[j,1],df_chrom.iloc[j,2]]
                
        r_df_chrom = df_chrom[df_chrom['flag']=='r']
        s_df_chrom = df_chrom[df_chrom['flag']=='s']
        r_clump = np.zeros((len(r_df_chrom),8)) # any pixel can only have a max of 8 adjacent pixels
        r_clump[:] = np.nan
                
        for r in range(len(r_df_chrom)):
            r_df_chrom['lat_dev'] = abs(r_df_chrom.iloc[:,1]-r_df_chrom.iloc[r,1])
            r_df_chrom['lon_dev'] = abs(r_df_chrom.iloc[:,2]-r_df_chrom.iloc[r,2])
            r_df_chrom.loc[r,('lat_dev','lon_dev')] = np.nan
            adj = r_df_chrom[r_df_chrom['lat_dev']<=1]
            if len(adj) == 0:
                pass
            else:
                adj = adj[adj['lon_dev']<=1].index.values
                glob_idx = r_df_chrom.iloc[r].name
                adj = np.append(glob_idx,adj) 
                r_clump[r,:len(adj)] = adj
            
            
        
#         for s in range(len(s_df_chrom)): ### do stream sites need to be a certain size?


        ### or instead of iterating again, count occurrences of the leading value
        ### maybe put a plug in this for now as may not be the best approach
        
        for r in range(r_clump.shape[0]):
            if np.isnan(r_clump[r,0]) == True:
                pass
            else:
                for q in range(r_clump.shape[0]): ### but dont double count the same array, ie r=/=q
                    if r_clump[r,0] in r_clump[q,:]:
                        
        
                
        aggregate = 0
        for j in range(pop.shape[1]):
            aggregate += chrom_pow[j,:]
        
        num_sites = 
        
        
        av_agg = aggregate/num_sites    
        p_mean = np.mean(av_agg)
        rmsd = np.sqrt(np.mean((av_agg-p_mean)**2))
        f_pop[i] = p_mean/rmsd
        
        
    return f_pop







def selection(pop,f_pop):
    '''
    tournament selection
    '''
    
    num_par = round(0.3*len(f_pop))
    f_idx = list(np.arange(len(f_pop)))
    parents = np.zeros((num_par,pop.shape[1]),dtype=int)
    
    for i in range(num_par):
        t_size = 5 # tournament size
        t_idx = np.array(sample(f_idx,t_size)) # random indices of chromosomes in tournament
        f_max = max(f_pop[t_idx]) # max fitness of chromosomes in tournament
        p_idx = np.where(f_pop==f_max)[0][0] # index of fittest chromosome
        parents[i,:] = pop[p_idx,:] # fittest chromosome
        
    return parents







def crossover(parents,pop):
    '''
    ring crossover
    '''
    
    pc = 0.7 # probability of crossover
    cross_num = round(pc*parents.shape[0]) # number of chromosomes used for crossover
    
    # want an even number offspring created from crossover
    # if num parents is even and cross num is odd, reprod is odd, num cross off is odd
    # if num parents is odd and cross num is even, reprod is odd, num cross off is odd
    if parents.shape[0]%2 == 0 and cross_num%2 == 1 or parents.shape[0]%2 == 1 and cross_num%2 == 0:
        cross_num -= 1
        
    
    par_idx = list(np.arange(parents.shape[0])) # indices of parents to select from 
    randp_idx = np.array(sample(par_idx,cross_num)) # random indices for pool of chromosomes used in crossover
    pool = parents[randp_idx,:] # pool of chromosomes used in crossover
    
    reprod_idx = []
    for idx in par_idx:
        if idx not in randp_idx:
            reprod_idx.append(idx)
    reprod_idx = np.array(reprod_idx)
    reprod = parents[reprod_idx,:] # remaining chromosomes not used in crossover to be passed unchanged to next gen
    
    
    cross_idx = list(np.arange(pool.shape[0])) # indices of pool to select from
    slice_idx = list(np.arange(1,pool.shape[1])) # indices of genes to select from, can't select zero which would be
                                                 # same as reproduction
    num_off = pop.shape[0]-reprod.shape[0] # number of chromosomes obtained from crossover (N-reprod)
    assert num_off%2 == 0
    
    children1 = np.zeros((int(num_off/2),pop.shape[1]),dtype=int)
    children2 = np.zeros((int(num_off/2),pop.shape[1]),dtype=int)
                    
    
    for i in range(int(num_off/2)):

        randc_idx = np.array(sample(cross_idx,2)) # pick 2 parents' indices randomly from pool for crossover
        parent1 = pool[randc_idx[0],:]
        parent2 = pool[randc_idx[1],:]
                     
        rands_idx = choice(slice_idx) # pick an index to slice genes
        
        child1 = np.concatenate((parent1[rands_idx:],parent2[-rands_idx:]))
        child2 = np.concatenate((parent1[:rands_idx],parent2[:-rands_idx]))
        
        children1[i,:] = child1
        children2[i,:] = child2
                     
    offspring = np.concatenate((children1,children2,reprod),axis=0)
        
    return offspring
    
    

    
    
    
    
def mutation(g,offspring):
    '''
    sequential mutation method 
    '''
#     pm = 0.01 # probability of mutation

    # constraint such that all genes in one chromosome are unique
    for i in range(offspring.shape[0]):
        chrom = offspring[i,:]
        u_idx = [x for x in df_idx if x not in chrom]
        mut = choice(u_idx)
        if g > (len(chrom)-1):
            g = g%len(chrom)
        chrom[g] = mut 
        
        if len(chrom) != len(np.unique(chrom)):
            new_chrom = np.unique(chrom)
            num_replace = len(chrom)-len(new_chrom)
            u_idx = [x for x in df_idx if x not in new_chrom]
            if num_replace == 1:
                rep_idx = np.array([choice(u_idx)])
            elif num_replace > 1:
                rep_idx = np.array(sample(u_idx,num_replace))
            chrom = np.concatenate((new_chrom,rep_idx))
            
        offspring[i,:] = chrom
           
            
    return offspring







def ga(G,N,k):
    '''
    '''
    start = time()
    
    fit_gen = np.zeros((G+1))

    pop = initialise(N,k)
    for gen in range(G):
        f = fitness(pop)
        print(f'generation {gen}: done')
        fit_gen[gen] = np.mean(f)
        parents = selection(pop,f)
        offspring = crossover(parents,pop)
        pop = mutation(gen,offspring)
    
    f = fitness(pop)
    print(f'generation {G}: done')
    fit_gen[-1] = np.mean(f)
    f_max = max(f)
    idx = np.where(f==f_max)[0][0]
    sol = pop[idx,:]
    
    
    end = time()
    print('\nrun time: ',end-start)
    
    
    
    return sol,f_max,fit_gen      
### implement termination condition: average fitness value becomes constant