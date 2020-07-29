"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import rba
import re
import copy
import numpy as np
import pandas as pd


def main():
    
    # set input and output paths
    xml_dir = 'model/'
    output_dir = 'simulation/mixotrophy/'
    
    # load model, build matrices
    model = rba.RbaModel.from_xml(xml_dir)
    
    # optionally modify medium
    orig_medium = model.medium
    orig_medium['M_fru'] = 0.0
    
    # optionally add additional flux constraints
    #set_flux_boundary(model, 'R_SUCCt2_2', 0.0)
    #set_flux_boundary(model, 'R_SUCCt_3', 0.0)
    #set_flux_boundary(model, 'R_SUCCtr', 0.0)
    
    # A) simulation for different substrates
    substrate = pd.read_csv('simulation/substrate_mixotrophy.csv')
    simulate_substrate(model, substrate, orig_medium, output_dir)
    
    # B) simulation for different k_apps
    #iterations = 200
    #simulate_variability(model, iterations, orig_medium, output_dir)


def simulate_substrate(model, substrate, orig_medium, output_dir):
    
    # run several simulations in a loop
    for index, row in substrate.iterrows():
        
        # add desired substrate concentration to minimal medium
        model2 = copy.deepcopy(model)
        new_medium = orig_medium.copy()
        new_medium[row['carbon_source']] = row['carbon_conc']
        new_medium[row['nitrogen_source']] = row['nitrogen_conc']
        model2.medium = new_medium
        
        # optionally set flux boundary for reactions
        if 'substrate_uptake' in row.index:
            set_flux_boundary(model2, row['substrate_TR'], row['substrate_uptake'])
        # force flux through Rubisco, from 0 to 5 mmol/gDCW
        set_flux_boundary(model2, 'R_RBPC', float(index))
        
        # solve model
        try:
            result = model2.solve()
            # report results; for yield calculation supply transport
            # reaction and MW of substrate
            report_results(result,
                output_dir = output_dir,
                output_suffix = '_{}_{}_{}_{}_{}.tsv'.format(*row.to_list()[0:4], index),
                substrate_TR = row['substrate_TR'],  
                substrate_MW = row['substrate_MW']
                )
        except TypeError:
            print('model not solvable due to matrix inconsistency')


def simulate_variability(model, iterations, orig_medium, output_dir):
    
    # assign medium, and iterator
    model.medium = orig_medium
    completed_cycles = 1
    
    # run several simulations in a loop
    while completed_cycles <= iterations:
        
        # randomly sample kapp values
        model2 = copy.deepcopy(model)
        randomize_efficiency(model2, log10_mean = 4, log10_sd = 1.06)
        
        # solve model
        try:
            result = model2.solve()
            # report results; for yield calculation supply transport
            # reaction and MW of substrate
            if result.mu_opt > 0:                
                report_results(result,
                    output_dir = output_dir,
                    output_suffix = '_iteration_{:0>3}.tsv'.format(completed_cycles),
                    substrate_TR = 'R_FORt',
                    substrate_MW = 0.04603
                    )
                completed_cycles = completed_cycles + 1
            else:
                print('growth rate is zero, discarding result')
        except TypeError:
            print('model not solvable due to matrix inconsistency')


def report_results(
    result, output_dir, output_suffix,
    substrate_TR = None, substrate_MW = None):
    
    # calculate yield
    # flux in mmol g_bm^-1 h^-1 needs to be converted to g substrate
    # using transporter TR and molecular weight MW
    if substrate_TR:
        yield_subs = result.mu_opt / (result.reaction_fluxes()[substrate_TR] * substrate_MW)
    
    # calculate compartment occupancy ('density status')
    ds_mem = result.density_status("Cell_membrane")
    ds_cyt = result.density_status("Cytoplasm")
    
    # export summary fluxes per reaction
    result.write_fluxes(
        output_dir + 'fluxes' + output_suffix,
        file_type = 'tsv',
        remove_prefix = True)
    
    # optionally re-export fluxes to comma separated values as well
    fluxes = pd.read_csv(output_dir + 'fluxes' + output_suffix, 
        sep = '\t', 
        index_col = 0, 
        header = None)
    fluxes.to_csv(output_dir + 'fluxes' + re.sub('tsv', 'csv', output_suffix))
    
    # export enzyme concentrations
    result.write_proteins(
        output_dir + 'proteins' + output_suffix,
        file_type = 'csv')
    
    # export growth rate, yield, and process machinery concentrations
    ma = result.process_machinery_concentrations()
    ma['P_ENZ'] = sum(result.enzyme_concentrations().values())
    ma['mu'] = result.mu_opt
    if substrate_TR:
        ma['yield'] = yield_subs
        ma['qS'] = result.reaction_fluxes()[substrate_TR]
    with open(output_dir + 'macroprocesses' + output_suffix, 'w') as fout:
        fout.write('\n'.join(['{}\t{}'.format(k, v) for k, v in ma.items()]))
    
    # print µ_max and yield to terminal
    print('\n----- SUMMARY -----\n')
    print('Optimal growth rate is {}.'.format(result.mu_opt))
    print('Yield on substrate is {}.'.format(yield_subs))
    print('Cell membrane occupancy is {} %.'.format(round(100*ds_mem[1]/ds_mem[0], 3)))
    print('Cytoplasm occupancy is {} %.'.format(round(100*ds_cyt[1]/ds_cyt[0], 3)))
    print('\n----- BOUNDARY FLUXES -----\n')
    for r in result.sorted_boundary_fluxes():
        print(r)


def set_flux_boundary(model, reaction, value):
    
    # construct target from input reaction ID and flux boundary
    boundary_id = reaction + '_flux_boundary'
    new_target = rba.xml.targets.TargetReaction(reaction)
    new_target.value = boundary_id
    
    # add target to model
    mp = model.targets.target_groups.get_by_id('metabolite_production')
    if reaction not in [i.reaction for i in mp.reaction_fluxes]:
        mp.reaction_fluxes.append(new_target)
        
        # add a new parameter with the actual value to model
        model.parameters.functions.append(
            rba.xml.Function(boundary_id, 'constant', {'CONSTANT': value})
        )
    else:
        # or update existing parameter
        fn = model.parameters.functions.get_by_id(boundary_id)
        fn.parameters.get_by_id('CONSTANT').value = value


def randomize_efficiency(
    model, log10_mean = 4, 
    log10_sd = 1):
    
    # get default efficiency
    fn = model.parameters.functions.get_by_id('default_efficiency')
    def_efficiency = fn.parameters.get_by_id('CONSTANT').value
    
    for e in model.enzymes.enzymes:
        
        # generate a random enzyme efficiency from a log normal distribution
        # (see O'Brien et al., PLOS Comp Bio, 2016)
        rand_eff = 10**((np.random.randn(1, )[0]*log10_sd)+log10_mean)
        
        for eff in ['forward_efficiency', 'backward_efficiency']:
            # 2 scenarios: either the enzyme has a defined enzyme efficiency
            # or it has the default efficiency
            # A) default efficiency: we add a new custom enzyme efficiency 
            # randomly sampled from default efficiency
            if getattr(e, eff) == 'default_efficiency':
                id_eff = e.id + '_efficiency'
                setattr(e, eff, id_eff)
                model.parameters.functions.append(
                    rba.xml.Function(
                        id_eff, 
                        'constant', 
                        {'CONSTANT': rand_eff}
                    )
                )
            
            # B) we permutate the existing efficiency (excluding transporters)
            else:
                # in case the efficiency is an aggregate, first need to identify eff ID
                agg = model.parameters.aggregates.get_by_id(getattr(e, eff))
                if agg:
                    is_transporter = ('default_transporter_efficiency' in 
                        [i.function for i in agg.function_references])
                    id_eff = e.id + '_' + eff
                else:
                    is_transporter = getattr(e, eff) == 'default_transporter_efficiency'
                    id_eff = getattr(e, eff)
                if not is_transporter:
                    fn = model.parameters.functions.get_by_id(id_eff)
                    #orig_eff = fn.parameters.get_by_id('CONSTANT').value
                    fn.parameters.get_by_id('CONSTANT').value = rand_eff


if __name__ == '__main__':
    main()
