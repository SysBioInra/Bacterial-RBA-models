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
    output_dir = 'simulation/variability_analysis/fructose/'
    
    # load model, build matrices
    model = rba.RbaModel.from_xml(xml_dir)
    
    # optionally modify medium
    orig_medium = model.medium
    orig_medium['M_fru'] = 0
    substrate = pd.read_csv('simulation/substrate_variability.csv')
    
    # run several simulations in a loop
    for index, row in substrate.iterrows():
        
        # optionally randomly sample kapp values
        model2 = copy.deepcopy(model)
        randomize_efficiency(model2, log10_boundary = 2)
        
        # add target substrate concentration to minimal medium
        new_medium = orig_medium.copy()
        new_medium[row['carbon_source']] = row['carbon_conc']
        new_medium[row['nitrogen_source']] = row['nitrogen_conc']
        model2.medium = new_medium
        
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


def report_results(
    result, output_dir, output_suffix,
    substrate_TR = None, substrate_MW = None):
    
    # calculate yield
    # flux in mmol g_bm^-1 h^-1 needs to be converted to g substrate;
    # for fructose: MW = 180.16 g/mol = 0.18 g/mmol
    if substrate_TR:
        yield_subs = result.mu_opt / (result.reaction_fluxes()[substrate_TR] * substrate_MW)
    
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
    with open(output_dir + 'macroprocesses' + output_suffix, 'w') as fout:
        fout.write('\n'.join(['{}\t{}'.format(k, v) for k, v in ma.items()]))
    
    # print µ_max and yield to terminal
    print('\n----- SUMMARY -----\n')
    print('Optimal growth rate is {}.'.format(result.mu_opt))
    print('Yield on substrate is {}.'.format(yield_subs))
    print('\n----- BOUNDARY FLUXES -----\n')
    for r in result.sorted_boundary_fluxes():
        print(r)


def randomize_efficiency(model, log10_boundary = 1):
    
    # get default efficiency
    fn = model.parameters.functions.get_by_id('default_efficiency')
    def_efficiency = fn.parameters.get_by_id('CONSTANT').value
    
    for e in model.enzymes.enzymes:
        
        # generate a random factor from a log continuous interval of 10^-1 to 10^1
        rand_factor = 10**(log10_boundary*2*np.random.random_sample((1, ))[0]-log10_boundary)
        
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
                        {'CONSTANT': def_efficiency * rand_factor}
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
                    orig_eff = fn.parameters.get_by_id('CONSTANT').value
                    fn.parameters.get_by_id('CONSTANT').value = orig_eff * rand_factor


if __name__ == '__main__':
    main()
