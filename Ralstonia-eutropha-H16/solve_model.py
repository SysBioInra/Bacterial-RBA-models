"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import rba
import numpy as np


def main():
    
    # set input and output paths
    xml_dir = 'model/'
    output_dir = 'simulation/'
    
    # load model, build matrices
    model = rba.RbaModel.from_xml(xml_dir)
    
    # optionally modify medium
    c_fruc = [round(10**i, 5) for i in np.arange(-3, -0.875, 0.125)]#[0.001, 0.01, 0.1, 1]
    new_medium = model.medium
    
    # loop through a set of conditions
    for conc in c_fruc:
        new_medium['M_fru'] = conc
        model.medium = new_medium
        
        # solve model
        result = model.solve()
        
        # report results; for yield calculation supply transport
        # reaction and MW of substrate
        report_results(result,
            output_dir = output_dir,
            output_suffix = '_fru_' + str(conc) + '.tsv',
            substrate_TR ='R_FRUabc',  
            substrate_MW = 0.18016
            )


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


if __name__ == '__main__':
    main()
