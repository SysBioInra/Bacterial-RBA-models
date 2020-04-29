"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import rba


def main():
    
    # 'model/' or '../Bacillus-subtilis-168-WT/' or
    # '../Escherichia-coli-K12-WT'
    xml_dir = 'model/'
    output_dir = 'simulation/'
    
    # load model, build matrices and solve
    model = rba.RbaModel.from_xml(xml_dir)
    results = model.solve()
    
    # print µ_max and yield 
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    # flux in mmol g_bm^-1 h^-1 needs to be converted to g substrate
    # fructose: MW = 180.16 g/mol = 0.18 g/mmol
    yield_frc = results.mu_opt / (results.reaction_fluxes()['R_FRUpts2'] * 0.18)
    print('Yield on substrate is {}.'.format(yield_frc))
    
    # export results
    results.write(output_dir)
    
    
    # also export
    results.write_fluxes(
        output_dir + 'fluxes.csv',
        file_type = 'csv',
        merge_isozyme_reactions = True,
        only_nonzero = True,
        remove_prefix = True)


if __name__ == '__main__':
    main()
