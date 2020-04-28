"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import rba


def main():
    
    xml_dir = "model/"
    output_dir = "simulation/"
    
    # load model, build matrices and solve
    model = rba.RbaModel.from_xml(xml_dir)
    results = model.solve()
    
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(output_dir)


if __name__ == '__main__':
    main()
