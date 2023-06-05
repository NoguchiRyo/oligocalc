import argparse,sys
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    args = sys.argv[1:]
    parser = argparse.ArgumentParser()

    parser.add_argument('input_excel_path',type=str)
    parser.add_argument('-e','--ec-excel-path',type=str,default=Path(__file__).resolve().parent.parent.joinpath("coefficients.xlsx"))
    parser.add_argument('-w','--mw-excel-path',type=str,default=Path(__file__).resolve().parent.parent.joinpath("mass.xlsx"))
    parser.add_argument('-o','--output-path',type=str,default='output.xlsx')
    args = parser.parse_args()

    df = pd.read_excel(args.ec_excel_path)
    stack_of_monomer = list(df['stack of monomer'])
    extinction_coefficient = list(df['extinction coefficient'])
    coefficient_dict = dict(zip(stack_of_monomer,extinction_coefficient))

    def calc_ec(sequence):
        neighbors = [coefficient_dict[str(sequence[i:i+2])] for i in range(0,len(sequence)-1)]
        individuals = [coefficient_dict[i] for i in sequence[1:len(sequence)-1]]
        lpermolpercm = np.sum(neighbors) - np.sum(individuals)
        return lpermolpercm

    df = pd.read_excel(args.input_excel_path)
    names = list(df['Name'])
    sequences = list(df['Sequence'])
    extinction_coefficients = [calc_ec(s) for s in sequences]
    nmol_per_od = [1/e*1e+6 for e in extinction_coefficients]

    mass = pd.read_excel(args.mw_excel_path)
    mass_dict = dict(zip(list(mass['nucleotide']),list(mass['molecular weight'])))
    def calc_mass(sequence,five_phos,three_phos):
        weight = np.sum([mass_dict[b]-18 for b in sequence])+17
        if not five_phos:
            weight -= 80
        if three_phos:
            weight += 80
        return weight
    
    mass_list =  [calc_mass(s,True,False) for s in sequences]

    df['nmol/OD'] = nmol_per_od
    df['M.W.'] = mass_list
    df['μg/OD'] = np.array(mass_list)*np.array(nmol_per_od)/1000 #ug/OD

    df.to_excel(args.output_path)

if __name__ == '__main__':
    main()
