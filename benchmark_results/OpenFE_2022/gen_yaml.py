import yaml
import pandas as pd
from openfe_benchmarks import tyk2, p38, mcl1, hif2a, thrombin, tnsk2


def write_network_with_results(network, calc_df, name):
    output = {}
    
    output['info'] = "OpenFE 2022 PLB results, LOMAP + MST graph"
    output['pdb'] = f"../../openfe_benchmarks/data/{name}_protein.pdb"
    output['sdf'] = f"../../openfe_benchmarsk/data/{name}_ligands.sdf"
    
    edges_dict = dict()
    
    for edge in network.edges:
        nameA = edge.componentA.name
        nameB = edge.componentB.name
        edge_name = f'edge_{edge.componentA.name}_{edge.componentB.name}'
        col_keys = ['Ligand1', 'Ligand2']
        col_values = [nameA, nameB]
        results_df = calc_df[(calc_df[col_keys] == col_values).all(1)]
        print(results_df, col_keys, col_values, results_df.empty)
        if results_df.empty:
            continue
        
        edges_dict[edge_name] = {
            'ligand_a': nameA,
            'ligand_b': nameB,
            'atom mapping': edge.componentA_to_componentB,
            'DDG': float(results_df['calc_DDG']),
            'dDDG': float(results_df['calc_dDDG(additional)']),}
        
    output['edges'] = edges_dict
    
    with open(f"{name}_results.yaml", 'w') as yaml_file:
        yaml.dump(output, yaml_file, default_flow_style=None, sort_keys=False)


systems = {
    'tyk2': tyk2.get_system(),
    'p38': p38.get_system(),
    'mcl1': mcl1.get_system(),
    'hif2a': hif2a.get_system(),
    'thrombin': thrombin.get_system(),
    'tnsk2': tnsk2.get_system(),
}

calc_dfs = {
    'tyk2': pd.read_csv('csv-files/tyk2-calc.csv', dtype=str),
    'p38': pd.read_csv('csv-files/p38-calc.csv', dtype=str),
    'mcl1': pd.read_csv('csv-files/mcl1-calc.csv', dtype=str),
    'hif2a': pd.read_csv('csv-files/hif2a-calc.csv', dtype=str),
    'thrombin': pd.read_csv('csv-files/thrombin-calc.csv', dtype=str),
    'tnsk2': pd.read_csv('csv-files/tnsk2-calc.csv', dtype=str),
}


for var in systems.keys():
    write_network_with_results(systems[var].ligand_network, calc_dfs[var], name=var)
