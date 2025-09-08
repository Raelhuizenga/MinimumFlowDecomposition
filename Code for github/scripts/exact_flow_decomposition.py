import networkx as nx
import argparse
import flowpaths as fp
import numpy as np
import csv
import os
from copy import deepcopy
import time
import json


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vgflow_constraints', type=str2bool, nargs='?', const=True, default=False)
    parser.add_argument('--strain', type=int, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--length', type=int, nargs='?', const=True, default=2000)
    parser.add_argument('--graph_directory', type=str, nargs='?', const=True, default='/mnt/data/simulation_data/gridsearch_haplotypes')
    parser.add_argument('--output_directory', type=str, nargs='?', const=True, default='/mnt/run_flowpaths/output/haplotype_gridsearch')
    args = parser.parse_args()

    constraints = args.vgflow_constraints
    strain = args.strain
    seed = args.seed
    graph_directory = args.graph_directory
    output_folder = args.output_directory
    length = args.length
    os.makedirs(output_folder, exist_ok=True)

    fp.utils.configure_logging(
        level=fp.utils.logging.WARNING,
        log_to_console=False,
    )

    print(f"Running for strain={strain}, seed={seed}")
    graph_file = f'{graph_directory}/graph_{length}_{strain}_{seed}.gfa'
    abundance_file = f'{graph_directory}/abundances_{length}_{strain}_{seed}.txt'

    try:
        runtime, nodes, edges, k = process_graph(graph_file, abundance_file, strain, seed, constraints)
        with open(f"{output_folder}/solution_{strain}_{seed}{length}.json", "w") as f:
            json.dump({ "time": runtime, "num_vertices": nodes, "num_edges": edges, "paths": k }, f)
    except Exception as e:
        print(f"Failed for strain={strain}, seed={seed}: {e}")



def convert_gfa_to_networkx(gfa_file, abundance_file):
    network = nx.DiGraph()
    network.graph["id"] = gfa_file
    subpaths = []
    with open(abundance_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(':')
            network.add_node(line[0], flow=float(line[1]))

    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'L':
                network.add_edge(line[1], line[3])
            if line[0] == 'P':
                path_nodes = [node.strip('+') for node in line[2].split(',')]
                if path_nodes[0] != path_nodes[0].strip('-'):
                    path_nodes = [node.strip('-') for node in line[2].split(',')][::-1]
                subpaths.append(path_nodes)
            if line[0] == 'S':
                network.nodes[line[1]]['sequence'] = line[2]
    return network, subpaths


def process_graph(graph_file, abundance_file, strain, seed, constraints):
    graph, subpath_constraints_nodes = convert_gfa_to_networkx(graph_file, abundance_file)

    neGraph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")
    neGraph.graph["id"] = f'graph_strain_{strain}_seed_{seed}'
    efd_model = fp.MinFlowDecomp(
            neGraph,
            flow_attr="flow",
            subpath_constraints=neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes),
            elements_to_ignore=neGraph.edges_to_ignore,
            strict_subpath_constraints=constraints,
            # solver_options = {"time_limit": 14400}
        )
    efd_model.solve()
    solution = efd_model.get_solution()
    time = efd_model.solve_statistics["mfd_solve_time"]
    k = solution["paths"]
    return time, neGraph.number_of_nodes(), neGraph.number_of_edges(), k


if __name__ == "__main__":
    main()
