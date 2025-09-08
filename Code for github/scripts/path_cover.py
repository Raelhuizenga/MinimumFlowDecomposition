import networkx as nx
import argparse
import flowpaths as fp
import numpy as np
import os
from copy import deepcopy
import time
import json

## In flow_paths container:
# apptainer shell \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flowpaths:/mnt/flowpaths \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data:/mnt/data \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths:/mnt/run_flowpaths \
#   containers/flow_paths/flow-paths.sif


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
    parser.add_argument('--strain', type=int, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--vgflow_constraints', type=str2bool, nargs='?', const=True, default=False)
    parser.add_argument('--subpaths', type=str2bool, nargs='?', const=True, default=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--input_dir', type=str, required=False)
    args = parser.parse_args()

    strain = args.strain
    seed = args.seed
    constraints = args.vgflow_constraints
    subpaths = args.subpaths
    output_dir = args.output_dir

    fp.utils.configure_logging(
        level=fp.utils.logging.WARNING,
        log_to_console=False,
    )

    graph_directory = '/mnt/data/simulation_data/simulation_data/seqwish_graphs/HCV_frequency'
    if args.input_dir:
        graph_directory = args.input_dir

    graph_file = f'{graph_directory}/{strain}-strain_seed_{seed}/final_graph.gfa'
    abundance_file = f'{graph_directory}/{strain}-strain_seed_{seed}/final_node_abundances.txt'

    if not os.path.exists(graph_file):
        graph_file = f'{graph_directory}/{strain}-strain_seed_{seed}/mod_graph.gfa'
        abundance_file = f'{graph_directory}/{strain}-strain_seed_{seed}/node_abundances.txt'

    result = {
        "strain": strain,
        "seed": seed,
    }

    try:
        start = time.time()
        found_paths = process_graph(graph_file, abundance_file, strain, seed, constraints, subpaths)
        elapsed = time.time() - start

        result.update({
            "found_paths": found_paths,
            "runtime_seconds": round(elapsed, 3)
        })

    except Exception as e:
        result.update({
            "error": str(e),
            "found_paths": None,
            "runtime_seconds": None
        })

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'result_strain_{strain}_seed_{seed}.json')
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"Result saved to {output_file}")



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


def round_flow(flow, allowed_values):
    return min(allowed_values, key=lambda w: abs(flow - w))


def process_graph(graph_file, abundance_file, strain, seed, constraints, subpaths):
    print('Start processing graph')
    graph, subpath_constraints_nodes = convert_gfa_to_networkx(graph_file, abundance_file)
    print('GFA to network conversion successfull')
    neGraph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")
    neGraph.graph["id"] = f'graph_strain_{strain}_seed_{seed}'
    if subpaths:
        print("Number of subpath constraints:", len(neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes)))
        mpc_model = fp.MinPathCover(
            neGraph,
            subpath_constraints=neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes),
            elements_to_ignore=neGraph.edges_to_ignore,
            strict_subpath_constraints=constraints,
        )
    else:
        mpc_model = fp.MinPathCover(
            neGraph,
            elements_to_ignore=neGraph.edges_to_ignore)
    mpc_model.solve()
    solution = mpc_model.get_solution()
    print("solution:", solution)
    k = solution["paths"]

    print("Ground truth paths:", strain)
    print("Paths found:", len(k))
    return len(k)


if __name__ == "__main__":
    main()
