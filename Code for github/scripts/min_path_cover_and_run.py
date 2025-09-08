import networkx as nx
import argparse
import flowpaths as fp
import numpy as np
from copy import deepcopy
import json
import time
import os


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
    parser.add_argument('--graph_dir', type=str, required=True)
    parser.add_argument('--k_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--vgflow_constraints', type=str2bool, nargs='?', const=True, default=False)
    parser.add_argument('--subpaths', type=str2bool, nargs='?', const=True, default=True)
    parser.add_argument('--set_weights_dir', type=str, nargs='?', const=True, default=None)

    args = parser.parse_args()
    graph_directory = args.graph_dir # '/mnt/data/simulation_data/simulation_data/seqwish_graphs/HCV_frequency'
    estimate_k_directory = args.k_dir # '/mnt/run_flowpaths/output/estimate_k/HCV_frequency'
    output_dir = args.output_dir # /mnt/run_flowpaths/output/min_cover_HCV_frequency
    strain = args.strain
    seed = args.seed
    constraints = args.vgflow_constraints
    subpaths = args.subpaths
    set_weights_dir = args.set_weights_dir
    
    fp.utils.configure_logging(
        level=fp.utils.logging.WARNING,  # Or ERROR
        log_to_console=True,            # Disable console output entirely
    )

    graph = f'{graph_directory}/{strain}-strain_seed_{seed}/final_graph.gfa'
    abundance = f'{graph_directory}/{strain}-strain_seed_{seed}/final_node_abundances.txt'
    k_file = f'{estimate_k_directory}/result_strain_{strain}_seed_{seed}.json'

    if not os.path.exists(graph):
        graph = f'{graph_directory}/{strain}-strain_seed_{seed}/mod_graph.gfa'
        abundance = f'{graph_directory}/{strain}-strain_seed_{seed}/node_abundances.txt'

    _, k = get_k(k_file)
    print(f"Found number of paths: {k}")
    if set_weights_dir:
        set_weights = get_set_of_allowed_weights(set_weights_dir, k)
    else:
        set_weights = None
    print(f"Set of allowed weights: {set_weights}")
    process_graph(graph, abundance, k_file, strain, seed, constraints, subpaths, output_dir, set_weights)



def get_set_of_allowed_weights(set_weights_dir, k):
    weights = []
    with open(f'{set_weights_dir}/ground_truth_abundances.txt', 'r') as f:
        for line in f:
            if line.startswith('>') and 'weight:' in line:
                try:
                    # Extract weight part (e.g., 'weight: 624.0x')
                    weight_str = line.split('weight:')[1].strip()
                    # Remove trailing 'x' and convert to float
                    weight = float(weight_str.rstrip('x').strip())
                    weights.append(weight)
                except (IndexError, ValueError):
                    print(f"Warning: Could not parse weight from line: {line.strip()}")
    
    # Pad with zeros if fewer than k elements
    if len(weights) < k:
        weights += [0.0] * (k - len(weights))
    
    return weights



def get_k(k_file):
    if not os.path.exists(k_file):
        raise FileNotFoundError(f"File not found: {k_file}")

    with open(k_file) as f:
        try:
            data = json.load(f)
            if "found_paths" not in data:
                raise KeyError(f"'found_paths' key not found in {k_file}")
            k = int(data["found_paths"])
        except json.JSONDecodeError as e:
            raise ValueError(f"Failed to parse JSON in {k_file}: {e}")
    return data, k


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



def process_graph(graph_file, abundance_file, k_file, strain, seed, constraints, subpaths, output_dir, set_weights):
    graph, subpath_constraints_nodes = convert_gfa_to_networkx(graph_file, abundance_file)
    data, k = get_k(k_file)
    

    # We create a node expanded graph, where the weights are taken from the attribute "flow"
    neGraph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")
    neGraph.graph["id"] = f'graph_strain_{strain}_seed_{seed}'


    if subpaths:

        start = time.time()

        if set_weights:
            print('running with predefined set of weights')

            klae_model = fp.kLeastAbsErrors(
                neGraph,
                k=k,
                flow_attr="flow",
                elements_to_ignore=neGraph.edges_to_ignore,
                subpath_constraints=neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes),
                strict_subpath_constraints=constraints,
                solver_options={"external_solver": "gurobi", "threads": 8, "log_to_console": "true"},
                solution_weights_superset=set_weights,
                optimization_options = {"allow_empty_paths": True},
                )
            klae_model.solve()
        else:
            klae_model = fp.kLeastAbsErrors(
                neGraph,
                k=k,
                flow_attr="flow",
                elements_to_ignore=neGraph.edges_to_ignore,
                subpath_constraints=neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes),
                strict_subpath_constraints=constraints,
                solver_options={"external_solver": "gurobi", "threads": 8, "log_to_console": "true"},
                )
            klae_model.solve()
        elapsed = start - time.time()
    else:
        start = time.time()

        klae_model = fp.kLeastAbsErrors(
            neGraph,
            k=k,
            flow_attr="flow",
            elements_to_ignore=neGraph.edges_to_ignore,
            solver_options={"external_solver": "gurobi", "threads": 8, "log_to_console": "true"}
            )
        klae_model.solve()

        elapsed = start - time.time()

        # Modify the data
    data["time_klae"] = round(elapsed, 3)

    # Write it back to the file
    with open(k_file, "w") as f:
        json.dump(data, f, indent=2)

    output = f'{output_dir}/solution_{strain}_{seed}'
    
    process_expanded_solution(graph, neGraph, klae_model, output)
    
    

def process_expanded_solution(graph: nx.DiGraph, neGraph: fp.NodeExpandedDiGraph, model: fp.MinFlowDecomp, output_fasta: str):
    if model.is_solved():
        solution = model.get_solution()
        expanded_paths = solution["paths"]
        original_paths = neGraph.get_condensed_paths(expanded_paths)
        print("Weights:", solution["weights"])

        # Build sequences for each path
        path_sequences = []
        for path in original_paths:
            seq = "".join(graph.nodes[node]["sequence"] for node in path if "sequence" in graph.nodes[node])
            path_sequences.append(seq)

        # Save to a FASTA file
        with open(output_fasta, "w") as fasta_file:
            for idx, (seq, weight) in enumerate(zip(path_sequences, solution["weights"]), start=1):
                fasta_file.write(f">Genome {idx}, weight: {weight}x\n")
                fasta_file.write(f"{seq}\n")
        print(f"Saved {len(path_sequences)} sequences to {output_fasta}")


        # fp.utils.draw_solution(
        #     G=graph, 
        #     filename="expanded_graph_17_clusters.pdf", 
        #     flow_attr="flow", 
        #     paths=original_paths,
        #     weights=solution["weights"],
        #     draw_options={
        #         "show_node_weights": True,
        #     })
    else:
        print("Model could not be solved.")


if __name__ == "__main__":
    main()