import networkx as nx
import argparse
import flowpaths as fp
import numpy as np
import os
import time
from sklearn.cluster import KMeans
from kneed import KneeLocator  
from copy import deepcopy
import json

## In flow_paths container:
# apptainer shell \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flowpaths:/mnt/flowpaths \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data:/mnt/data \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths:/mnt/run_flowpaths \
#   containers/flow_paths/flow-paths.sif

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strain', type=int, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--graph_dir', type=str, default='/mnt/data/simulation_data/simulation_data/seqwish_graphs/random_abundances/wrong_nodes_removed')
    parser.add_argument('--k_dir', type=str, default='/mnt/run_flowpaths/output/estimate_k/HCV_nodes_removed')
    parser.add_argument('--out_dir', type=str, default='/mnt/run_flowpaths/output/k-means_clustering/zero_weight')
    args = parser.parse_args()

    # graph_directory = '/mnt/data/simulation_data/simulation_data/seqwish_graphs/HIV_8000'
    strain = args.strain
    seed = args.seed
    graph_directory = args.graph_dir
    k_dir = args.k_dir
    out_dir = args.out_dir
    fp.utils.configure_logging(
        level=fp.utils.logging.DEBUG,  # Or ERROR
        log_to_console=True,            # Disable console output entirely
    )

    graph_file = (f'{graph_directory}/{strain}-strain_seed_{seed}/mod_graph.gfa', f'{graph_directory}/{strain}-strain_seed_{seed}/node_abundances.txt')
    k_file = f'{k_dir}/result_strain_{strain}_seed_{seed}.json'
    out_file = f'{out_dir}/solution_{strain}_{seed}'

    process_graph(graph_file[0], graph_file[1], k_file, out_file)



def preprocess_weights_k_means(weights_array, min_k: int = 4, max_k: int = 50):
    distortions = []
    K_range = range(min_k, max_k + 1)
    
    for k in K_range:
        kmeans = KMeans(n_clusters=k, random_state=0).fit(weights_array)
        distortions.append(kmeans.inertia_)

    try:
        kneedle = KneeLocator(K_range, distortions, curve="convex", direction="decreasing")
        best_k = kneedle.elbow
        if best_k is None:
            raise ValueError("KneeLocator did not find a valid elbow.")
    except Exception as e:
        print(f"[Warning] Skipping clustering: {e}")
        # Use raw weights (rounded to int)
        raw_weights = np.unique(np.round(weights_array.flatten()).astype(int)).tolist()
        print(f"Using all unique raw weights: {raw_weights}")
        return raw_weights

    print(f"Optimal number of clusters (k): {best_k}")
    kmeans = KMeans(n_clusters=best_k, random_state=0).fit(weights_array)
    centers = np.round(kmeans.cluster_centers_.flatten()).astype(int)
    print(f"Cluster centers: {centers}")

    # # Remove the highest weight
    # centers = centers[centers != centers.max()]

    # Add a zero weights
    centers = np.append(centers, 0.0)  
    
    return centers.tolist()


def create_weights_array_from_abudance_file(abundance_file: str):
    node_weights = {}
    with open(abundance_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(':')
            node_weights[line[0]] = float(line[1])

    weights_array = np.array(list(node_weights.values())).reshape(-1, 1)
    return weights_array


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



def create_weights_array_from_graph(graph):
    weights = []
    for _, _, data in graph.edges(data=True):
        if 'flow' in data:
            weights.append(data['flow'])
    return np.array(weights).reshape(-1, 1)



def process_graph(graph_file, abundance_file, k_file, outfile):
    # w = create_weights_array_from_abudance_file(abundance_file)
    # new_weights = preprocess_weights_k_means(w, min_k=4, max_k=100)
    # print("Weights:", new_weights)
    graph, subpath_constraints_nodes = convert_gfa_to_networkx(graph_file, abundance_file)

    # We create a node expanded graph, where the weights are taken from the attribute "flow"
    neGraph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")

    correction_model = fp.MinErrorFlow(
        neGraph,
        flow_attr="flow",
        weight_type=int,
        elements_to_ignore=neGraph.edges_to_ignore,
    )
    correction_model.solve()
    corrected_graph = correction_model.get_corrected_graph()

    w = create_weights_array_from_graph(corrected_graph)
    new_weights = preprocess_weights_k_means(w, min_k=4, max_k=100)
    print("Weights:", new_weights)



    with open(k_file) as f:
        data = json.load(f)
        k = int(data["found_paths"])

  
    print("_______ Starting optimization _________")


    start_time = time.time()
    klae_model = fp.kLeastAbsErrors(
            neGraph,
            k=k,
            flow_attr="flow",
            elements_to_ignore=neGraph.edges_to_ignore,
            subpath_constraints=neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes),
            strict_subpath_constraints=False,
            solution_weights_superset=new_weights,
            optimization_options = {"allow_empty_paths": True},
            solver_options={"external_solver": "gurobi", "threads": 8, "log_to_console": "true"}
            )
    klae_model.solve()
    elapsed = start_time - time.time()

        # Modify the data
    data["time_klae_superset"] = round(elapsed, 3)

    # Write it back to the file
    with open(k_file, "w") as f:
        json.dump(data, f, indent=2)

    process_expanded_solution(graph, neGraph, klae_model, outfile)
    
    

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

    else:
        print("Model could not be solved.")


if __name__ == "__main__":
    main()