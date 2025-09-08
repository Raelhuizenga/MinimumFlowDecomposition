import networkx as nx
import argparse
import flowpaths as fp
import numpy as np
import os
import time
from copy import deepcopy

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
    args = parser.parse_args()
    # graph_directory = '/mnt/data/simulation_data/simulation_data/seqwish_graphs/random_abundances'
    graph_directory = '/mnt/data/simulation_data/simulation_data/seqwish_graphs/HIV_8000_dagify'
    strain = args.strain
    seed = args.seed
    fp.utils.configure_logging(
        level=fp.utils.logging.DEBUG,  # Or ERROR
        log_to_console=True,            # Disable console output entirely
    )

    # graph_file = (f'{graph_directory}/{strain}-strain_seed_{seed}/mod_graph.gfa', f'{graph_directory}/{strain}-strain_seed_{seed}/node_abundances.txt',)
    graph_file = (f'{graph_directory}/{strain}-strain_seed_{seed}/final_graph.gfa', f'{graph_directory}/{strain}-strain_seed_{seed}/final_node_abundances.txt',)

    process_graph(graph_file[0], graph_file[1], strain, seed)


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



def process_graph(graph_file, abundance_file, strain, seed):
    
    graph, subpath_constraints_nodes = convert_gfa_to_networkx(graph_file, abundance_file)

    # We create a node expanded graph, where the weights are taken from the attribute "flow"
    neGraph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")

    fp.utils.draw(
            G=neGraph, 
            filename=f"/mnt/run_flowpaths/output/graphs/{strain}_strain_seed_{seed}_graph_with_errors_version_29_07.pdf", 
            flow_attr="flow", 
            subpath_constraints = neGraph.get_expanded_subpath_constraints(subpath_constraints_nodes),
            draw_options={
            "show_graph_edges": True,
            "show_edge_weights": True,
            "show_path_weights": False,
            "show_path_weight_on_first_edge": False,
            "pathwidth": 2,
        })


if __name__ == "__main__":
    main()