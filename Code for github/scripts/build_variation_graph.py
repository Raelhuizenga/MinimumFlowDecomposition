import argparse
import subprocess
import json
import os
import sys
from collections import defaultdict
from graph_tool.all import Graph


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('contigfile', type=str, help='The file containing the contigs, format is .fasta')
    parser.add_argument('forward_reads', type=str, help='The file containing the forward reads, format is .fastq')
    parser.add_argument('reverse_reads', type=str, help='The file containing the reverse reads, format is .fastq')
    parser.add_argument('output_directory', type=str, help='Directory where the output should be saved')
    parser.add_argument('-g', '--gfa_file', dest='gfa_file', type=str, default=None, help="Location of the GFA file, if not provided a new graph is created")
    parser.add_argument('-vg', '--vg_directory', type=str, help="Location of the VG executable", default='vg')
    args = parser.parse_args()

    # Set the VG executable path
    vg_path = args.vg_directory 

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    if args.gfa_file is None:
        all_to_all_alignment = align_all_contigs(args.contigfile, args.output_directory)
        graph = build_variation_graph(args.contigfile, all_to_all_alignment, args.output_directory, vg_path)
        simple_graph = simplify_graph(graph, args.output_directory, vg_path)
        gfa_file_final = convert_vg_to_gfa(simple_graph, vg_path)
        read_mapping = map_reads_to_graph(simple_graph, args.forward_reads, args.reverse_reads, args.output_directory, vg_path)
    else:
        gfa_file = args.gfa_file
        read_mapping, gfa_file_final = map_reads_to_gfa_graph(gfa_file, args.forward_reads, args.reverse_reads, args.output_directory, vg_path)

    read_mapping_json = convert_gam_to_json(read_mapping, args.output_directory, vg_path)
    node_lengths = get_sequence_length_nodes(gfa_file_final)
    log_path = os.path.join(args.output_directory, 'stats.log')
    print("saving abundances")
    with open(log_path, 'w') as f:
        sys.stdout = f
        save_abundances(read_mapping_json, node_lengths, args.output_directory)
        # Check if the read mapping was successful
        vg_stats_result = subprocess.run([vg_path, "stats", "-a", read_mapping],  capture_output=True, text=True, check=True)
        print("=== vg stats ===")
        print(vg_stats_result.stdout)
        contig_graph = read_gfa(gfa_file_final)
        print(f'HAS CYCLES: {cyclic(contig_graph)}')
        sys.stdout = sys.__stdout__ 
        print("Filtering low abundance nodes...")
        filter_low_abundance_nodes(
            gfa_file_final,
            os.path.join(args.output_directory, 'node_abundances.txt'),
            args.output_directory
        )
        print("Saved final_graph.gfa and final_node_abundances.txt in output directory.")


def read_gfa(gfa_file):
    """
    Reads a graph from a GFA-file and returns graph in gt-format.
    """
    print("Reading graph from {}".format(gfa_file))
    # Define a graph with its vertex properties
    g = Graph(directed=True)
    g.vp.seq = g.new_vertex_property('string')
    g.vp.contigs = g.new_vertex_property('vector<string>')
    g.ep.ori = g.new_edge_property('string')

    # read gfa and add vertices to graph
    node_dict_old2new = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'S':
                # vertex line
                node_id = line[1]
                v = g.add_vertex()
                assert node_id not in node_dict_old2new
                node_dict_old2new[node_id] = int(v)
                seq = line[2].upper()
                if len(seq) == 0:
                    print("WARNING: empty sequence in GFA")
                    print(line)
                g.vp.seq[v] = seq

    # parse through gfa again to add edges and contig paths to graph
    paths = {}
    path_count = 0
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'L':
                # edge line
                v1 = node_dict_old2new[line[1]]
                v2 = node_dict_old2new[line[3]]
                e = g.add_edge(v1, v2)
                g.ep.ori[e] = "{}{}".format(line[2], line[4])
            elif line[0] == 'P':
                # path line
                path_count += 1
                path_name = line[1]
                path_line = line[2].split(',')
                path = [(node_dict_old2new[node_inf[:-1]],
                         node_inf[-1]) for node_inf in path_line]
                if path_name not in paths.keys():
                    paths[path_name] = path
                else:
                    print("WARNING: sequence {} has multiple "
                          "alignments".format(path_name))
                    while path_name in paths.keys():
                        path_name = path_name + "*"
                    paths[path_name] = path
                for (node, ori) in path:
                    g.vp.contigs[node].append(path_name)

    print("Vertex count: {}".format(len(list(g.vertices()))))
    print("Edge count: {}".format(len(list(g.edges()))))
    print("Path count: {}".format(path_count))
    return g

def cyclic(graph):
    """Return True if the directed input graph has a cycle."""
    visited = set()
    path = [object()]
    path_set = set(path)
    stack = [graph.vertices()]
    while stack:
        for v in stack[-1]:
            if v in path_set:
                print(list(path_set))
                print(v)
                return True
            elif v not in visited:
                visited.add(v)
                path.append(v)
                path_set.add(v)
                stack.append(iter(v.out_neighbors()))
                break
        else:
            path_set.remove(path.pop())
            stack.pop()
    return False

def align_all_contigs(contigfile, output_directory):
    """Align all contigs against each other using minimap2."""
    output_file = f'{output_directory}/contigs_overlap.paf'
    cmd = ["minimap2", "-c", "-X", contigfile, contigfile, "--for-only"] # -c for CIGAR strings, 
    with open(output_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    return output_file  # Return the output file for further processing


def build_variation_graph(fasta_file, paf_file, output_directory, vg_path):
    """Use seqwish to construct a variation graph from the PAF alignment."""
    output_gfa = f'{output_directory}/initial_seqwish.gfa'
    cmd = ["seqwish", "-s", fasta_file, "-p", paf_file, "-g", output_gfa, "-k", "20", "-r",  "1", "-l", "1000"]
    subprocess.run(cmd, check=True)

    # Convert GFA to VG format
    graph_vg = f'{output_directory}/seqwish_graph.vg'
    with open(graph_vg, "w") as vg_file:
        subprocess.run([vg_path, "convert", output_gfa], stdout=vg_file, check=True)

    return graph_vg  # Return the VG file


def simplify_graph(graph_vg, output_directory, vg_path):
    """Simplify the variation graph using vg mod."""
    mod_graph = f'{output_directory}/mod_graph.vg'
    with open(mod_graph, "w") as out_file:
        subprocess.run([vg_path, "mod", "-X", "1024", "-n", "-c", "-s", "-u", "-U", "100", graph_vg], stdout=out_file, check=True)
    return mod_graph  # Return the pruned graph file


# def simplify_graph(graph_vg, output_directory, vg_path):
#     """Simplify a variation graph by DAGifying it with `vg unfold`, retaining paths."""

#     os.makedirs(output_directory, exist_ok=True)

#     # Output files
#     xg_index = os.path.join(output_directory, "graph.xg")
#     gcsa_index = os.path.join(output_directory, "graph.gcsa")
#     path_names_file = os.path.join(output_directory, "paths.txt")
#     dagified_graph = os.path.join(output_directory, "dag_with_paths.vg")

#     # Step 1: Create XG and GCSA indexes
#     subprocess.run([
#         vg_path, "index", "-x", xg_index, "-g", gcsa_index, "-k", "16", graph_vg
#     ], check=True)

#     # Step 2: Extract path names to a file
#     with open(path_names_file, "w") as f:
#         subprocess.run([
#             vg_path, "paths", "-x", xg_index, "-L"
#         ], stdout=f, check=True)

#     # Step 3: Use `vg unfold` to DAGify and retain paths
#     with open(dagified_graph, "w") as out_file:
#         subprocess.run([
#             vg_path, "unfold", "-x", xg_index, "-g", gcsa_index, "-p", "-P", path_names_file
#         ], stdout=out_file, check=True)

#     return dagified_graph  # Return path to the final simplified graph




def map_reads_to_graph(graph_vg, forward_reads, reverse_reads, output_directory, vg_path):
    """Map reads to the variation graph using vg map."""

    # Index the graph (XG and GCSA indices required for mapping)
    xg_index = f'{output_directory}/mod_graph.xg'
    gcsa_index = f'{output_directory}/mod_graph.gcsa'
    subprocess.run([vg_path, "index", graph_vg, "-x", xg_index, "-g", gcsa_index], check=True)

    # Map reads using vg map
    gam_output = f'{output_directory}/mapped_reads.gam'
    with open(gam_output, "w") as gam_file:
        subprocess.run(
            [vg_path, "map", "-d",  f'{output_directory}/mod_graph', "-f", forward_reads, "-f", reverse_reads, '--exclude-unaligned'], # -d basename
            stdout=gam_file,
            check=True
        )

    return gam_output  # Return the mapped reads file


def map_reads_to_gfa_graph(graph_gfa, forward_reads, reverse_reads, output_directory, vg_path):
    """Convert GFA to VG, chop, and map reads using vg map."""

    graph_vg = f"{output_directory}/mod_graph.vg"
    
    # Convert GFA to VG
    with open(graph_vg, "w") as vg_file:
        subprocess.run([vg_path, "convert", "-g", graph_gfa], stdout=vg_file, check=True)
    
    # Chop long nodes
    chopped_graph = f"{output_directory}/mod_graph_chopped.vg"
    with open(chopped_graph, "w") as chopped_file:
        subprocess.run([vg_path, "mod", "-X", "1024", graph_vg], stdout=chopped_file, check=True)

    # Save chopped graph as GFA
    chopped_gfa = f"{output_directory}/mod_graph_chopped.gfa"
    with open(chopped_gfa, "w") as gfa_file:
        subprocess.run([vg_path, "view", chopped_graph], stdout=gfa_file, check=True)

    # Index
    xg_index = f'{output_directory}/mod_graph.xg'
    gcsa_index = f'{output_directory}/mod_graph.gcsa'
    subprocess.run([vg_path, "index", chopped_graph, "-x", xg_index, "-g", gcsa_index], check=True)

    # Map reads
    gam_output = f'{output_directory}/mapped_reads.gam'
    with open(gam_output, "w") as gam_file:
        subprocess.run(
            [vg_path, "map", "-d", f'{output_directory}/mod_graph', "-f", forward_reads, "-f", reverse_reads, '--exclude-unaligned'],
            stdout=gam_file,
            check=True
        )

    return gam_output, chopped_gfa  # Return the mapped reads file and chopped graph



def convert_gam_to_json(gam_file, output_directory, vg_path):
    """Convert a GAM file to JSON format using vg view."""
    json_file = f'{output_directory}/mapped_reads.json'
    with open(json_file, "w") as out:
        subprocess.run([vg_path, "view", "-a", gam_file], stdout=out, check=True)
    return json_file  # Return the JSON file


def save_abundances(mapped_json, node_lengths, output_directory):
     """Save read counts per node to a text file."""
     output_file = f'{output_directory}/node_abundances.txt'
     node_basepairs = defaultdict(int)
     non_proper_paired = 0
     with open(mapped_json) as f:
        for line in f:
            alignment = json.loads(line)
            try:
                for mapping in alignment["path"]["mapping"]:
                    node_id = int(mapping["position"]["node_id"]) 
                    for edit in mapping["edit"]:
                        from_len = edit.get("from_length", 0) 
                        to_len = edit.get("to_length", 0)
                        aligned_bp = min(from_len, to_len)
                        node_basepairs[node_id] += aligned_bp  # Add to node count
            except KeyError:
                non_proper_paired += 1  # Count non-properly mapped reads
                continue
     # Print the number of non-properly mapped reads
     print(f"Number of non-properly mapped reads: {non_proper_paired}")

    # Compute normalized abundance (aligned basepairs per basepair in node)
     node_abundance = {
        node_id: aligned_bp / node_lengths[node_id]
        for node_id, aligned_bp in node_basepairs.items()
    }
    

    # Write output to a text file (tab-separated)
     with open(output_file, "w") as f:
        for node_id, abundance in sorted(node_abundance.items()):
            f.write(f"{node_id}:{abundance:.6f}\n")


def convert_vg_to_gfa(vg_file, vg_path):
    """Convert a VG file to GFA format using vg view."""
    gfa_file = f'{vg_file[:-2]}gfa'
    with open(gfa_file, "w") as out:
        subprocess.run([vg_path, "view", vg_file], stdout=out, check=True)
    return gfa_file  # Return the GFA file


def get_sequence_length_nodes(graph_gfa):
    """Get the sequence length of each node in the graph."""
    node_lengths = {}
    with open(graph_gfa, "r") as f:
        for line in f:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                node_id = parts[1]
                sequence = parts[2]
                node_lengths[int(node_id)] = len(sequence)
    return node_lengths  # Return a dictionary of node_id -> length

def filter_low_abundance_nodes(gfa_file, abundance_file, output_dir):
    """Filter out low-abundance nodes from a GFA graph."""
    os.makedirs(output_dir, exist_ok=True)

    def load_abundance(abundance_file):
        abundances = {}
        with open(abundance_file) as f:
            for line in f:
                node_id, val = line.strip().split(":")
                abundances[node_id] = float(val)
        return abundances

    def parse_gfa(gfa_file):
        headers = []
        segments = {}
        links = []
        paths = []

        with open(gfa_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts[0] == "H":
                    headers.append(line.strip())
                elif parts[0] == "S":
                    segments[parts[1]] = parts
                elif parts[0] == "L":
                    links.append(parts)
                elif parts[0] == "P":
                    paths.append(parts)
        return headers, segments, links, paths

    def write_gfa(outfile, headers, segments, links, paths):
        with open(outfile, "w") as f:
            for line in headers:
                f.write(line + "\n")
            for seg in segments.values():
                f.write("\t".join(seg) + "\n")
            for link in links:
                f.write("\t".join(link) + "\n")
            for path in paths:
                f.write("\t".join(path) + "\n")

    def write_abundances(outfile, abundances, retained_nodes):
        with open(outfile, "w") as f:
            for node in retained_nodes:
                f.write(f"{node}:{abundances[node]:.6f}\n")

    abundances = load_abundance(abundance_file)
    max_abundance = max(abundances.values())
    threshold = 0.005 * max_abundance

    headers, segments, links, paths = parse_gfa(gfa_file)

    all_gfa_nodes = set(segments.keys())
    abundant_nodes = set(abundances.keys())
    missing_nodes = all_gfa_nodes - abundant_nodes

    to_remove = {
        node for node in abundant_nodes
        if abundances[node] < threshold
    }.union(missing_nodes)


    print(f"Removing {len(to_remove)} nodes below threshold ({threshold:.2f})")
    print(f"Removing nodes: {sorted(to_remove)}")

    for node in to_remove:
        segments.pop(node, None)

    filtered_links = [link for link in links if link[1] not in to_remove and link[3] not in to_remove]

    updated_paths = []
    for path in paths:
        path_id = path[1]
        node_str = path[2]
        nodes = node_str.strip().split(",")
        new_nodes = [n for n in nodes if n.rstrip("+-") not in to_remove]
        if new_nodes:
            updated_path = path[:]
            updated_path[2] = ",".join(new_nodes)
            updated_paths.append(updated_path)
        else:
            print(f"Skipping empty path {path_id}")

    filtered_gfa_file = os.path.join(output_dir, "final_graph.gfa")
    filtered_abundance_file = os.path.join(output_dir, "final_node_abundances.txt")


    write_gfa(filtered_gfa_file, headers, segments, filtered_links, updated_paths)
    retained_nodes = [n for n in abundances if n not in to_remove]
    write_abundances(filtered_abundance_file, abundances, retained_nodes)


if __name__ == "__main__":
    main()