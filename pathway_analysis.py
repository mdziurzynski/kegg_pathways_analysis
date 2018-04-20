import argparse
import urllib.request
from multiprocessing import Pool, cpu_count
import subprocess as sp
import networkx as nx


# 1. Parse KO text output from KAAS
def parse_text_output(file_path):
    """
    requirements: all protein ids must be unique

    @return: dict protein_id: KO number
    """
    res = {}
    with open(file_path, 'r') as kaas_output:
        for line in kaas_output:
            if len(line.split()) == 2:
                key, value = line.split()
                res[key] = value
    return res


# 2. Using rest kegg api, get all pathways that the KO's are involved in
def  get_kegg_pathways(enzyme_id):
    """
    Get pathways with more than one found enzyme.
    @return: dict of pathway names in 'ko' standard as keys and 
        list of enzyme ids present in the pathway 
    """
    res = {}
    pathway_data = urllib.request.urlopen("http://rest.kegg.jp/link/pathway/{0}".format(enzyme_id)).read()
    pathway_data = str(pathway_data, 'utf-8')
    for pathway in pathway_data.split('\n'):
        if pathway:
            ko_id, pathway_id = pathway.split()
            if pathway_id.startswith('ko'):
                if pathway_id.split(':')[1] in res:
                    res[pathway_id.split(':')[1]] = res[pathway_id.split(':')[1]].append(enzyme_id)
                else:
                    res[pathway_id.split(':')[1]] = [enzyme_id]

    return res  


# def  get_kegg_pathways(pid_ko_dict):
#     """
#     Get pathways with more than one found enzyme.
#     @return: dict of pathway names in 'ko' standard as keys and 
#         list of enzyme ids present in the pathway 
#     """
#     unique_enzyme_kos = set(pid_ko_dict.values())
#     res = {}
#     for enzyme in unique_enzyme_kos:
#         pathway_data = urllib.request.urlopen("http://rest.kegg.jp/link/pathway/{0}".format(enzyme)).read()
#         pathway_data = str(pathway_data, 'utf-8')
#         for pathway in pathway_data.split('\n'):
#             if pathway:
#                 ko_id, pathway_id = pathway.split()
#                 if pathway_id.startswith('ko'):
#                     if pathway_id.split(':')[1] in res:
#                         res[pathway_id.split(':')[1]] = res[pathway_id.split(':')[1]].append(enzyme)
#                     else:
#                         res[pathway_id.split(':')[1]] = [enzyme]

#     res = {k: v for k, v in res.items() if len(v) > 1} 
#     return res

# 3. Download kegg xml files for pathways from point 2.
# 4. Extract KO links from kegg xml files and build graphs from them
def build_pathway_graph(pathway_id, enzyme_list):
    G = nx.Graph(name=pathway_id)
    pathway_kgml = urllib.request.urlopen("http://rest.kegg.jp/get/{0}/kgml".format(pathway_id)).read()
    root = ET.fromstring(pathway_kgml)
    for relation in root.findall('relation'):
        entry1_elem = root.find("./entry[@id='{0}']".format(relation.attrib['entry1']))
        entry1_elem = root.find("./entry[@id='{0}']".format(relation.attrib['entry2']))
        for record1 in entry1_elem.attrib['name'].split():
            if not record1.startswith('ko:'):
                continue
            for record2 in entry2_elem.attrib['name'].split():
                if not record2.startswith('ko:'):
                    continue
                G.add_node(record1.split(':')[-1], present=record1.split(':')[-1] in enzyme_list)
                G.add_node(record2.split(':')[-1], present=record2.split(':')[-1] in enzyme_list)
                G.add_edge(record1.split(':')[-1], record2.split(':')[-1])

    return G



# 
# ...? output csv with the following columns:
# - pathway name
# - percentage full map
# - number of present KO's/ number of all KO's in the pathway
# - longest link
# - single KO gaps number

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='KEGG Pathways Analysis\nThis script extends KAAS output by adding aditional information to per pathway data.')
    parser.add_argument('-i', type=str, help='path to KAAS output text file')
    parser.add_argument('-o', type=str, help='path to output file')

    args = parser.parse_args()
    if not args.i or not args.o:
        parser.print_help()
        sys.exit(1)

    numer_of_workers = cpu_count() - 1

    print("Parsing KAAS output...\n")
    pid_ko_dict = parse_text_output(args.i)
    print('done.\n')
    print("Building pathways dict...\n")

    unique_enzyme_kos = set(pid_ko_dict.values())
    print("Number of enzymes to process: {0}".format(len(unique_enzyme_kos)))
    with Pool(processes=numer_of_workers) as pool:
        multiple_results = [pool.apply_async(get_kegg_pathways, (enzyme_id,))
            for enzyme_id in unique_enzyme_kos]
        pathways = [res.get() for res in multiple_results]

    #pathway_ids_dict = get_kegg_pathways(pid_ko_dict)
    print("done.")
    print(pathways)
    

    # print("Building pathway graphs...\n")
    # with Pool(processes=numer_of_workers) as pool:
    #     multiple_results = [pool.apply_async(build_pathway_graph, (pathway_id, enzyme_list))
    #         for (pathway_id, enzyme_list) in pathway_ids_dict]

    #     graphs_list = [res.get() for res in multiple_results]
    # print("done.")
    # print(graphs_list)