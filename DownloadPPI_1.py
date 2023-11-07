import MySQLdb
import networkx as nx
from matplotlib import pylab as plt
import numpy as np


def get_ppi(lcc):
    '''
    Main function to extract the PPI from our local database.
    Connect to GenesGO Database and extract PPI edges from

    '''
    # Open database connection
    db = MySQLdb.connect("<menchelab_server>", "readonly", "<MencheLabPW>", "GenesGO")
    # prepare a cursor object using cursor() method
    cursor = db.cursor()

    sql = """
            SELECT
            e.entrez_1,
            e.entrez_2,
            g1.Locus_Type,
            g1.Locus_Group,
            g2.Locus_Type,
            g2.Locus_Group
            FROM networks.PPI_hippie2017 e
            INNER JOIN GenesGO.hgnc_complete g1 ON e.entrez_1 = g1.Entrez_Gene_ID_NCBI
            INNER JOIN GenesGO.hgnc_complete g2 ON e.entrez_2 = g2.Entrez_Gene_ID_NCBI
            WHERE
                (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'T cell receptor gene' AND g2.Locus_Type = 'T cell receptor gene')                     # 0 links   
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'immunoglobulin gene' AND g2.Locus_Type = 'immunoglobulin gene')              # 4 links   
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'immunoglobulin gene' AND g2.Locus_Type = 'T cell receptor gene')             # 0 links
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'T cell receptor gene' AND g2.Locus_Type = 'immunoglobulin gene')                 # 0 links        
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'T cell receptor gene' AND g2.Locus_Type = 'gene with protein product')       # 17 links        
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'gene with protein product' AND g2.Locus_Type = 'T cell receptor gene')       # 1 links        
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'immunoglobulin gene' AND g2.Locus_Type = 'gene with protein product')        # 115 links        
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'gene with protein product' AND g2.Locus_Type = 'immunoglobulin gene')        # 295 links   
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'gene with protein product' AND g2.Locus_Type = 'gene with protein product')  # 309602 links 


        """
    try:
        # execute SQL query
        cursor.execute(sql)
        data = cursor.fetchall()

    except:
        print('SQL error')

    db.close()

    l_nodes = []
    for x in data:
        l_nodes.append(x[0])
        l_nodes.append(x[1])
    l_nodes = list(set(l_nodes))

    G = nx.Graph()
    G.add_nodes_from(l_nodes)

    for x in data:
        G.add_edge(x[0], x[1])

    print('PPI All:')
    print('Number of genes found: %d' % len(G.nodes()))
    print('Number of interactions found: %d' % len(G.edges()))

    if lcc == 1:
        Nl_l = sorted(nx.connected_components(G))  # generates list of components node lists
        l = [len(x) for x in Nl_l]  # returns list of length of node lists
        idx = l.index(max(l))  # find the index of the maximal length i.e. lcc
        Nlcc = Nl_l[idx]  # pin down lcc
        G_lcc = G.subgraph(Nlcc)  # extract lcc graph
        G = G_lcc.copy()

        print('PPI Only Largest Connected Component:')
        print('Number of genes found: %d' % len(G.nodes()))
        print('Number of interactions found: %d' % len(G.edges()))

    else:
        pass

    return G