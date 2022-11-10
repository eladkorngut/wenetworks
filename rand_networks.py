import csv

import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
from collections import defaultdict
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from itertools import combinations
from itertools import zip_longest
from itertools import chain
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import uniform

import json


def draw_basic_nx_g(G):
    fig = plt.gcf()
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos,
                           node_color='r',
                           node_size=500,
                           alpha=1)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=1.0,arrowsize=20)
    plt.tight_layout()
    labels = nx.draw_networkx_labels(G, pos=pos,font_color='w')
    fig.savefig('network_strcture.png', dpi=100)
    plt.show()
    return pos


def draw_infections_nx_g(G,pos=None,headline=None,highest_rate=None):
    fig = plt.gcf()
    infected_list = []
    disinfected_list=[]
    for key, value in nx.get_node_attributes(G, 'infected').items():
        if value == True:
            infected_list.append(key)
        else:
            disinfected_list.append(key)
    if pos==None:
        pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos,
                           nodelist=infected_list,
                           node_color='r',
                           node_size=500,
                           alpha=1)
    nx.draw_networkx_nodes(G, pos,
                           nodelist=disinfected_list,
                           node_color='b',
                           node_size=500,
                           alpha=1)
    if highest_rate is not None:
        nx.draw_networkx_nodes(G, pos,
                               nodelist=highest_rate,
                               node_color='g',
                               node_size=500,
                               alpha=1)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=1.0)
    plt.tight_layout()
    labels = nx.draw_networkx_labels(G, pos=pos,font_color='w')
    fig.savefig('network_infection'+' '+str(headline)+'.png', dpi=100)
    plt.show()


# def draw_rate_nx_g(G):
#     Rate=[]
#     for n in range(G.number_of_nodes()):
#         Rate.append(G.nodes[n]['rate'])
#     Rate=np.dot(Rate,1/np.max(Rate))
#     Rate=np.array(Rate)
#     fig = plt.gcf()
#     pos = nx.spring_layout(G)
#     nx.draw_networkx_nodes(G, pos,
#                            cmap=plt.cm.Blues,
#                            node_size=500,
#                            alpha=1)
#     nx.draw_networkx_edges(G, pos, width=1.0, alpha=1.0)
#     plt.tight_layout()
#     labels = nx.draw_networkx_labels(G, pos=pos,font_color='w')
#     fig.savefig('/home/elad/Michael Assaf/binomal/photos/network_strcture.png', dpi=100)
#     plt.show()
#     return pos

def export_network(G):
    f = open('network_data.txt', "w+")
    for n in range(G.order()):
        f.write(str(G.degree(n)) + ' ')
        for neighbor in G[n]:
            f.write(str(neighbor) + ' ')
        f.write('\n')
    f.close()


def export_adj_network_cpp(G,infile):
    f = open(infile, "w+")
    for n in range(G.order()):
        f.write(str(G.degree(n)) + ',')
        for neighbor in G[n]:
            f.write(str(neighbor) + ',')
        f.write('\n')
    f.close()


def random_bimodal_graph(d1,d2, n, seed=None):
    """Return a random bimodal graph of n nodes each with degree d1 and d2.
        the degree distubition is half for each one of the two degrees

    The resulting graph G has no self-loops or parallel edges.

    Parameters
    ----------
    d1 : int
      Degre
    d2 : int
      Degree
    n : integer
      Number of nodes. The value of n*d must be even.
    seed : hashable object
        The seed for random number generator.

    Notes
    -----
    The nodes are numbered form 0 to n-1.

    Kim and Vu's paper [2]_ shows that this algorithm samples in an
    asymptotically uniform way from the space of random graphs when
    d = O(n**(1/3-epsilon)).

    References
    ----------
    .. [1] A. Steger and N. Wormald,
       Generating random regular graphs quickly,
       Probability and Computing 8 (1999), 377-396, 1999.
       http://citeseer.ist.psu.edu/steger99generating.html

    .. [2] Jeong Han Kim and Van H. Vu,
       Generating random regular graphs,
       Proceedings of the thirty-fifth ACM symposium on Theory of computing,
       San Diego, CA, USA, pp 213--222, 2003.
       http://portal.acm.org/citation.cfm?id=780542.780576

       The multiplcation of n*d1 and n*d2 needs to be even. For the safe side choose even number for both nodes and degree
    """
    if (n * d1) % 2 != 0 or (n * d2) % 2!=0:
        raise nx.NetworkXError("n * d must be even")

    if not (0 <= d1 < n or 0 <= d2 < n):
        raise nx.NetworkXError("the 0 <= d < n inequality must be satisfied")

    if d1 == 0 or d2 == 0:
        return nx.empty_graph(n)

    if seed is not None:
        random.seed(seed)

    def _suitable(edges, potential_edges):
    # Helper subroutine to check if there are suitable edges remaining
    # If False, the generation of the graph has failed
        if not potential_edges:
            return True
        for s1 in potential_edges:
            for s2 in potential_edges:
                # Two iterators on the same dictionary are guaranteed
                # to visit it in the same order if there are no
                # intervening modifications.
                if s1 == s2:
                    # Only need to consider s1-s2 pair one time
                    break
                if s1 > s2:
                    s1, s2 = s2, s1
                if (s1, s2) not in edges:
                    return True
        return False

    def _try_creation(n):
        # Attempt to create an edge set

        edges = set()
        stub_d1=random.sample(list(range(n)),int(n/2))
        stub_d2= [item for item in list(range(n)) if item not in stub_d1]
        stub_d1=stub_d1 * d1
        stub_d2=stub_d2 * d2
        stubs = stub_d1+stub_d2

        while stubs:
            potential_edges = defaultdict(lambda: 0)
            random.shuffle(stubs)
            stubiter = iter(stubs)
            for s1, s2 in zip(stubiter, stubiter):
                if s1 > s2:
                    s1, s2 = s2, s1
                if s1 != s2 and ((s1, s2) not in edges):
                    edges.add((s1, s2))
                else:
                    potential_edges[s1] += 1
                    potential_edges[s2] += 1

            if not _suitable(edges, potential_edges):
                return None # failed to find suitable edge set

            stubs = [node for node, potential in potential_edges.items()
                     for _ in range(potential)]
        return edges

    # Even though a suitable edge set exists,
    # the generation of such a set is not guaranteed.
    # Try repeatedly to find one.
    edges = _try_creation(n)
    while edges is None:
        edges = _try_creation(n)

    G = nx.Graph()
    G.name = "random_bimodal_graph(%s, %s)" % (d1, n)
    G.add_edges_from(edges)

    return G



def random_bimodal_directed_graph(d1_in,d1_out,d2_in,d2_out, n, seed=None):
    """Return a random bimodal directed graph of n nodes each with degree d1 and d2.
        the degree distubition is half for each one of the two degrees

    The resulting graph G has no self-loops or parallel edges.

    Parameters
    ----------
    d1 : int
      Degre
    d2 : int
      Degree
    n : integer
      Number of nodes. The value of n*d must be even.
    seed : hashable object
        The seed for random number generator.

    Notes
    -----
    The nodes are numbered form 0 to n-1.

    Kim and Vu's paper [2]_ shows that this algorithm samples in an
    asymptotically uniform way from the space of random graphs when
    d = O(n**(1/3-epsilon)).

    References
    ----------
    .. [1] A. Steger and N. Wormald,
       Generating random regular graphs quickly,
       Probability and Computing 8 (1999), 377-396, 1999.
       http://citeseer.ist.psu.edu/steger99generating.html

    .. [2] Jeong Han Kim and Van H. Vu,
       Generating random regular graphs,
       Proceedings of the thirty-fifth ACM symposium on Theory of computing,
       San Diego, CA, USA, pp 213--222, 2003.
       http://portal.acm.org/citation.cfm?id=780542.780576

       The multiplcation of n*d1 and n*d2 needs to be even. For the safe side choose even number for both nodes and degree
    """
    if (n * d1_in) % 2 != 0 or (n * d2_in) % 2!=0:
        raise nx.NetworkXError("n * d must be even")

    if not (0 <= d1_in < n or 0 <= d2_in < n):
        raise nx.NetworkXError("the 0 <= d < n inequality must be satisfied")

    if not (d1_in+d2_in==d1_out+d2_out):
        raise nx.NetworkXError("the d1_in+d2_in==d1_out+d2_out equality must be satisfied")

    if seed is not None:
        random.seed(seed)

    def _suitable(edges, potential_edges_in,potential_edges_out):
    # Helper subroutine to check if there are suitable edges remaining
    # If False, the generation of the graph has failed
        if not (potential_edges_in and potential_edges_out):
            return True
        for s1 in potential_edges_in:
            for s2 in potential_edges_out:
                # Two iterators on the same dictionary are guaranteed
                # to visit it in the same order if there are no
                # intervening modifications.
                if s1 == s2:
                    # Only need to consider s1-s2 pair one time
                    break
                if (s1, s2) not in edges:
                    return True
        return False

    def _try_creation(n):
        # Attempt to create an edge set

        edges = set()
        stub_d1_in = random.sample(list(range(n)),int(n/2))
        stub_d2_in = [item for item in list(range(n)) if item not in stub_d1_in]
        stub_d1_in =stub_d1_in * d1_in
        stub_d2_in = stub_d2_in * d2_in
        stubs_in = stub_d1_in+stub_d2_in

        stub_d1_out = [item for item in list(range(n)) if item in stub_d1_in]
        stub_d2_out = [item for item in list(range(n)) if item in stub_d2_in]
        stub_d1_out = stub_d1_out * d1_out
        stub_d2_out = stub_d2_out * d2_out
        stubs_out = stub_d1_out + stub_d2_out


        while stubs_in and stubs_out:
            potential_edges_in = defaultdict(lambda: 0)
            potential_edges_out = defaultdict(lambda: 0)
            random.shuffle(stubs_in)
            stubiter_in = iter(stubs_in)
            random.shuffle(stubs_out)
            stubiter_out = iter(stubs_out)
            for s1, s2 in zip(stubiter_in, stubiter_out):
                if s1 != s2 and ((s1, s2) not in edges):
                    edges.add((s1, s2))
                else:
                    potential_edges_in[s1] += 1
                    potential_edges_out[s2] += 1

            if not _suitable(edges, potential_edges_in,potential_edges_out):
                return None # failed to find suitable edge set

            stubs_in = [node for node, potential in potential_edges_in.items()
                     for _ in range(potential)]
            stubs_out = [node for node, potential in potential_edges_out.items()
                     for _ in range(potential)]
        return edges

    # Even though a suitable edge set exists,
    # the generation of such a set is not guaranteed.
    # Try repeatedly to find one.
    edges = _try_creation(n)
    while edges is None:
        edges = _try_creation(n)

    G = nx.DiGraph()
    G.name = "random_bimodal_graph(%s, %s)" % (d1_in, n)
    G.add_edges_from(edges)
    return G

def create_exact_gauss_dis(epsilon,avg_degree,N):
    low, high,normalization = (avg_degree * (1 - 4 * epsilon)).astype(int), (avg_degree * (1 + 4 * epsilon)).astype(int),N
    possible_degrees = np.arange(low, high)
    deg_dist = (norm.pdf(possible_degrees, loc=avg_degree, scale=avg_degree * epsilon) * normalization).astype(int)
    deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    while np.size(deg_for_algo)<N:
        normalization=normalization+1
        deg_dist = (norm.pdf(np.arange(low, high), loc=avg_degree, scale=avg_degree * epsilon) * normalization).astype(int)
        deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    deg_dist, possible_degrees = deg_dist[deg_dist != 0], possible_degrees[deg_dist != 0]
    for i in range(np.size(deg_for_algo)-N): deg_dist[i] =deg_dist[i] - 1
    # deg_dist[0] = deg_dist[0] - 1 if np.size(deg_for_algo) != N else deg_dist[0]
    deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    return np.array(deg_for_algo)


def create_deg_distubtion_conf_model(epsilon,avg_degree,N,dist_type):

    def create_stubs_degree(epsilon,avg_degree,dist_type,normalization):
        if dist_type=='gauss_c':
            low, high = (avg_degree * (1 - 4 * epsilon)).astype(int), (
                        avg_degree * (1 + 4 * epsilon)).astype(int)
            possible_degrees = np.arange(low, high)
            deg_dist = (norm.pdf(possible_degrees, loc=avg_degree, scale=avg_degree * epsilon) * normalization).astype(int)
        elif dist_type == 'bimodal_c':
            possible_degrees = np.array([(avg_degree * (1 - epsilon)).astype(int), (avg_degree * (1 + epsilon)).astype(int)])
            deg_dist = np.array([N/2,N/2]).astype(int)
        elif dist_type == 'gamma_c':
            low, high = (avg_degree * (1 - 100 * epsilon)).astype(int), (
                        avg_degree * (1 + 100 * epsilon)).astype(int)
            possible_degrees = np.arange(low, high)
            deg_dist = (gamma.pdf(possible_degrees, a=1/epsilon**2,scale=epsilon**2*avg_degree)*normalization).astype(int)
        elif dist_type == 'uniform_c':
            low,high = (avg_degree * (1 - 100 * epsilon)).astype(int), (
                        avg_degree * (1 + 100 * epsilon)).astype(int)
            possible_degrees = np.arange(low, high)
            deg_dist = (uniform.pdf(possible_degrees, loc=int(avg_degree*(1-epsilon)),scale=2*np.sqrt(3)*epsilon*avg_degree)*normalization).astype(int)
        return list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist))),deg_dist,possible_degrees

    def trim_edges_to_fit_N_nodes(deg_for_algo,deg_dist,possible_degrees,N):
        deg_dist, possible_degrees = deg_dist[deg_dist != 0], possible_degrees[deg_dist != 0]
        reduce_by_one =  np.random.uniform(0, np.size(deg_dist)-1, np.size(deg_for_algo) - N).astype(int)
        for i in reduce_by_one: deg_dist[i]=deg_dist[i]-1
        # for i in range(np.size(deg_for_algo) - N): deg_dist[i] = deg_dist[i] - 1
        return deg_dist,possible_degrees

    normalization=N
    deg_for_algo,deg_dist,possible_degrees=create_stubs_degree(epsilon,avg_degree,dist_type,normalization)
    while np.size(deg_for_algo)<=N:
        normalization=normalization+1
        deg_for_algo,deg_dist,possible_degrees = create_stubs_degree(epsilon,avg_degree,dist_type,normalization)
    deg_dist, possible_degrees = trim_edges_to_fit_N_nodes(deg_for_algo,deg_dist,possible_degrees,N)
    deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    return np.array(deg_for_algo)




def create_random_degree_sequence(name,epsilon,avg_degree,N):
    if epsilon==0.0:
        return np.ones(N).astype(int)*int(avg_degree)
    if name=='gauss':
        return np.random.normal(avg_degree, epsilon * avg_degree, N)
    elif name=='uniform':
        return np.random.uniform(avg_degree*(1-np.sqrt(3)*epsilon),avg_degree*(1+np.sqrt(3)*epsilon),N).astype(int)
    elif name=='gamma':
        return np.random.gamma(1/epsilon**2,avg_degree*epsilon**2,N)
    elif name=='gauss_e':
        return create_exact_gauss_dis(epsilon,avg_degree,N)
    elif name=='gauss_c' or name=='bimodal_c' or name=='gamma_c' or name=='uniform_c':
        return create_deg_distubtion_conf_model(epsilon,avg_degree,N,name)


def equate_degree_sequence(din,dout,N,avg_degree,eps_in=0.0,eps_out=0.0):
    din[din < 0] = 0
    dout[dout < 0] = 0
    return_order_din_dout=True
    if eps_in==0.0:
        for i in range(N*avg_degree-np.sum(dout)):
            node = int(N * np.random.random())
            dout[node] = dout[node] + 1
        return din,dout
    if eps_out==0.0:
        for i in range(N*avg_degree-np.sum(din)):
            node = int(N * np.random.random())
            din[node] = din[node] + 1
        return din,dout
    if np.sum(din)>np.sum(dout):
        high, low = din, dout
    else:
        return_order_din_dout = False
        high, low = dout, din
    prob=np.random.uniform(0,1,np.sum(high)-np.sum(low))
    for p in prob:
        if p>0.5:
            node=int(N*np.random.random())
            while high[node]==0:
                node = int(N * np.random.random())
            high[node]=high[node]-1
        else:
            node=int(N*np.random.random())
            while low[node]>=N:
                node = int(N * np.random.random())
            low[node]=low[node]+1
    if return_order_din_dout: return high,low
    return low,high


chaini = chain.from_iterable

def _to_stublist_simple_graph(degree_sequence):
    """Returns a list of degree-repeated node numbers.

    ``degree_sequence`` is a list of nonnegative integers representing
    the degrees of nodes in a graph.

    This function returns a list of node numbers with multiplicities
    according to the given degree sequence. For example, if the first
    element of ``degree_sequence`` is ``3``, then the first node number,
    ``0``, will appear at the head of the returned list three times. The
    node numbers are assumed to be the numbers zero through
    ``len(degree_sequence) - 1``.

    Examples
    --------

    >>> degree_sequence = [1, 2, 3]
    >>> _to_stublist(degree_sequence)
    [0, 1, 1, 2, 2, 2]

    If a zero appears in the sequence, that means the node exists but
    has degree zero, so that number will be skipped in the returned
    list::

    >>> degree_sequence = [2, 0, 1]
    >>> _to_stublist(degree_sequence)
    [0, 0, 2]

    """
    return list(chaini([n] * d for n, d in enumerate(degree_sequence)))


def _configuration_model_simple_graph(
    create_using,name,epsilon_in,epsilon_out,avg_degree,n):
    """Helper function for generating either undirected or directed
    configuration model graphs.

    ``deg_sequence`` is a list of nonnegative integers representing the
    degree of the node whose label is the index of the list element.

    ``create_using`` see :func:`~networkx.empty_graph`.

    ``directed`` and ``in_deg_sequence`` are required if you want the
    returned graph to be generated using the directed configuration
    model algorithm. If ``directed`` is ``False``, then ``deg_sequence``
    is interpreted as the degree sequence of an undirected graph and
    ``in_deg_sequence`` is ignored. Otherwise, if ``directed`` is
    ``True``, then ``deg_sequence`` is interpreted as the out-degree
    sequence and ``in_deg_sequence`` as the in-degree sequence of a
    directed graph.

    .. note::

       ``deg_sequence`` and ``in_deg_sequence`` need not be the same
       length.

    ``seed`` is a random.Random or numpy.random.RandomState instance

    This function returns a graph, directed if and only if ``directed``
    is ``True``, generated according to the configuration model
    algorithm. For more information on the algorithm, see the
    :func:`configuration_model` or :func:`directed_configuration_model`
    functions.

    """
    G = nx.empty_graph(n, create_using)
    # If empty, return the null graph immediately.
    if n == 0:
        return G
    # Build a list of available degree-repeated nodes.  For example,
    # for degree sequence [3, 2, 1, 1, 1], the "stub list" is
    # initially [0, 0, 0, 1, 1, 2, 3, 4], that is, node 0 has degree
    # 3 and thus is repeated 3 times, etc.
    #
    # Also, shuffle the stub list in order to get a random sequence of
    # node pairs.
    def _suitalbe(edges,potential_edges_in,potential_edges_out):
        # Helper suroutine to check if there are suitable edges remaining
        # If false generation of graph failed
        if not (potential_edges_in and potential_edges_out):
            return True
        for s1 in potential_edges_in:
            for s2 in potential_edges_out:
                if s1==s2:
                    break
                if (s2,s1) not in edges:
                    return True
        return False
    def _try_creation(n,epsilon_in,epsilon_out,avg_degree):
        edges = set()
        din, dout = equate_degree_sequence(create_random_degree_sequence(name, np.abs(epsilon_in), avg_degree, n),
                                           create_random_degree_sequence(name, np.abs(epsilon_out), avg_degree, n), n,avg_degree,epsilon_in,epsilon_out)
        din = np.sort(din)
        dout = np.sort(dout) if epsilon_out * epsilon_in > 0 else np.sort(dout)[::-1]
        pairs = zip_longest(dout, din, fillvalue=0)
        # Unzip the list of pairs into a pair of lists.
        out_deg, in_deg = zip(*pairs)
        out_stublist = _to_stublist_simple_graph(out_deg)
        in_stublist = _to_stublist_simple_graph(in_deg)
        while in_stublist and out_stublist:
            potential_edge_in = defaultdict(lambda: 0)
            potential_edge_out = defaultdict(lambda: 0)
            np.random.shuffle(in_stublist)
            subiter_in = iter(in_stublist)
            np.random.shuffle(out_stublist)
            subiter_out = iter(out_stublist)
            for s1,s2 in zip(subiter_in,subiter_out):
                if s1!=s2 and ((s2,s1) not in edges):
                    edges.add((s2,s1))
                else:
                    potential_edge_in[s1] += 1
                    potential_edge_out[s2] += 1
            if not _suitalbe(edges,potential_edge_in,potential_edge_out):
                return None # failed to find suitable edge set
            in_stublist = [node for node, potential in potential_edge_in.items() for _ in range(potential)]
            out_stublist = [node for node,potential in potential_edge_out.items() for _ in range(potential)]
        return edges
    edges = _try_creation(n,epsilon_in,epsilon_out,avg_degree)
    while edges is None:
        edges = _try_creation(n,epsilon_in,epsilon_out,avg_degree)
    G =nx.DiGraph()
    G.name = "simple_directed_config_graph"
    G.add_edges_from(edges)
    return G



def configuration_model_directed_graph(
    name, epsilon_in,epsilon_out,avg_degree,n, create_using=None):
    """Returns a directed_random graph with the given degree sequences.

    The configuration model generates a random directed pseudograph
    by randomly assigning edges to match the given degree sequences.
    Self edges and loops are excluded by randomly choosing
    Parameters
    ----------
    in_degree_sequence :  list of nonnegative integers
       Each list entry corresponds to the in-degree of a node.
    out_degree_sequence :  list of nonnegative integers
       Each list entry corresponds to the out-degree of a node.
    create_using : NetworkX graph constructor, optional (default MultiDiGraph)
        Graph type to create. If graph instance, then cleared before populated.
    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.

    Returns
    -------
    G : MultiDiGraph
        A graph with the specified degree sequences.
        Nodes are labeled starting at 0 with an index
        corresponding to the position in deg_sequence.

    Raises
    ------
    NetworkXError
        If the degree sequences do not have the same sum.

    See Also
    --------
    configuration_model

    Notes
    -----
    Algorithm as described by Newman [1]_.

    References
    ----------
    .. [1] Newman, M. E. J. and Strogatz, S. H. and Watts, D. J.
       Random graphs with arbitrary degree distributions and their applications
       Phys. Rev. E, 64, 026118 (2001)
    """

    if create_using is None:
        create_using = nx.MultiDiGraph

    G = _configuration_model_simple_graph(
        create_using,
        name,
        epsilon_in,
        epsilon_out,
        avg_degree,
        n
    )

    name = "directed configuration_model {} nodes {} edges"
    return G


def configuration_model_undirected_graph(epsilon,avg_degree,N):
    d=np.random.normal(avg_degree, epsilon * avg_degree, N).astype(int)
    if np.sum(d)%2!=0:
        d[int(N*np.random.random())]+=1
    G = nx.configuration_model(d.astype(int))
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


def import_jason_network(file):
    a=[]
    with open(file,"r") as f:
        reader = csv.reader(f,delimiter=",")
        for i, line in enumerate(reader):
            a.append(line)
    return np.array(a)

def create_jason_edges(jason_edge):
    edges=[]
    for i in range(np.size(jason_edge,0)):
        x=list(filter(('').__ne__, jason_edge[i,:]))
        for j in x:
            if j is '':
                pass
            edges.append((i, int(j)))
    return edges

def jason_graph(file_name):
    csv_edges_jason=import_jason_network(file_name)
    edges=create_jason_edges(csv_edges_jason)
    G = nx.Graph()
    G.name = "myFile.csv"
    G.add_edges_from(edges)
    return G


if __name__ == '__main__':
    N,k=2000,50
    interaction_strength = 100.5
    # d1_in, d1_out, d2_in, d2_out = 45,45,55,55
    # beta_inf,beta_sus = 0.1,0.1
    # susceptibility_avg = 1.0
    # infectability_avg = 1.0
    #
    # # G = nx.complete_graph(N)
    # G = nx.random_regular_graph(k, N)
    # susceptibility = 'bimodal'
    # infectability = 'bimodal'
    # choose_beta = lambda net_dist, avg, epsilon: np.random.normal(avg, epsilon * avg, N) \
    #     if net_dist == 'gauss' else np.random.gamma((avg / epsilon) ** 2, epsilon ** 2 / avg, N) \
    #     if net_dist == 'gamma' else np.zeros(N) if net_dist == 'z' else np.ones(
    #     N) if net_dist == 'ones' else netinithomo.bi_beta(N, epsilon, avg)
    #
    # epsilon_sus,epsilon_inf=0.1,0.1
    # beta_sus = choose_beta(susceptibility, susceptibility_avg, epsilon_sus)
    # beta_inf = choose_beta(infectability, infectability_avg, epsilon_inf)
    #
    # # random_bimodal_directed_graph(d1_in, d1_out, d2_in, d2_out, N)
    # G = netinithomo.intalize_lam_graph(G, N, beta_sus, beta_inf)
    # temp=nx.adjacency_data(G)
    # with open('json_data.json', 'w') as outfile:
    #     json.dump(temp, outfile)
    # export_network(G)
    # export_adj_network_cpp(G,'Gnull.csv')
    # print('This no love song')

    fig=plt.figure()
    # G = rand_weighted_networks.bulid_weighted_network(N, k, 0.0, 0.0, interaction_strength)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # degree_in = [G.in_degree(n,'weight') for n in G.nodes()]
    # degree_out = [G.out_degree(n,'weight') for n in G.nodes()]
    # # plt.hist((degree_in, degree_out), bins=100)
    # plt.xlabel('Degree')
    # plt.ylabel('count')
    # plt.title('N='+str(N)+' k='+str(k)+' eta='+str(interaction_strength))
    # plt.hist(degree_out, bins=100)
    # fig.savefig('hist.png', dpi=400)
    # plt.show()

    eta=np.linspace(0.1,2.5,5)
    std_in,std_out=[],[]
    for i in eta:
        G = rand_weighted_networks.bulid_weighted_network(N, k, 0.0, 0.0, i)
        G = netinithomo.set_graph_attriubute_DiGraph(G)
        degree_in = [G.in_degree(n, 'weight') for n in G.nodes()]
        degree_out = [G.out_degree(n, 'weight') for n in G.nodes()]
        std_in.append(np.std(degree_in))
    plt.plot(eta,std_in)
    plt.show()
    # G=random_bimodal_directed_graph(20, 50, 80, 50, 100, seed=None)
    # draw_basic_nx_g(G)
    # np.histogram(create_random_degree_sequence('gauss', 0.1, 108.5, 600),100)
    # plt.hist(create_random_degree_sequence('gauss', 0.1, 108.5, 600), bins='auto')
    # plt.show()
    # temp=configuration_model_directed_graph('gauss',0.1,0.1,108.5,600)
    # din,dout=create_random_degree_sequence('gauss', 0.1, 108.5, 600),create_random_degree_sequence('gauss', 0.1, 108.5, 600)
    # plt.hist((din,dout),bins='auto')
    # plt.show()
    # temp=equate_degree_sequence(din,dout,600)
    # plt.hist(temp, bins='auto')
    # plt.show()
    # degree=create_random_degree_sequence('gauss', 0.1, 108.5, 600)
    # temp=nx.utils.discrete_sequence(1000,degree)
    # v_in,v_out=[],[]
    # # eps_in,eps_out=[0.1,0.1,0.1,0.1,0.1,0.1],[0.05,0.1,0.15,0.2,0.25,0.3]
    # eps_in,eps_out=[0.0],[0.0]
    # for ein,eout in zip(eps_in,eps_out):
    #     G=configuration_model_directed_graph('gauss_c', ein, eout, 200, 1000)
    #     degree_in=[G.in_degree(n) for n in G.nodes()]
    #     degree_out=[G.out_degree(n) for n in G.nodes()]
    #     # v_in.append(np.std(degree_in)/np.mean(degree_in))
    #     # v_out.append(np.std(degree_out)/np.mean(degree_out))
    #     plt.hist((degree_in,degree_out),bins=500)
    #     # plt.hist((degree_out),bins=400)
    #     plt.show()
    # print(*v_in, sep = ", ")
    # print(*v_out, sep = ", ")
    # degrees_in = [G.in_degree(n) for n in G.nodes()]
    # degrees_out = [G.out_degree(n) for n in G.nodes()]
    # plt.hist((degrees_in, degrees_out), bins='auto')
    # plt.show()
    # temp=import_jason_network('jason_trans_file_no_degree.csv')
    # temp_edge=create_jason_edges(temp)
    # G=jason_graph('jason_trans_file_no_degree.csv')
    # print('This no love song')
# k=4
# d1=78
# d2=66
# N=796

#This set of commands is used in order to create a random homogenous graph and draw it

# G=nx.random_regular_graph(k,N)
# draw_basic_nx_g(G)
# export_network(G)


#This set of commands is used in order to create a bimodal graph and draw it

# G=random_bimodal_graph(d1 , d2, N, seed=None)
# export_network(G)
# draw_basic_nx_g(G)



