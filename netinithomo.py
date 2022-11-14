import networkx as nx
import numpy as np
import random as rand
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import gamma
from itertools import chain
import rand_networks


def intialize_graph(G, N, Num_inf, Beta, Alpha):
    ######################################################
    # Initaialize infection, infected network and history
    #####################################################
    nx.set_node_attributes(G, False, 'infected')
    nx.set_node_attributes(G, 0, 'contact_rate')
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    for i in range(N):
        G.nodes[i]['contact_rate'] = Beta[i]
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G[l]:
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + G.nodes[i]['contact_rate']
    R_tot = np.sum(Rates)
    return R_tot, Rates

def intalize_lam_graph(G, N, beta_sus,beta_inf):
    nx.set_node_attributes(G, False, 'infected')
    nx.set_node_attributes(G, 0, 'contact_rate')
    nx.set_node_attributes(G, 0, 'spread_rate')
    for i in range(N):
        G.nodes[i]['contact_rate'] = beta_sus[i]
        G.nodes[i]['spread_rate'] = beta_inf[i]
    return G


def intalize_group_graph(G, N, beta_sus,beta_inf):
    nx.set_node_attributes(G, False, 'infected')
    nx.set_node_attributes(G, 0, 'contact_rate')
    nx.set_node_attributes(G, 0, 'spread_rate')
    nx.set_node_attributes(G, 0, 'group')
    for i in range(N):
        G.nodes[i]['contact_rate'] = beta_sus[i]
        G.nodes[i]['spread_rate'] = beta_inf[i]
    return G


def set_graph_attriubute_DiGraph(G):
    nx.set_node_attributes(G, False, 'infected')
    nx.set_node_attributes(G, 0, 'contact_rate')
    nx.set_node_attributes(G, 0, 'spread_rate')
    return G


def set_graph_attriubute_markov(G):
    nx.set_node_attributes(G, False, 'infected')
    nx.set_node_attributes(G, 0, 'contact_rate')
    nx.set_node_attributes(G, 0, 'spread_rate')
    nx.set_node_attributes(G, 0, 'infected_neighbors')
    return G


def inatlize_inf_graph(G,Num_inf,N,Alpha,Beta):
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G[l]:
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] +Beta*(G.nodes[i]['contact_rate'] * G.nodes[l]['spread_rate'])
    R_tot = np.sum(Rates)
    return R_tot, Rates


def inatlize_group_inf_graph(G,Num_inf,N,Alpha,Beta):
    groups =[]
    count=0
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    group_num_inf = []
    for l in range(N):
        exist_before = False
        for j in range(len(groups)):
            if G.nodes[l]['spread_rate'] == groups[j]:
                exist_before = True
                G.nodes[l]['group'] = j
                # group_num_inf[j] = group_num_inf[j] + 1
                break
        if exist_before == False:
            groups.append(G.nodes[l]['spread_rate'])
            group_num_inf.append(0)
            G.nodes[l]['group'] = count
            count = count + 1
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            group_num_inf[G.nodes[l]['group']] = group_num_inf[G.nodes[l]['group']] +1
            for i in G[l]:
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] +Beta*(G.nodes[i]['contact_rate'] * G.nodes[l]['spread_rate'])
    return group_num_inf, Rates

# def inatlize_quarntine_graph(G,N,Alpha,Beta):
#     Rates = np.zeros(N)
#     for l in range(N):
#         if G.nodes[l]['infected'] == True:
#             Rates[l] = Alpha
#             for i in G[l]:
#                 if (G.nodes[i]['infected'] == False):
#                     Rates[i] = Rates[i] +Beta*(G.nodes[i]['contact_rate'] * G.nodes[l]['spread_rate'])
#     R_tot = np.sum(Rates)
#     return R_tot, Rates


def inatlize_quarntine_graph(G,N,Alpha,Beta):
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G.successors(l):
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + Beta
    R_tot = np.sum(Rates)
    return R_tot, Rates


def inatlize_inf_DiGraph(G,Num_inf,N,Alpha,Beta):
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G.successors(l):
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + Beta
    R_tot = np.sum(Rates)
    return R_tot, Rates

def inatlize_inf_weighted_graph(G,Num_inf,N,Alpha,Beta):
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G.successors(l):
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + Beta*G[l][i]['weight']
    R_tot = np.sum(Rates)
    return R_tot, Rates


def inatlize_markov_graph(G,Num_inf,N,Alpha,Beta):
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G.successors(l):
                # if (G.nodes[i]['infected'] == False):
                G.nodes[i]['infected_neighbors']=G.nodes[i]['infected_neighbors']+1
                Rates[i] = Rates[i] + Beta
    R_tot = np.sum(Rates)
    return R_tot, Rates


def inatlize_inf_weighted_graph(G,Num_inf,N,Alpha,Beta):
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G.successors(l):
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + Beta*G[l][i]['weight']
    R_tot = np.sum(Rates)
    return R_tot, Rates




def inatlize_inf_undirected(G,Num_inf,N,Alpha,Beta):
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G[l]:
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + Beta
    R_tot = np.sum(Rates)
    return R_tot, Rates


def intialize_graph_track_bimodal(G, N, Num_inf, Beta, Alpha,lam_plus,lam_minus):
    ######################################################
    # Initaialize infection, infected network and history
    #####################################################
    Num_inf_plus,Num_inf_minus = 0,0
    nx.set_node_attributes(G, False, 'infected')
    nx.set_node_attributes(G, 0, 'contact_rate')
    for i in range(N):
        G.nodes[i]['contact_rate'] = Beta[i]
    for i in rand.sample(range(0, N - 1), Num_inf):
        G.nodes[i]['infected'] = True
        if round(G.nodes[i]['contact_rate'],6) == lam_plus:
            Num_inf_plus = Num_inf_plus+1
        else:
            Num_inf_minus = Num_inf_minus +1
    Rates = np.zeros(N)
    for l in range(N):
        if G.nodes[l]['infected'] == True:
            Rates[l] = Alpha
            for i in G[l]:
                if (G.nodes[i]['infected'] == False):
                    Rates[i] = Rates[i] + G.nodes[i]['contact_rate']
    R_tot = np.sum(Rates)
    return R_tot, Rates, Num_inf_plus, Num_inf_minus


def bi_beta(N,epsilon,l):
    b=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon)),np.full((1,int(N/2)),l*(1-epsilon))),axis=1)
    np.random.shuffle(b[0])
    return b[0]


def bi_beta_correlated(N,epsilon_lam,epsilon_mu,l):
    if epsilon_lam==0.0:
        a=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon_mu)),np.full((1,int(N/2)),l*(1-epsilon_mu))),axis=1)
        np.random.shuffle(a[0])
        return l*np.ones(N),a[0] 
    b=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon_lam)),np.full((1,int(N/2)),l*(1-epsilon_lam))),axis=1)
    np.random.shuffle(b[0])
    a=[]
    for x in b[0]:
        a.append(1+epsilon_mu) if x>1.0 else a.append(1-epsilon_mu)
    a=np.array(a)
    return b[0],np.array(a)


def bi_beta_anti_correlated(N,epsilon_lam,epsilon_mu,l):
    if epsilon_lam==0.0:
        a=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon_mu)),np.full((1,int(N/2)),l*(1-epsilon_mu))),axis=1)
        np.random.shuffle(a[0])
        return l*np.ones(N),a[0] 
    b=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon_lam)),np.full((1,int(N/2)),l*(1-epsilon_lam))),axis=1)
    np.random.shuffle(b[0])
    a=[]
    for x in b[0]:
        a.append(1+epsilon_mu) if x<1.0 else a.append(1-epsilon_mu)
    a=np.array(a)
    return b[0],np.array(a)


def bi_beta_any(N,epsilon_lam,epsilon_mu,l):
    if epsilon_lam==0.0:
        a=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon_mu)),np.full((1,int(N/2)),l*(1-epsilon_mu))),axis=1)
        np.random.shuffle(a[0])
        return l*np.ones(N),a[0]
    b=np.concatenate((np.full((1,int(N/2)),l*(1+epsilon_lam)),np.full((1,int(N/2)),l*(1-epsilon_lam))),axis=1)
    np.random.shuffle(b[0])
    a=[]
    if epsilon_lam*epsilon_mu<0:
        for x in b[0]:
            a.append(1+epsilon_mu) if x<1.0 else a.append(1-epsilon_mu)
    else:
        for x in b[0]:
            a.append(1 + epsilon_mu) if x > 1.0 else a.append(1 - epsilon_mu)
    a=np.array(a)
    return b[0],np.array(a)


def triangular_beta(N,epsilon_lam,epsilon_mu,correlation):
    a_triangular=lambda s: 1-np.sqrt(6)*s
    b_triangular=lambda s: 1+np.sqrt(6)*s
    beta_lam=np.random.triangular(a_triangular(epsilon_lam),1.0,b_triangular(epsilon_lam),N)
    beta_mu=np.random.triangular(a_triangular(epsilon_mu),1.0,b_triangular(epsilon_mu),N)
    beta_lam=sorted(beta_lam)
    beta_mu=sorted(beta_mu,reverse=True) if correlation is 'a' else sorted(beta_mu,reverse=False)
    return beta_lam,beta_mu


def uniform_beta(N,epsilon_lam,epsilon_mu,correlation):
    a_unifrom=lambda s: 1-np.sqrt(3)*s
    b_uniform=lambda s: 1+np.sqrt(3)*s
    beta_lam = np.linspace(a_unifrom(epsilon_lam),b_uniform(epsilon_lam),N)
    beta_mu = np.linspace(b_uniform(epsilon_mu),a_unifrom(epsilon_mu),N)
    beta_mu=sorted(beta_mu,reverse=True) if correlation is 'a' else sorted(beta_mu,reverse=False)
    return beta_lam,beta_mu

def general_beta(N,epsilon_lam,epsilon_mu,name,avg_degree):
    beta_lam = np.array(rand_networks.create_random_degree_sequence(name, np.abs(epsilon_lam), avg_degree, N)/avg_degree)
    beta_mu = np.array(rand_networks.create_random_degree_sequence(name, np.abs(epsilon_mu), avg_degree, N)/avg_degree)
    beta_lam=np.sort(beta_lam)
    # beta_mu=sorted(beta_mu,reverse=True) if epsilon_lam * epsilon_mu > 0 else sorted(beta_mu,reverse=False)
    beta_mu=np.sort(beta_mu) if epsilon_lam * epsilon_mu > 0 else np.sort(beta_mu)[::-1]
    return beta_lam, beta_mu

def fraction_beta(N,epsilon_lam,epsilon_mu,l,f):
    if epsilon_lam==0.0:
        a=np.concatenate((np.full((1,int(N*f)),l*(1+epsilon_mu)),np.full((1,int(N*(1-f))),l*(1-epsilon_mu))),axis=1)
        np.random.shuffle(a[0])
        return l*np.ones(N),a[0]
    b=np.concatenate((np.full((1,int(N*f)),l*(1+epsilon_lam)),np.full((1,int(N*(1-f))),l*(1-epsilon_lam))),axis=1)
    np.random.shuffle(b[0])
    a=[]
    for x in b[0]:
        if epsilon_lam*epsilon_mu>0:
            a.append(1+epsilon_mu) if x>1.0 else a.append(1-epsilon_mu)
        else:
            a.append(1+epsilon_mu) if x<1.0 else a.append(1-epsilon_mu)
    a=np.array(a)
    return b[0],np.array(a)

if __name__ == '__main__':
    general_beta(1000, 0.1, 0.1, 'gauss_c',1000)