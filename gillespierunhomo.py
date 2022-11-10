import numpy as np
import bisect
import netinithomo
import networkx as nx
import csv
import pickle
import sys
import rand_networks
import rand_weighted_networks

import matplotlib.pyplot as plt



def fluctuation_init(epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename,type,ending=''):
    choose_beta = lambda type: np.random.normal(avg_beta, epsilon * avg_beta, N) \
        if type == 'gauss' else np.random.gamma((avg_beta / epsilon) ** 2, epsilon ** 2 / avg_beta, N) \
        if type == 'gamma' else netinithomo.bi_beta(N, epsilon, avg_beta)
    Beta = choose_beta(type)
    G = nx.read_gpickle(G_name+ending)
    T = []
    Num_inf = int(x * N)
    R_tot, Rates = netinithomo.intialize_graph(G, N, Num_inf, Beta, Alpha)
    Total_time = 0.0
    T.append(Total_time)

    count = 0
    r = np.random.uniform(0, 1, (bank, 2))
    ######################
    # Main Gillespie Loop
    ######################
    while Num_inf > 0 and Total_time<Time_limit:
        R_norm = np.cumsum(Rates)
        r_pos = R_tot * r[count, 1]
        person = bisect.bisect_left(R_norm, r_pos)

        tau= np.log(1 / r[count, 0]) / R_tot
        Total_time = Total_time + tau

        # rand_networks.draw_infections_nx_g(G,pos_nodes_plt,'frame'+str(count),[person])

        if G.nodes[person]['infected'] == True:
            Num_inf = Num_inf - 1
            Rates[person] = 0.0
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] - G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot - G.nodes[Neighbor]['contact_rate']
                else:
                    Rates[person] = Rates[person] + G.nodes[person]['contact_rate']
            R_tot = R_tot + Rates[person] - Alpha
            G.nodes[person]['infected'] = False
        else:
            Num_inf = Num_inf + 1
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] + G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot + G.nodes[Neighbor]['contact_rate']
            R_tot = R_tot - Rates[person] + Alpha
            Rates[person] = Alpha
            G.nodes[person]['infected'] = True
        count = count + 1
        if count >= bank:
            r = np.random.uniform(0, 1, (bank, 2))
            count = 0
    outfile = open (outfilename + ending + '_Rates.pickle','wb')
    pickle.dump(Rates,outfile)
    outfile.close()
    outfile = open(outfilename + ending + '_R_tot.pickle','wb')
    pickle.dump(R_tot, outfile)
    outfile.close()
    outfile = open(outfilename + ending + '_Num_inf.pickle','wb')
    pickle.dump(Num_inf, outfile)
    outfile.close()
    nx.write_gpickle(G,outfilename + ending + '_G.pickle')
    return 0

def fluctuation_run(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    for run_loop_counter in range(runs):
        T = []
        I = []
        runs_csv=[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)

        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time-T[-1]>=0.1:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
        f = open(outfile + '.csv', "a+")
        l = [T, I, runs_csv,net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0


def fluctuation_run_catastrophe(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,network_number,Beta,factor,duration,time_q):
    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    Beta_org=Beta
    for run_loop_counter in range(runs):
        T = []
        I = []
        runs_csv=[]
        Beta=Beta_org
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        nx.set_node_attributes(G,False,'infected')
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)
        quarntine = False

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:


            if Total_time>time_q and Total_time<=time_q+duration and quarntine == False:
                # start a quarntine if there isn't one already and if the time is right
                Beta = Beta*factor
                R_tot, Rates = netinithomo.inatlize_quarntine_graph(G, G.number_of_nodes(), Alpha, Beta)
                quarntine = True
            elif Total_time>time_q and Total_time>time_q+duration and quarntine==True:
                # stop the quarantine and resume the previous infection rate
                Beta = Beta/factor
                R_tot, Rates = netinithomo.inatlize_quarntine_graph(G, G.number_of_nodes(), Alpha, Beta)
                quarntine = False


            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number,Num_inf]])
        f.close()
    return 0


def fluctuation_run_no_decay(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,network_number,Beta,start_recording_time):
    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    for run_loop_counter in range(runs):
        T = []
        I = []
        I_type_above_avg,I_type_below_avg= [],[]
        runs_csv=[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))

        for l in range(G.number_of_nodes()):
            G.nodes[l]['infected'] = False
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)

        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        num_inf_above_avg, num_inf_below_avg = 0, 0
        for node_type_number in range(G.number_of_nodes()):
            if G.nodes[node_type_number]['infected'] == True:
                if G.nodes[node_type_number]['contact_rate'] > 1.0:
                    num_inf_above_avg = num_inf_above_avg + 1
                else:
                    num_inf_below_avg = num_inf_below_avg + 1
        I_type_above_avg.append(num_inf_above_avg)
        I_type_below_avg.append(num_inf_below_avg)

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True and Num_inf>1:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            elif G.nodes[person]['infected'] == False :
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time-T[-1]>=0.1 and Total_time>=start_recording_time:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
                num_inf_above_avg,num_inf_below_avg=0,0
                for node_type_number in range(G.number_of_nodes()):
                    if G.nodes[node_type_number]['infected'] == True:
                        if G.nodes[node_type_number]['contact_rate']>1.0:
                            num_inf_above_avg = num_inf_above_avg + 1
                        else:
                            num_inf_below_avg = num_inf_below_avg + 1
                I_type_above_avg.append(num_inf_above_avg)
                I_type_below_avg.append(num_inf_below_avg)
        f = open(outfile + '.csv', "a+")
        l = [T, I,I_type_above_avg,I_type_below_avg ,runs_csv,net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0

def fluctuation_run_autocorrelation_hetro_rates(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,time_diff,Beta,start_recording_time):
    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    t0_array = np.linspace(start_recording_time,Time_limit,runs)
    for t0_start_recording in t0_array:
        T = []
        autocorrelation_sus_above_1,autocorrelation_sus_below_1 = [],[]
        runs_csv=[]
        runs_csv.append(t0_start_recording)
        Total_time = 0.0
        # T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))

        for l in range(G.number_of_nodes()):
            G.nodes[l]['infected'] = False
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)


        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit+(t0_start_recording-start_recording_time):
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True and Num_inf>1:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            elif G.nodes[person]['infected'] == False :
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time>=t0_start_recording:
                if len(T)==0:
                    infected_at_t0 = np.array(list(nx.get_node_attributes(G, "infected").values()))
                    T.append(Total_time)
                    autocorrelation_sus_above_1.append(1.0)
                    autocorrelation_sus_below_1.append(1.0)
                    correlation_above_avg,correlation_below_avg = 0.0,0.0
                    for node_type_number in range(G.number_of_nodes()):
                        if G.nodes[node_type_number]['contact_rate'] > 1.0:
                            correlation_above_avg = correlation_above_avg+infected_at_t0[node_type_number] * G.nodes[node_type_number]['infected']
                        else:
                            correlation_below_avg = correlation_below_avg+infected_at_t0[node_type_number] * G.nodes[node_type_number]['infected']

                elif Total_time-T[-1]>=time_diff:
                    correlation_above_avg,correlation_below_avg = 0.0,0.0
                    for node_type_number in range(G.number_of_nodes()):
                        if G.nodes[node_type_number]['contact_rate'] > 1.0:
                            correlation_above_avg = correlation_above_avg+infected_at_t0[node_type_number] * G.nodes[node_type_number]['infected']
                        else:
                            correlation_below_avg = correlation_below_avg+infected_at_t0[node_type_number] * G.nodes[node_type_number]['infected']
                    autocorrelation_sus_above_1.append(correlation_above_avg/len(infected_at_t0))
                    autocorrelation_sus_below_1.append(correlation_below_avg/len(infected_at_t0))
                    runs_csv.append(t0_start_recording)
                    T.append(Total_time)
        f = open(outfile + '.csv', "a+")
        l = [T, autocorrelation_sus_above_1,autocorrelation_sus_below_1,runs_csv]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0

def fluctuation_run_extinction(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1


            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0



def fluctuation_run_extinction_undirected_graph(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_undirected(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1


            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                    else:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0


def fluctuation_run_extinction_DiGraph(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0



def fluctuation_run_extinction_DiGraph_record(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta,start_recording_time,Time_limit):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        T,I,I_type_above_avg,I_type_below_avg,runs_csv=[],[],[],[],[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        for l in range(G.number_of_nodes()):
            G.nodes[l]['infected'] = False
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        num_inf_above_avg, num_inf_below_avg = 0, 0
        for node_type_number in range(G.number_of_nodes()):
            if G.nodes[node_type_number]['infected'] == True:
                if G.nodes[node_type_number]['contact_rate'] > 1.0:
                    num_inf_above_avg = num_inf_above_avg + 1
                else:
                    num_inf_below_avg = num_inf_below_avg + 1
        I_type_above_avg.append(num_inf_above_avg)
        I_type_below_avg.append(num_inf_below_avg)


        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True and Num_inf>1:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            elif G.nodes[person]['infected'] == False:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time - T[-1] >= 0.1 and Total_time >= start_recording_time:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
                num_inf_above_avg, num_inf_below_avg = 0, 0
                for node_type_number in range(G.number_of_nodes()):
                    if G.nodes[node_type_number]['infected'] == True:
                        if G.nodes[node_type_number]['contact_rate'] > 1.0:
                            num_inf_above_avg = num_inf_above_avg + 1
                        else:
                            num_inf_below_avg = num_inf_below_avg + 1
                I_type_above_avg.append(num_inf_above_avg)
                I_type_below_avg.append(num_inf_below_avg)
        f = open(outfile + '.csv', "a+")
        l = [T, I, I_type_above_avg, I_type_below_avg, runs_csv, net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0



# def run_markov_chain_bimodal_network(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta,start_recording_time,Time_limit):
#     G = nx.read_gpickle(infile)
#     seed_nodes = Num_inf
#     for run_loop_counter in range(runs):
#         T,I,runs_csv,tau_list=[],[],[],[]
#         runs_csv.append(run_loop_counter)
#         Total_time = 0.0
#         T.append(Total_time)
#         count = 0
#         Num_inf = seed_nodes
#         r = np.random.uniform(0, 1, (bank, 2))
#         for l in range(G.number_of_nodes()):
#             G.nodes[l]['infected'] = False
#         R_tot, Rates = netinithomo.inatlize_markov_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
#         net_num = []
#         I.append(Num_inf)
#         net_num.append(network_number)
#
#
#         ######################
#         # Main Gillespie Loop
#         ######################
#         while Num_inf > 0 and Total_time<Time_limit:
#             R_norm = np.cumsum(Rates)
#             r_pos = R_norm[-1] * r[count, 1]
#             person = bisect.bisect_left(R_norm, r_pos)
#             tau= np.log(1 / r[count, 0]) / R_norm[-1]
#             Total_time = Total_time + tau
#
#             try:
#                 if G.nodes[person]['infected'] == True:
#                   pass
#             except:
#                   print('Accessing G.node[person][infected] failed value of person is ',person)
#                   if person == G.number_of_nodes():
#                       person =G.number_of_nodes()-1
#             if G.nodes[person]['infected'] == True and Num_inf>1:
#                 Num_inf = Num_inf - 1
#                 # Rates[person] = 0.0
#                 for Neighbor in G.successors(person):
#                     if G.nodes[Neighbor]['infected'] == False:
#                         G.nodes[Neighbor]['infected_neighbors'] = G.nodes[Neighbor]['infected_neighbors'] - 1
#                 for Neighbor in G.predecessors(person):
#                     if G.nodes[Neighbor]['infected'] == True:
#                         G.nodes[person]['infected_neighbors'] = G.nodes[person]['infected_neighbors'] + 1
#                 G.nodes[person]['infected'] = False
#             elif G.nodes[person]['infected'] == False:
#                 Num_inf = Num_inf + 1
#                 G.nodes[person]['infected_neighbors'] = 0
#                 for Neighbor in G.successors(person):
#                     if G.nodes[Neighbor]['infected'] == False:
#                         G.nodes[Neighbor]['infected_neighbors'] = G.nodes[Neighbor]['infected_neighbors'] + 1
#                 # Rates[person] = Alpha
#                 G.nodes[person]['infected'] = True
#             for i in G.nodes():
#                 if G.nodes[i]['infected'] == False:
#                     Rates[i] = 1-(1-Beta*tau)**G.nodes[i]['infected_neighbors']
#                 else:
#                     Rates[i] = Alpha * tau
#             count = count + 1
#             if count >= bank:
#                 r = np.random.uniform(0, 1, (bank, 2))
#                 count = 0
#             if Total_time - T[-1] >= 0.1 and Total_time >= start_recording_time:
#                 I.append(Num_inf)
#                 T.append(Total_time)
#                 tau_list.append(tau)
#                 net_num.append(network_number)
#                 runs_csv.append(run_loop_counter)
#         f = open(outfile + '.csv', "a+")
#         l = [T, I, runs_csv, net_num,tau_list]
#         l = zip(*l)
#         with f:
#             writer = csv.writer(f)
#             writer.writerows(l)
#         f.close()
#         with open('time_'+ outfile + '.npy', 'wb') as f:
#             np.save(f, np.array(T))
#         f.close()
#         with open('infected_'+ outfile + '.npy', 'wb') as f:
#             np.save(f, np.array(I))
#         f.close()
#     return 0



def run_markov_chain_bimodal_network(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta,start_recording_time,Time_limit,dt):
    pmu = Alpha * dt
    pbeta = Beta * dt
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        T,I,runs_csv=[],[],[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, bank)
        for l in range(G.number_of_nodes()):
            G.nodes[l]['infected'] = False
        R_tot, Rates = netinithomo.inatlize_markov_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        ######################
        # Synchronous updating
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            Total_time = Total_time + dt
            infected_neighbs_old = [G.nodes[i]['infected_neighbors'] for i in range(G.number_of_nodes())]
            for person in range(G.number_of_nodes()):
                if G.nodes[person]['infected'] == True:
                    if r[count]<pmu and Num_inf>1:
                        G.nodes[person]['infected'] = False
                        Num_inf = Num_inf - 1
                        for Neighbor in G.successors(person):
                            # if G.nodes[Neighbor]['infected'] == False:
                            G.nodes[Neighbor]['infected_neighbors'] = G.nodes[Neighbor]['infected_neighbors'] - 1
                else:
                    if r[count]< 1-(1-pbeta)**infected_neighbs_old[person]:
                        Num_inf = Num_inf + 1
                        G.nodes[person]['infected'] = True
                        for Neighbor in G.successors(person):
                            # if G.nodes[Neighbor]['infected'] == False:
                            G.nodes[Neighbor]['infected_neighbors'] = G.nodes[Neighbor]['infected_neighbors'] + 1
                count = count + 1
                if count >= bank:
                    r = np.random.uniform(0, 1, bank)
                    count = 0
                if Total_time - T[-1] >= 0.1 and Total_time >= start_recording_time:
                    I.append(Num_inf)
                    T.append(Total_time)
                    net_num.append(network_number)
                    runs_csv.append(run_loop_counter)
        f = open(outfile + '.csv', "a+")
        l = [T, I, runs_csv, net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
        with open('time_'+ outfile + '.npy', 'wb') as f:
            np.save(f, np.array(T))
        f.close()
        with open('infected_'+ outfile + '.npy', 'wb') as f:
            np.save(f, np.array(I))
        f.close()
    return 0



def markov_extinction_bimodal_network(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_markov_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                # Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        G.nodes[Neighbor]['infected_neighbors'] = G.nodes[Neighbor]['infected_neighbors'] - 1
                        R_tot = R_tot - Rates[Neighbor]
                        Rates[Neighbor] = 1-(1-Beta*tau)**G.nodes[Neighbor]['infected_neighbors']
                        R_tot = R_tot + Rates[Neighbor]
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        G.nodes[person]['infected_neighbors'] = G.nodes[person]['infected_neighbors'] + 1
                Rates[person] = 1-(1-Beta*tau)**G.nodes[person]['infected_neighbors']
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            elif G.nodes[person]['infected'] == False:
                Num_inf = Num_inf + 1
                G.nodes[person]['infected_neighbors'] = 0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        G.nodes[Neighbor]['infected_neighbors'] = G.nodes[Neighbor]['infected_neighbors'] + 1
                        R_tot = R_tot - Rates[Neighbor]
                        Rates[Neighbor] = 1-(1-Beta*tau)**G.nodes[Neighbor]['infected_neighbors']
                        R_tot = R_tot + Rates[Neighbor]
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile + '.csv', "a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
        with open('ext_'+ outfile + '.npy', 'wb') as f:
            np.save(f, np.array(Total_time))
        f.close()
        with open('net_'+ outfile + '.npy', 'wb') as f:
            np.save(f, np.array(network_number))
        # with open('tau_'+ outfile + '.npy', 'wb') as f:
        #     np.save(f, np.array(tau))
        f.close()
    return 0

def fluctuation_run_extinction_weighted_graph(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_weighted_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*G[person][Neighbor]['weight']
                        R_tot = R_tot - Beta*G[person][Neighbor]['weight']
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta*G[Neighbor][person]['weight']
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*G[person][Neighbor]['weight']
                        R_tot = R_tot + Beta*G[person][Neighbor]['weight']
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0


def fluctuation_weighted_run_infected_no_decay(Alpha,bank,outfile,infile,runs,Num_inf,network_number,
                                               Beta,start_recording_time,Time_limit):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        T,I=[],[]
        runs_csv=[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        nx.set_node_attributes(G, False, 'infected')
        R_tot, Rates = netinithomo.inatlize_inf_weighted_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True and Num_inf>1:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*G[person][Neighbor]['weight']
                        R_tot = R_tot - Beta*G[person][Neighbor]['weight']
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta*G[Neighbor][person]['weight']
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*G[person][Neighbor]['weight']
                        R_tot = R_tot + Beta*G[person][Neighbor]['weight']
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time-T[-1]>=0.1 and Total_time>=start_recording_time:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
        f = open(outfile + '.csv', "a+")
        l = [T, I,runs_csv,net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0


def fluctuation_init_track_lam(epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename):
    lam_plus, lam_minus = avg_beta*(1+epsilon), avg_beta*(1-epsilon)
    Beta = netinithomo.bi_beta(N, epsilon, avg_beta)
    G=nx.read_gpickle(G_name)
    T = []
    Num_inf = int(x * N)
    R_tot, Rates,Num_inf_plus, Num_inf_minus = netinithomo.intialize_graph_track_bimodal(G, N, Num_inf, Beta, Alpha,lam_plus,lam_minus)
    Total_time = 0.0
    T.append(Total_time)

    count = 0
    r = np.random.uniform(0, 1, (bank, 2))
    ######################
    # Main Gillespie Loop
    ######################
    while Num_inf > 0 and Total_time<Time_limit:
        R_norm = np.cumsum(Rates)
        r_pos = R_tot * r[count, 1]
        person = bisect.bisect_left(R_norm, r_pos)

        tau= np.log(1 / r[count, 0]) / R_tot
        Total_time = Total_time + tau

        # rand_networks.draw_infections_nx_g(G,pos_nodes_plt,'frame'+str(count),[person])

        if G.nodes[person]['infected'] == True:
            Num_inf = Num_inf - 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus -1
            else:
                Num_inf_minus = Num_inf_minus - 1
            Rates[person] = 0.0
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] - G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot - G.nodes[Neighbor]['contact_rate']
                else:
                    Rates[person] = Rates[person] + G.nodes[person]['contact_rate']
            R_tot = R_tot + Rates[person] - Alpha
            G.nodes[person]['infected'] = False
        else:
            Num_inf = Num_inf + 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus + 1
            else:
                Num_inf_minus = Num_inf_minus + 1
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] + G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot + G.nodes[Neighbor]['contact_rate']
            R_tot = R_tot - Rates[person] + Alpha
            Rates[person] = Alpha
            G.nodes[person]['infected'] = True
        count = count + 1
        if count >= bank:
            r = np.random.uniform(0, 1, (bank, 2))
            count = 0
    outfile = open (outfilename+'_Rates.pickle','wb')
    pickle.dump(Rates,outfile)
    outfile.close()
    outfile = open(outfilename + '_R_tot.pickle','wb')
    pickle.dump(R_tot, outfile)
    outfile.close()
    outfile = open(outfilename+'_Num_inf.pickle','wb')
    pickle.dump(Num_inf, outfile)
    outfile.close()
    outfile = open(outfilename+'_Num_inf_plus.pickle','wb')
    pickle.dump(Num_inf_plus, outfile)
    outfile.close()
    outfile = open(outfilename+'_Num_inf_minus.pickle','wb')
    pickle.dump(Num_inf_minus, outfile)
    outfile.close()
    nx.write_gpickle(G,outfilename+'_G.pickle')
    return 0


def fluctuation_run_track_lam(Alpha,Time_limit,bank,outfile,infile,epsilon,avg_beta):
    lam_plus, lam_minus = avg_beta*(1+epsilon), avg_beta*(1-epsilon)
    G=nx.read_gpickle(infile+'_G.pickle')
    with open(infile+'_Num_inf.pickle','rb') as pickle_file:
        Num_inf = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile+'_Num_inf_plus.pickle','rb') as pickle_file:
        Num_inf_plus = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile+'_Num_inf_minus.pickle','rb') as pickle_file:
        Num_inf_minus = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile + '_R_tot.pickle','rb') as pickle_file:
        R_tot = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile + '_Rates.pickle','rb') as pickle_file:
        Rates = pickle.load(pickle_file)
    pickle_file.close()
    T, I, I_plus, I_minus = [], [], [], []
    I.append(Num_inf)
    I_plus.append(Num_inf_plus)
    I_minus.append(Num_inf_minus)
    Total_time = 0.0
    T.append(Total_time)

    count = 0
    r = np.random.uniform(0, 1, (bank, 2))
    ######################
    # Main Gillespie Loop
    ######################
    while Num_inf > 0 and Total_time<Time_limit:
        R_norm = np.cumsum(Rates)
        r_pos = R_tot * r[count, 1]
        person = bisect.bisect_left(R_norm, r_pos)

        tau= np.log(1 / r[count, 0]) / R_tot
        Total_time = Total_time + tau

        # rand_networks.draw_infections_nx_g(G,pos_nodes_plt,'frame'+str(count),[person])

        if G.nodes[person]['infected'] == True:
            Num_inf = Num_inf - 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus - 1
            else:
                Num_inf_minus = Num_inf_minus -1
            Rates[person] = 0.0
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] - G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot - G.nodes[Neighbor]['contact_rate']
                else:
                    Rates[person] = Rates[person] + G.nodes[person]['contact_rate']
            R_tot = R_tot + Rates[person] - Alpha
            G.nodes[person]['infected'] = False
        else:
            Num_inf = Num_inf + 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus + 1
            else:
                Num_inf_minus = Num_inf_minus +1
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] + G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot + G.nodes[Neighbor]['contact_rate']
            R_tot = R_tot - Rates[person] + Alpha
            Rates[person] = Alpha
            G.nodes[person]['infected'] = True
        count = count + 1
        if count >= bank:
            r = np.random.uniform(0, 1, (bank, 2))
            count = 0
        if Total_time-T[-1]>=1:
            I.append(Num_inf)
            I_plus.append(Num_inf_plus)
            I_minus.append(Num_inf_minus)
            T.append(round(Total_time))
    f = open(outfile+'_total.csv',"a+")
    with f:
        writer = csv.writer(f)
        writer.writerows([I])
    f.close()
    f = open(outfile+'infected_plus.csv',"a+")
    with f:
        writer = csv.writer(f)
        writer.writerows([I_plus])
    f.close()
    f = open(outfile+'infected_minus.csv',"a+")
    with f:
        writer = csv.writer(f)
        writer.writerows([I_minus])
    f.close()
    return 0


def actasmain():
    Epsilon_sus = [0.0]
    Epsilon_inf = [0.0]
    Epsilon=[0.0]
    N = 500
    k = 200
    x = 0.2
    eps_din,eps_dout = 0.0,0.0
    eps_sus,eps_lam = 0.0,0.0
    Num_inf = int(x * N)
    Alpha = 1.0
    susceptibility = 'bimodal'
    infectability = 'bimodal'
    directed_model='gauss_c'
    prog = 'q' #can be either 'i' for the inatilization and reaching eq state or 'r' for running and recording fluc
    Lam = 1.5
    Time_limit = 200
    Start_recording_time = 0.0
    Beta_avg =Alpha* Lam / k
    Num_different_networks= 1
    Num_inital_conditions= 1
    bank = 1000000
    parts = 1
    graphname  = 'GNull'
    count = 0
    susceptibility_avg = 1.0
    infectability_avg = 1.0
    foldername ='base'
    graphname  = 'GNull'
    outfile ='o'
    sus_inf_correlation = 'a'
    # Beta = Beta_avg / (1 + Epsilon_sus[0] * Epsilon_inf[0])
    # Beta = Beta_avg / (1 - Epsilon_sus[0] * Epsilon_inf[0]) if sus_inf_correlation is 'a' else Beta_avg / (
    #             1 + Epsilon_sus[0] * Epsilon_inf[0])
    Beta = Beta_avg / (1 + eps_lam * eps_sus)
    factor, duration, time_q = 0.5, 15.0, 100.0
    interaction_strength = 2.0
    auto_correlation_time_diff = 0.1
    dt_discrite = 0.01


    # G = nx.random_regular_graph(k, N)
    # G = nx.complete_graph(N)
    # beta_inf, beta_sus = netinithomo.general_beta(N, eps_lam, eps_sus, directed_model, k)
    # beta_inf, beta_sus = netinithomo.bi_beta_correlated(N, eps_lam, eps_sus, 1.0)
    # G = netinithomo.intalize_lam_graph(G, N, beta_sus, beta_inf)
    d1_in, d1_out, d2_in, d2_out = int(k * (1 - eps_din)), int(k * (1 - eps_dout)), int(k * (1 + eps_din)), int(
        k * (1 + eps_dout))
    G = rand_networks.random_bimodal_directed_graph(d1_in, d1_out, d2_in, d2_out, N)
    G = netinithomo.set_graph_attriubute_markov(G)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)

    # choose_beta = lambda net_dist, avg, epsilon: np.random.normal(avg, epsilon * avg, N) \
    #     if net_dist == 'gauss' else np.random.gamma((avg / epsilon) ** 2, epsilon ** 2 / avg, N) \
    #     if net_dist == 'gamma' else np.zeros(N) if net_dist == 'z' else np.ones(
    #     N) if net_dist == 'ones' else netinithomo.bi_beta(N, epsilon, avg)

    # beta_sus = choose_beta(susceptibility, susceptibility_avg,Epsilon_sus[0])
    # beta_inf = choose_beta(infectability, infectability_avg,Epsilon_inf[0])
    # G = netinithomo.intalize_lam_graph(G, N, beta_sus, beta_inf)
    # beta_inf, beta_sus = netinithomo.triangular_beta(N, Epsilon_inf[0], Epsilon_sus[0], sus_inf_correlation)
    # G = rand_networks.configuration_model_directed_graph(directed_model,eps_dout,eps_din,k,N)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # G = rand_networks.configuration_model_undirected_graph(directed_model,eps_dout,eps_din,k,N)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # G = rand_networks.configuration_model_undirected_graph(0.1,k,N)
    # G = rand_networks.jason_graph('jason_trans_file_no_degree.csv')
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # G = netinithomo.intalize_lam_graph(G, N, np.ones(N), np.ones(N))
    # Beta = Beta_avg / (1 + eps_din * eps_dout)
    # Beta = Beta/np.mean([G.in_degree(n) for n in G.nodes()])
    # Beta = Lam/np.mean([G.in_degree(n) for n in G.nodes()])
    n=0

    # G = rand_weighted_networks.bulid_weighted_network(N, k, eps_din, eps_dout, interaction_strength)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n) + '.pickle'
    nx.write_gpickle(G, infile)
    outfile = 'o_eps_in' + str(np.abs(eps_din)).replace('.', '') + 'eps_dout' + str(np.abs(eps_dout)).replace('.', '')
    # k_avg_graph, eps_in_graph, eps_out_graph = rand_weighted_networks.weighted_epsilon(G)
    # Beta_graph = Lam / k_avg_graph
    # Beta = Beta_graph / (1 + np.sign(eps_din) * eps_in_graph * np.sign(eps_dout) * eps_out_graph)

    # fluctuation_run_extinction_weighted_graph(Alpha, bank, outfile, infile,
    #                                           Num_inital_conditions, Num_inf, n,
    #                                           Beta)
    # fluctuation_weighted_run_infected_no_decay(Alpha,bank,outfile,infile,Num_inital_conditions,Num_inf,n,Beta,Start_recording_time,Time_limit)
    # fluctuation_run_autocorrelation_hetro_rates(Alpha,Time_limit,bank,outfile,infile,Num_inital_conditions,Num_inf,auto_correlation_time_diff,Beta,Start_recording_time)
    run_markov_chain_bimodal_network(Alpha,bank,outfile,infile,Num_inital_conditions,Num_inf,n,Beta,Start_recording_time,Time_limit,dt_discrite)
    # markov_extinction_bimodal_network(Alpha,bank,outfile,infile,Num_inital_conditions,Num_inf,n,Beta)
    # nx.write_gpickle(G, graphname)
    # infile = graphname + '_' + str(epsilon).replace('.', '') + '_' + str(n)+'.pickle'
    # infile=graphname
    # fluctuation_run_extinction(Alpha,bank,outfile,infile,Num_inital_conditions,Num_inf,1,Beta)
    # fluctuation_run_catastrophe(Alpha,Time_limit,bank,outfile,infile,Num_inital_conditions,Num_inf,n,Beta,factor,duration,time_q)
    # fluctuation_run_no_decay(Alpha, Time_limit, bank, outfile, infile, Num_inital_conditions,
    #                 Num_inf, 1, Beta,Start_recording_time)
    # fluctuation_run_extinction_undirected_graph(Alpha, bank, outfile, infile, Num_inital_conditions,
    #                            Num_inf, 1, Beta)
    # fluctuation_run_extinction_DiGraph(Alpha, bank, outfile, infile, Num_inital_conditions,
    #                                    Num_inf, 1, Beta)

if __name__ == '__main__':
    submit = True
    if submit==False:
        actasmain()
    else:
         if sys.argv[1]=='i':
             # Parmeters order: epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename,type
             fluctuation_init(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]),sys.argv[10],sys.argv[11])
         elif sys.argv[1]=='r':
             # Parameters order: Alpha,Time_limit,bank,outfile,infile
             fluctuation_run(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),float(sys.argv[10]))
         elif sys.argv[1]=='b':
             #Parameters: Alpha,Time_limit,bank,outfile,infile,epsilon,avg_beta
             fluctuation_run_track_lam(float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6],float(sys.argv[7]),float(sys.argv[8]))
         elif sys.argv[1]=='bi':
             fluctuation_init_track_lam(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]),sys.argv[6],float(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),sys.argv[10])
         elif sys.argv[1]=='si':
             # Parmeters order: epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename,type,ending
             fluctuation_init(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]),sys.argv[10],sys.argv[11],sys.argv[12])
         elif sys.argv[1]=='sr':
             # Parameters order: Alpha,Time_limit,bank,outfile,infile,ending
             fluctuation_run(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],sys.argv[7])
         elif sys.argv[1]=='e' or sys.argv[1]=='ec' or sys.argv[1]=='ac' or sys.argv[1] == 'g' or sys.argv[1] == 'cr' or sys.argv[1] == 'f':
             # Parameters order: Alpha,bank,outfile,infile
             fluctuation_run_extinction(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),float(sys.argv[9]))
         elif sys.argv[1] == 'ri' or sys.argv[1] == 'rg':
             # Parameters order: Alpha,bank,outfile,infile
             fluctuation_run_no_decay(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]))
         elif sys.argv[1] == 'bd' or sys.argv[1] == 'co':
             # Parameters order: Alpha,bank,outfile,infile
             fluctuation_run_extinction_DiGraph(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),float(sys.argv[9]))
         elif sys.argv[1] == 'q':
             # activate quarntine at a given time
             fluctuation_run_catastrophe(float(sys.argv[2]), float(sys.argv[3]) ,int(sys.argv[4]), sys.argv[5],
                                         sys.argv[6],int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]),float(sys.argv[12]),float(sys.argv[13]))
         elif sys.argv[1] == 'wn':
             # Simulate an extinction process on a graph where there is connection between next nearest neighbors
             fluctuation_run_extinction_weighted_graph(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],
                                                       int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]),
                                                       float(sys.argv[9]))
         elif sys.argv[1] == 'wnr' or sys.argv[1]=='an' :
             # Mesaure the number of infected nodes (endemic state) when the network is with next nearest interactions
             fluctuation_weighted_run_infected_no_decay(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],
                                                int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]),
                                                float(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]))
         elif sys.argv[1] == 'au':
             # Finds autocorrelation for two groups sus below 1 and sus above 1
             fluctuation_run_autocorrelation_hetro_rates(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],
                                      int(sys.argv[7]), int(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10]),
                                      float(sys.argv[11]))
         elif sys.argv[1] == 'br':
            # Run and record number of infected for bimodal degree graph
             fluctuation_run_extinction_DiGraph_record(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],
                                                int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]),
                                                float(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]))
         elif sys.argv[1] == 'mcr' or sys.argv[1]== 'mcer':
            # Run and record number of infected for bimodal degree graph with infection process 1-(1-lam)^m
             run_markov_chain_bimodal_network(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],
                                                int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]),
                                                float(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]),
                                              float(sys.argv[12]))
         elif sys.argv[1] == 'mc':
             # Extinction times for bimodal degree graph with infection process 1-(1-lam)^m
             markov_extinction_bimodal_network(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],
                                              int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]),
                                              float(sys.argv[9]))

