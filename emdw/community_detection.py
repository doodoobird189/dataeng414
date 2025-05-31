
# Author: Corné van Daalen
# Date: 18 April 2025

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

#========================= function generateGraph =========================
# This function generates a graph from the specified stochastic blockmodel (SBM)
# Arguments:
# R: the number of communities
# N: the number of nodes
# c: categorical distribution over the communities of each node (leave out for uniform)
# w_d: the probability of generating an edge between nodes of different communities
# w_s: the probability of generating an edge between nodes of the same community
# node_filename: name of the file in which the community ground truth should be stored
#   in the form of one-hot vectors (i.e., categorical distributions) for each node
# edge_filename: name of the file in which the adjacency matrix should be stored
def generateGraph(R, N, w_s, w_d, node_filename, edge_filename, c=None):
    # Set up parameters:
    if (c == None):
        c = np.ones(R)/R # uniform community distribution
    W = np.ones((R,R))*w_d + np.diag(np.ones(R))*(w_s-w_d) # edge probability matrix

    # Generate graph:
    V = np.random.choice(R,N,p=c) # sample a community for each node
    # Sample edges and create adjacency matrix:
    E = np.zeros((N,N))
    for i in range(N):
        for j in range(i+1,N):
            E[i,j] = np.random.choice(2,1,p=[1.-W[V[i],V[j]], W[V[i],V[j]]])
            E[j,i] = E[i,j]
    # Convert node communities to one-hot vectors (i.e., categorical distributions):
    V_onehot = np.zeros((N,R))
    V_onehot[np.arange(N), V] = 1
    # Save the node community distributions:
    np.savetxt(node_filename, V_onehot, fmt='%.6f')
    # Save the sampled edges as an adjacency matrix:
    np.savetxt(edge_filename, E, fmt='%d')

    
#========================= function visualiseRawData =========================
# This function visualises the raw data; i.e., it draws the nodes and the observed edges
# between these nodes, but without any clustering of the nodes
# Arguments:
# edge_filename: the name of the file that contains the adjacency matrix with observed edges
def visualiseRawData(edge_filename):
    E = np.loadtxt(edge_filename, dtype=int) # read adjacency matrix
    N = E.shape[0] # number of nodes (i.e., no. of rows or columns)
    devs = np.random.randn(2,N) # random positions for nodes drawn from Gaussian distribution
    plt.figure()
    # Plot an edges between nodes according to adjacency matrix
    for i in range(N):
        for j in range(i+1,N):
            if (E[i,j] == 1):
                plt.plot([devs[0,i],devs[0,j]], [devs[1,i],devs[1,j]], 'lightgrey')
    # Plot nodes:
    plt.plot(devs[0],devs[1],'o',color='grey')
    for i in range(N):
        plt.annotate(i, (devs[0,i],devs[1,i]+0.05), ha='center', va='bottom')
    # Format figure:
    plt.title("Raw data")
    plt.gca().set_axis_off()
    plt.gca().set_aspect('equal')
    plt.show()
    
    
#========================= function visualiseClusteredNodes =========================
# This function visualises clustered nodes: it draws nodes that belong to the same community close
# together, and colour them similarly; it also draws the observed edges between the nodes
# Arguments:
# node_distribution_filename: the name of the file containing the categorical distributions
#    over the communities for each node
# edge_filename: the name of the file that contains the adjacency matrix with observed edges
def visualiseClusteredNodes(node_distribution_filename, edge_filename):
    V_distribs = np.loadtxt(node_distribution_filename) # read the node community distributions
    N = V_distribs.shape[0] # number of nodes (i.e., number of rows)
    R = V_distribs.shape[1] # number of communities (i.e., number of colums)
    E = np.loadtxt(edge_filename, dtype=int) # read adjacency matrix
    # The positions of the centres of each community:
    centres = np.array([np.cos(np.arange(R)/R*2.*np.pi),np.sin(np.arange(R)/R*2.*np.pi)])
    # Determine position of each node according to its community distribution:
    positions = np.zeros((2,N))
    for i in range(N):
        positions[:,i] = 5.*R/np.pi*np.dot(centres,V_distribs[i])
    positions = positions + np.random.randn(2,N)
    # Determine colour of each node according to its community distribution:
    tab_colours = np.array([mcolors.to_rgb(key) for key, value in mcolors.TABLEAU_COLORS.items()])
    node_colours = []
    for i in range(N):
        node_colour = np.zeros(3)
        for j in range(R):
            node_colour = node_colour + V_distribs[i,j]*tab_colours[j%10]
        node_colours.append(tuple(node_colour))
    plt.figure()
    # Plot an edges between nodes according to adjacency matrix
    for i in range(N):
        for j in range(i+1,N):
            if (E[i,j] == 1):
                plt.plot([positions[0,i],positions[0,j]], [positions[1,i],positions[1,j]], 'lightgrey')
    # Plot nodes
    for i in range(N):
        plt.plot(positions[0,i],positions[1,i],'o',color=node_colours[i])
        plt.annotate(i, (positions[0,i],positions[1,i]+0.05), ha='center', va='bottom')
    # Format figure:
    plt.title("Clustered data")
    plt.gca().set_axis_off()
    plt.gca().set_aspect('equal')
    plt.show()


generateGraph(R=3, N=50, w_s=0.7, w_d=0.2, node_filename="test2truth.txt", edge_filename="test2.txt", c=None)