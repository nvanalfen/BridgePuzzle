import numpy as np
import pandas as pd
import itertools

# Class representing the node to be connected in the bridge connection puzzle
class Node:
    node_id = 0
    max_single_edges = 2
    
    def __init__(self, max_degree, location):
        self.id = Node.node_id
        Node.node_id += 1
        self.max_degree = max_degree
        self.loc = location
        self.neighbors = []                     # list of connected node ids
    
    # Returns the degree of the node
    def degree(self):
        return len(self.neighbors)
    
    # Returns the number of links to a designated node
    def links_to(self, node_id):
        return sum( [ n in self.neighbors for n in self.neighbors if n == node_id ] )
    
    # Returns the number of links that can still be made to the node
    def remaining_degree(self):
        return self.max_degree - self.degree()
    
    # True if no more links can be made
    def is_filled(self):
        return self.degree() == self.max_degree
    
    # Adds a node id to its neighbors
    def add_edge(self, new_id):
        if sum( [ el == new_id for el in self.neighbors ] ) >= Node.max_single_edges:
            return False
        self.neighbors.append( new_id )
        return True

# Class representing an edge between two nodes. Mostly for administrative and visual purposes
class Edge:
    max_weight = 2                          # Typically the maximum number of links between two nodes
    def __init__(self, id_A, id_B, nodes):
        self.end_points = ( min([id_A, id_B]), max([id_A, id_B]) )
        self.weight = 1
        self.horizontal = False
        self.set_direction(nodes)
    
    # Returns true if the id given is one of the endpoints
    def contains(self, node_id):
        return node_id in self.end_points
    
    # If the edge is already present, but we want to add another link, increase the weight
    def increase_weight(self):
        if self.weight >= Edge.max_weight:
            return False
        self.weight += 1
        return True
    
    # Purely for visual purposes when printing. Which way is the link oriented
    def set_direction(self, nodes):
        loc_A = nodes[ self.end_points[0] ].loc
        loc_B = nodes[ self.end_points[1] ].loc
        if loc_A[0] == loc_B[0]:
            self.horizontal = True
        else:
            self.horizontal = False
    
    # Purely for visual purposes when printing. Gives the symbol closest to what the link would look like
    def get_symbol(self):
        if self.horizontal:
            if self.weight == 1:
                return "-"
            elif self.weight == 2:
                return "="
            else:
                return str(self.weight)+"H"
        else:
            if self.weight == 1:
                return "|"
            elif self.weight == 2:
                return "||"
            else:
                return str(self.weight)+"V"

# The class representing a bridge connection puzle
class BridgePuzzle:
    def __init__(self, file_name=None, autosolve=True):
        self.heaviest_link = 2
        self.nodes = {}
        self.grid = None
        self.initialized = False
        # Loads a grid if a file name has been given
        if not file_name is None:
            self.read_file(file_name)
        # If the grid has been initialized and we want to solve the grid right away, do it
        if autosolve and self.initialized:
            self.solve()
    
    # Reads the csv file containing the representation of the puzzle
    def read_file(self, f_name):
        layout = np.array( pd.read_csv(f_name, header=None) )
        self.grid = np.repeat(None, layout.shape[0]*layout.shape [1] ).reshape( layout.shape )
        for row in range(len(layout)):
            for col in range(len(layout[row])):
                if layout[row, col] > 0:
                    n = Node( layout[row, col], (row,col) )
                    self.nodes[ n.id ] = n
                    self.grid[row, col] = n
        
        self.initialized = True
    
    # Adds an edge between two nodes and updates the gridto show this edge
    def add_edge(self, id_A, id_B):
        self.nodes[id_A].add_edge(id_B)
        self.nodes[id_B].add_edge(id_A)
        x_step = ( self.nodes[id_B].loc[0] - self.nodes[id_A].loc[0] )
        y_step = ( self.nodes[id_B].loc[1] - self.nodes[id_A].loc[1] )
        if x_step != 0:
            x_step = int(x_step/abs(x_step))
        if y_step != 0:
            y_step = int(y_step/abs(y_step))
        x = int(self.nodes[id_A].loc[0] + x_step)
        y = int(self.nodes[id_A].loc[1] + y_step)
        while x != self.nodes[id_B].loc[0] or y != self.nodes[id_B].loc[1]:
            if self.grid[x,y] is None:
                self.grid[x,y] = Edge(id_A, id_B, self.nodes)
            else:
                self.grid[x,y].increase_weight()
            x += x_step
            y += y_step
    
    # Calls sweep (which solves with all available information) iteratively
    # until no new links have been added
    def solve(self, of_name="Output.csv"):
        pre_degree = -1
        post_degree = sum( [ self.nodes[n].degree() for n in self.nodes ] )
        while pre_degree != post_degree:
            pre_degree = post_degree
            self.sweep()
            post_degree = sum( [ self.nodes[n].degree() for n in self.nodes ] )
        self.write_grid(of_name)
    
    # Performs a single iteration through the grid, trying to connect all possible links
    def sweep(self):
        for node_id in self.nodes:
            neighbors = self.find_neighbors(node_id)
            potential_links = self.get_potential_links(node_id, neighbors)
            self.choose_links(node_id, potential_links)
    
    # Given a list of all the tchnically possible sets of links, eliminate those that would isolate
    # an incomplete subgraph
    def eliminate_impossible_link_sets(self, central_id, links):
        potentials = []
        for link_set in links:
            if not self.isolates_incomplete_graph(central_id, link_set):
                potentials.append( link_set )
        return potentials
    
    # Get the subgraph made of all nodes connected to the chosen node
    def get_subgraph(self, node_id):
        nodes = set()
        nodes.add(node_id)
        check = [node_id]
        while len(check) > 0:
            new_check = []
            for n_id in check:
                new_check += [ s for s in self.nodes[ n_id ].neighbors if not s in nodes ]
            for n_id in new_check:
                nodes.add(n_id)
            check = new_check
        return nodes
    
    # Given a subgraph, dinds how many links the subgraph can currently make
    def open_links_in_subgraph(self, graph):
        count = 0
        for node_id in graph:
            count += self.nodes[ node_id ].remaining_degree()
        return count
    
    # Checks to see if connecting the central_id node to the others in link_set will create an isolated subgraph
    # Context will be checked in another function (an isolated subgraph is allowed if it contains all nodes)
    def creates_isolated_subgraph(self, central_id, link_set):
        remaining_degrees = 0           # Current possible links in the subgraphs
        added_degree = 2*len(link_set)  # Multiply by 2 because a link from 0->1 is also a link from 1->0
        subgraphs = self.generate_unique_subgraphs([central_id]+list(link_set))
        for graph in subgraphs:
            remaining_degrees += self.open_links_in_subgraph(graph)
        
        # If we add as many links as we have left, we have created an isolated graph
        return added_degree >= remaining_degrees        # Should never be more, but added in case
    
    # Returns true if connecting the set of links given will create an isolated subgraph that 
    # does not contain all of the nodes
    def isolates_incomplete_graph(self, central_id, link_set):
        subgraphs = self.generate_unique_subgraphs([central_id]+list(link_set))
        N = sum( [ len(g) for g in subgraphs ] )            # Number of nodes in the resulting subgraph
        complete = N == len(self.nodes)                     # True if all nodes have been included
        isolated = self.creates_isolated_subgraph(central_id, link_set)
        return (not complete) and isolated
    
    # Given a list of node ids, generate all unique subgraphs represented
    def generate_unique_subgraphs(self, node_ids):
        subgraphs = []
        for node_id in node_ids:
            if not any( [ node_id in g for g in subgraphs ] ):
                subgraphs.append( self.get_subgraph(node_id) )
        return subgraphs
    
    # Find the neighbors of a goven node. Up to four possible (one for each direction)
    def find_neighbors(self, node_id):
        if not self.initialized:
            return None
        neighbors = []
        start = self.nodes[ node_id ]
        neighbors.append( self.find_next_node( node_id, start.loc[0]-1, start.loc[1], -1, 0 ) )
        neighbors.append( self.find_next_node( node_id, start.loc[0]+1, start.loc[1], 1, 0 ) )
        neighbors.append( self.find_next_node( node_id, start.loc[0], start.loc[1]-1, 0, -1 ) )
        neighbors.append( self.find_next_node( node_id, start.loc[0], start.loc[1]+1, 0, 1 ) )
        neighbors = [ n for n in neighbors if not n is None ]
        return neighbors
    
    # Move one space farther to one direction until you either find a node or nothing
    def find_next_node(self, node_id, x, y, x_change, y_change):
        if not self.initialized or x < 0 or x >= self.grid.shape[0] or y < 0 or y >= self.grid.shape[1]:
            return None
        
        obj = self.grid[x,y]
        
        if isinstance( obj, Edge ) and ( obj.weight >= self.heaviest_link or not obj.contains(node_id) ):
            return None
        elif isinstance( obj, Node ):
            return obj.id
        return self.find_next_node(node_id, x+x_change, y+y_change, x_change, y_change)
    
    # Get list of all possible sets of links for a central node to its neighbor nodes
    def get_potential_links(self, central_id, neighbors):
        potential = []
        central = self.nodes[ central_id ]
        for node in neighbors:
            count = min( [ central.remaining_degree(), self.nodes[node].remaining_degree(), self.heaviest_link - central.links_to(node) ] )
            for i in range(count):
                potential.append(node)
        
        return list( itertools.combinations(potential, central.remaining_degree()) )
    
    # Of all the potential sets of links, choose those that must be present (if any)
    def choose_links(self, central_id, links):
        # First, only consider those links which would not create an incomplete subgraph
        potential_links = self.eliminate_impossible_link_sets(central_id, links)
        if len(potential_links) == 0:
            # Only one possible set of links, connect them
            self.link_set(central_id, potential_links[0])
        else:
            # Find if there are any links that are present in every possibility
            common_links = self.find_common_links(potential_links)
            self.link_set(central_id, common_links)
            
    # Create links from the node at central_id to every node listed in links
    def link_set(self, central_id, links):
        for link in links:
            self.add_edge(central_id, link)
    
    # Finds links that are present in all sets
    def find_common_links(self, links):
        seed_set = set( links[0] )
        for link_set in links:
            seed_set = seed_set & set(link_set)
        
        common_links = list(seed_set)
        # There may be multi links, but because we're using a set, those multiplicities will be ignored
        # Here, we call the same 
        if len(common_links) > 0:
            new_links = [ list(el) for el in links ]
            for common in common_links:
                for link_set in new_links:
                    link_set.remove(common)             # As everything in common should be in all of these sets, we should get no issues
            common_links += self.find_common_links(new_links)
        return common_links
    
    # Output the solved grid in a human-readable format
    def write_grid(self, of_name="Output.csv"):
        symbols = np.repeat( None, self.grid.shape[0]*self.grid.shape[1] ).reshape( self.grid.shape )
        for row in range(len(self.grid)):
            for col in range(len(self.grid[row])):
                if isinstance( self.grid[row,col], Node ):
                    symbols[row,col] = str( self.grid[row,col].max_degree )
                elif isinstance( self.grid[row,col], Edge ):
                    symbols[row,col] = self.grid[row,col].get_symbol()
                else:
                    symbols[row,col] = ""
        df = pd.DataFrame(symbols)
        df.to_csv(of_name, header=None, index=False)
    