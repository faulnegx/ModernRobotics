# Python 2.7 script
from pprint import pprint
import csv

class Node(object):
    def __init__(self, num, coord, h):
        # Node knows its badge number
        self.num = num
        # Coord in (x,y)
        self.parent = None
        self.coord = coord
        self.heuristic = h
        self.__past_cost = float("inf")
        self.__est_total_cost = float("inf") + self.heuristic
    
    def __repr__(self):
        return ("Node #:" + str(self.num)
             + ", Parent: " + str(self.parent)
             + ", Coords: "+ str(self.coord)
             + ", Heuristic: " + str(self.heuristic)
             + ", Past_cost: " + str(self.__past_cost)
             + ", Est_total_cost: " + str(self.__est_total_cost)
        )

    @property
    def past_cost(self):
        return self.__past_cost

    # When past_cost is updated, est_total_cost also has to be updated.
    @past_cost.setter
    def past_cost(self,value):
        self.__past_cost = value
        self.__est_total_cost = value + self.heuristic

    @property
    def est_total_cost(self):
        return self.__est_total_cost

class Edge(object):
    def __init__(self,node1,node2,cost):
        # Nodes are represented in their number
        self.nodes = frozenset([node1,node2])
        self.cost = cost
    
    def __repr__(self):
        return ("Nodes: " + str(self.nodes)
             + ", cost: " + str(self.cost)
        )

def Index_of_node(l, desired_num):
    """
    Given a list of nodes, return index of node matching desired_num
    If not found, return None
    """
    for i in range(len(l)):
        if l[i].num == desired_num:
            return i
    return None

def Edge_cost(edges_list, node1, node2):
    # Look up edge cost in edge_list, if it is not in edge_list, return -1.0
    for edge in edges_list:
        if frozenset([node1, node2]) == edge.nodes:
            return edge.cost
    return -1.0

def A_star(nodes_list, edges_list, start_node, goal_node):
    assert (start_node != goal_node), "Start_node == goal_node"
    assert (Index_of_node(nodes_list, start_node) != None), "Start_node not found in nodes_list"
    assert (Index_of_node(nodes_list, goal_node) != None), "Goal_node not found in nodes_list"

    matrix_cost = {}
    nodes_count = len(nodes_list)
    # Fill out matrix_cost, where the key is the pair of nodes, and the value is provided
    # by edge list or is -1 when not given by edge list.
    for i in range(1,nodes_count+1):
        for j in range(1,nodes_count+1):
            if i == j:
                pass
            elif i > j:
                pass
            else:
                matrix_cost[frozenset([i,j])] = Edge_cost(edges_list, i, j)

    open_nodes_sorted_list = []
    closed_nodes_list = []
    for i in range(nodes_count):
        print "Iteration " + str(i) + ":"
        # First iteration set the start node's past_cost to 0 and put the start_node 
        # into closed_nodes_list 
        if i == 0:
            current_node = nodes_list.pop(Index_of_node(nodes_list, start_node))
            current_node.past_cost = 0.0
            print "Current_nodes:" 
            print current_node
            closed_nodes_list.append(current_node)
        # Any other node: pop the last element(lowest est_total_cost), and put into closed_nodes_list
        else:
            if len(open_nodes_sorted_list) == 0:
                print "No more available nodes to explore. No solution found."
            else:
                current_node = open_nodes_sorted_list.pop()
                print "Current_nodes:"
                print current_node
                closed_nodes_list.append(current_node)
        # Loop is complete when current_node is goal_node
        if current_node.num == goal_node:
            print "Found goal node!"
            break
        else:
            # Find all the neighbors from nodes_list
            neighbors_list = []
            for j in range(len(nodes_list)):
                # If cost < 0, that means it is not a neighbor
                if matrix_cost[frozenset([current_node.num,nodes_list[j].num])] >= 0.0:
                    neighbors_list.append(nodes_list[j].num)
                    print (str(current_node.num) + " to " 
                         + str(nodes_list[j].num) + " edge cost: "
                         + str(matrix_cost[frozenset([current_node.num,nodes_list[j].num])])
                    )
            # Add neighbors into open_nodes_sorted_list
            print "neighbors are: " + str(neighbors_list)
            for x in neighbors_list:
                open_nodes_sorted_list.append(nodes_list.pop(Index_of_node(nodes_list,x)))
            # Update the past cost of all nodes in open_nodes_sorted_list
            # if they are connected and the new past_cost < existing past cost
            for node in open_nodes_sorted_list:
               if (matrix_cost[frozenset([current_node.num,node.num])] >= 0.0 and
                  node.past_cost > matrix_cost[frozenset([current_node.num,node.num])]):
                    node.past_cost = current_node.past_cost + matrix_cost[frozenset([current_node.num,node.num])]
                    node.parent = current_node.num
            # After adding all the neighbors and updating their past cost, 
            # sort open_nodes_sorted_list based on est_total_cost, lowest cost at the end.
            open_nodes_sorted_list.sort(key=lambda node: node.est_total_cost, reverse=True)
            print "Open_nodes_sorted_list: "
            pprint(open_nodes_sorted_list)
            print "Closed_nodes_list: "
            pprint(closed_nodes_list)

    # Look up and add node parent 1 by 1, starting with the goal node, until start node
    node_sequence_int = [goal_node]
    while node_sequence_int[-1] != start_node:
        node_sequence_int.append(
            closed_nodes_list[
                Index_of_node(closed_nodes_list,node_sequence_int[-1])
            ].parent
        )
    node_sequence_int.reverse()
    return node_sequence_int


if __name__=="__main__":
    nodes_list = []
    edges_list = []
    start_node = 1
    goal_node = 12

    with open('../results/nodes.csv') as nodefile:
        nodescsvobj = csv.reader(nodefile)
        for line in nodescsvobj:
            if line[0].startswith('#'):
                pass
            else:
                nodes_list.append(
                    Node(int(line[0]),(float(line[1]),float(line[2])), float(line[3]))
                )

    with open('../results/edges.csv') as edgefile:
        edgescsvobj = csv.reader(edgefile)
        for line in edgescsvobj:
            if line[0].startswith('#'):
                pass
            else:
                edges_list.append(
                    Edge(int(line[0]), int(line[1]), float(line[2]))
                )

    node_sequence_int = A_star(nodes_list, edges_list, start_node, goal_node)
    print "Lowest cost path is :" + ",".join([str(item) for item in node_sequence_int])
    with open('../results/path.csv', 'wb') as pathfile:
        pathcsvwriter = csv.writer(pathfile)
        pathcsvwriter.writerow(node_sequence_int)