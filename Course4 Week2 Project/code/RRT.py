# Python 2.7 script
from pprint import pprint
import csv
from random import uniform

class Node(object):
    def __init__(self, num, coord): #, h):
        # Node knows its badge number
        self.num = num
        # Coord in (x,y)
        self.parent = None
        self.coord = coord
    
    def __repr__(self):
        return ("Node #:" + str(self.num)
             + ", Parent: " + str(self.parent)
             + ", Coords: "+ str(self.coord)
        )

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


def collisions_checker(point1, point2, circle_obstacles):
    """
    point1 = (x1,y1)
    point2 = (x2,y2)
    circle_obstacle = [((x0,y0), diameter0)), ...etc ]

    point1 and point 2 can form a line in the form ax+by=c
    a = y1-y2
    b = x2-x1
    c = (y2-y1)*x1 + (x1-x2)*y1
    Shortest distance = abs(a*x0 + b*y0 + c) / sqrt(a^2 + b^2)

    if the shortest distance is less than diameter/2(radius), then there is collision.

    Shortcut: When we found 1 collision, we stop checking the rest. 

    returns True if collision, returns False when no collision. 
    """
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[0]
    y2 = point2[1]
    a = y1 - y2
    b = x2 - x1
    c = (x1-x2)*y1 + (y2-y1)*x1

    for circle_obstacle in circle_obstacles:
        x0 = circle_obstacle[0][0]
        y0 = circle_obstacle[0][1]
        radius = circle_obstacle[1]*0.5
        shortest_dist = abs(a*x0 + b*y0 + c)/ ((a**2 + b**2) ** 0.5)
        if shortest_dist < radius:
            # Collision for this circle_obstacle, done and get out
            return True
    # No collision after checking all the distances
    return False

def nearest_node(nodes_l, sample_node_coord):
    """
    Check distance between sample_node_coord and each node in nodes_l
    returns nearest_node_id
    """
    shortest_dist = float('inf')
    nearest_node_id = 0
    for node in nodes_l:
        dist = distance_func(node.coord, sample_node_coord)
        if dist < shortest_dist:
            shortest_dist = dist
            nearest_node_id = node.num

    return nearest_node_id

def gen_potentially_new_node_coord(nearest_node_coord, sample_node_coord, small_distance):
    """
    if sample_node_coord is closer than small_distance away, then just return sample_node
    if the sample_node_coord is further than small_distance away,
        generate potentially_new_node_coord
        and return the potentially_new_node_coord
    """
    distance = distance_func(nearest_node_coord, sample_node_coord)
    if distance < small_distance:
        return sample_node_coord
    else:
        # dir/distance normalizes the (x_dir, y_dir) vector
        x_dir = sample_node_coord[0] - nearest_node_coord[0]
        y_dir = sample_node_coord[1] - nearest_node_coord[1]
        new_x = small_distance * x_dir/distance + nearest_node_coord[0]
        new_y = small_distance * y_dir/distance + nearest_node_coord[1]
        return (new_x, new_y)

def distance_func(node1_coord, node2_coord):
    """
    Return regular old euclidean distance
    """
    x_dir = node1_coord[0] - node2_coord[0]
    y_dir = node1_coord[1] - node2_coord[1]
    return (x_dir**2 + y_dir**2) ** 0.5

def within_goal_space(new_node, goal_node, tolerance=0.1):
    """
    Returns True when new_node is within tolerance distance away from goal_node
    """
    if distance_func(new_node.coord, goal_node.coord) < tolerance:
        return True
    else:
        return False

if __name__=="__main__":
    nodes_list = []
    edges_list = []
    obstacles_list = []
    start_node = Node(1, (-0.5, -0.5))
    goal_node = Node(-1, (0.5, 0.5))
    small_distance = 0.1

    with open('../results/obstacles.csv') as obstaclesfile:
        obstaclescsvobj = csv.reader(obstaclesfile)
        for circle in obstaclescsvobj:
            if circle[0].startswith('#'):
                pass
            else:
                obstacles_list.append(
                    ((float(circle[0]), float(circle[1])), float(circle[2]))
                )

    max_iteration = 1000
    max_nodes_count = 200
    print "Max iteration set to {}. Max nodes count set to {}".format(max_iteration, max_nodes_count)

    nodes_list.append(start_node)
    next_node_num = 1
    iteration = 0
    while len(nodes_list) < max_nodes_count and iteration < max_iteration:
        print "iteration {}".format(iteration)
        # sample some random coord
        sample_node_coord = (uniform(-0.5, 0.5), uniform(-0.5, 0.5))
        # print "sample_node_coord"
        # print sample_node_coord
        # find the closest node
        nearest_node_num = nearest_node(nodes_list, sample_node_coord)
        # print nearest_node_num
        # potentially_new_node_coord starts from nearest node and move small_distance toward sample_node_coord
        potentially_new_node_coord = gen_potentially_new_node_coord(
            nodes_list[Index_of_node(nodes_list, nearest_node_num)].coord, 
            sample_node_coord, 
            small_distance
        )
        # print "potentially_new_node"
        # print potentially_new_node_coord
        collided = collisions_checker(
            nodes_list[Index_of_node(nodes_list, nearest_node_num)].coord,
            potentially_new_node_coord, 
            obstacles_list
        )
        if collided:
            # we will not save this potentially_new_node since there is a collision
            print "collided"
        else:
            next_node_num += 1
            new_node = Node(next_node_num, sample_node_coord)
            nodes_list.append(new_node)
            new_node.parent = nearest_node_num
            edges_list.append(
                Edge(
                    nodes_list[Index_of_node(nodes_list, nearest_node_num)], 
                    new_node, 
                    distance_func(
                        nodes_list[Index_of_node(nodes_list, nearest_node_num)].coord, 
                        new_node.coord,
                    ),
                )
            )
            print "len(nodes_list) = {}".format(len(nodes_list))
            if within_goal_space(new_node, goal_node):
                print "success"
                break
        iteration += 1 


    node_sequence = [new_node.num]
    while new_node.parent != None:
        node_sequence.append(new_node.parent)
        new_node = nodes_list[Index_of_node(nodes_list, new_node.parent)]
    node_sequence.reverse()

    with open('../results/nodes.csv', 'wb') as nodefile:
        nodecsvwriter = csv.writer(nodefile)
        for node in nodes_list:
            nodecsvwriter.writerow([node.num, node.coord[0], node.coord[1]])
    with open('../results/edges.csv', 'wb') as edgefile:
        edgecsvwriter = csv.writer(edgefile)
        edgecsvwriter.writerow(["# node1 num", "node2 num", "cost"])
        for edge in edges_list:
            items_to_write = [node.num for node in edge.nodes]
            items_to_write.append(edge.cost)
            edgecsvwriter.writerow(items_to_write)
    with open('../results/path.csv', 'wb') as pathfile:
        pathcsvwriter = csv.writer(pathfile)
        pathcsvwriter.writerow(node_sequence)