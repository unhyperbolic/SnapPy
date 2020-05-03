"""

"""

from .knot import *

def two_three_simplify(MC, num_blowups, num_attempts):
    """
    Attempt to simplify the Mcomplex triangulation using only
    2-3 and 3-2 moves in a very naive way.
    """
    if len(MC)<=2:
        return MC
    for i in range(num_attempts):
        print('T:{},E:{},V:{}'.format( len(MC.Tetrahedra), len(MC.Edges), len(MC.Vertices) ))
        for j in range(num_blowups):
            print('going up')
            apply_best_two_to_three_move(MC)
            check_arcs_connected(MC)
        print('going_down')
        eliminate_random_valence_three(MC)
        check_arcs_connected(MC)


def search_for_best_two_to_three(MC):
    num_valence_three = MC.EdgeValences.count(3)
    best_face = None
    best_diff = 0
    for face in MC.Faces:
        corner = face.Corners[0]
        tet = corner.Tetrahedron
        copy = MC.copy()
        copy_tet = copy.Tetrahedra[tet.Index]
        copy.two_to_three(corner.Subsimplex, copy_tet)
        copy.rebuild()
        diff = copy.EdgeValences.count(3)-num_valence_three
        if diff > best_diff:
            best_face = face
            best_diff = diff

    return best_face, best_diff
        
def apply_best_two_to_three_move(MC):
    face, diff = search_for_best_two_to_three(MC)
    if face == None:
        face = random.choice(MC.Faces)

    corner = face.Corners[0]
    MC.two_to_three(corner.Subsimplex, corner.Tetrahedron)
    MC.rebuild()


def apply_many_two_to_three(MC, n):
    for i in range(n):
        for tet in MC.Tetrahedra:
            if tet.arcs:
                f = random.choice([F0,F1,F2,F3])
                MC.two_to_three(f, tet)
#                MC.rebuild()
                break


def apply_one_three_to_two(MC):
    did_simplify = 0
    for edge in MC.Edges:
        if edge.valence() == 3:
            if MC.three_to_two(edge):
                did_simplify = 1
                MC.rebuild()
                break
            
    return did_simplify    
        




def eliminate_random_valence_three(MC):
    did_simplify = 0
    progress = 1
    while progress:
        progress = 0
        valence_three_edges = [edge for edge in MC.Edges if edge.valence()==3]
        random.shuffle(valence_three_edges)
        for edge in valence_three_edges:
            if MC.three_to_two(edge):
                progress, did_simplify = 1, 1
                break
        MC.rebuild()
    return did_simplify


def check_arcs_connected(MC):
    for tet in MC.Tetrahedra:
        pts = [arc.start for arc in tet.arcs]
        pts.extend([arc.end for arc in tet.arcs])        
        for pt in pts:
            is_on_face = False
            for x in pt.vector:
                if x == 0:
                    is_on_face = True
                    assert pts.count(pt) == 1
            if not is_on_face:
                assert pts.count(pt) == 2


def trace_knot(MC):
    return trace_knot_starting_at(*get_arc(MC))
                
def trace_knot_starting_at(arc, tet):
    """
    Follow the arcs around from one tetrahedron to the next.
    """
    
    start_point, start_tet = arc.start, tet
    knot = [(arc,tet)]
    next_arc, next_tet = continue_arc(arc,tet)

    while next_arc.start != start_point or next_tet != start_tet:
        print(next_arc, next_tet)
        knot.append((next_arc,next_tet))
        next_arc, next_tet = continue_arc(next_arc,next_tet)
        
    return knot


FacesByIndex = [F0,F1,F2,F3]

def continue_arc(arc, tet):

    endpoint = arc.end
    if arc.end.is_on_boundary_face(): #move to neighboring tetrahedron
        face_bitmap = FacesByIndex[arc.end.zero_coordinates()[0]]
        perm = tet.Gluing[face_bitmap]
        endpoint = endpoint.permute(perm)
        tet = tet.Neighbor[face_bitmap]
        other_arcs = tet.arcs[:] 
    else: 
        other_arcs = tet.arcs[:]
        other_arcs.remove(arc)
        
    for other_arc in other_arcs:
        if other_arc.start == endpoint:
            return other_arc, tet
        elif other_arc.end == endpoint:
            return other_arc.reversed(), tet
        else:
            continue
    raise Exception("No matching arc found")

def get_arc(MC):
    for tet in MC.Tetrahedra:
        for arc in tet.arcs:
            return arc, tet


def low_valence_simplify(MC, n):
    if len(MC)<=2:
        return MC
    for i in range(n):        
        lower_valence = low_valence_reducing_moves(MC)
        while lower_valence:
            tet_index, subsimplex = random.choice(lower_valence)
            MC.two_to_three(subsimplex, MC.Tetrahedra[tet_index])
            lower_valence = low_valence_reducing_moves(MC)

        MC.eliminate_valence_three()
        MC.rebuild()

    return MC
            

def combined_simplify(MC, n, lv, tt, up):
    isosigs = set([MC.isosig()])
    for i in range(n):

        low_valence_simplify(MC, lv)
        two_three_simplify(MC, tt, up)
        isosig = MC.isosig()
        if isosig in isosigs:
            print('\nREPEAT\n')
        isosigs.add(isosig)
        if len(MC)<=2:
            break
        
    

def all_two_to_three(MC):
    adjacent_MC = []
    for face in MC.Faces:
        corner = face.Corners[0]
        tet = corner.Tetrahedron
        copy = MC.copy()
        copy_tet = copy.Tetrahedra[tet.Index]
        if copy.two_to_three(corner.Subsimplex, copy_tet):
            copy.rebuild()
            adjacent_MC.append(copy)
    return adjacent_MC

def low_valence_reducing_moves(MC):
    num_low_valence = MC.EdgeValences.count(1)+MC.EdgeValences.count(2)
    num_valence_three = MC.EdgeValences.count(3)
    adjacent_MC = []
    print('num low valence:{}'.format(num_low_valence))

    locations = []
    for face in MC.Faces:
        corner = face.Corners[0]
        tet = corner.Tetrahedron
        copy = MC.copy()
        copy_tet = copy.Tetrahedra[tet.Index]

        if copy.two_to_three(corner.Subsimplex, copy_tet):
            copy.rebuild()
            new_num_low_valence = copy.EdgeValences.count(1)+copy.EdgeValences.count(2)
            if new_num_low_valence < num_low_valence:
                print('new num low valence:{}'.format(new_num_low_valence))
                adjacent_MC.append(copy)
                locations.append((tet.Index,corner.Subsimplex))
    return locations



def all_good_two_to_three(MC):
    num_valence_three = MC.EdgeValences.count(3)
    adjacent_MC = []
    for face in MC.Faces:
        corner = face.Corners[0]
        tet = corner.Tetrahedron
        copy = MC.copy()
        copy_tet = copy.Tetrahedra[tet.Index]
        if copy.two_to_three(corner.Subsimplex, copy_tet):
            copy.rebuild()
            if copy.EdgeValences.count(3)>num_valence_three:
                adjacent_MC.append(copy)
    return adjacent_MC

def simplify_three_to_two_avoiding(MC, isosigs):
    progress = True
    while progress:
        progress = False
        

def random_simplify(MC, n, max_up=4):
    for i in range(n):
        if len(MC) == 2:
            return MC
        isosigs = set([MC.isosig()])
        for i in range(max_up):
            moves_up = all_good_two_to_three(MC)
            if moves_up:
                MC = random.choice(moves_up)
            else:
                moves_up = all_two_to_three(MC)
                MC = random.choice(moves_up)
            isosigs.add(MC.isosig())
        moves_down = all_three_to_two(MC)
        random.shuffle(moves_down)
        new_move = True
        while moves_down and new_move:

            new_move = False
            for other_MC in moves_down:
                if other_MC.isosig() not in isosigs:
                    MC = other_MC
                    new_move = True
                    break
            moves_down = all_three_to_two(MC)
            random.shuffle(moves_down)
                    
        MC.rebuild()
    return MC
        

def all_three_to_two(MC):
    adjacent_MC = []
    for i, edge in enumerate(MC.Edges):
        if edge.valence() == 3:
            copy = MC.copy()
            copy_edge = copy.Edges[i]
            assert copy_edge.valence() == 3
            if copy.three_to_two(copy_edge):
                copy.rebuild()
                adjacent_MC.append(copy)            
    return adjacent_MC
    
from collections import Counter
def greedy_simplify(MC):
    seen_isosigs = set()
    search_stack = [MC.copy()]
    while search_stack:
        print(Counter([len(M) for M in search_stack]))
        MC = search_stack.pop()
        if len(MC) == 1:
            return MC
        isosig = MC.isosig()
        if isosig not in seen_isosigs:
            seen_isosigs.add(isosig)
            for next_MC in all_three_to_two(MC):
                search_stack.append(next_MC)
            for next_MC in all_two_to_three(MC):
                search_stack.insert(0,next_MC)            


def remove_duplicates(MC_list):
    isosigs = set()
    MC_list_no_duplicates = []
    for MC in MC_list:
        isosig = MC.isosig()
        if isosig not in isosigs:
            isosigs.add(isosig)
            MC_list_no_duplicates.append(MC)
    return MC_list_no_duplicates


def bounded_up_simplify(MC, max_search_height = 5):

    while len(MC)>2:
        MCs = [[MC]]
        print('Best num tets:{}'.format(len(MC)))
        for i in range(max_search_height):
            
            print('Increasing size by {}'.format(i+1))
            found_better_MC = False
            last_layer = MCs[-1]
            current_layer = []
            isosigs = set()
            for C in last_layer:
                new_MCs = all_two_to_three(C)
                for new_MC in new_MCs:
                    isosig = new_MC.isosig()
                    if isosig not in isosigs:
                        isosigs.add(isosig)
                        current_layer.append(new_MC)
                        copy = new_MC.copy()
                        copy.eliminate_valence_three()
                        if len(copy)<len(MC):
                            MC = new_MC
                            MC.eliminate_valence_three()
                            MC.rebuild()
                            found_better_MC = True
                            break
                                                
            print('{} new triangulations'.format(len(current_layer)))
            if found_better_MC:
                break
            else:    
                MCs.append(current_layer)
            
    return MC

def valence_two_edges(MC):
    return [edge for edge in MC.Edges if edge.valence() == 2]

def valence_one_edges(MC):
    return [edge for edge in MC.Edges if edge.valence() == 1]


def all_valence_one_moves(MC):
    for i, edge in enumerate(MC.Edges):
        if edge.valence() == 1:
            copy = MC.copy()
            edge_copy = MC.Edges[i]
            corner = edge_copy.Corners[0]
            for j, face in enumerate(FacesByIndex):
                print('F{}'.format(j))
                
            op_face = RightFace[comp(corner.Subsimplex)]
            tet = corner.Tetrahedron
            if MC.two_to_three(op_face, tet):
                progress = True
                MC.rebuild()
                break
            else:
                print('found valence 1, could not divide')
                break
    

def subdivide_valence_one(MC):
    print('num valence one before: {}'.format(MC.EdgeValences.count(1)))
    print('num valence two before: {}'.format(MC.EdgeValences.count(2)))
    print('num valence three before: {}'.format(MC.EdgeValences.count(3)))
    for edge in MC.Edges:
        if edge.valence() == 1:
            corner = edge.Corners[0]
            op_face = RightFace[comp(corner.Subsimplex)]
            tet = corner.Tetrahedron
            if MC.two_to_three(op_face, tet):
                progress = True
                MC.rebuild()
                break
            else:
                print('found valence 1, could not divide')
                break
    print('')
    print('num valence one after: {}'.format(MC.EdgeValences.count(1)))
    print('num valence two after: {}'.format(MC.EdgeValences.count(2)))
    print('num valence three after: {}'.format(MC.EdgeValences.count(3)))
    print('\n')

def subdivide_valence_two(MC):
    print('num valence one before: {}'.format(MC.EdgeValences.count(1)))
    print('num valence two before: {}'.format(MC.EdgeValences.count(2)))
    print('num valence three before: {}'.format(MC.EdgeValences.count(3)))    
    for edge in MC.Edges:
        if edge.valence() == 2:
            corner = edge.Corners[0]
            op_face = RightFace[corner.Subsimplex]
            tet = corner.Tetrahedron
            if MC.two_to_three(op_face, tet):
                progress = True
                MC.rebuild()
                break
            else:
                print('found valence 2, could not divide')
                break
    print('')
    print('num valence one after: {}'.format(MC.EdgeValences.count(1)))
    print('num valence two after: {}'.format(MC.EdgeValences.count(2)))
    print('num valence three after: {}'.format(MC.EdgeValences.count(3)))    
    print('\n')
    
def subdivide_low_valence(MC):
    progress = True    
    while progress:
        progress = False
        for edge in MC.Edges:
            if edge.valence() == 2:
                corner = edge.Corners[0]
                op_face = RightFace[corner.Subsimplex]
                tet = corner.Tetrahedron
                if MC.two_to_three(op_face, tet):
                    progress = True
                    MC.rebuild()
                    continue
                else:
                    print('found valence 2, could not divide')
            elif edge.valence() == 1:
                corner = edge.Corners[0]
                op_face = LeftFace[comp(corner.Subsimplex)]
                tet = corner.Tetrahedron
                if MC.two_to_three(op_face, tet):
                    progress = True
                    MC.rebuild()
                    continue
                else:
                    print('found valence 1, could not divide')
                    
        MC.eliminate_valence_three()
        MC.rebuild()
        print('num valence one: {}'.format(MC.EdgeValences.count(1)))
        print('num valence two: {}'.format(MC.EdgeValences.count(2)))




import heapq

class PriorityQueue:
    def __init__(self):
        self.elements = []
    
    def empty(self):
        return len(self.elements) == 0
    
    def put(self, item, priority):
        heapq.heappush(self.elements, (priority, item))
    
    def get(self):
        return heapq.heappop(self.elements)[1]

    
def heuristic(a, b):
    (x1, y1) = a
    (x2, y2) = b
    return abs(x1 - x2) + abs(y1 - y2)

def a_star_search(manifold, goal_manifold):
    goal = goal_manifold.isosig()
    goal_length = len(goal_manifold)
    start = manifold.isosig()

    frontier = PriorityQueue()
    frontier.put(start, 0)
    came_from = {}
    cost_so_far = {}
    came_from[start] = None
    cost_so_far[start] = 0
    
    while not frontier.empty():
        print(len(came_from))
        current = frontier.get()

        if current == goal:
            break
        current_manifold = Mcomplex(current)        
        for next_manifold in neighbors(current_manifold):
            new_cost = cost_so_far[current] + 1
            next = next_manifold.isosig()
            if next not in cost_so_far or new_cost < cost_so_far[next]:
                cost_so_far[next] = new_cost
                priority = new_cost + abs(len(next_manifold)-goal_length)
                frontier.put(next, priority)
                came_from[next] = current
    
    return came_from, cost_so_far



def neighbors(manifold):
    return all_three_to_two(manifold)+all_two_to_three(manifold)
