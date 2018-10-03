# unfolder.py
# Started 2018-10-01 by Tevis Tsai
# Adapted from parabolic_reflector.py
# for the general problem of unfolding a 3d (.stl) model

import math
import cairo
import sys
import itertools
from stl import mesh
from copy import copy

PX_WIDTH, PX_HEIGHT = 2550, 3300
MM_WIDTH, MM_HEIGHT = 215.9, 279.4

class Vertex:
    def __init__(self, xi, yi, zi):
        self.x = float(xi)
        self.y = float(yi)
        self.z = float(zi)
    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
    def dot_product(self, other):
        return self.x*other.x + self.y*other.y + self.z*other.z
    def distance_from(self, other):
        return ((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)**.5
    def minus(self, other):
        return Vertex(self.x - other.x, self.y - other.y, self.z - other.z)
    def equals(self, other):
        return (self.x == other.x and self.y == other.y and self.z == other.z)
    def magnitude(self):
        return (self.x**2 + self.y**2 + self.z**2)**.5

class Triangle:
    def __init__(self, v0, v1, v2):
        self.vertices = [v0, v1, v2]
        self.neighbors = [None, None, None]
        self.visited = False
        self.id = ""
    def __str__(self):
        vertex_string = "vertices[" + str(self.vertices[0]) + ", " + str(self.vertices[1]) + ", " + str(self.vertices[2]) + "]"
        neighbor_string = "neighbors["
        if self.neighbors[0] == None:
            neighbor_string = neighbor_string + "None\n"
        else:
            neighbor_string = neighbor_string + "[" + str(self.neighbors[0].vertices[0]) + ", " + str(self.neighbors[0].vertices[1]) + ", " + str(self.neighbors[0].vertices[2]) + "],\n"
        if self.neighbors[1] == None:
            neighbor_string = neighbor_string + "None\n"
        else:
            neighbor_string = neighbor_string + "[" + str(self.neighbors[1].vertices[0]) + ", " + str(self.neighbors[1].vertices[1]) + ", " + str(self.neighbors[1].vertices[2]) + "]\n"
        if self.neighbors[2] == None:
            neighbor_string = neighbor_string + "None\n"
        else:
            neighbor_string = neighbor_string + "[" + str(self.neighbors[2].vertices[0]) + ", " + str(self.neighbors[2].vertices[1]) + ", " + str(self.neighbors[2].vertices[2]) + "]]"
        return "Triangle:\n" + vertex_string + "\n" + neighbor_string
    def get_normal(self):
        "on the assumption that v0 through v2 are counterclockwise, return normal vector"
        a = self.vertices[1].minus(self.vertices[0])
        b = self.vertices[2].minus(self.vertices[0])
        normal = Vertex(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)
        return normal
    def translate_to_origin(self, n):
        x = self.vertices[n].x
        y = self.vertices[n].y
        z = self.vertices[n].z
        for i in range(3):
            self.vertices[i].x -= x
            self.vertices[i].y -= y
            self.vertices[i].z -= z
    def get_neighbor_index(self, t):
        "Returns -1 if not neighbors, which neigbhor of t1 t2 is otherwise."
        for i in range(3):
            if self.neighbors[i] != None and self.neighbors[i].vertices == t.vertices:
                return i
        return -1
    def get_sidelengths(self):
        sidelengths = []
        for i in range(3):
            sidelengths.append(self.vertices[i%3].distance_from(self.vertices[(i + 1)%3]))
        return sidelengths
    def center(self):
        x = (self.vertices[0].x + self.vertices[1].x + self.vertices[2].x)/3.0
        y = (self.vertices[0].y + self.vertices[1].y + self.vertices[2].y)/3.0
        z = (self.vertices[0].z + self.vertices[1].z + self.vertices[2].z)/3.0
        return Vertex(x, y, z)
    def apply_2d_transform(self, a, b, c, d):
        for v in self.vertices:
            new_x = a*v.x + b*v.y
            new_y = c*v.x + d*v.y
            v.x = new_x
            v.y = new_y
        return
    def align_to_neighbor(self, n):
        "Rotate a flattened triangle so that it aligns with its nth neighbor."
        initial_sidelengths = self.get_sidelengths()
        # translate self so A is at origin
        self.translate_to_origin((n + 1)%3)

        # set up some helpful shorthand
        # We go counterclockwise within self, counterclockwise within neighbor
        vertexA = self.vertices[(n + 1)%3]
        vertexB = self.vertices[(n + 2)%3]
        vertexC = self.vertices[n]
        neighbor = self.neighbors[n]
        self_index = neighbor.get_neighbor_index(self)
        neighborA = neighbor.vertices[(self_index + 2)%3]
        neighborB = neighbor.vertices[(self_index + 1)%3]
        relative_neighborB = neighborB.minus(neighborA)

        # rotate to parallel neighbor
        denom = (vertexB.dot_product(vertexB)*relative_neighborB.dot_product(relative_neighborB))**.5
        # dot product for cosine
        cos_theta = vertexB.dot_product(relative_neighborB)/denom
        # cross product for sine
        sin_theta = (vertexB.x*relative_neighborB.y - vertexB.y*relative_neighborB.x)/denom

        self.apply_2d_transform(cos_theta, -sin_theta, sin_theta, cos_theta)

        # translate it to neighbor
        vertexA.x = neighborA.x
        vertexA.y = neighborA.y
        vertexB.x += vertexA.x
        vertexB.y += vertexA.y
        vertexC.x += vertexA.x
        vertexC.y += vertexA.y

        #if self.get_sidelengths() != initial_sidelengths:
        #    print "Sidelengths changed!"
        return

    def flatten(self):
        vertices = self.vertices
        side_0_1 = vertices[0].distance_from(vertices[1])
        side_0_2 = vertices[0].distance_from(vertices[2])
        self.translate_to_origin(0)
        dp = vertices[1].dot_product(vertices[2])
        vertices[1].x = side_0_1
        vertices[1].y = 0.0
        vertices[1].z = 0.0
        vertices[2].x = dp/side_0_1
        vertices[2].y = (side_0_2**2 - vertices[2].x**2)**.5
        vertices[2].z = 0.0
        return

def midpoint(v0, v1, coords):
    "Gives the 2d midpoint of two vertices based on current coordinate system."
    if coords == "rectangular":
        mid = Vertex((v0.x + v1.x)/2.0, (v0.y + v1.y)/2.0, 0)
        return mid
    else:
        # midpointing for spherical
        # make two unit vectors, add them
        u0 = (math.cos(v0.x)*math.cos(v0.y), math.sin(v0.x)*math.cos(v0.y), math.sin(v0.y))
        u1 = (math.cos(v1.x)*math.cos(v1.y), math.sin(v1.x)*math.cos(v1.y), math.sin(v1.y))
        u_sum = (u0[0] + u1[0], u0[1] + u1[1], u0[2] + u1[2])
        # get the azimuth and altitude of the result
        theta = math.atan2(u_sum[1], u_sum[0])
        phi = math.atan2(u_sum[2], (u_sum[0]**2 + u_sum[1]**2)**.5)
        mid = Vertex(theta, phi, 1.0)
        return mid

def rotate_triangles_to_plane(triangles):
    # starting with triangles[0] as seed, use breadth-first web traversal
    # and flatten into plane

    # place the seed triangle on the plane
    for t in triangles:
        t.flatten()

    # recur on its unvisited neighbors
    triangles[0].visited = True
    frontier = [triangles[0]]
    while len(frontier) > 0:
        new_frontier = []
        for triangle in frontier:
            for neighbor in triangle.neighbors:
                if neighbor != None and neighbor.visited == False:
                    self_index = neighbor.get_neighbor_index(triangle)
                    neighbor.align_to_neighbor(self_index)
                    neighbor.visited = True
                    new_frontier.append(neighbor)
        frontier = new_frontier


    return

def draw_triangles(cr, center, triangles):
        cr.set_source_rgb(0.0, 0.0, 0.0)

        for triangle in triangles:
            draw_triangle(cr, center, triangle)
        return

def draw_triangle(cr, center, triangle):
    vertices = triangle.vertices
    x0 = center[0] + vertices[0].x
    y0 = center[1] + vertices[0].y
    x1 = center[0] + vertices[1].x
    y1 = center[1] + vertices[1].y
    x2 = center[0] + vertices[2].x
    y2 = center[1] + vertices[2].y
    cr.move_to(x0, y0)
    cr.line_to(x1, y1)
    cr.line_to(x2, y2)
    cr.line_to(x0, y0)
    cr.stroke()
    return

def draw_template(cr, triangles):
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.rectangle(0, 0, PX_WIDTH, PX_HEIGHT)
    cr.paint()
    cr.stroke()

    xc, yc = find_template_center(triangles)
    center = (PX_WIDTH*.5 - xc, PX_HEIGHT*.5 - yc)

    draw_triangles(cr, center, triangles)
    return

def draw_neighborhoods(cr, triangles):

    xc, yc = find_template_center(triangles)
    center = (PX_WIDTH*.5 - xc, PX_HEIGHT*.5 - yc)

    cr.set_source_rgb(0.0, 1.0, 0.0)
    for t in triangles:
        tc = t.center()
        tx, ty, tz = tc.x, tc.y, tc.z
        for n in [q for q in t.neighbors if q != None]:
            cr.move_to(tx + center[0], ty + center[1])
            nc = n.center()
            nx, ny, nz = nc.x, nc.y, nc.z
            cr.line_to((nx + tx)/2.0 + center[0], (ny + ty)/2.0 + center[1])
            cr.stroke()
    return

def find_template_center(triangles):
    x_vals = []
    y_vals = []
    for t in triangles:
        for vertex in t.vertices:
            x_vals.append(vertex.x)
            y_vals.append(vertex.y)
    xmin = min(x_vals)
    xmax = max(x_vals)
    ymin = min(y_vals)
    ymax = max(y_vals)
    return (xmax + xmin)/2.0, (ymax + ymin)/2.0

def link_all_triangles(triangles):
    for triangle in triangles:
        for i in range(3):
            v0 = triangle.vertices[(i + 1)%3]
            v1 = triangle.vertices[(i + 2)%3]
            for candidate in triangles:
                if candidate != triangle:
                    v0_match = False
                    v1_match = False
                    for v in candidate.vertices:
                        if v0.equals(v):
                            v0_match = True
                        if v1.equals(v):
                            v1_match = True
                    if v0_match and v1_match:
                        triangle.neighbors[i] = candidate
    return

def link_and_split_triangles(triangles):
    "Returns several disconnected polygons/tilings that are guaranteed not to self-overlap."
    # first, link everybody.  We'll split neighbors that create saddle points later.
    link_all_triangles(triangles)
    polygons = []
    while len(triangles) > 0:
        print("Starting new polygon " + str(len(polygons)) + ".")
        print(str(len(triangles)) + " triangles remain.")
        polygon = []
        # use the first one in line as a seed
        triangle = triangles.pop()
        if triangle == None:
            break
        else:
            polygon.append(triangle)
        finished = False
        while not finished:
            print("Searching for frontier . . .")
            print(str(len(triangles)) + " triangles remain.")
            finished = True
            for triangle in polygon:
                print("Searching polygon . . .")
                for j in range(len(triangle.neighbors)):
                    print("Searching neighbor . . .")
                    print(str(len(triangles)) + " triangles remain.")
                    if triangle.neighbors[j] != None and triangle.neighbors[j] in triangles:
                        print("Neighbor not null or already in polygon . . .")
                        if not forms_saddle_point(triangle, j, polygon + [triangle.neighbors[j]]):
                            print("Found new triangle for polygon.")
                            polygon.append(triangle.neighbors[j])
                            triangles.pop(triangles.index(triangle.neighbors[j]))
                            finished = False
        polygons.append(polygon)
        print("Appending new polygon (" + str(len(polygons)) + " total).")
    # unlink everyone, then link them only to their fellows in their polygon
    for polygon in polygons:
        for triangle in polygon:
            triangle.neighbors = [None, None, None]
        link_all_triangles(polygon)
    return polygons

def forms_saddle_point(triangle, j, polygon):
    "Returns true if triangle (element of polygon) would gain saddle point by adding its jth neighbor"
    if not (triangle in polygon):
        print("Undefined precondition for forms_saddle_point().")
        exit(1)
    # counterclockwise case
    angle_sum = 0.0
    current_triangle = triangle.neighbors[j]
    neighbor_index = current_triangle.neighbors.index(triangle) # we have to start with the proposed neighbor and go backward
    finished = False
    while not finished:
        # establish some shorthand
        pivot_vertex = current_triangle.vertices[(neighbor_index + 2)%3]
        forward_leg = current_triangle.vertices[(neighbor_index + 1)%3].minus(pivot_vertex)
        rear_leg = current_triangle.vertices[neighbor_index].minus(pivot_vertex)
        # add the angle of our current triangle
        angle_sum += math.acos(forward_leg.dot_product(rear_leg)/forward_leg.magnitude()/rear_leg.magnitude())
        # move current triangle to the next neighbor around this vertex
        next_triangle = current_triangle.neighbors[neighbor_index]
        if next_triangle == None or not (next_triangle in polygon) or next_triangle == triangle.neighbors[j]:
            finished = True
        else:
            next_neighbor_index = (next_triangle.neighbors.index(current_triangle) + 2)%3
            current_triangle = next_triangle
            neighbor_index = next_neighbor_index
    if angle_sum > 2*math.pi:
        return True

    #clockwise case
    angle_sum = 0.0
    current_triangle = triangle.neighbors[j]
    neighbor_index = current_triangle.neighbors.index(triangle) # we have to start with the proposed neighbor and go backward
    finished = False
    while not finished:
        # establish some shorthand
        pivot_vertex = current_triangle.vertices[(neighbor_index + 1)%3]
        forward_leg = current_triangle.vertices[(neighbor_index + 2)%3].minus(pivot_vertex)
        rear_leg = current_triangle.vertices[neighbor_index].minus(pivot_vertex)
        # add the angle of our current triangle
        angle_sum += math.acos(forward_leg.dot_product(rear_leg)/forward_leg.magnitude()/rear_leg.magnitude())
        # move current triangle to the next neighbor around this vertex
        next_triangle = current_triangle.neighbors[neighbor_index]
        if next_triangle == None or not (next_triangle in polygon) or next_triangle == triangle.neighbors[j]:
            finished = True
        else:
            next_neighbor_index = (next_triangle.neighbors.index(current_triangle) + 1)%3
            current_triangle = next_triangle
            neighbor_index = next_neighbor_index
    if angle_sum > 2*math.pi:
        return True

    return False

def get_triangles_from_file(path, scale):
    my_mesh = mesh.Mesh.from_file(path)
    triangles = []
    for i in range(len(my_mesh.points)):
        point0 = my_mesh.points[i][0:3]
        point1 = my_mesh.points[i][3:6]
        point2 = my_mesh.points[i][6:9]
        v0 = Vertex(point0[0]*scale, point0[1]*scale, point0[2]*scale)
        v1 = Vertex(point1[0]*scale, point1[1]*scale, point1[2]*scale)
        v2 = Vertex(point2[0]*scale, point2[1]*scale, point2[2]*scale)
        triangle = Triangle(v0, v1, v2)
        # check that this agrees with the STL-specified normal vector
        computed_normal = triangle.get_normal()
        normal = my_mesh.normals[i]
        if normal[0]*computed_normal.x + normal[1]*computed_normal.y + normal[2]*computed_normal.z < 0:
            # if not, change the order so it does
            v0 = Vertex[point0[0], point0[2], point0[1]]
            print("Reversing triangle " + str(i) + " orientation.")
        triangles.append(triangle)
    return triangles

def main(argv):
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, PX_WIDTH, PX_HEIGHT)
    cr = cairo.Context(surface)

    if len(argv) < 2:
        print("Need .stl argument.")
        exit(1)

    scale = 1.0
    if len(argv) >= 3:
        scale = float(argv[2])
    triangles = get_triangles_from_file(argv[1], scale)

    polygons = link_and_split_triangles(triangles)

    print(polygons[0][0].neighbors)
    for triangles in polygons[0:1]:
        path = "preelevation.png"
        draw_template(cr, triangles)
        draw_neighborhoods(cr, triangles)
        surface.write_to_png(path)
    
        rotate_triangles_to_plane(triangles)
        draw_template(cr, triangles)


        path = "unfolded.png"
        surface.write_to_png(path)

if __name__ == "__main__":
        main(sys.argv)