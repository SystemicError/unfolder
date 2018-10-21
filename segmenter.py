# segmenter.py
# Started 2018-10-21 by Tevis Tsai
# To make sure that STL files are papercraft friendly by disallowing
# multiple neighbors for a single edge or edges longer than a certain threshold

import numpy
from stl import mesh
from unfolder import *

def point_splits_edge(point, endpoint0, endpoint1, tolerance = .001):
    "Returns true if point is collinear to endpoints and in between them."
    vector = endpoint1.minus(endpoint0)
    rel_point = point.minus(endpoint0)
    # project rel_point onto the edge's unit vector
    t = vector.dot_product(rel_point)/vector.magnitude()
    if t != rel_point.magnitude():
        return False
    if t > tolerance and t < 1.0 - tolerance:
        return True
    return False

def split_triangle(triangle, edge, point):
    "Given that point is on edge, break triangle into two smaller triangles so that point is vertex of each."
    t0 = Triangle(point, triangle.vertices[edge], triangle.vertices[(edge + 1)%3])
    t1 = Triangle(point, triangle.vertices[(edge + 2)%3], triangle.vertices[edge])
    return [t0, t1]

def split_edge_points(triangles):
    "Finds triangles whose edges contain vertices somewhere other than their endpoints and splits them at those vertices."
    # adding new triangles won't create new vertices, so we don't have to re-check old ones if we're smart about it
    i = 0
    while i < len(triangles):
        j = i + 1
        while j < len(triangles):
            m = 0
            while m < len(triangles[j].vertices):
                vertex = triangles[j].vertices[m]
                k = 0
                while k < 3:
                    if point_splits_edge(vertex, triangles[i].vertices[k], triangles[i].vertices[(k + 1)%3]):
                        print("Split!")
                        triangles = triangles + split_triangle(triangles[i], (k + 2)%3, vertex)
                        triangles.pop(i)
                        i = -1
                        j = len(triangles)
                        k = 4
                        m = 4
                    k += 1
                m += 1
            j += 1
        i += 1
    return


def save_triangles_to_file(triangles, path):
    "Saves triangles to mesh at filepath."
    data = numpy.zeros(len(triangles)*3, dtype=mesh.Mesh.dtype)
    my_mesh = mesh.Mesh(data)

    vertices = numpy.zeros((len(triangles)*3, 3))
    for i, triangle in enumerate(triangles):
        for j, vtx in enumerate(triangle.vertices):
            vertices[i*3 + j][0] = vtx.x
            vertices[i*3 + j][1] = vtx.y
            vertices[i*3 + j][2] = vtx.z

    faces = numpy.array([[3*i, 3*i + 1, 3*i + 2] for i in range(len(triangles))])

    my_mesh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            my_mesh.vectors[i][j] = vertices[f[j],:]
    my_mesh.save(path)
    return

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Need .stl argument.")
        sys.stderr.write("segmenter.py file.stl [maximum_length]")
        exit(1)
    path = argv[1]
    max_length = None
    if len(argv) >= 3:
        max_length = float(argv[2])

    triangles = get_triangles_from_file(argv[1], 1.0)

    split_edge_points(triangles)

    save_triangles_to_file(triangles, argv[1])

    return


if __name__ == "__main__":
    main(sys.argv)
