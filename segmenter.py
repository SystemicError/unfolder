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
    edge_length = vector.magnitude()
    rel_point = point.minus(endpoint0)
    # project rel_point onto the edge's unit vector
    t = vector.dot_product(rel_point)/edge_length
    if abs(t - rel_point.magnitude()) > tolerance:
        return False
    if t > tolerance*edge_length and t < edge_length*(1.0 - tolerance):
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
    for splitter in triangles:
        for i, split in enumerate(triangles):
            for corner in splitter.vertices:
                for edge in range(len(split.vertices)):
                    if point_splits_edge(corner, split.vertices[(edge + 1)%3], split.vertices[(edge + 2)%3]):
                        print("Split!")
                        children = split_triangle(split, edge, corner)
                        triangles.append(children[0])
                        triangles.append(children[1])
                        triangles.pop(i)
                        return split_edge_points(triangles)
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
