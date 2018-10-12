import numpy
from stl import mesh
import cairo
import sys

def get_triangles_from_file(path):
    my_mesh = mesh.Mesh.from_file(path)
    return my_mesh.points, my_mesh.normals

def find_offsets(triangles, scale):
    count = len(triangles)*3
    xs = [x[0] for x in triangles] + [x[3] for x in triangles] + [x[6] for x in triangles]
    ys = [x[1] for x in triangles] + [x[4] for x in triangles] + [x[7] for x in triangles]
    zs = [x[2] for x in triangles] + [x[5] for x in triangles] + [x[8] for x in triangles]
    x_offset = -(max(xs) + min(xs))/2.0*scale
    y_offset = -(max(ys) + min(ys))/2.0*scale
    z_offset = -(max(zs) + min(zs))/2.0*scale
    return [x_offset, y_offset, z_offset]

def draw_projections(path, triangles, normals, WIDTH = 1024, HEIGHT = 768, scale = 1.0):
    offsets = find_offsets(triangles, scale)
    draw_projection(path, triangles, normals, WIDTH, HEIGHT, scale, offsets, 0, 1, 2)
    draw_projection(path, triangles, normals, WIDTH, HEIGHT, scale, offsets, 0, 2, 1)
    draw_projection(path, triangles, normals, WIDTH, HEIGHT, scale, offsets, 1, 2, 0)
    return

def draw_projection(path, triangles, normals, WIDTH, HEIGHT, scale, offsets, x_dim, y_dim, zdim):
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    cr = cairo.Context(surface)
    for i in range(len(triangles)):
        tri = triangles[i]
        cr.move_to(tri[6 + x_dim]*scale + offsets[x_dim] + WIDTH/2.0, tri[6 + y_dim]*scale + offsets[y_dim] + HEIGHT/2.0)
        for vertex_index in range(0, 9, 3):
            cr.line_to(tri[x_dim + vertex_index]*scale + offsets[x_dim] + WIDTH/2.0, tri[y_dim + vertex_index]*scale + offsets[y_dim] + HEIGHT/2.0)

        shade = normals[i][zdim]/(normals[i][0]**2 + normals[i][1]**2 + normals[i][2]**2)
        alpha = 0.7
        if shade < 0.0:
            shade = 0.0
            alpha = 0.0
        cr.set_source_rgba(255*shade, 205*shade, 205*shade, alpha)
        cr.fill()

    surface.write_to_png(path + "_" + str(x_dim) + str(y_dim) + ".png")
    return


triangles, normals = get_triangles_from_file(sys.argv[1] + ".stl")
draw_projections(sys.argv[1], triangles, normals, 1024, 768, 6.0)
