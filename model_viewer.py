import numpy
import sys
from stl import mesh
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *

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

def get_triangle_color(normal):
    magnitude = (normal[0]**2 + normal[1]**2 + normal[2]**2)**.5
    rgb = [normal[j]/magnitude/4.0 + .75 for j in range(3)]
    for j in range(3):
        if rgb[j] < .5:
            rgb[j] =.5 
    return rgb

def draw_projections(triangles, normals, offsets, scale = 1.0):
    glBegin(GL_TRIANGLES)
    for i in range(len(triangles)):
        tri = triangles[i]
        # decide color
        glColor3fv(get_triangle_color(normals[i]))
        # draw vertices
        for j in range(0, 9, 3):
            glVertex3fv([tri[j]*scale + offsets[0], tri[j + 1]*scale + offsets[1], tri[j + 2]*scale + offsets[2]])

    glEnd()

    return

def on_mouse_move():
    dx, dy = pygame.mouse.get_rel()
    buttons = pygame.mouse.get_pressed()
    if buttons[0] == 1:
        glRotate(dx, 0, 0, 1)
        glRotate(dy, 1, 0, 0)
    return

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Argument expected.")
        exit(1)
    triangles, normals = get_triangles_from_file(argv[1])
    scale = 1.0
    offsets = find_offsets(triangles, scale)

    pygame.init()
    display = (1024, 768)
    pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
    glEnable(GL_DEPTH_TEST)

    gluPerspective(35, display[0]/display[1], .01, 500.0)
    glTranslate(0,0,-250)


    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()
            elif event.type == pygame.MOUSEMOTION:
                on_mouse_move()

        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        draw_projections(triangles, normals, offsets, scale)
        pygame.display.flip()
        pygame.time.wait(1)
    return

if __name__ == "__main__":
    main(sys.argv)
