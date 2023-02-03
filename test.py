from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import mpl_toolkits.mplot3d.art3d as art3d

def distance(p1, p2):
    return sqrt(((p2[0] - p1[0]) ** 2) + ((p2[1] - p1[1]) ** 2))

def triangle(p1, p2, p3):
    points = (p1, p2, p3)
    return Polygon(points, closed=True)

def eq_triangle(center, size):
    p1 = (center[0], center[1] + sqrt(3)/3 * size)
    p2 = (center[0] - size/2, center[1] - sqrt(3)/6 * size)
    p3 = (center[0] + size/2, center[1] - sqrt(3)/6 * size)
    if (distance(p1, p2) != distance (p1, p3) or distance(p2, p3) != distance(p1, p2)):
        print("Uh oh")
    return triangle(p1, p2, p3)
    

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

p = eq_triangle((5, 5), 10)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_zlim(0, 10)

plt.show()
