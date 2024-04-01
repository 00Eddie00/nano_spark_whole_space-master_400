# -*- coding: UTF-8 -*-
import turtle as t
import numpy as np

coordinates = np.loadtxt("open_grid_coordinates.csv", delimiter=',')
x = (coordinates[:, 1] - 2500) / 4
y = (coordinates[:, 0] - 1250) / 4

line = np.loadtxt("open_neighbor.csv", int, delimiter=',')
line1 = line[:, 0]
line2 = line[:, 1]
line3 = line[:, 2]
line4 = line[:, 3]

t.screensize(1000, 1000)

t.tracer(False)
for i in range(len(line)):
    if line1[i] != -1:
        t.penup()
        t.goto(x[i], y[i])
        t.pendown()
        t.goto(x[line1[i]], y[line1[i]])
    if line2[i] != -1:
        t.penup()
        t.goto(x[i], y[i])
        t.pendown()
        t.goto(x[line2[i]], y[line2[i]])
    if line3[i] != -1:
        t.penup()
        t.goto(x[i], y[i])
        t.pendown()
        t.goto(x[line3[i]], y[line3[i]])
    if line4[i] != -1:
        t.penup()
        t.goto(x[i], y[i])
        t.pendown()
        t.goto(x[line4[i]], y[line4[i]])
t.hideturtle()
ts = t.getscreen()
ts.getcanvas().postscript(file="open_grid.eps")
t.done()
