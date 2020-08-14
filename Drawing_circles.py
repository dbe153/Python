# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:40:25 2020

Draw Circles 
    Acer 1920 x 1080 - 46 dpi
"""
import turtle

turtle.pen(pen=None,speed=10)

""" Draw Tubing 5.5" OD  =(OD*DPI-pensize/2) pensize=(wall thickness*dpi)"""
turtle.pen(fillcolor="black", pensize=28)
turtle.right(90)
turtle.penup()
turtle.forward(239)
turtle.left(90)
turtle.pendown()
turtle.circle(239)

"""Draw Tubing ID
turtle.pen(fillcolor="black", pensize=1)
turtle.left(90)
turtle.penup()
turtle.forward(28)
turtle.right(90)
turtle.pendown()
turtle.circle(225)"""

"""Draw Pump Housing OD 4.628" 211-pensize/2 pensize=wall thickness"""
turtle.pen(pencolor="red",fillcolor="red", pensize=15.5)
turtle.left(90)
turtle.penup()
turtle.forward(35)
turtle.right(90)
turtle.pendown()
turtle.circle(203.25)

"""End Turtle"""
turtle.done()
