# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 01:16:52 2018

@author: deniz
"""
from numpy import arccos, pi
class point:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x  =   x
        self.y  =   y
        self.z  =   z
    def __repr__(self):
        x   = "{}i".format(self.x) if self.x != 0 else ""
        y   = "{}{}j".format(("-",("+","")[self.x == 0])[self.y > 0], abs(self.y)) if self.y != 0 else ""
        z   = "{}{}k".format(("-",("+","")[self.x == 0 and self.y == 0])[self.z > 0], abs(self.z)) if self.z != 0 else ""
        pos = x+y+z if x+y+z != "" else "origin"
        return "[point ]\n\t"+pos
    def __add__(self, other):
        return point(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other):
        return point(self.x - other.x, self.y - other.y, self.z - other.z)
    def __mul__(self,other):
        return point(self.x * other, self.y * other, self.z * other)
    def __truediv__(self,other):
        return point(self.x / other, self.y / other, self.z / other)
class vector:
    def __init__(self, tail = point(), head = point()):
        self.tail = tail
        self.head = head
        self.vec  = head - tail
        self.mag  = (self.vec.x**2 + self.vec.y**2 + self.vec.z**2)**.5
    def __repr__(self):
        x   = "{}i".format(self.vec.x) if self.vec.x != 0 else ""
        y   = "{}{}j".format(("-",("+","")[self.vec.x == 0])[self.vec.y > 0], abs(self.vec.y)) if self.vec.y != 0 else ""
        z   = "{}{}k".format(("-",("+","")[self.vec.x == 0 and self.vec.y == 0])[self.vec.z > 0], abs(self.vec.z)) if self.vec.z != 0 else ""
        pos = x+y+z if x+y+z != "" else "origin"
        return "[vector]\n\t"+pos
    def __add__(self, other):
        return vector(self.tail, self.head + other.vec)
    def __sub__(self, other):
        return vector(self.tail, self.head - other.vec)
    def __mul__(self, other):
        "dot product"
        if type(other) != vector: return vector(head = self.vec * other)
        else: return self.vec.x * other.vec.x + self.vec.y * other.vec.y + self.vec.z * other.vec.z
    def __truediv__(self, other):
        return vector(head = self.vec / other)
    def __matmul__(self, other):
        x = self.vec.y * other.vec.z - self.vec.z * other.vec.y
        y = self.vec.z * other.vec.x - self.vec.x * other.vec.z
        z = self.vec.x * other.vec.y - self.vec.y * other.vec.x
        return vector(head = point(x = x, y = y, z = z))
    def unit(self):
        return vector(head = self.vec / self.mag)
    def phi(self,other):
        return arccos(self*other/(self.mag * other.mag))*180/pi
"""
a = point(-2,3,1)
b = point(0,4,0)
c = point(-1,3,3)

A = vector(head = a)
B = vector(head = b)
C = vector(head = c)
"""
