#!/usr/bin/env python
# geometry functions
import sys
import numpy as np

# Geometry operations are through a 4x4 operation matrix.
# Operations can be stacked in the matrix.
# The matrix then is applied to a series of points.

class LINE:
    def __init__(self):
        self.p0 = (0, 0, 0)   # line pass this  point
        self.t = (1, 0, 0)    # direction cosine
        return

    def from2p(self, p0, p1):
        self.p0 = np.array(p0)
        self.p1 = np.array(p1)
        v = np.subtract(p1, p0)
        v_norm = np.linalg.norm(v)
        self.t = v/v_norm
        return

class OPERATION:
    def __init__(self):
        self.operation = np.zeros((4, 4))
        self.operation[0][0] = 1.0
        self.operation[1][1] = 1.0
        self.operation[2][2] = 1.0
        self.operation[3][3] = 1.0
        return

    def reset(self):
        self.operation.fill(0)
        self.operation[0][0] = 1.0
        self.operation[1][1] = 1.0
        self.operation[2][2] = 1.0
        self.operation[3][3] = 1.0
        return

    def printme(self):
        print(self.operation)
        return

    def move(self, v):
        "Move by vector"
        op = np.zeros((4, 4))
        op[0][0] = 1.0
        op[1][1] = 1.0
        op[2][2] = 1.0
        op[3][3] = 1.0
        op[0][3] = v[0]
        op[1][3] = v[1]
        op[2][3] = v[2]
        self.operation = np.matmul(op, self.operation)

    def clone(self):
        clone_op = OPERATION()
        clone_op.operation = self.operation.copy()
        return clone_op

    def inverse(self):
        op = OPERATION()
        op.operation = np.linalg.inv(self.operation)
        return op

    def transpose(self):
        self.operation = self.operation.transpose()

    def roll(self, phi, axis):
        # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

        # validate the axis
        v = np.copy(axis.t)
        if np.linalg.norm(v) > 0.0000001:  # validate direction cosines
            v = v / np.linalg.norm(v)

            # translate to origin
            rotate_op = np.zeros((4, 4))
            rotate_op[0][0] = 1.0
            rotate_op[1][1] = 1.0
            rotate_op[2][2] = 1.0
            rotate_op[3][3] = 1.0
            rotate_op[0][3] = -axis.p0[0]
            rotate_op[1][3] = -axis.p0[1]
            rotate_op[2][3] = -axis.p0[2]
            self.operation = np.matmul(rotate_op, self.operation)

            # rotate
            SIN = np.sin(phi)
            COS = np.cos(phi)
            rotate_op = np.zeros((4, 4))
            rotate_op[3][3] = 1.0

            rotate_op[0][0] = v[0]*v[0]*(1-COS) + COS
            rotate_op[0][1] = v[0]*v[1]*(1-COS) - v[2]*SIN
            rotate_op[0][2] = v[0]*v[2]*(1-COS) + v[1]*SIN

            rotate_op[1][0] = v[0]*v[1]*(1-COS) + v[2]*SIN
            rotate_op[1][1] = v[1]*v[1]*(1-COS) + COS
            rotate_op[1][2] = v[1]*v[2]*(1-COS) - v[0]*SIN

            rotate_op[2][0] = v[0]*v[2]*(1-COS) - v[1]*SIN
            rotate_op[2][1] = v[1]*v[2]*(1-COS) + v[0]*SIN
            rotate_op[2][2] = v[2]*v[2]*(1-COS) + COS

            self.operation = np.matmul(rotate_op, self.operation)

            # translate back
            rotate_op = np.zeros((4, 4))
            rotate_op[0][0] = 1.0
            rotate_op[1][1] = 1.0
            rotate_op[2][2] = 1.0
            rotate_op[3][3] = 1.0
            rotate_op[0][3] = axis.p0[0]
            rotate_op[1][3] = axis.p0[1]
            rotate_op[2][3] = axis.p0[2]
            self.operation = np.matmul(rotate_op, self.operation)

        return

def geom_apply(operation, v):
    v1 = np.append(np.array(v), 1)[np.newaxis].transpose()
    v2 = np.matmul(operation.operation, v1)
    v3 = v2[:-1]
    return v3.T[0]


def avv(v1, v2):
    "Angle between v1 to v2 in [0, pi], around axis with right-hand rule."
    v1_n = np.linalg.norm(v1)
    v2_n = np.linalg.norm(v2)
    COS = np.dot(v1, v2) / v1_n / v2_n

    return np.arccos(COS)

def geom_2v_onto_2v(v1, v2, t1, t2):
    "superimpose v1 to t1, then align v2 to t2"
    op = OPERATION()

    # superimpose v1 to t1
    op.move(np.array(t1) - np.array(v1))

    # align v1-v2 to t1-t2
    line1 = LINE()
    line1.from2p(v1, v2)
    line2 = LINE()
    line2.from2p(t1, t2)
    angle = avv(line1.t, line2.t)
    axis = LINE()
    axis.p0 = np.array(t1)
    axis.t = np.cross(line1.t, line2.t)
    op.roll(angle, axis)
    return op

if __name__ == "__main__":
    p1 = (1, 1, 1)
    p2 = (9, 3, 3)

    t1 = (5, 1, 0)
    t2 = (0, 2, 0)

    op = geom_2v_onto_2v(p1, p2, t1, t2)
    op2 = op.inverse()

    new_p2 = geom_apply(op, p2)
    recovered = geom_apply(op2, new_p2)
    print(new_p2)
    print(recovered)
