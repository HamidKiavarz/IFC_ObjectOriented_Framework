import numpy

def norm2(X):
	return numpy.sqrt(numpy.sum(X ** 2))

def normalized(X):
	return X / norm2(X)

'''
Returns the intersection, a line, between the plane A and B
  - A and B are planes equations, such as A0 * x + A1 * y + A2 * z + A3 = 0
  - The line is returned as (U, V), where any point of the line is t * U + C, for all values of t
  - U is a normalized vector
  - C is the line origin, with the triangle (Ao, Bo, C) is orthogonal to the plane A and B,
    with Ao and Bo being the origin of plane A an B
  - If A and B are parallel, a numpy.linalg.LinAlgError exception is raised
'''
def get_plane_plane_intersection(A, B):
	U = normalized(numpy.cross(A[:-1], B[:-1]))
	M = numpy.array((A[:-1], B[:-1], U))
	X = numpy.array((-A[-1], -B[-1], 0.))
	return U, numpy.linalg.solve(M, X)