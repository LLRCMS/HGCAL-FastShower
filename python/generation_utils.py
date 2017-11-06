import math as m
from ROOT import TVector2

# FIXME: the cell side can be used as a global multiplicative factor

def vertices(type, i, a):
  if type is 'Triangles':
    if i%2: 
      return [(0., -a/m.sqrt(3.)),
              (a/2., a/(2.*m.sqrt(3.))),
              (-a/2., a/(2.*m.sqrt(3.)))
             ]
    else:
      return [(0., a/m.sqrt(3.)),
              (a/2.,-a/(2.*m.sqrt(3.))),
              (-a/2.,-a/(2.*m.sqrt(3.)))
             ]

  elif type is 'Hexagons':
      return [(a*m.sqrt(3)/2., -a/2.),
              (a*m.sqrt(3)/2., a/2.),
              (0., a),
              (-a*m.sqrt(3)/2., a/2.),
              (-a*m.sqrt(3)/2., -a/2.),
              (0., -a)
             ]
  return []

def edge_centers(type, i, a):
  if type is 'Triangles':
    if i%2: 
      return [(0., a/(2*m.sqrt(3.))),
              (a/4., -a/(4*m.sqrt(3.))),
              (-a/4., -a/(4*m.sqrt(3.)))
             ]
    else:
      return [(0., -a/(2*m.sqrt(3.))),
              (a/4., a/(4*m.sqrt(3.))),
              (-a/4., a/(4*m.sqrt(3.)))
             ]

  elif type is 'Hexagons':
      return [(a*m.sqrt(3)/2., 0.),
              (a*m.sqrt(3)/4., a*3./4.),
              (-a*m.sqrt(3)/4., a*3./4.),
              (-a*m.sqrt(3)/2., 0.),
              (-a*m.sqrt(3)/4., -a*3./4.),
              (a*m.sqrt(3)/4., -a*3./4.)
             ]
  return []

def shoot_in_cell(eta, phi, eta_min, eta_max, phi_min, phi_max, z, shoot_type, small_cell_side,\
                      large_cell_side, limit_first_zone, type, vertex_number, edge_number):

  if type not in ['Triangles', 'Hexagons']:
      raise Exception('shoot_cell_center() not implemented for geometry type '+type)
  if eta < eta_min or eta > eta_max or\
     TVector2.Phi_mpi_pi(phi-phi_min) < 0 or TVector2.Phi_mpi_pi(phi-phi_max) > 0:
      raise Exception("Error: trying to shoot particle outside geometry window")

  # Find center
  theta0 = 2.*m.atan(m.exp(-eta_max))
  theta = 2.*m.atan(m.exp(-eta))
  r0 = z*m.tan(theta0)
  r = z*m.tan(theta)
  x0 = r0*m.cos(phi_min)
  y0 = r0*m.sin(phi_min)
  x = r*m.cos(phi)
  y = r*m.sin(phi)

  if x*x+y*y <= limit_first_zone*limit_first_zone:
    cell_side = small_cell_side
  else:
    cell_side = large_cell_side

  if type is 'Triangles':
    dxdi = cell_side/2.
    dxdj = cell_side/2.
    dydj = cell_side*m.sqrt(3.)/2.
  else:
    dxdi = cell_side*m.sqrt(3.)
    dxdj = cell_side*m.sqrt(3.)/2.
    dydj = cell_side*3./2.

  j = round((y-y0)/dydj)
  i = round((x-x0 - j*dxdj)/dxdi)
  x_center = x0 + i*dxdi + j*dxdj
  y_center = y0 + j*dydj

  if type is 'Triangles' and int(i)%2:
    y_center += cell_side*m.sqrt(3)/6.

  if shoot_type == "cell_center":
    r_center = m.sqrt(x_center**2 + y_center**2)
    theta_center = m.atan(r_center/z)
    eta_shoot = -m.log(m.tan(theta_center/2.))
    phi_cell = m.copysign(m.acos(x_center/r_center), y_center)

  elif shoot_type == "cell_vertex":
    # Get position of vertices
    vertex = vertices(type, int(i), cell_side)
    if vertex_number >= len(vertex):
        raise Exception('Vertex {0} doesn\'t exist for type {1}'.format(vertex_number, type))
    x_vertex, y_vertex = vertex[vertex_number]
    x_vertex += x_center
    y_vertex += y_center
    r_vertex = m.sqrt(x_vertex**2 + y_vertex**2)
    theta_vertex = m.atan(r_vertex/z)
    eta_shoot = -m.log(m.tan(theta_vertex/2.))
    phi_cell = m.copysign(m.acos(x_vertex/r_vertex), y_vertex)

  elif shoot_type == "cell_edge":
    edges = edge_centers(type, int(i), cell_side)
    if edge_number >= len(edges):
        raise Exception('Edge {0} doesn\'t exist for type {1}'.format(edge_number, type))
    x_edge, y_edge = edges[edge_number]
    x_edge += x_center
    y_edge += y_center
    r_edge = m.sqrt(x_edge**2 + y_edge**2)
    theta_edge = m.atan(r_edge/z)
    eta_shoot = -m.log(m.tan(theta_edge/2.))
    phi_cell = m.copysign(m.acos(x_edge/r_edge), y_edge)

  else:
    raise Exception('The shoot type should be cell_center, cell_edge or cell_vertex')

  return eta_shoot, phi_cell
