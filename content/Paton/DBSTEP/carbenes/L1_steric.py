from pymol.cgo import *
from pymol import cmd

def cgoCircle(r=3.5, cr=1.0, cg=0.4, cb=0.8, w=8.0, z=0.0, name=None):
  """
  Create a CGO circle
  from https://pymolwiki.org/index.php/CgoCircle
  PARAMS
        x, y, z
          X, Y and Z coordinates of the origin
        r
          Radius of the circle
        cr, cg, cb
          Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].
        w
          Line width of the circle
  RETURNS
        the CGO object (it also loads it into PyMOL, too).

  """
  x, y = 0.0, 0.0

  r = abs(float(r))
  cr = abs(float(cr))
  cg = abs(float(cg))
  cb = abs(float(cb))
  w = float(w)

  obj = [ BEGIN, LINES, COLOR, cr, cg, cb ]
  for i in range(180):
        obj.append( VERTEX )
        obj.append(r*math.cos(i) + x )
        obj.append(r*math.sin(i) + y )
        obj.append(z)
        obj.append( VERTEX )
        obj.append(r*math.cos(i+0.1) + x )
        obj.append(r*math.sin(i+0.1) + y )
        obj.append(z)
  obj.append(END)
  if cb == 0.0:
        if cr == 0.0:
              cName = "BminC_"+str(z)
        else:
              cName = "BmaxC_"+str(z)
  else:
        if name is not None:
            cName = "circle_"+str(r)+name
        else:
            cName = "circle_"+str(r)
  cmd.load_cgo( obj, cName )
  cmd.set("cgo_line_width", w, cName )
  return obj

cmd.extend( "cgoCircle", cgoCircle )
cgoCircle(r=3.5,name="a")
cgoCircle(r=3.5,name="b")
cgoCircle(r=3.5,name="c")
cmd.rotate([1,0,0],90,object="circle_3.5b",origin=[0,0,0])
cmd.rotate([1,0,0],90,object="circle_3.5c",origin=[0,0,0])
cmd.rotate([0,0,1],90,object="circle_3.5c",origin=[0,0,0])

cylinder = [
]
cmd.load_cgo(cylinder, "axes")


cmd.load("/Users/rpaton/DBSTEP_example/carbenes/L1_transform.xyz")
cmd.show_as("spheres", "L1_transform")
cmd.set("sphere_transparency", 0.5)
cmd.set("orthoscopic", "on")
