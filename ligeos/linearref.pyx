"""
This package is designed to offer FAST linear referencing tools without outside dependencies
on other libraries (except cython).  

Much of the code is a direct translation from the PostGIS repository's liblwgeom
abstraction. One important note: Points maintain a "geo" flag, while indicates the point 
is (or is not) a geographic coordinate.  Geographic coordinates use a spherical distance 
calculation rather than a planar distance calculation.
"""


cdef extern from "math.h":
    double asin(double d)
    double sin(double d)
    double cos(double d)
    double fabs(double d)
    double sqrt(double d)


# Constants
DEF __PI = 3.14159265
DEF EARTH_RADIUS = 6370986.884258304


cpdef double distance_earth(float x1, float y1, float x2, float y2):
    """ Calculates the spherical distance between two geographic points."""
    cdef double long1, lat1, long2, lat2
    cdef double longdiff
    cdef double sino
    
    long1 = -2 * (x1 / 360.0) * __PI
    lat1 = 2 * (y1 / 360.0) * __PI

    long2 = -2 * (x2 / 360.0) * __PI
    lat2 = 2 * (y2 / 360.0) * __PI

    # compute difference in longitudes - want < 180 degrees
    longdiff = fabs(long1 - long2)
    if longdiff > __PI:
        longdiff = (2 * __PI) - longdiff;

    sino = sqrt(sin(fabs(lat1 - lat2) / 2.) * sin(fabs(lat1 - lat2) / 2.) + \
            cos(lat1) * cos(lat2) * sin(longdiff / 2.) * sin(longdiff / 2.))
    if sino > 1.0: sino = 1.0
    return 2.0 * EARTH_RADIUS * asin(sino);

#def distance_earth(float x1, float y1, float x2, float y2):
#    return cdistance_earth(x1,y1,x2,y2)

    
cdef class CPoint:
    """ A coordinate pair.  If geo is set, distance_earth will be used in distance calculations."""
    cdef double _x, _y
    cdef bool geo

    def __init__(self, double x, double y, bool geo):
        self._x = x; self._y = y; self.geo = geo;
        
    def __str__(self):
        return "(%f %f)" % (self._x, self._y)
    
    def clone(self):
        """ Return a copy of this point."""
        return CPoint(self._x, self._y, self.geo)
    
    property x:
        def __set__(self, double v):
            self._x = v
        def __get__(self):
            return self._x
    
    property y:
        def __set__(self, double v):
            self._y = v
        def __get__(self):
            return self._y
    
    def is_geo(self):
        """Is this a geographic point."""
        return self.geo == True
        
        
    
    cpdef double distance_pt(self, CPoint p2):
        if self.geo:
            return distance_earth(self._x, self._y, p2._x, p2._y)
        # else
        cdef double hside = p2._x - self._x
        cdef double vside = p2._y - self._y
        return ( hside*hside + vside*vside )**0.5
    

def Point(float x, float y):
    """Helper function to create a simple coordinate."""
    return CPoint(x, y, False)

def GPoint(float x, float y):
    """ Helper function to return a geographic coordinate."""
    return CPoint(x, y, True)

cdef class Line:
    """ A line as defined by two points."""
    cdef CPoint p1, p2
    
    def __init__(self, CPoint p1, CPoint p2):
        self.p1 = p1; self.p2 = p2;
        
    def __str__(self):
        return "%f,%f->%f,%f" %(self.p1._x, self.p1._y, self.p2._x, self.p2._y)

    cpdef double length(self):
        """ Length of the line."""
        return self.p1.distance_pt(self.p2)

    cpdef CPoint closest_pt(self, CPoint p):
        """Return the point on the line nearest the given point."""
        return seg_closest_pt(self.p1, self.p2, p)
    
    cpdef double distance_pt(self, CPoint p):
        """The distance from the nearest point on this line to the given point."""
        return seg_distance_pt(self.p1, self.p2, p)

    cpdef double distance_seg(self, Line s):
        """The closest distance from this line to another line."""
        return seg_distance_seg(self.p1, self.p2, s.p1, s.p2)

cdef struct Vector2D:
    double x, y

cpdef CPoint seg_intersect_seg(CPoint a, CPoint b, CPoint c, CPoint d):
    if ( a._x == b._x) and (a._y == b._y):
        if seg_distance_pt(c, d, a) == 0:
            return a.clone()
        return None

    # U and V are the same point
    if ( c._x == d._x) and (c._y == d._y):
        if seg_distance_pt(a, b, c) == 0:
            return c.clone()
        return None
    
    if (( a._x == c._x ) and ( a._y == c._y )) or (( a._x == d._x ) and ( a._y == d._y )):
        return a.clone()
    
    if ( b._x == c._x ) and ( b._y == c._y ) or (( b._x == d._x ) and ( b._y == d._y )):
        return b.clone()

    cdef float r_bot = (b._x-a._x)*(d._y-c._y) - (b._y-a._y)*(d._x-c._x)
    # colinear or parallel lines
    if r_bot == 0:
        return None    
    
    cdef double x, y
    x = ((a._x*b._y - a._y*b._x)*(c._x - d._x) - (a._x - b._x)*(c._x*d._y - c._y*d._x)) / ((a._x - b._x)*(c._y - d._y) - (a._y - b._y)*(c._x - d._x))
    y = ((a._x*b._y - a._y*b._x)*(c._y - d._y) - (a._y - b._y)*(c._x*d._y - c._y*d._x)) / ((a._x - b._x)*(c._y - d._y) - (a._y - b._y)*(c._x - d._x))
    return CPoint(x, y, a.geo)
    

cpdef double seg_distance_seg(CPoint a, CPoint b, CPoint c, CPoint d):
    """The closest distance from one line to another line."""
    # A and B are the same point 
    if ( a._x == b._x) and (a._y == b._y):
        return seg_distance_pt(c, d, a)

    # U and V are the same point
    if ( c._x == d._x) and (c._y == d._y):
        return seg_distance_pt(a, b, c)
    
    """
    AB and CD are line segments
    from comp.graphics.algo

    Solving the above for r and s yields
                (Ay-Cy)(Dx-Cx)-(Ax-Cx)(Dy-Cy)
               r = ----------------------------- (eqn 1)
                (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)

             (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = ----------------------------- (eqn 2)
            (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
    Let P be the position vector of the intersection point, then
        P=A+r(B-A) or
        Px=Ax+r(Bx-Ax)
        Py=Ay+r(By-Ay)
    By examining the values of r & s, you can also determine some other limiting conditions:
        If 0<=r<=1 & 0<=s<=1, intersection exists
        r<0 or r>1 or s<0 or s>1 line segments do not intersect
        If the denominator in eqn 1 is zero, AB & CD are parallel
        If the numerator in eqn 1 is also zero, AB & CD are collinear.
    """
    cdef float r_top = (a._y-c._y)*(d._x-c._x) - (a._x-c._x)*(d._y-c._y)
    cdef float r_bot = (b._x-a._x)*(d._y-c._y) - (b._y-a._y)*(d._x-c._x)

    cdef float s_top = (a._y-c._y)*(b._x-a._x) - (a._x-c._x)*(b._y-a._y)
    cdef float s_bot = r_bot
    
    if  r_bot==0 or s_bot == 0:
        return min(seg_distance_pt(c, d, a),
                   min(seg_distance_pt(c, d, b),   
                       min(seg_distance_pt(a, b, c), seg_distance_pt(a, b, d))))
    cdef float sv = s_top/s_bot
    cdef float r = r_top/r_bot

    if r<0 or r>1 or sv<0 or sv>1:
        #no intersection 
        #print "R/S", a, b, c, d, r, sv, s_top, s_bot, r_top, r_bot
        return min(seg_distance_pt(c, d, a),
                   min(seg_distance_pt(c, d, b),   
                       min(seg_distance_pt(a, b, c), 
                           seg_distance_pt(a, b, d))))

    else:
        return 0#; intersection exists 
        
cpdef double seg_distance_pt(CPoint a, CPoint b, CPoint p):
    """The closest distance from a line to a point."""
    cdef CPoint c = seg_closest_pt(a, b, p)
    #print "Closest", c
    return c.distance_pt(p)

cpdef CPoint seg_closest_pt(CPoint a, CPoint b, CPoint p):
    """The point from a line to another line."""
    if a._x == b._x and b._y == a._y:
        #print "Case 1"
        return a
    
    """
     * otherwise, we use comp.graphics.algorithms
     * Frequently Asked Questions method
     *
     *  (1)               AC dot AB
         *         r = ---------
         *               ||AB||^2
     *    r has the following meaning:
     *    r=0 P = A
     *    r=1 P = B
     *    r<0 P is on the backward extension of AB
     *    r>1 P is on the forward extension of AB
     *    0<r<1 P is interior to AB
     *
     * Further ref: http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
    """

    r = ((p._x-a._x) * (b._x-a._x) + (p._y-a._y) * (b._y-a._y)) / ((b._x-a._x)*(b._x-a._x) +(b._y-a._y)*(b._y-a._y))
    #print "r", r
    if r<0: 
        #print "Case 2"
        return a.clone()
    if r>1: 
        #print "Case 3"
        return b.clone()

    return CPoint(a._x + r*(b._x - a._x), a._y + r*(b._y - a._y), a.geo)

cdef interpolate_point2d(CPoint A, CPoint B, float F):
        return CPoint(A.x+((B.x-A.x)*F), A.y+((B.y-A.y)*F), A.geo)

cdef CPoint _to_cpoint(object arg, object geographic):
    if type(arg) == tuple or type(arg) == list:
        return CPoint(float(arg[0]), float(arg[1]), geographic)
    if hasattr(arg,"x") and hasattr(arg,"y"):
        return CPoint(float(arg.x), float(arg.y), geographic)
    raise TypeError("unable to convert %s to a coordinate" % type(arg))

cdef class LineString:
    """A sequence of points that define a linestring."""
    cdef object _pts
    cdef object geographic # cant be a bool, because cython chokes on it.
    def __init__(self, *args, geographic=False):
        self._pts = []
        self.geographic = bool(geographic) == True
        pts = []
        for a in args:
            #print type(a), type(a[0])
            if type(a) == list or type(a) == tuple and type(a[0]) == tuple:
                self._pts.extend([_to_cpoint(_a, geographic) for _a in a])
            else:
                self._pts.append(_to_cpoint(a, geographic))
    
    property pts:
        def __get__(self):
            return self._pts
    
    def __len__(self):
        return len(self._pts)
    
    def __str__(self):
        return "(" + ",".join(map(str, self._pts)) + ")"
    
    property coords:
        def __get__(self):
            return [(c.x, c.y) for c in self._pts]
    
    def length(self):
        """Length of all segments comprising this line."""
        return sum([self._pts[i-1].distance_pt(self._pts[i]) for i in range(1,len(self._pts))])
                
    cpdef CPoint closest_pt(self, CPoint p):
        """The nearest point on this linestring to a given point."""
        cdef int min_dist = -1
        closest = None
        for i in range(1,len(self._pts)):
            a = self._pts[i-1]
            b = self._pts[i]
            c = seg_closest_pt(a, b, p)
            d = c.distance_pt(p)
            if d == 0:
                return c
            if not closest or min_dist > d:
                min_dist = d
                closest = c
        return closest
    
    cpdef double distance_pt(self, CPoint p):
        """The distance from the nearest point on this linestring to the given point."""
        return self.closest_pt(p).distance_pt(p)
    
    cpdef CPoint closest_pt_to_line(self, LineString l2):
        """The distance from this line to another."""
        cdef int min_dist = -1
        cdef CPoint closest = None
        cdef double closest_u, closest_v
        for i in range(1,len(self._pts)):
            a = self._pts[i-1]
            b = self._pts[i]
            for j in range(1,len(l2._pts)):
                u = l2._pts[j-1]
                v = l2._pts[j]
                d = seg_distance_seg(a, b, u, v)
                if d == 0:
                    return seg_intersect_seg(a, b, u, v)
                if min_dist == -1 or min_dist > d:
                    min_dist = d
                    closest_u = seg_distance_pt(a, b, u)
                    closest_v = seg_distance_pt(a, b, v)
                    if seg_distance_pt(a, b, u) < seg_distance_pt(a, b, v):
                        closest = seg_closest_pt(a, b, u)
                    else:
                        closest = seg_closest_pt(a, b, v)
        return closest
    
    cpdef double distance_to_line(self, LineString l2):
        """The distance from this line to another."""
        cdef int min_dist = -1
        for i in range(1,len(self._pts)):
            a = self._pts[i-1]
            b = self._pts[i]
            for j in range(1,len(l2._pts)):
                u = l2._pts[j-1]
                v = l2._pts[j]
                d = seg_distance_seg(a, b, u, v)
                #print "SDS",d
                if d == 0:
                    return 0
                if min_dist == -1 or min_dist > d:
                    min_dist = d
        return d

    
    cpdef double locate_point(self, CPoint p):
        """Find the measure [0,1] on this linestring to which the given point falls nearest."""
        cdef CPoint nearest = None
        cdef double mindist = -1 # min dist found
        cdef CPoint start = None # current 
        cdef CPoint end = None # current + 1
        cdef int seg = -1 # closest segment index
        cdef double dist

        start = self._pts[0]
        for i in range(1, len(self)):
            end = self._pts[i]
            dist = seg_distance_pt(start, end, p)
            if i == 1 or dist < mindist:
                mindist = dist
                seg = i-1
            if mindist <= 0.0: 
                break
            
            start = end
        if mindist > 0:
            nearest = seg_closest_pt(self._pts[seg], self._pts[seg+1], p)
        else:
            nearest = p
    
        #print "Min dist", mindist, nearest, seg_closest_pt(self._pts[seg+1], self._pts[seg+2], p)
        
        
        cdef double tlen = self.length() # total length
        cdef double plen = 0 # length so far
        
        start = self._pts[0]
        for i in range(0, seg):
            end = self._pts[i+1]
            plen = plen + start.distance_pt(end)
            start = end
        plen += start.distance_pt(nearest)
        #print plen, tlen, start.distance_pt(nearest), plen/tlen
        return plen / tlen
    
    cpdef LineString substring(self, double frm, double to):
        """Returns the substring of this linestring between the two measures."""
        cdef int state = 0 # 0=before, 1=inside        
        cdef object dpa = [] # the output array
        cdef object ipa = self._pts # the input array
        cdef double length = self.length()
        cdef double tlength = 0 # traversed length
        cdef CPoint p1 = ipa[0]
        cdef CPoint p2 = None

        if frm > to:
            raise ArgumentError("from is less than to")
        # Get 'from' and 'to' lengths 
        frm = length*frm
        to = length*to

        for i in range(0, len(ipa)-1):
            p2 = ipa[i+1]
            #LWDEBUGF(3 ,"Segment %d: (%g,%g,%g,%g)-(%g,%g,%g,%g)",
            #    i, p1.x, p1.y, p1.z, p1.m, p2.x, p2.y, p2.z, p2.m);
            # Find the length of this segment 
            slength = p1.distance_pt(p2)
    
            """
             * We are before requested start.
            """
            if state == 0: #  before 
                """
                 * Didn't reach the 'from' point,
                 * nothing to do
                """
                if frm > tlength + slength:
                    tlength += slength
                    p1 = p2
                    continue
                    # goto END
    
                elif frm == tlength + slength:
                    #LWDEBUG(3, "  Second point is our start");
                    """
                     * Second point is our start
                    """
                    dpa.append(p2)
                    state=1 #;  we're inside now 
                    
                elif frm == tlength:
                    #LWDEBUG(3, "  First point is our start");
                    """
                     * First point is our start
                    """
                    dpa.append(p1)
                    """
                     * We're inside now, but will check
                     * 'to' point as well
                    """
                    state=1
    
                else:  # tlength < from < tlength+slength 
                    """
                     * Our start is between first and
                     * second point
                    """
                    dseg = (frm- tlength) / slength;
                    pt = interpolate_point2d(p1, p2,dseg)
                    dpa.append(pt)
                    """
                     * We're inside now, but will check
                     * 'to' point as well
                     """
                    state=1
    
            if state == 1: # inside     
                #LWDEBUG(3, " Inside");
                """
                 * Didn't reach the 'end' point,
                 * just copy second point
                 """
                if to > tlength + slength:
                    dpa.append(p2)
                    tlength += slength
                    p1 = p2
                    continue
                    #goto END;
                    
                elif to == tlength + slength:
                    """
                     * 'to' point is our second point.
                    """
                    dpa.append(p2)
                    break #  substring complete 
    
                elif to == tlength:
                    """
                     * 'to' point is our first point.
                     * (should only happen if 'to' is 0)
                    """
                    dpa.append(p1)
                    break#;  substring complete 
    
               
                elif ( to < tlength + slength ):
                    """
                     * 'to' point falls on this segment
                     * Interpolate and break.
                    """
                    dseg = (to - tlength) / slength
                    pt = interpolate_point2d(p1, p2, dseg)
                    dpa.append(pt)
                    break
                else:
                    print "Unhandled case"
            
            tlength += slength
            p1 = p2
    
        return LineString(dpa, geographic=ipa[0].is_geo())
    
    cpdef LineString concatenate(self, LineString l2, float merge_tolerance=0):
        """Returns true if l2 was joined or None, if the min distance is above the max distance allowed."""
        cdef CPoint a, b, c, d
        a = self._pts[0]
        b = self._pts[-1]
        c = l2._pts[0]
        d = l2._pts[-1]
        cdef double ac = a.distance_pt(c)
        cdef double ad = a.distance_pt(d)
        cdef double bc = b.distance_pt(c)
        cdef double bd = b.distance_pt(d)
        
        if ac < min(ad, min(bc, bd)):
            #print "Case 1", a, c, ac
            if ac <= merge_tolerance:
                return LineString(self._pts[::-1], l2._pts[1:], geographic=self.geographic)
            else:
                return LineString(self._pts[::-1], l2._pts, geographic=self.geographic)

        elif ad < min(bc, bd):
            #print "Case 2", a, d, ad
            if ad <= merge_tolerance:
                return LineString(l2._pts, self._pts[1:], geographic=self.geographic)
            else:
                return LineString(l2._pts, self._pts, geographic=self.geographic)
        elif bc < bd:
            #print "Case 3", b, c, bc
            if bc <= merge_tolerance:
                return LineString(self._pts, l2._pts[1:], geographic=self.geographic)
            else:
                return LineString(self._pts, l2._pts, geographic=self.geographic)
        else:
            #print "Case 4", b, d, bd
            if bd <= merge_tolerance:
                return LineString(self._pts, l2._pts[1::-1], geographic=self.geographic)
            else:
                return LineString(self._pts, l2._pts[::-1], geographic=self.geographic)
        
    
class ArgumentError(Exception):
    pass
