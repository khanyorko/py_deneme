class angle:
    def __init__(self, arg = 0, base = "rad"):
        consts          = {"deg": 180, "rad": pi, "gon": 200}
        self.base       = base
        
        if base == "dms":
            arg         = self.dms2deg(arg)
            base        = "deg"
            
        if not isinstance(arg, (int, float)): 
            raise ValueError("Given data type in not valid for an angle class!")
            
        self.val        = {i:consts[i]/consts[base]*arg for i in consts.keys()}
        self.val["dms"] = self.deg2dms(self.val["deg"])
        
    def __repr__(self, precision = 5):
        if self.base == "dms": return self.val[self.base] + "\tdms"
        return format(getattr(self,"val")[self.base], ".{}f".format(precision)) +"  "+self.base
    
    def sin(self):
        return sin(self.val["rad"])
    
    def cos(self):
        return cos(self.val["rad"])
   
    def __add__(self, other):
        return angle(self.val[self.base]+other.val[self.base], self.base)
    
    def __sub__(self, other):
        return angle(self.val[self.base]+other.val[self.base], self.base)
   
    def __mul__(self, const):
        if self.base == "dms": 
            return angle(self.deg2dms(self.dms2deg(self.val["dms"]) * const), "dms")
        return angle(self.val[self.base] * const, self.base)
    
    def __truediv__(self, const):
        if self.base == "dms": 
            return angle(self.deg2dms(self.dms2deg(self.val["dms"]) / const), "dms")
        return angle(self.val[self.base] / const, self.base)
    
    def __floordiv__(self, const):
        if self.base == "dms": 
            return angle(self.deg2dms(self.dms2deg(self.val["dms"]) // const), "dms")
        return angle(self.val[self.base] // const, self.base)
    
    def deg2dms(self, deg):
        d = int(deg)
        m = int((deg -d) * 60)
        s = (deg - d - m / 60) * 3600
        return "{} {} {:7.4f}".format(d,m,s)
    
    def dms2deg(self, dms):
        d, m, s = list(map(float, dms.split()))
        deg = d + m/60 + s/3600
        return deg
