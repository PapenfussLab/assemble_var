"""
diagram.Glyphs module
"""

from reportlab.graphics.widgetbase import Widget
from reportlab.lib.attrmap import *    
from reportlab.lib.validators import *
from reportlab.graphics import shapes
from reportlab.graphics.charts.textlabels import Label


class _Glyph(Widget):
    "Base class for Glyphs"
    
    _nodoc = 1
    _attrMap = AttrMap(
        x = AttrMapValue(isNumber, desc='symbol x coordinate'),
        y = AttrMapValue(isNumber, desc='symbol y coordinate'),
        dx = AttrMapValue(isNumber, desc='symbol x coordinate adjustment'),
        dy = AttrMapValue(isNumber, desc='symbol x coordinate adjustment'),
        length = AttrMapValue(isNumber, desc='symbol length'),
        height = AttrMapValue(isNumber, desc='symbol height'),
        aspectRatio = AttrMapValue(isNumber, desc='symbol aspect ratio'),
        fillColor = AttrMapValue(isColorOrNone, desc='symbol fill color'),
        strokeWidth = AttrMapValue(isNumber, desc='symbol outline stroke width'),
        strokeColor = AttrMapValue(isColorOrNone, desc='symbol outline stroke color'),
        label = AttrMapValue(isString, desc='label'),
        labeldx = AttrMapValue(isNumber, desc='label x coordinate adjustment'),
        labeldy = AttrMapValue(isNumber, desc='label x coordinate adjustment'),
        labelAngle = AttrMapValue(isNumber, desc='label angle'),
        fontSize = AttrMapValue(isNumber,desc='label font size'),
        boxAnchor = AttrMapValue(isString, desc='label box anchor (w,s,...)'),
        textAnchor = AttrMapValue(isString, desc='label text anchor (start,middle,end)'),
        wmin = AttrMapValue(isNumber, desc='minimum width (symbol is expanded to this)'),
        wNoTail = AttrMapValue(isNumber, desc='no tail on widths less than this'),
    )
    
    def __init__(self):
        assert self.__class__.__name__!='_Symbol', 'Abstract class _Symbol instantiated'
        self.x = self.y = self.dx = self.dy = 0
        self.length = 100
        self.height = 10
        self.aspectRatio = 0.5
        self.fillColor = colors.red
        self.strokeColor = None
        self.strokeWidth = 0
        self.fontSize = 12
        self.label = ''
        self.labeldx = 0
        self.labeldy = 4
        self.labelAngle = 0
        self.boxAnchor = 's'
        self.textAnchor = 'middle'
        self.wmin = 2.5
        self.wNoTail = 8


class Arrow(_Glyph):
    """This widget draws an arrow."""
    
    def __init__(self):        
        self.x = self.y = self.dx = self.dy = 0
        self.length = 100
        self.height = 20
        self.aspectRatio = 0.5
        self.fillColor = colors.red
        self.strokeColor = None
        self.strokeWidth = 0
        self.label = ''
        self.labeldx = self.labeldy = 0
        self.labelAngle = 0
        self.boxAnchor = 's'
        self.textAnchor = 'middle'
        self.fontSize = 14
        self.wmin = 3
        self.wNoTail = 6
    
    def draw(self):
        # general widget bits
        w = float(self.length)
        h = float(self.height)
        # print self.label,w,h
        
        # Set minimum size
        if abs(w)<self.wmin:
            xmid = self.x+0.5*w
            w = w/abs(w) * self.wmin
            self.x = xmid-0.5*w
        
        g = shapes.Group()
        if abs(w)>self.wNoTail:
            # arrow specific bits
            body = shapes.Rect(x=self.x, y=self.y-self.aspectRatio*h/2,
                width=2*(w/3),
                height=self.aspectRatio*h,
                fillColor=self.fillColor,
                strokeColor=self.strokeColor,
                strokeWidth=self.strokeWidth)
            g.add(body)
            
            head = shapes.Polygon(
                points=[self.x+w, self.y,
                    self.x+2*(w/3), self.y+h/2,
                    self.x+2*(w/3), self.y-h/2,
                    self.x+w, self.y],
                fillColor=self.fillColor,
                strokeColor=self.strokeColor,
                strokeWidth=self.strokeWidth)
            g.add(head)
        else:
            head = shapes.Polygon(
                points=[self.x+w, self.y,
                    self.x, self.y+h/2,
                    self.x, self.y-h/2,
                    self.x+w, self.y],
                fillColor=self.fillColor,
                strokeColor=self.strokeColor,
                strokeWidth=self.strokeWidth)
            g.add(head)
        
        if self.label:
            b = g.getBounds()
            s = Label()
            s.setText(self.label)
            s.setOrigin(self.x+0.5*w+self.labeldx, self.y-h/2+b[3]-b[1]+self.labeldy)
            s.boxAnchor = self.boxAnchor
            s.textAnchor = self.textAnchor
            s.fontName = 'Helvetica'
            s.fontSize = self.fontSize
            s.angle = self.labelAngle
            g.add(s)
        
        return g


class Block(_Glyph):
    """This widget draws a block."""
    
    def __init__(self):
        self.x = self.y = self.dx = self.dy = 0
        self.length = 100
        self.height = 20
        self.fillColor = colors.red
        self.strokeColor = None
        self.strokeWidth = 0
        self.label = ''
        
    def draw(self):
        # general widget bits
        w = float(self.length)
        h = float(self.height)
        g = shapes.Group()
        
        body = shapes.Rect(x=self.x, y=self.y-h/2,
            width=w, height=h,
            fillColor=self.fillColor,
            strokeColor=self.strokeColor,
            strokeWidth=self.strokeWidth)
        g.add(body)
        
        if self.label:
            b = g.getBounds()
            s = Label()
            s.setText(self.label)
            s.setOrigin(self.x+0.5*w, self.y-h/2+b[3]-b[1]+4)
            s.boxAnchor = self.boxAnchor
            s.textAnchor = self.textAnchor
            s.fontName = 'Helvetica'
            s.fontSize = self.fontSize
            s.angle = self.labelAngle
            g.add(s)
        
        return g


class Triangle(_Glyph):
    """This widget draws a triangle."""
    
    def __init__(self):
        self.x = self.y = self.dx = self.dy = 0
        self.length = 100
        self.height = 20
        self.fillColor = colors.red
        self.strokeColor = None
        self.strokeWidth = 0
        self.label = ''
        
    def draw(self):
        # general widget bits
        w = float(self.length)
        h = float(self.height)
        g = shapes.Group()
        
        body = shapes.Polygon(
            [self.x-0.5*w, self.y-0.5*w,
             self.x-0.5*w, self.y+0.5*w,
             self.x+0.5*w, self.y],
            fillColor=self.fillColor,
            strokeColor=self.strokeColor,
            strokeWidth=self.strokeWidth)
        g.add(body)
        
        if self.label:
            b = g.getBounds()
            s = Label()
            s.setText(self.label)
            s.setOrigin(self.x+0.5*w, self.y-h/2+b[3]-b[1]+4)
            s.boxAnchor = self.boxAnchor
            s.textAnchor = self.textAnchor
            s.fontName = 'Helvetica'
            s.fontSize = self.fontSize
            s.angle = self.labelAngle
            g.add(s)
        
        return g


def signum(x):
    if x>=0:
        return 1
    else:
        return -1


class PointyBlock(_Glyph):
    """This widget draws an arrow."""
    
    def __init__(self):        
        self.x = self.y = self.dx = self.dy = 0
        self.length = 100
        self.height = 20
        self.aspectRatio = 0.5
        self.fillColor = colors.red
        self.strokeColor = None
        self.strokeWidth = 0
    
    def draw(self):
        # general widget bits
        w = float(self.length)
        h = float(self.height)
        
        g = shapes.Group()
        block = shapes.Rect(x=self.x, y=self.y-h/2,
            width=self.length, height=h, fillColor=self.fillColor,
            strokeColor=self.strokeColor, strokeWidth=self.strokeWidth)
        g.add(block)
        
        point = shapes.Polygon(
            points=[self.x+w, self.y-h/2,
                self.x+w+signum(w)*h/4, self.y,
                self.x+w, self.y+h/2,
                self.x+w, self.y-h/2],
            fillColor=self.fillColor,
            strokeColor=self.strokeColor,
            strokeWidth=self.strokeWidth)
        g.add(point)
        return g
