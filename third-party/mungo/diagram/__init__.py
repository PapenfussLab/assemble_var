"""
diagram module
"""

from reportlab.graphics import renderPDF, renderSVG
from reportlab.lib import colors
from reportlab.graphics.shapes import *
from reportlab.graphics.charts.textlabels import Label
from Glyphs import *


DPI = 72

class A4Portrait:
    """A4 portrait dimensions"""
    margin = DPI
    height = 841.89 # 24*DPI/2.54
    width = 595.28 # 18*DPI/2.54


class A4Landscape:
    """A4 landscape dimensions"""
    margin = DPI
    height = 595.28 # 18/2.54*DPI
    width = 841.89 # 24/2.54*DPI


class Mapping:
    """Define a map between pixels and feature coordinates"""
    def __init__(self, x0, x1, start, end, flipped=False):
        """Constructor
        
        @param x0: Left-most pixel coord
        @param x1: Right-most pixel coord
        @param start: Left-most feature coord
        @param end: Right-most feature corrd
        @param flipped: Flip coords (default: False)
        """
        self.x0 = x0
        self.x1 = x1
        self.start = start
        self.end = end
        self.scale = float(self.x1-self.x0)/float(self.end-self.start)
        self.flipped = flipped
        if flipped:
            self.x0,self.x1 = self.x1,self.x0
            self.scale *= -1
    
    def flip(self):
        """Flip coords"""
        self.flipped = not self.flipped
        self.x0,self.x1 = self.x1,self.x0
        self.scale *= -1
    
    def __call__(self, i):
        """Return """
        return self.x0 + self.scale*(i-self.start)
    
    def shift(self, dx, newStart):
        """Shift page view; end-start is preserved; scale is not changed
        
        @param dx: Pixels to shift x0 & x1 by
        @param newStart: New feature start coord
        """
        self.x0 += dx
        self.x1 += dx
        L = self.end-self.start
        self.start = newStart
        self.end = newStart+L
    
    @staticmethod
    def fromScale(x0, start, end, scale, flipped=False):
        """Static method to construct a Mapping object with a fixed scale
        
        @param x0: Left-most pixel coord
        @param start: Left-most feature coord
        @param end: Right-most feature coord
        @param scale: Scale ie (x1-x0)/(end-start)
        @param flipped: Flip coords (default: False)
        """
        x1 = int(x0 + scale*(end-start))
        m = Mapping(x0, x1, start, end, flipped=flipped)
        return m


class UniformMapping:
    """Define a map for uniformly spaced genes"""
    def __init__(self, x0, x1, genes, flipped=False):
        """Constructor
        
        @param x0: Left-most pixel coord
        @param x1: Right-most pixel coord
        @param genes: List of genes. Each gene must be an object with start, end & strand attributes
        @param flipped: Flip coords (default: False)
        """
        self.x0 = x0
        self.x1 = x1
        self.flipped = flipped
        
        self.starts = [g.start for g in genes]
        self.starts.sort()
        if self.flipped: self.starts.reverse()
        self.positions = dict(zip(self.starts, range(len(self.starts))))
        self.start = self.starts[0]
        self.end = self.starts[-1]
    
    def flip(self):
        """Flip the coords"""
        self.flipped = not self.flipped
        self.starts = self.positions.keys()
        self.starts.sort()
        if self.flipped: self.starts.reverse()
        self.positions = dict(zip(self.starts, range(len(self.starts))))
    
    def __call__(self, start):
        """Return the pixel position of the gene at position start"""
        return self.x0 + self.positions[start]/float(len(self.positions))*(self.x1-self.x0)
    
    def __getitem__(self, i):
        """Return the pixel position of the i^th gene"""
        start = self.starts[i]
        return self.x0 + self.positions[start]/float(len(self.positions))*(self.x1-self.x0)


class UniformMapping2:
    def __init__(self, x0, x1, genes, flipped=False):
        self.x0 = x0
        self.x1 = x1
        self.flipped = flipped
        
        n = 2*len(genes)
        dx = (x1-x0)/float(n-1)
        self.pos = []
        self.x = []
        for i,g in enumerate(genes):
            self.pos.append(g.start)
            self.pos.append(g.end)
            self.x.append(x0 + i*dx)
            self.x.append(x0 + (i+1)*dx)
        self.pos.sort()
        if self.flipped: self.pos.reverse()
        self.start = self.pos[0]
        self.end = self.pos[-1]
    
    def flip(self):
        self.flipped = not self.flipped
    
    def __call__(self, p):
        for i in xrange(len(self.pos)):
            if p<self.pos[i]: break
        return self.x[i-1] + (self.x[i]-self.x[i-1])/(self.pos[i]-self.pos[i-1])*(p-self.pos[i-1])


def addHRule(drawing, x1, x2, y, strokeColor=colors.black, strokeWidth=0.5):
    """Add a horizontal line to the drawing.
    
    @param drawing:
    @param x1:
    @param x2:
    @param y:
    @param strokeColor:
    @param strokeWidth:
    """
    line = Line(x1, y, x2, y, strokeColor=strokeColor, strokeWidth=strokeWidth)
    drawing.add(line)


def addScale(drawing, xmap, y, start, end, tickLen=10, dx=3, dy=6,
  textAnchor='middle', boxAnchor='s', fontSize=12,
  strokeWidth=1, strokeColor=colors.black, scale=1.0, format='%ibp'):
    x1 = xmap(start)
    x2 = xmap(end)
    line = Line(x1+dx,y,x2-dx,y,
        strokeWidth=strokeWidth, strokeColor=strokeColor)
    drawing.add(line)
    
    leftTick = Line(x1+dx,y-0.5*tickLen,x1+dx,y+0.5*tickLen,
        strokeWidth=strokeWidth, strokeColor=strokeColor)
    drawing.add(leftTick)
    
    rightTick = Line(x2-dx,y-0.5*tickLen,x2-dx,y+0.5*tickLen,
        strokeWidth=strokeWidth, strokeColor=strokeColor)
    drawing.add(rightTick)
    
    label = Label()
    label.setOrigin(0.5*(x1+x2), y+dy)
    
    distance = float(end-start)/scale
    label.setText(format % (distance/scale))
    label.fontSize = fontSize
    label.textAnchor = textAnchor
    label.boxAnchor = boxAnchor
    drawing.add(label)

addRuler = addScale


def tick_generator(start, end, n=10, convert=None):
    """Generate tick positions for addAxis"""
    dp = float(end-start)/n
    for i in xrange(n+1):
        t = start + i*dp
        if convert:
            t = convert(t)
        yield t


def addAxis(drawing, xmap, y, strokeWidth=1, minorStrokeWidth=0.5, 
  tickDir='down', autoTicks=False, nTicks=20, tickLen=5, fontSize=10, nMinorTicks=80, 
  minorTickLen=2, angle=0, dx=0, dy=-2, textAnchor='middle', boxAnchor=None, 
  scale=1.0, format='%i'):
    """Add a horizontal axis to the drawing.
    
    To do: Round tick positions
    """
    line = Line(xmap.x0, y, xmap.x1, y, strokeWidth=strokeWidth)
    drawing.add(line)
    
    if not boxAnchor:
        if tickDir=='down':
            boxAnchor = 'n'
        else:
            boxAnchor = 's'
    signum = {'up': -1, 'down': 1}[tickDir]
    
    if nTicks>0:
        ticks = tick_generator(xmap.start, xmap.end, n=nTicks, convert=int)
    
    for p in ticks:
        x = xmap(p)
        line = Line(x, y, x, y-signum*tickLen, strokeWidth=strokeWidth)
        drawing.add(line)
        s = Label()
        s.setOrigin(x, y-signum*tickLen)
        s.setText(format % (p/scale))
        s.dx = dx
        s.dy = signum*dy
        s.fontName = 'Helvetica'
        s.fontSize = fontSize
        s.textAnchor = textAnchor
        s.boxAnchor = boxAnchor
        s.angle = angle
        drawing.add(s)
    
    minorticks = tick_generator(xmap.start, xmap.end, n=nMinorTicks, convert=int)
    for p in minorticks:
        x = xmap(p)
        line = Line(x, y, x, y-signum*minorTickLen, strokeWidth=minorStrokeWidth)
        drawing.add(line)


def addUniformAxis(self, xmap, y):
    """Add a horizontal axis suitable for uniformly a spaced gene map to the drawing.
    
    Not yet finished.
    """
    line = Line(xmap.x0, y, xmap.x1, y, strokeWidth=strokeWidth)
    drawing.add(line)


def addLabel(drawing, x, y, text, fontName='Helvetica', fontSize=11, dy=0,
             angle=0, boxAnchor='sw', textAnchor='start'):
    """Add a label to the drawing. 
    This interface here is inconsistent in that it requires pixel coords. FIX
    This just sets convenient defaults for Label."""
    label = Label()
    label.setText(text)
    label.setOrigin(x, y)
    label.fontName = fontName
    label.fontSize = fontSize
    label.boxAnchor = boxAnchor
    label.textAnchor = textAnchor
    label.dy = dy
    label.angle = angle
    drawing.add(label)


def addBlock(drawing, xmap, y, block, height=10, width=None, fillColor=colors.red, strokeColor=None):
    """Add a colored block to the drawing.
    This just sets convenient defaults for Rect."""
    x = xmap(block.start)
    if not width:
        width = xmap(block.end)-xmap(block.start)
    r = Rect(x,y,width,height,strokeColor=strokeColor,fillColor=fillColor)
    drawing.add(r)


def addFixedLengthFeature(drawing, xmap, y, feature, length, glyph=Arrow, 
  fillColor=colors.red, strokeColor=None,
  height=20, fontSize=14, labeldy=3, labelAngle=90, textAnchor='start', 
  boxAnchor='w', aspectRatio=0.5, wmin=3, wNoTail=6):
    x1 = xmap(feature.start)
    x2 = x1 + length
    if feature.strand=='-':
        x1,x2 = x2,x1
    
    g = glyph()
    g.x = x1
    g.y = y+height/2
    g.height = height
    g.length = x2-x1
    g.fillColor = fillColor
    if strokeColor:
        g.strokeColor = strokeColor
    else:
        g.strokeColor = fillColor
    g.fontSize = fontSize
    g.label = feature.name
    g.labeldy = labeldy
    g.labelAngle = labelAngle
    g.textAnchor = textAnchor
    g.boxAnchor = boxAnchor
    g.aspectRatio = aspectRatio
    g.wmin = wmin
    g.wNoTail = wNoTail
    drawing.add(g)


def addFeature(drawing, xmap, y, feature, glyph=Arrow, 
  fillColor=colors.red, strokeColor=None, strokeWidth=0,
  height=20, fontSize=14, labeldx=0, labeldy=3, labelAngle=90, textAnchor='start', 
  boxAnchor='w', aspectRatio=0.5, wmin=3, wNoTail=6):
    """Adds a feature (typically an arrow) with label to the drawing"""
    if feature.strand=='+':
        x1,x2 = xmap(feature.start), xmap(feature.end)
    else:
        x2,x1 = xmap(feature.start), xmap(feature.end)
    
    g = glyph()
    g.x = x1
    g.y = y+height/2
    g.height = height
    g.length = x2-x1
    g.fillColor = fillColor
    g.strokeWidth = strokeWidth
    if strokeColor:
        g.strokeColor = strokeColor
    else:
        g.strokeColor = fillColor
    
    g.fontSize = fontSize
    g.label = feature.name
    g.labeldx = labeldx
    g.labeldy = labeldy
    g.labelAngle = labelAngle
    g.textAnchor = textAnchor
    g.boxAnchor = boxAnchor
    g.aspectRatio = aspectRatio
    g.wmin = wmin
    g.wNoTail = wNoTail
    drawing.add(g)


def addCompoundFeature(drawing, xmap, y, gene, 
  strokeColor=None, fillColor=colors.blue, 
  intronColor=colors.blue, intronWidth=0.5,
  glyph=Block, height=12, utrHeight=6,
  labeldy=10, fontSize=10, textAnchor='middle', boxAnchor='s'):
    """Adds a compund feature to the drawing.
    A compound feature is typically several exons joined by zig-zag lines."""
    rise = height + utrHeight
    
    intronStarts = [None]
    intronEnds = []
    heights = []
    for exon in gene:
        x1,x2 = xmap(exon.start), xmap(exon.end)
        
        kind = exon.kind.lower()
        if kind in ['exon', 'utr']:
            intronStarts.append(exon.end)
            intronEnds.append(exon.start)
        
        g = glyph()
        g.x = x1
        g.y = y+height/2
        if exon.kind.lower()=='exon':
            g.height = height
            heights.append(height)
        else:
            g.height = utrHeight
            heights.append(utrHeight)
        
        g.length = x2-x1
        g.fillColor = fillColor
        if strokeColor:
            g.strokeColor = strokeColor
        else:
            g.strokeColor = fillColor
        g.fontSize = fontSize
        drawing.add(g)
    
    for i,(intronStart,intronEnd) in enumerate(zip(intronStarts[1:], intronEnds[1:])):
        x1 = xmap(intronStart)
        x2 = xmap(0.5*(intronStart+intronEnd))
        x3 = xmap(intronEnd)
        # if abs(x3-x1)<3: continue
        # print intronStart,intronEnd,heights[i],heights[i+1]
        
        y1 = y+heights[i]/2+height/2
        y2 = y+rise
        y3 = y+heights[i+1]/2+height/2
        
        line1 = Line(x1,y1,x2,y2,strokeColor=intronColor,strokeWidth=intronWidth)
        line2 = Line(x2,y2,x3,y3,strokeColor=intronColor,strokeWidth=intronWidth)
        drawing.add(line1)
        drawing.add(line2)
    
    # Draw arrows
    if xmap.flipped:
        signum = -1
    else:
        signum = 1
    
    if gene.strand=='+':
        x1 = xmap(gene.end)
        x2 = x1 + signum*15
        x3 = x1 + signum*10
        y1 = y + 0.5*height
        y2 = y + 0.75*height
        y3 = y + 0.25*height
        line1 = Line(x1,y1,x2,y1,strokeColor=intronColor,strokeWidth=intronWidth)
        line2 = Line(x2,y1,x3,y2,strokeColor=intronColor,strokeWidth=intronWidth)
        line3 = Line(x2,y1,x3,y3,strokeColor=intronColor,strokeWidth=intronWidth)
        drawing.add(line1)
        drawing.add(line2)
        drawing.add(line3)
    else:
        x1 = xmap(gene.start)
        x2 = x1 - signum*15
        x3 = x1 - signum*10
        y1 = y + 0.5*height
        y2 = y + 0.75*height
        y3 = y + 0.25*height
        line1 = Line(x1,y1,x2,y1,strokeColor=intronColor,strokeWidth=intronWidth)
        line2 = Line(x2,y1,x3,y2,strokeColor=intronColor,strokeWidth=intronWidth)
        line3 = Line(x2,y1,x3,y3,strokeColor=intronColor,strokeWidth=intronWidth)
        drawing.add(line1)
        drawing.add(line2)
        drawing.add(line3)
    
    # if gene has attribute name...
    label = Label()
    label.setText(gene.name)
    pos = 0.5*(gene.start+gene.end)
    x = xmap(pos)
    label.setOrigin(x,y)
    label.dy = labeldy
    label.textAnchor = textAnchor
    label.boxAnchor = boxAnchor
    drawing.add(label)


def addPointyCompoundFeature(drawing, xmap, y, gene, 
  strokeColor=None, fillColor=colors.blue, intronColor=colors.blue,
  glyph=PointyBlock, height=12, utrHeight=6, rise=8, 
  labeldy=10, fontSize=10, textAnchor='middle', boxAnchor='s'):
    """Adds a pointy compound feature to the drawing. This is typically
    several exons joined by zig-zag lines with an arrow showing strand."""
    if gene.strand=='+':
        x1,x2 = xmap(gene.start), xmap(gene.end)
    else:
        x2,x1 = xmap(gene.start), xmap(gene.end)
    y = y+height/2
    y1 = y
    line = Line(x1,y1,x2,y1,strokeColor=intronColor)
    drawing.add(line)
    
    for exon in gene:
        if exon.strand=='+':
            x1,x2 = xmap(exon.start), xmap(exon.end)
        else:
            x2,x1 = xmap(exon.start), xmap(exon.end)
        
        g = glyph()
        g.x = x1
        g.y = y
        if exon.kind.lower()=='utr':
            g.height = utrHeight
        else:
            g.height = height
        g.length = x2-x1
        g.fillColor = fillColor
        if strokeColor:
            g.strokeColor = strokeColor
        else:
            g.strokeColor = fillColor
        g.fontSize = fontSize
        drawing.add(g)
    
    label = Label()
    label.setText(gene.name)
    x = 0.5*(gene.start+gene.end)
    label.setOrigin(x,y)
    label.dy = labeldy
    label.textAnchor = textAnchor
    label.boxAnchor = boxAnchor
    drawing.add(label)


def addCropMarks(drawing, xsize, ysize, margin, L=10, verbose=False):
    """Add crop marks to the drawing. This is helpful for lining up several pdfs in Illustrator."""
    x1 = margin
    x2 = xsize-margin
    y1 = margin
    y2 = ysize-margin
    if verbose:
        print x1,y1
        print x2,y2
    
    drawing.add(Line(x1,y1,x1+L,y1,strokeColor=colors.black))
    drawing.add(Line(x1,y1,x1,y1+L,strokeColor=colors.black))
    drawing.add(Line(x2,y2,x2-L,y2,strokeColor=colors.black))
    drawing.add(Line(x2,y2,x2,y2-L,strokeColor=colors.black))
    drawing.add(Line(x1,y2,x1+L,y2,strokeColor=colors.black))
    drawing.add(Line(x1,y2,x1,y2-L,strokeColor=colors.black))
    drawing.add(Line(x2,y1,x2-L,y1,strokeColor=colors.black))
    drawing.add(Line(x2,y1,x2,y1+L,strokeColor=colors.black))

