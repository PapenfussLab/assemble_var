"""
WIG module
"""

from useful import smartopen


def load(iFileHandle):
    """Load WIG file.
    
    @param iFileHandle: Input file or filename
    @return: (header, data)
    """
    iFile = smartopen(iFileHandle)
    header = iFile.readline()
    data = []
    for line in iFile:
        score = float(line.strip())
        data.append(score)
    return header,data
