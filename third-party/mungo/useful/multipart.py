"""
Multipart form module
"""

import sys
import httplib
import mimetypes
import urlparse


class FormSubmissionException(Exception):
    pass


def post(url, fields, files, proxy=None, proxyPort=None):
    """
    Post fields and files to an http host as multipart/form-data.
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded as files
    Return the server's response page.
    """
    content_type, body = encode_formdata(fields, files)
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; U; PPC Mac OS X Mach-O; en-GB; rv:1.7.10) Gecko/20050717 Firefox/1.0.6',
        'Content-type': content_type
    }
    
    if proxy and proxyPort:
        h = httplib.HTTPConnection(proxy, proxyPort)
        h.request('POST', url, body, headers)
    else:
        urlparts = urlparse.urlsplit(url)
        h = httplib.HTTPConnection(urlparts[1])
        h.request('POST', urlparts[2], body, headers)
    
    try:
        resp = h.getresponse()
    except httplib.BadStatusLine:
        print >> sys.stderr, "Weird status error encountered."
        raise FormSubmissionException
    
    if resp.status==501:
        print >> sys.stderr, "Error 501: This doesn't seem to be a valid GenomeScan server."
        raise FormSubmissionException
    elif resp.status!=200:
        print >> sys.stderr, "Error %d: %s"  % (resp.status, resp.reason)
        raise FormSubmissionException
    else:
        data = resp.read()
        return data


def encode_formdata(fields, files):
    """
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded as files
    Return (content_type, body) ready for httplib.HTTP instance
    """
    BOUNDARY = '----------ThIs_Is_tHe_bouNdaRY_$'
    CRLF = '\r\n'
    L = []
    for (key, value) in fields:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % key)
        L.append('')
        L.append(value)
    for (key, filename, value) in files:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename))
        L.append('Content-Type: %s' % get_content_type(filename))
        L.append('')
        L.append(value)
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
    return content_type, body


def get_content_type(filename):
    return mimetypes.guess_type(filename)[0] or 'application/octet-stream'
