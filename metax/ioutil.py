import gzip
import bz2

def compressed_open(fname, *args, **kwargs):
    if fname.endswith('.gz'):
        return gzip.open(fname, *args, **kwargs)
    elif fname.endswith('.bz2'):
        return bz2.open(fname, *args, **kwargs)
    else:
        return open(fname, *args, **kwargs)
