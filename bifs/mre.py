from functools import wraps
import inspect

def mywrap(f):
    "wrap a function f"
    @wraps(f)
    def mywrapper(*args, **kwds):
        "just a silly example"
        return f(*args, **kwds)
    return mywrapper

class Bar:

    @mywrap
    def foo(self):
        return 4

def quickinfo(f):
    "quick info about signature of f"
    s = inspect.signature(f, follow_wrapped=False)
    p = s.parameters  # does not respond to usual dictionary protocol
    for k in p:
        v = p[k]
        print("{}: {} {}".format(k, v.kind, v))

b=Bar()
#print(b.foo(False)) # fails: foo takes 1 arg, 2 given
print(b.foo())
quickinfo(Bar.foo)
quickinfo(b.foo)
