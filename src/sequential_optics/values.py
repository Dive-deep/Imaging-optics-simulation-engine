"""Small immutable mapping that retains dataclass/JSON dictionary compatibility."""
class FrozenDict(dict):
    def _immutable(self,*args,**kwargs):
        raise TypeError('Configuration mappings are immutable; create a replacement value')
    __setitem__=__delitem__=clear=pop=popitem=setdefault=update=__ior__=_immutable
    def __copy__(self):return self
    def __deepcopy__(self,memo):
        from copy import deepcopy
        return FrozenDict((deepcopy(k,memo),deepcopy(v,memo)) for k,v in self.items())
