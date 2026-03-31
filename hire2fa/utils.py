# ------------------------------------------------------------------------------
def argsort(arr: tuple | list) -> tuple[int]:
    """Helper function to get the indices that would sort a tuple/list"""
    if not len(arr): return ()
    idxs,_ = zip(*sorted(
        enumerate(arr), key = lambda t: int(t[1])
    ))
    return idxs

# ------------------------------------------------------------------------------
def slice_tuple(arr: tuple | list, idxs: list[int] | tuple[int]) -> tuple:
    """Helper function to slice tuples/lists by using a list of indices (similar to how it behaves in numpy)"""
    return tuple(arr[i] for i in idxs)

# ------------------------------------------------------------------------------
