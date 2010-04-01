
import sets
set = sets.Set



def power_set(A):
    n = len(A)
    S = set()
    
    for k in range(2**n):
        next = set()
        for j,a in enumerate(A):
            if bool(k & 2**j):
                next.add(a)
        S.add(next)
    
    return S

def isset(A):
    return isinstance(A, sets.BaseSet )

def issetofsets(A):
    if not isset(A):
        return false
    
    for a in A:
        if not isset(a):
            return False
    return True

def set_fmt(X):
    s = ""
    if isset(X):
        s += "{"
        for k,x in enumerate(X):
            s += set_fmt(x)
            if k < len(X)-1:
                s +=  ", "
        s +=  "}"
    
    if isinstance(X, list ):
        s += "["
        for k,x in enumerate(X):
            s += set_fmt(x)
            if k < len(X)-1:
                s +=  ", "
        s +=  "]"
    
    if s == '':
        s = str(X)
    
    return s


def gen_sig_alg(C):
    
    if issetofsets(C):
        U = set()
        for c in C:
            U = U.union(c)
        
        E = list(C)
        
        A = set()
        A.add(set())
        A.add(U.copy())
        A.add(E[0].copy())
        A.add(U - E[0])
        
        for n in range(len(E)-1):
            S = set()
            S.add( E[n+1].copy() )
            S.add( U - E[n+1] )
            for B in A:
                S.add( B.union(E[n+1]) )
                S.add( U - B.union(E[n+1]) )
                S.add( B.union(E[n+1]) )
                S.add( U - B.union(U-E[n+1]) )
            
            A = A.union(S)
            
            A = list(A)
            A = [a for a in A]
            A.sort( lambda x, y: len(x) - len(y) )
            A = set(A)
            
        return A
        
    return set()


def is_sig_alg(A):
    
    U = set()
    for B in A:
        U = U.union(B)
    
    E = list(A)
    
    for B in A:
        if not (U - B) in A:
            return false
    
    n = len(A)
    
    for k in range(2**n):
        S = set()
        for j,a in enumerate(A):
            if bool(k & 2**j):
                S = S.union(a)
        
        if not S in A:
            return false
    
    return True
    
    




if __name__ == "__main__":
    
    N = 8
    U = range(1,N+1 )
    
    A = set([ set([1]), set([2]), set([3, 4]), set([3, 4,5,6,7]), set([8]) ])
    
    A = gen_sig_alg(A)
    P = power_set(U)
    
    print '\n',len(P), len(A)
    print A.issubset(P)
    
    print '\n', set_fmt(P), '\n\n', set_fmt(A)
    
    print is_sig_alg(A)







