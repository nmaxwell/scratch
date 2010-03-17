
import sets
set = sets.Set



def gen_sig_alg(A):
    n = len(A)
    S = []
    for k in range(2**n):
        S.append(set())
    
    for k in range(2**n):
        for j in range(n):
            if bool(k & 2*j):
                S[k].add(A[j])
    return S




if __name__ == "__main__":
    
    
    A = ['a','b','c']
    
    print gen_sig_alg(A)




