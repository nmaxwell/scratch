
from sigma_algebra import *




if __name__ == "__main__":
    
    A = set([ set(['AxA']), set(['AxB']), set(['BxA']), set(['BxB']) ])
    
    A = gen_sig_alg(A)
    
    print '\n', len(A)
    
    print '\n\n', set_fmt(A)
    
    print is_sig_alg(A)







