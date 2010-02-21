




#include <mathlib/math/std_math.h>
#include <mathlib/math/transforms/fft.h>
#include <mathlib/math/random/ml_random.h>


#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>





/*

# FFTW_R2HC computes a real-input DFT with output in “halfcomplex” format, i.e. real and imaginary parts for a transform of size n stored as:

r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
(Logical N=n, inverse is FFTW_HC2R.)
# FFTW_HC2R computes the reverse of FFTW_R2HC, above. (Logical N=n, inverse is FFTW_R2HC.)
# FFTW_DHT computes a discrete Hartley transform. (Logical N=n, inverse is FFTW_DHT.)
# FFTW_REDFT00 computes an REDFT00 transform, i.e. a DCT-I. (Logical N=2*(n-1), inverse is FFTW_REDFT00.)
# FFTW_REDFT10 computes an REDFT10 transform, i.e. a DCT-II (sometimes called “the” DCT). (Logical N=2*n, inverse is FFTW_REDFT01.)
# FFTW_REDFT01 computes an REDFT01 transform, i.e. a DCT-III (sometimes called “the” IDCT, being the inverse of DCT-II). (Logical N=2*n, inverse is FFTW_REDFT=10.)
# FFTW_REDFT11 computes an REDFT11 transform, i.e. a DCT-IV. (Logical N=2*n, inverse is FFTW_REDFT11.)
# FFTW_RODFT00 computes an RODFT00 transform, i.e. a DST-I. (Logical N=2*(n+1), inverse is FFTW_RODFT00.)
# FFTW_RODFT10 computes an RODFT10 transform, i.e. a DST-II. (Logical N=2*n, inverse is FFTW_RODFT01.)
# FFTW_RODFT01 computes an RODFT01 transform, i.e. a DST-III. (Logical N=2*n, inverse is FFTW_RODFT=10.)
# FFTW_RODFT11 computes an RODFT11 transform, i.e. a DST-IV. (Logical N=2*n, inverse is FFTW_RODFT11.) 

*/

double err( double approx, double exact )
{
    if ( fabs(approx-exact) > 0.0 and exact != 0.0 )
        return log10(fabs(approx-exact)) - log10(fabs(exact));
    if ( fabs(approx-exact) > 0.0 and exact == 0.0 )
        return log10(fabs(approx-exact));
    else return -17;
}



int main()
{
    std_setup();
    
    ml_random rng;
    
    int n=32;
    
    fftw_r2r_mz_1d mem;
    
    double * in = ml_alloc<double> (n);
    double * out = ml_alloc<double> (n);
    double * temp = ml_alloc<double> (n);
    
    for (int k=0; k<n; k++ )
        in[k] = out[k] = temp[k] = 0;
    
    for (int k=0; k<n; k++ )
        in[k] = rng.gen_double();
    
    int plan = mem.plan( n, FFTW_REDFT00  );
    assert(!mem.execute(plan, in, temp ));
    
    plan = mem.plan( n, FFTW_REDFT00  );
    assert(!mem.execute(plan, temp, out ));
    
    for (int k=0; k<n; k++ )
        out[k] /= 2*n-2;
    
    for (int k=0; k<n; k++ )
        cout << k << "\t" << err( out[k], in[ k] ) << endl;
        
    //sprintf(fname, "/workspace/output/temp/%d", j);
    //output( out, n, fname );
    
    
    
    //int plan = mem.plan( n, FFTW_REDFT01  );
    
    
    
    
    std_exit();
}
















