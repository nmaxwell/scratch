
#include <mathlib/math/std_math.h>


#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>



int main()
{
    std_setup();
    
    int x = 10;
    
    cout << "x  : " << x << endl;
    cout << "x++: " << (x++) << endl;
    
    cout << "\nx  : " << x << endl;
    cout << "++x: " << (++x) << endl;    
    
    cout << "\nx  : " << x << endl;
    cout << "x--: " << (x--) << endl;
    
    cout << "\nx  : " << x << endl;
    cout << "--x: " << (--x) << endl;
    
    cout << endl;
    
    cout << "4/5:  " << (4/5) << endl;
    cout << "5/5:  " << (5/5) << endl;
    cout << "6/5:  " << (6/5) << endl;
    
    cout << endl;
    
    for (int k=1000; k>=0; k-- )
    {
        int *p = ml_alloc<int > (k);
        
        for (int j=0; j<k; j++ )
            p[j] = 0;
        
        ml_free(p);
    }
    
    std_exit();
}




